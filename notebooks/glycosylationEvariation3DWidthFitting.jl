#%%
# flux_νₑ = (diffusive_flux_ν + advective_flux_ν)
# flux_νₑ = K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*cₑ    where cᵥ is concentration over vertices, cₑ is concentration over edges 
# cₑ = Aᵤₚ*cᵥ
# flux_νₑ = (K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ)*cᵥ
# flux_xyₑ = Dₑ*hₑ*diffusive_flux_xy
# flux_xyₑ = Dₑ*hₑ*Pxy*∇ₑ*cᵥ
# ċ = aE∇⋅flux_νₑ + a∇⋅flux_xyₑ
# ċ = a*E*∇⋅(K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*Aᵤₚ*cᵥ) + a∇⋅(Dₑ*hₑ*Pxy*∇ₑ*cᵥ)
# Dₑ constant over edges 
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝓓.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ

# L = -W⁻¹*Aᵀ*𝓓*l⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 



#


using DifferentialEquations
using SparseArrays
using UnPack
using CairoMakie 
using HCubature
using FromFile
using DrWatson
using Printf
using LaTeXStrings
using SciMLOperators
using DataFrames
using XLSX
using Interpolations
using Dates

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices

simpsonsRule(xs) = sum(xs[1:end-1].+xs[2:end])./2.0

function productionTotalM(u, W, ghostVertexMask, dims, ϕ)
    uInternal = reshape((W*u)[ghostVertexMask], dims)
    return sum(selectdim(uInternal, 1, round(Int, ϕ*dims[1])))
end

df = DataFrame(XLSX.readtable(datadir("exp_raw", "ModellingData.xlsx"), 1))

cisternaSeriesID = 1

hs = filter(x->x.series_ID==cisternaSeriesID, df)[1,5:end]
hs = collect(skipmissing(hs))


# PDE discretisation parameters 
Nx = 101             # Number of discretisation points in space
Ny = 101             # Number of discretisation points in space
Nν = 101             # Number of discretisation points in polymerisation space
Nghost = 1           # Number of ghost points on each side of the domain 
Nνplus = Nν+2*Nghost # Number of discretised points including ghost points 
Nxplus = Nx+2*Nghost # Number of discretised points including ghost points 
Nyplus = Ny+2*Nghost # Number of discretised points including ghost points 

xMax = (length(hs)-1)   
xs   = collect(range(0.0..xMax, Nxplus)) # Positions of discretised vertices in space
dx   = xs[2]-xs[1]
νMax = 1.0
νs   = collect(range(0.0..νMax, Nνplus)) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]


# Basic parameters: geometry
Ω = 1.0         # Lumen volume
Ωperp = 1000.0  # Lumen footprint area
N = 100         # Maximum polymer length 

# Basic parameters: rate constants
k_Cd = 1.0 # Complex desorption rate
k_Ca = 1.0 # Complex adsorption rate
k_Sd = 1.0 # Substrate desorption rate
k_Sa = 1.0 # Substrate adsorption rate
k₁ = 1.0   # Complex formation forward reaction rate 
k₂ = 1.0   # Complex dissociation reverse reaction rate 
k₃ = 1.0   # Product formation
k₄ = 1.0   # Product dissociation 

# Basic parameters: concentrations 
C_b = 1.0  # Initial bulk monomer concentration
S_b = 1.0  # Initial bulk substrate concentration
S_0 = 1.0  # Early surface substrate concentration 
E_0 = 0.01 # Mean enzyme concentration

# Basic parameters: diffusivities
D_C = 1.0  # Monomer/polymer diffusivity
D_S = 1.0  # Substrate diffusivity

# Basic parameters: Timescale 
Tᵣ⁰ = 1.0  # Release time

# Derived quantities: rates
α_C = (k_Cd*Ω)/(2*k_Ca*Ωperp) # Balance of complex in bulk to complex on membrane
α_S = (k_Sd*Ω)/(2*k_Sa*Ωperp) # Balance of substrate in bulk to substrate on membrane
K₂  = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω)) # Non-dimensionalised complex formation net reaction rate
K₃  = k₃/k₁    # Non-dimensionalised product formation rate
K₄  = k₄/k₁    # Non-dimensionalised prodict dissociation rate

# Derived quantities: masses and concentrations 
h₀  = Ω/Ωperp                   # Mean thickness 
L₀  = sqrt(π)*Ω / (Ωperp)^(1.5) # Mean radius 
C_0 = C_b*h₀/(2*(1+α_C))        # Early surface monomer concentration
𝓒   = C_b*Ω                     # Initial monomer mass
σ   = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
𝓢   = S_b*Ω                     # Initial substrate mass
𝓔   = 2*E_0*Ωperp               # Total enzyme mass
ε   = 𝓔*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ωperp)

β = N*(σ*K₃ - K₂*K₄)

# Derived quantities: diffusivities
δ_C = π*D_C/(k₁*𝓔)
δ_S = π*D_S/(k₁*𝓔)
𝓓   = α_C*δ_C*N^2*(K₂ + σ*K₃)

# Derived quantities: non-dimensionalised time
Tᵣ  = k₁*𝓔*Tᵣ⁰/(2*Ωperp)

# T̃ᵣ $ (\ref{eq:ttr})

# K₂ = 1.0
# K₃ = 1.0
# K₄ = 1.0  
# α_C = 100.0
# δ_C = 1.0
# σ = 10.0
# N = 100
# β = N*(σ*K₃ - K₂*K₄)
# 𝓓 = α_C*δ_C*N^2*(K₂+σ*K₃)
# tMax = 60.0

# Create directory for run data labelled with current time.
paramsName = @savename cisternaSeriesID K₂ K₃ K₄ α_C δ_C σ N Tᵣ
folderName = "$(paramsName)_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
# Create frames subdirectory to store system state at each output time
subFolder = "hFitting"
mkpath(datadir("sims",subFolder,folderName))

#%%

itp_cubic = cubic_spline_interpolation(0:xMax, hs)
hFun(x) = itp_cubic(x)

u0fun(x, μx, σx, y, μy, σy) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2)
μνu0 = 0.0; σνu0 = νMax/10.0
μxu0 = xMax/2.0; σxu0 = 10.0*xMax

fFun(x, μx, σx) = 0.1 #+ exp(-(x-μx)^2/σx^2)
μxF = xMax/2.0; σxF=xMax/10.0

#%%

A = makeIncidenceMatrix3D(Nνplus, Nxplus,1)
Ā = abs.(A)
Aᵀ = transpose(A)
Aᵤₚ = dropzeros((Ā-A).÷2)

nVerts = Nνplus*Nxplus*Nyplus       # Total number of vertices 
nEdgesi = (Nνplus-1)*Nxplus*Nyplus  # Number of i-directed edges (ν, in this case)
nEdgesj = Nνplus*(Nxplus-1)*Nyplus  # Number of j-directed edges (x, in this case)
nEdgesk = Nνplus*Nxplus*(Nyplus-1)  # Number of k-directed edges (y, in this case)
nEdges = nEdgesi+nEdgesj+nEdgesk    # Total number of edges over all dimensions 

# Ghost point masks
ghostVertexMask = makeGhostVertexMask((Nνplus, Nxplus))
ghostVertexMaskSparse = spdiagm(ghostVertexMask)
ghostEdgeMask = makeGhostEdgeMask((Nνplus, Nxplus))
ghostEdgeMaskSparse = spdiagm(ghostEdgeMask)

# Weights
W = vertexVolumeWeightsMatrix((Nνplus, Nxplus), (dν, dx))
W⁻¹ =  vertexVolumeWeightsInverseMatrix((Nνplus, Nxplus), (dν, dx))
l⁻¹ = edgeLengthInverseMatrix((Nνplus, Nxplus), (dν, dx))

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
mat_h = zeros(Nνplus, Nxplus, Nyplus)
for j=1:Nxplus
    # mat_h[:, j] .= hFun(xs[j])
    selectdim(mat_h, 2, j) .= hFun(xs[j])
end
hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
aᵥ = spdiagm(1.0./(1.0 .+ α_c.*hᵥ_vec)) # Prefactor 1/(1+α_c*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience
aₑ = spdiagm(1.0./(1.0 .+ α_c.*hₑ_vec)) # Prefactor 1/(1+α_c*hₑ(x)) evaluated over edges, packaged into a sparse diagonal matrix for convenience

# Velocity field 
V_i = fill(β, (Nνplus-1, Nxplus, Nyplus))
V_j = fill(0.0, (Nνplus, Nxplus-1, Nyplus))
V_k = fill(0.0, (Nνplus, Nxplus, Nyplus-1))
Vvec = vcat(reshape(V_i, nEdgesi), reshape(V_j, nEdgesj), reshape(V_k, nEdgesk))
V = ghostEdgeMaskSparse*spdiagm(Vvec)*aₑ   # Diagonal matrix of advection velocities at each edge

# Diffusivity field over edges 
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
D_i = fill(dx*dy*K₂*K₄, (Nνplus-1, Nxplus, Nyplus))
D_j = fill(dν*dy*K₂*K₄, (Nνplus, Nxplus-1, Nyplus))
D_k = fill(dν*dx*K₂*K₄, (Nνplus, Nxplus, Nyplus-1))
Dvec = vcat(reshape(D_i, nEdgesi), reshape(D_j, nEdgesj), reshape(D_k, nEdgesk))
𝒟 = ghostEdgeMaskSparse*spdiagm(Dvec)*aₑ # Diagonal matrix of advection velocities at each edge

# Matrices for picking out ν and xy directions in derivatives 
# P = ghostEdgeMaskSparse*spdiagm(vcat(reshape(P_i_mat, nEdgesi), reshape(P_j_mat, nEdgesk))) # Diagonal sparse matrix to exclude all edges adjacent to ghost points  
P = ghostEdgeMaskSparse*spdiagm(vcat(ones(Int64, nEdgesi), ones(Int64, nEdgesj), ones(Int64, nEdgesk)))     # Diagonal sparse matrix to exclude all edges adjacent to ghost points  
Pν = ghostEdgeMaskSparse*spdiagm(vcat(ones(Int64, nEdgesi), zeros(Int64, nEdgesj), zeros(Int64, nEdgesk)))   # Diagonal sparse matrix to exclude all xy edges and ν edges adjacent to ghost points  
# Px = ghostEdgeMaskSparse*spdiagm(vcat(zeros(Int64, nEdgesi), ones(Int64, nEdgesj), ones(Int64, nEdgesk)))           # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 
Pxy = ghostEdgeMaskSparse*spdiagm(vcat(zeros(Int64, nEdgesi), ones(Int64, nEdgesj), ones(Int64, nEdgesk)))           # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 

# Diagonal matrix of edge lengths
l_i = fill(dx, (Nνplus-1, Nxplus, Nyplus))
l_j = fill(dy, (Nνplus, Nxplus-1, Nyplus))
l_k = fill(dν, (Nνplus, Nxplus, Nyplus-1))
lvec = vcat(reshape(l_i, nEdgesi), reshape(l_j, nEdgesj), reshape(l_k, nEdgesk))
l = spdiagm(lvec)
l⁻¹ = spdiagm(1.0./lvec)

# Initial conditions using Gaussian
uMat = zeros(Float64, Nνplus, Nxplus, Nyplus)
for yy=1:Nyplus, xx=1:Nxplus, νν=1:Nνplus
    uMat[νν, xx, yy] = u0fun(νs[νν], μνu0, σνu0, xs[xx], μxu0, σxu0)            
end
u0 = reshape(uMat, nVerts)
u0[ghostVertexMask.!=true] .= 0.0
integ = sum(W*u0)
u0 ./= integ

∇ₑ = l⁻¹*A       # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

matFₑ = zeros(Nνplus, Nxplus, Nyplus)
for j=1:Nxplus
    # matFₑ[:, j] .= fFun(xs[j], μxF, σxF)
    selectdim(matFₑ, 2, j) .= fFun(xs[j], μxF, σxF)
    # matFₑ[i] = 1.0
end
matE = zeros(Nνplus, Nxplus, Nyplus)
vecE = reshape(matE, nVerts)
E = spdiagm(vecE)


# Cνν = W⁻¹*Aᵀ*Pν*l⁻¹*A
# Cν = Aᵀ*l⁻¹*Pν*Aᵤₚ
# flux_νₑ = (diffusive_flux_ν + advective_flux_ν)
# flux_νₑ = K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*cₑ    where cᵥ is concentration over vertices, cₑ is concentration over edges 
# cₑ = Aᵤₚ*cᵥ
# flux_νₑ = (K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ)*cᵥ
# flux_xyₑ = Dₑ*hₑ*diffusive_flux_xy
# flux_xyₑ = Dₑ*hₑ*Pxy*∇ₑ*cᵥ
# ċ = aE∇⋅flux_νₑ + a∇⋅flux_xyₑ
# ċ = a*E*∇⋅(K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*Aᵤₚ*cᵥ) + a∇⋅(Dₑ*hₑ*Pxy*∇ₑ*cᵥ)
# Dₑ constant over edges 
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝒟.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ

L1 = aᵥ*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
L2 = aᵥ*∇cdot*(𝒟*hₑ*Pxy*∇ₑ)

p = (L1 = L1,
    L2 = L2,
    Nνplus = Nνplus,
    Nxplus = Nxplus,
    Nyplus = Nyplus,
    K₂ = K₂,
    matE = matE,
    E = E,
    matFₑ = matFₑ)

function update_func!(L, u, p, t)
    @unpack L1,
        L2,
        Nνplus,
        Nxplus,
        Nyplus,
        K₂,
        matE,
        E,
        matFₑ = p

    cs = reshape(u, (Nνplus, Nxplus, Nyplus))     
    for k = 1:Nyplus
        for j= 1:Nxplus
            integrationFactor = K₂/(K₂ + simpsonsRule(cs[:,j,k]))
            matE[:,j,k] .= matFₑ[:,j,k].*integrationFactor
        end
    end
    E .= spdiagm(reshape(matE, nVerts)) 
    L .= E*L1 .+ L2
end

L = MatrixOperator(E*L1.+L2, update_func! = update_func!)
prob = ODEProblem(L, u0, (0.0, tMax), p)
sol = solve(prob, Vern9(), saveat=tMax/100.0)

#%%

isdir(datadir("sims", "hFitting")) ? nothing : mkdir(datadir("sims", "hFitting"))

fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax = CairoMakie.Axis3(fig[1, 1])
ax.xlabel = "ν"
ax.ylabel = "M, ∱cdxdy"
uInternal2D = reshape((W*sol.u[end])[ghostVertexMask], (Nν, Nx))
M = sum(uInternal2D, dims=1)[1,:]
minima = Float64[]
maxima = Float64[]
for i=1:length(sol.t)
    uInternal2D .= reshape(sol.u[i][ghostVertexMask], (Nν, Nx))
    M .= sum(uInternal2D, dims=2)[:,1]
    push!(minima, minimum(M))
    push!(maxima, maximum(M))
end
globalmin = minimum(minima)
globalmax = maximum(maxima)

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "M, ∱cdxdy"
ax.title = "Integral of Cₛ over x against ν"
M = Observable(zeros(Nν))
lines!(ax, νs[1:Nghost:end-2*Nghost], M)
ylims!(ax, (globalmin, globalmax))
record(fig, datadir("sims",subFolder, folderName, "Mvsν.mp4"), 1:length(sol.t); framerate=10) do i
    uInternal2D .= reshape(sol.u[i][ghostVertexMask], (Nν, Nx))
    M[] .= sum(uInternal2D, dims=2)[:,1]
    M[] = M[]
end
save(datadir("sims",subFolder,folderName,"finalνVsM.png"), fig)

function productionTotalM(u, W, ghostVertexMask, dims, ϕ)
    uInternal = reshape((W*u)[ghostVertexMask], dims)
    return sum(selectdim(uInternal, 1, round(Int, ϕ*dims[1])))
end

# # Matrices for picking out ν and xy directions in derivatives 
# P_i_mat = ones(Int64, Nνplus-1, Nxplus, Nyplus) # Matrix of i-directed edge accessibility
# P_j_mat = ones(Int64, Nνplus, Nxplus-1, Nyplus) # Matrix of j-directed edge accessibility
# P_k_mat = ones(Int64, Nνplus, Nxplus, Nyplus-1) # Matrix of k-directed edge accessibility
# for m in [P_i_mat, P_j_mat, P_k_mat]
#     for i=1:ndims(m)
#         if size(m, i)>=1
#             selectdim(m, i, 1) .= 0
#             selectdim(m, i, size(m, i)) .= 0
#         end
#     end 
# end
# # Matrices for picking out ν and xy directions in derivatives 
# # Pν = spdiagm(vcat(zeros(nEdgesi), zeros(nEdgesj), ones(nEdgesk)))
# Pν = spdiagm(vcat(reshape(P_i_mat, nEdgesi), zeros(nEdgesj), ones(nEdgesk)))
# # Pxy = spdiagm(vcat(ones(nEdgesi), ones(nEdgesj), zeros(nEdgesk)))
# Pxy = spdiagm(vcat(ones(nEdgesi), reshape(P_j_mat, nEdgesj), reshape(P_k_mat, nEdgesk)))



# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "x"
# ax.ylabel = "y"
# uInternal3D = reshape(sol.u[1][ghostVertexMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nx,Ny))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# clims = (globalmin,globalmax)
# clims = (minimum(uInternal3D), maximum(uInternal3D))
# heatmap!(ax, xs[Nghost+1:end-Nghost], ys[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
# record(fig, datadir("sims",folderName, "LargeNuTimeScan.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal3D .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))
#     uInternal[] .= uInternal3D[:,:,end]
#     uInternal[] = uInternal[]
#     ax.title = "xy profile of ν=1.0 at t=$(sol.t[i])"
# end

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "c"
# uInternal3D = reshape(sol.u[1][ghostVertexMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nν))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# ylims!(ax, (globalmin, globalmax))
# lines!(ax, νs[Nghost+1:end-Nghost], uInternal)
# ax.title = "c against ν at x=0.5, y=0.5 at t=0.0"
# record(fig, datadir("sims",folderName, "NuProfileAtxyOverTime.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal3D .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))
#     uInternal[] .= uInternal3D[Nx÷2,Nx÷2,:]
#     uInternal[] = uInternal[]
#     ax.title = "c against ν at x=0.5, y=0.5 at t=$(sol.t[i])"
# end

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "x"
# ax.ylabel = "ν"
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# uInternal = Observable(zeros(Nx,Nν))
# clims = (globalmin,globalmax)
# ax.title = "c against x and ν at y=0.5 at t=0.0"
# heatmap!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
# record(fig, datadir("sims",folderName, "xνOverTimeAty.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal[] .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))[:,Nyplus÷2,:]
#     uInternal[] = uInternal[]
#     ax.title = "c against x and ν at y=0.5 at t=$(sol.t[i])"
# end

# #%%
# p = L
# function model!(du, u, p, t)
#     du .= p*u
# end
# prob = ODEProblem(model!, u0, (0.0,tMax), p)