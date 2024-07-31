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
# using InvertedIndices

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices

Nx = 31
Ny = 31
Nν = 31
Nghost = 1 # number of ghost points on each side of the domain 
Nxplus = Nx+2*Nghost # number of discretised points including ghost points 
Nyplus = Ny+2*Nghost # number of discretised points including ghost points 
Nνplus = Nν+2*Nghost # number of discretised points including ghost points 
xMax = 100.0
yMax = 100.0
νMax = 1.0
xs = collect(range(0.0..xMax, Nxplus))
dx = xs[2]-xs[1]
ys = collect(range(0.0..yMax, Nyplus))
dy = ys[2]-ys[1]
νs = collect(range(0.0..νMax, Nνplus))
dν = νs[2]-νs[1]

# K₂ = 0.05
# K₄ = 0.05  
K₁ = 1.0
K₂ = 1.0
K₃ = 2.0
K₄ = 1.0  
α = 10.0
δ = 1.0
σ = 10.0
N = 100
β = N*(σ*K₃ - K₂*K₄)
D = α*δ*N^2*(K₂+σ*K₃)
tMax = 2.0

#%%

hFun(x, μx, σx, y, μy, σy) = 1.0 + exp(-(x-μx)^2/σx^2 -(y-μy)^2/σy)
μxh = xMax/2.0; σxh = xMax/10.0
u0fun(x, μx, σx, y, μy, σy, z, μz, σz) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2 - (z-μz)^2/σz^2)
μx = 0.5*xMax; σx=10.0*xMax; μy = 0.5*yMax; σy=10.0*yMax; μν=0.0; σν=0.1*νMax
fFun(x, μx, σx, y, μy, σy) = 1.0 + exp(-(x-μx)^2/σx^2 -(y-μy)^2/σy)

#%%

A = makeIncidenceMatrix3D(Nxplus,Nyplus,Nνplus)
Ā = abs.(A)
Aᵀ = transpose(A)
Aᵤₚ = dropzeros((Ā-A).÷2)

nVerts = Nxplus*Nyplus*Nνplus       # Total number of vertices 
nEdgesi = (Nxplus-1)*Nyplus*Nνplus  # Number of i-directed edges (x, in this case)
nEdgesj = Nxplus*(Nyplus-1)*Nνplus  # Number of j-directed edges (y, in this case)
nEdgesk = Nxplus*Nyplus*(Nνplus-1)  # Number of k-directed edges (ν, in this case)
nEdges = nEdgesi+nEdgesj+nEdgesk    # Total number of edges over all dimensions 

# Ghost point masks
ghostVertexMask = makeGhostVertexMask((Nxplus, Nyplus, Nνplus))
ghostVertexMaskSparse = spdiagm(ghostVertexMask)
ghostEdgeMask = makeGhostEdgeMask((Nxplus, Nyplus, Nνplus))
ghostEdgeMaskSparse = spdiagm(ghostEdgeMask)

# Weights
W = vertexVolumeWeightsMatrix((Nxplus, Nyplus, Nνplus), (dx, dy, dν))
W⁻¹ =  vertexVolumeWeightsInverseMatrix((Nxplus, Nyplus, Nνplus), (dx, dy, dν))
l⁻¹ = edgeLengthInverseMatrix((Nxplus, Nyplus, Nνplus), (dx, dy, dν))

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
mat_h = zeros(Nxplus, Nyplus, Nνplus)
for j=1:Nyplus, i=1:Nxplus
    mat_h[i, j, :] .= hFun(xs[i], μxh, σxh, ys[j], μyh, σyh)
end
hᵥ_vec = reshape(mat_h, nVerts)
hₑ_vec = 0.5.*Ā*hᵥ_vec
hᵥ = spdiagm(hᵥ_vec)
hₑ = spdiagm(hₑ_vec)
hₑ = spdiagm(hₑ_vec)
aᵥ = spdiagm(1.0./(1.0 .+ α.*hᵥ_vec))
aₑ = spdiagm(1.0./(1.0 .+ α.*hₑ_vec))

# Velocity field 
V_i = fill(0.0, (Nxplus-1, Nyplus, Nνplus))
V_j = fill(0.0, (Nxplus, Nyplus-1, Nνplus))
V_k = fill(β, (Nxplus, Nyplus, Nνplus-1))
Vvec = vcat(reshape(V_i, nEdgesi), reshape(V_j, nEdgesj), reshape(V_k, nEdgesk))
V = ghostEdgeMaskSparse*spdiagm(Vvec)*aₑ   # Diagonal matrix of advection velocities at each edge

# Diffusivity field over edges 
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
D_i = fill(dy*dν*K₂*K₄, (Nxplus-1, Nyplus, Nνplus))
D_j = fill(dx*dν*K₂*K₄, (Nxplus, Nyplus-1, Nνplus))
D_k = fill(dx*dy*K₂*K₄, (Nxplus, Nyplus, Nνplus-1))
Dvec = vcat(reshape(D_i, nEdgesi), reshape(D_j, nEdgesj), reshape(D_k, nEdgesk))
D = ghostEdgeMaskSparse*spdiagm(Dvec)*aₑ # Diagonal matrix of advection velocities at each edge

# Matrices for picking out ν and xy directions in derivatives 
# Matrix of i-directed edge accessibility
P_i_mat = ones(Nxplus-1, Nyplus, Nνplus)
# Matrix of j-directed edge accessibility  
P_j_mat = ones(Nxplus, Nyplus-1, Nνplus)
# Matrix of k-directed edge accessibility  
# P_k_mat = fill(dx*dy*dν, (Nxplus, Nyplus, Nνplus-1))
P_k_mat = ones(Nxplus, Nyplus, Nνplus-1)
P = ghostEdgeMaskSparse*spdiagm(vcat(reshape(P_i_mat, nEdgesi), reshape(P_j_mat, nEdgesj), reshape(P_k_mat, nEdgesk))) # Diagonal sparse matrix to exclude all edges adjacent to ghost points  
Pν = ghostEdgeMaskSparse*spdiagm(vcat(zeros(nEdgesi), zeros(nEdgesj), reshape(P_k_mat, nEdgesk))) # Diagonal sparse matrix to exclude all xy edges and ν edges adjacent to ghost points  
Pxy = ghostEdgeMaskSparse*spdiagm(vcat(reshape(P_i_mat, nEdgesi), reshape(P_j_mat, nEdgesj), zeros(nEdgesk))) # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 


# Initial conditions using Gaussian
uMat = zeros(Float64, Nxplus, Nyplus, Nνplus)
for νν=1:Nνplus, yy=1:Nyplus, xx=1:Nxplus
    uMat[xx,yy,νν] = u0fun(xs[xx], μx, σx, ys[νν], μy, σy, νs[νν], μν, σν)            
end
u0 = reshape(uMat, nVerts)
u0[ghostVertexMask.!=true] .= 0.0
integ = sum(W*u0)
u0 ./= integ


∇ₑ = l⁻¹*A      # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

matFₑ = zeros(Nxplus, Nyplus, Nνplus)
for j=1:Nyplus, i=1:Nxplus
    matFₑ[i, j, :] .= fFun(xs[i], μxh, σxh, ys[j], μyh, σyh)
    # matFₑ[i] = 1.0
end
matE = zeros(Nxplus, Nyplus, Nνplus)
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
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + D.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ


p = (L1 = L1,
     L2 = L2,
     Nxplus = Nxplus,
     Nyplus = Nyplus,
     K₂ = K₂,
     matE = matE,
     vecE = vecE,
     E = E,
     matFₑ = matFₑ,
     strides = strides(matE))

function update_func!(L, u, p, t)
    @unpack L1,
        L2,
        Nxplus,
        Nyplus,
        K₂,
        matE,
        vecE,
        E,
        matFₑ,
        strides = p

    for j=1:Nyplus, i=1:Nxplus
        cs = @view u[(1+(i-1)*strides[1] + (j-1)*strides[2]):strides[3]:end]
        matE[i,j,:] .= matFₑ[i,j]*K₂/(K₂ + 0.5*sum(cs[1:end-1].+cs[2:end]))
    end
    vecE .= reshape(matE, nVerts)
    E .= spdiagm(vecE) 

    L .= E*L1 .+ L2
end

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
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + D.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ

# L = -W⁻¹*Aᵀ*D*l⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 
L1 = aᵥ*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
L2 = aᵥ*∇cdot*(D*hₑ*Pxy*∇ₑ)
L = MatrixOperator(E*L1.+L2, update_func! = update_func!)
prob = ODEProblem(L, u0, (0.0, tMax), p)
# sol = solve(prob, Trapezoid(), saveat=tMax/100.0)
sol = solve(prob, Vern9(), saveat=tMax/100.0)

#%%

isdir(datadir("glycosylationEvariation3D")) ? nothing : mkdir(datadir("glycosylationEvariation3D"))

fig = Figure(size=(1000,1000))
ax = Axis3(fig[1, 1], aspect=:equal, azimuth=2.275π)
ax.xlabel = "x"
ax.ylabel = "ν"
ax.zlabel = "c"
uInternal3D = reshape(sol.u[1][ghostVertexMask], (Nx, Ny, Nν))
uInternal = Observable(zeros(Nx, Nν))
globalmin = minimum([minimum(u[ghostVertexMask]) for u in sol.u])
globalmax = maximum([maximum(u[ghostVertexMask]) for u in sol.u])
zlims!(ax, (globalmin, globalmax))
clims = (globalmin,globalmax)
surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
record(fig, datadir("glycosylationEvariation3D","c_against_xν_at_halfy.mp4"), 1:length(sol.t); framerate=10) do i
    uInternal3D .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))
    uInternal[] .= uInternal3D[:,Ny÷2,:]
    uInternal[] = uInternal[]
end

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "M, ∱cdxdy"
uInternal3D = reshape((W*sol.u[end])[ghostVertexMask], (Nx, Ny, Nν))
M = sum(uInternal3D, dims=1)
M = sum(M, dims=2)
lines!(ax, νs[1:Nghost:end-2*Nghost], M[1,1,:])
ax.title = "Integral of c over x and y against ν at final time"
save(datadir("glycosylationEvariation3D", "finalνVsM.png"), fig)


fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "M, ∱cdxdy"
uInternal3D = reshape((W*sol.u[end])[ghostVertexMask], (Nx, Ny, Nν))
Mtmp = sum(uInternal3D, dims=1)
minima = Float64[]
maxima = Float64[]
for i=1:length(sol.t)
    uInternal3D .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))
    Mtmp .= sum(uInternal3D, dims=1)
    push!(minima, minimum(sum(Mtmp, dims=2)))
    push!(maxima, maximum(sum(Mtmp, dims=2)))
end
globalmin = minimum(minima)
globalmax = maximum(maxima)
M = Observable(sum(Mtmp, dims=2)[1,1,:])
lines!(ax, νs[1:Nghost:end-2*Nghost], M)
ax.title = "Integral of c over x and y against ν at final time"
ylims!(ax, (globalmin, globalmax))
record(fig, datadir("glycosylationEvariation3D", "Mvsν.mp4"), 1:length(sol.t); framerate=10) do i
    uInternal3D .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))
    uInternal[] .= uInternal3D[:,Ny÷2,:]
    uInternal[] = uInternal[]
    Mtmp .= sum(uInternal3D, dims=1)
    M[] .= sum(Mtmp, dims=2)[1,1,:]
    M[] = M[]
end

# #%%

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
# record(fig, datadir("glycosylationEvariation3D", "LargeNuTimeScan.mp4"), 1:length(sol.t); framerate=10) do i
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
# record(fig, datadir("glycosylationEvariation3D", "NuProfileAtxyOverTime.mp4"), 1:length(sol.t); framerate=10) do i
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
# record(fig, datadir("glycosylationEvariation3D", "xνOverTimeAty.mp4"), 1:length(sol.t); framerate=10) do i
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