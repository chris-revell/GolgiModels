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


using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
using DifferentialEquations
using SparseArrays
using UnPack
using CairoMakie 
using FromFile
using DrWatson
using SciMLOperators

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices

Nx = 51
Nν = 51
Nghost = 1           # number of ghost points on each side of the domain 
Nxplus = Nx+2*Nghost # number of discretised points including ghost points 
Nνplus = Nν+2*Nghost # number of discretised points including ghost points 
xMax = 100.0
νMax = 1.0
xs = collect(range(0.0, xMax, Nxplus))
dx = xs[2]-xs[1]
νs = collect(range(0.0, νMax, Nνplus))
dν = νs[2]-νs[1]

K₁ = 1.0
K₂ = 1.5
K₃ = 1.0
K₄ = 1.5  
α = 10.0
δ = 1.0
σ = 1.0
N = 1000
β = N*(σ*K₃ - K₂*K₄)
D = α*δ*N^2*(K₂+σ*K₃)
tMax = 0.001

#%%
fInitial(x, μx, σx, z, μz, σz) = exp(-(x-μx)^2/σx^2 - (z-μz)^2/σz^2)
μx = xMax/2.0; σx=xMax; μν=0.0; σν=0.1
hFun(x, μx, σx) = exp(-(x-μx)^2/σx^2)
μh = xMax/2.0; σh = xMax/10.0
#%%

A = makeIncidenceMatrix3D(Nxplus, 1, Nνplus)
Ā = abs.(A)
Aᵀ = transpose(A)
Aᵤₚ = dropzeros((Ā-A)./2)

nVerts = Nxplus*Nνplus       # Total number of vertices 
nEdgesi = (Nxplus-1)*Nνplus  # Number of i-directed edges (x, in this case)
nEdgesj = Nxplus*(Nνplus-1)  # Number of j-directed edges (ν, in this case)
nEdges = nEdgesi+nEdgesj     # Total number of edges over all dimensions 

# Ghost point masks
ghostVertexMask = makeGhostVertexMask((Nxplus, Nνplus))
ghostEdgeMask = makeGhostEdgeMask((Nxplus, Nνplus))

# Weights
W = vertexVolumeWeightsMatrix((Nxplus, Nνplus), (dx, dν))
W⁻¹ =  vertexVolumeWeightsInverseMatrix((Nxplus, Nνplus), (dx, dν))
l⁻¹ = edgeLengthInverseMatrix((Nxplus, Nνplus), (dx, dν))

# Matrices for picking out ν and xy directions in derivatives 
# Matrix of i-directed edge accessibility
P_i = ones(Nxplus-1, Nνplus); P_i[:, 1] .= 0.0; P_i[:, end] .= 0.0; P_i[1, :] .= 0.0; P_i[end, :] .= 0.0
# Matrix of j-directed edge accessibility  
P_j = ones(Nxplus, Nνplus-1); P_j[:, 1] .= 0.0; P_j[:, end] .= 0.0; P_j[1, :] .= 0.0; P_j[end, :] .= 0.0
P = dropzeros(spdiagm(vcat(reshape(P_i, nEdgesi), reshape(P_j, nEdgesj)))) # Diagonal sparse matrix to exclude all edges adjacent to ghost points  
Pν = dropzeros(spdiagm(vcat(zeros(nEdgesi), reshape(P_j, nEdgesj)))) # Diagonal sparse matrix to exclude all xy edges and ν edges adjacent to ghost points  
Pxy = dropzeros(spdiagm(vcat(reshape(P_i, nEdgesi), zeros(nEdgesj)))) # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
mat_h = zeros(Nxplus, Nνplus)
for i=1:Nxplus
    mat_h[i, :] .= hFun(xs[i], μh, σh).+1.0
    # mat_h[i, :] .= 1.0
end
h_vec = reshape(mat_h, nVerts)
hᵥ = spdiagm(h_vec)
hₑ = spdiagm(0.5.*Ā*h_vec); # hₑ = spdiagm(Aᵤₚ*h_vec)
# Prefactor diagonal matrix a to store 1/(1+αh) term 
a = spdiagm(1.0./(1.0 .+ α.*h_vec))
matFₑ = zeros(Nxplus)
for i=1:Nxplus
    matFₑ[i] = hFun(xs[i], xMax/2.0, xMax/10.0).+1.0
    # matFₑ[i] = 1.0
end
matE = zeros(Nxplus, Nνplus)
vecE = reshape(matE, nVerts)
E = spdiagm(vecE)

#%%

∇ₑ = l⁻¹*A      # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 
L1 = a*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ) 
L2 = D.*a*∇cdot*(hₑ*Pxy*∇ₑ)

#%%

p = (L1 = L1,
     L2 = L2,
     Nxplus = Nxplus,
     nVerts = nVerts,
     K₂ = K₂,
     matE = matE,
     vecE = vecE,
     E = E,
     matFₑ = matFₑ,
     strides = strides(matE))

#%%

# Initial conditions using Gaussian
uMat = zeros(Float64, Nxplus, Nνplus)
for νν=1:Nνplus
    for xx=1:Nxplus
        uMat[xx,νν] = fInitial(xs[xx], μx, σx, νs[νν], μν, σν)            
    end
end
u0 = reshape(uMat, nVerts)
u0[ghostVertexMask.!=true] .= 0.0
integ = sum((W*u0)[ghostVertexMask])
u0 ./= integ

#%%
function update_func!(A, u, p, t)
    for i=1:p.Nxplus
        cs = @view u[1+(i-1)*p.strides[1]:p.strides[2]:end]
        p.matE[i,:] .= p.matFₑ[i]*p.K₂/(p.K₂ + 0.5*sum(cs[1:end-1].+cs[2:end]))
    end
    p.vecE .= reshape(p.matE, p.nVerts)
    p.E .= spdiagm(p.vecE) 
    A .= p.E*L1 .+ L2
    return nothing 
end

L = MatrixOperator(E*L1.+L2, update_func! = update_func!)

#%%

prob = ODEProblem(L, u0, (0.0, tMax), p)
@time sol = solve(prob, Vern9(), saveat=tMax/100.0, progress=true); nothing

solChecks(sol, W, ghostVertexMask)

#%%

isdir(datadir("glycosylation2DSpatialEVariation4")) ? nothing : mkdir(datadir("glycosylation2DSpatialEVariation4"))

fig = Figure(size=(1000,1000))
ax = Axis3(fig[1, 1], aspect=:equal, azimuth=2.275π)
ax.xlabel = "x"
ax.ylabel = "ν"
ax.zlabel = "c"
uInternal = Observable(reshape(sol.u[1][ghostVertexMask], (Nx, Nν)))
globalmin = minimum([minimum(u[ghostVertexMask]) for u in sol.u])
globalmax = maximum([maximum(u[ghostVertexMask]) for u in sol.u])
zlims!(ax, (globalmin, globalmax))
clims = (globalmin,globalmax)
surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
record(fig, datadir("glycosylation2DSpatialEVariation4", "surface2D_4_sparse.mp4"), 1:length(sol.t); framerate=10) do i
    uInternal[] .= reshape(sol.u[i][ghostVertexMask], (Nx, Nν))
    uInternal[] = uInternal[]
end

# fig2 = Figure(size=(1000,1000))
# ax2 = Axis(fig2[1, 1])
# ax2.xlabel = "x"
# ax2.ylabel = "Fₑ"
# lines!(ax2, xs[Nghost+1:end-Nghost], matFₑ[Nghost+1:end-Nghost])
# display(fig2)
# save(datadir("matF_e.png"), fig2)