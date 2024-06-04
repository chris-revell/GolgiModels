using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using DifferentialEquations
using SparseArrays
using UnPack
using CairoMakie 
# using HCubature
using FromFile
using DrWatson
# using Printf
# using LaTeXStrings
using GaussianRandomFields
# using TerminalLoggers
# using ProgressLogging
using SciMLOperators

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix

Nx = 21
# Ny = 11
Nν = 21
Nghost = 1           # number of ghost points on each side of the domain 
Nxplus = Nx+2*Nghost # number of discretised points including ghost points 
# Nyplus = Ny+2*Nghost # number of discretised points including ghost points 
Nνplus = Nν+2*Nghost # number of discretised points including ghost points 
xMax = 100.0#; yMax = 100.0;
νMax = 1.0
xs = collect(range(0.0, xMax, Nxplus))
dx = xs[2]-xs[1]
# ys = collect(range(0.0, yMax, Nyplus))
# dy = ys[2]-ys[1]
νs = collect(range(0.0, νMax, Nνplus))
dν = νs[2]-νs[1]

K₁ = 1.0
K₂ = 0.9
K₃ = 1.0
K₄ = 0.9  
α = 1.0
δ = 1.0
σ = 1.0
N = 100
# β = 1.0
β = N*(σ*K₃ - K₂*K₄)
D = α*δ*N^2*(K₂+σ*K₃)
tMax = 1.0

# vx = 0.0            # Scalar advection speed in x
# vy = 0.0            # Scalar advection speed in y
# vν = β/(1+α*h)      # Scalar advection speed in ν
# Dx = K₂*K₄/(1+α*h)  # Diffusivity in x
# Dy = K₂*K₄/(1+α*h)  # Diffusivity in x
# Dν = K₂*K₄/(1+α*h)  # Diffusivity in ν


#%%

hFun(x, μx, σx) = exp(-(x-μx)^2/σx^2)
# hFun(x, xhz, y, yhz) = 1+0.5*sin(2π*(x*xhz+y*yhz))

#%%

A = makeIncidenceMatrix3D(Nxplus, 1, Nνplus)
Ā = abs.(A)
Aᵀ = transpose(A)
Aᵤₚ = dropzeros((Ā-A)./2)

nVerts = Nxplus*Nνplus       # Total number of vertices 
nEdgesi = (Nxplus-1)*Nνplus  # Number of i-directed edges (x, in this case)
nEdgesj = Nxplus*(Nνplus-1)  # Number of j-directed edges (ν, in this case)
nEdges = nEdgesi+nEdgesj     # Total number of edges over all dimensions 

# Ghost point mask is a 1D vector in which the value at component i 
# is true if vertex i in the flattened vector of vertices is an internal vertex
# but false if vertex i is a ghost vertex in the flattened vector of vertices 
# This can be used to exclude ghost points from calculations over the whole state vector
ghostMaskArray = fill(true, (Nxplus, Nνplus))
ghostMaskArray[:, 1] .= false
ghostMaskArray[:, end] .= false
ghostMaskArray[1, :] .= false
ghostMaskArray[end, :] .= false
ghostMask = reshape(ghostMaskArray, nVerts)

# Vertex weights
# Forming a diagonal matrix of volumes around each vertex, divided by 2 at the periphery
matW = fill(dx*dν, (Nxplus, Nνplus))
matW[:, 1] ./= 2.0
matW[:, end] ./= 2.0
matW[1, :] ./= 2.0
matW[end, :] ./= 2.0
vecW = reshape(matW, nVerts)
W = spdiagm(vecW)
W⁻¹ = spdiagm(1.0./vecW)

# Diagonal matrix of volumes around each edge, divided by 2 at the periphery
# Matrix of i-directed edge weights  
F_i = fill(dx*dν, (Nxplus-1, Nνplus))
F_i[:, 1] ./= 2.0
F_i[:, end] ./= 2.0
F_i[1, :] ./= 2.0
F_i[end, :] ./= 2.0
# Matrix of j-directed edge weights  
F_j = fill(dx*dν, (Nxplus, Nνplus-1))
F_j[:, :, 1] ./= 2.0
F_j[:, :, end] ./= 2.0
F_j[:, 1, :] ./= 2.0
F_j[:, end, :] ./= 2.0
F_j[1, :, :] ./= 2.0
F_j[end, :, :] ./= 2.0
# Matrix of k-directed edge weights  
Fvec = vcat(reshape(F_i, nEdgesi), reshape(F_j, nEdgesj))
F = spdiagm(Fvec)
F⁻¹ = spdiagm(1.0./Fvec)

# Diagonal matrix of edge lengths
l_i = fill(dx, (Nxplus-1, Nνplus))
l_j = fill(dν, (Nxplus, Nνplus-1))
lvec = vcat(reshape(l_i, nEdgesi), reshape(l_j, nEdgesj))
l = spdiagm(lvec)
l⁻¹ = spdiagm(1.0./lvec)

# Matrices for picking out ν and xy directions in derivatives 
# Matrix of i-directed edge accessibility
P_i = ones(Nxplus-1, Nνplus)
P_i[:, 1] .= 0.0
P_i[:, end] .= 0.0
P_i[1, :] .= 0.0
P_i[end, :] .= 0.0
# Matrix of j-directed edge accessibility  
P_j = ones(Nxplus, Nνplus-1)
P_j[:, 1] .= 0.0
P_j[:, end] .= 0.0
P_j[1, :] .= 0.0
P_j[end, :] .= 0.0
P = dropzeros(spdiagm(vcat(reshape(P_i, nEdgesi), reshape(P_j, nEdgesj)))) # Diagonal sparse matrix to exclude all edges adjacent to ghost points  
Pν = dropzeros(spdiagm(vcat(zeros(nEdgesi), reshape(P_j, nEdgesj)))) # Diagonal sparse matrix to exclude all xy edges and ν edges adjacent to ghost points  
Pxy = dropzeros(spdiagm(vcat(reshape(P_i, nEdgesi), zeros(nEdgesj)))) # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
mat_h = zeros(Nxplus, Nνplus)
for i=1:Nxplus
    # mat_h[i, :] .= hFun(xs[i], 0.5, 0.1)
    mat_h[i, :] .= 1.0
end
h_vec = reshape(mat_h, nVerts)
hᵥ = spdiagm(h_vec)
hₑ = spdiagm(0.5.*Ā*h_vec)
# hₑ = spdiagm(Aᵤₚ*h_vec)

# Prefactor diagonal matrix a to store 1/(1+αh) term 
a = spdiagm(1.0./(1.0 .+ α.*h_vec))

# Fₑ field over all vertices, for enxyme concentration calculation, derived from Gaussian random field 
# cov = CovarianceFunction(2, Exponential(.5))
# pts = range(0, stop=1, length=1001)
# grf = GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=2001)
# matFₑ = sample(grf)[10:9+Nxplus, 10:9+Nxplus][:,1]
# matFₑ .= matFₑ.-minimum(matFₑ)
matFₑ = zeros(Nxplus)
for i=1:Nxplus
    matFₑ[i] = hFun(xs[i], xMax/2.0, xMax/10.0).+1.0
    # mat_h[i, j, :] .= 1.0
end


# function to find enzyme concentration at each vertex by integrating over ν
# function enzymeE!(Niplus, Njplus, c, K₂, matE, vecE, E, Fₑ, stride)
#     # Fₑ(x)*K₂/(K₂ + ∱c(x,ν,t)dν)
#     for i=1:Niplus
#         for j=1:Njplus
#             cs = @view c[(1+(i-1)*stride[1] + (j-1)*stride[2]):stride[3]:end]
#             matE[i,j,:] .= Fₑ[i,j]*K₂/(K₂ + 0.5*sum(cs[1:end-1].+cs[2:end]))
#         end
#     end
#     vecE .= reshape(matE, nVerts)
#     E .= spdiagm(vecE) 
#     return nothing 
# end

matE = zeros(Nxplus, Nνplus)
vecE = reshape(matE, nVerts)
E = spdiagm(vecE)

∇ₑ = l⁻¹*A      # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

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

L1 = a*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ) 
L2 = D.*a*∇cdot*(hₑ*Pxy*∇ₑ)

p = (L1 = Matrix(L1),
     L2 = Matrix(L2),
     Nxplus = Nxplus,
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
        K₂,
        matE,
        vecE,
        E,
        matFₑ,
        strides = p

    for i=1:Nxplus
        cs = @view u[1+(i-1)*strides[1]:strides[2]:end]
        matE[i,:] .= matFₑ[i]*K₂/(K₂ + 0.5*sum(cs[1:end-1].+cs[2:end]))
    end
    vecE .= reshape(matE, nVerts)
    E .= spdiagm(vecE) 

    # L .= E*L1 .+ L2
    # L .= vecE.*L1 .+ L2
    L .= L1 .+ L2
end


#%%

# Initial conditions using Gaussian
fInitial(x, μx, σx, z, μz, σz) = exp(-(x-μx)^2/σx^2 - (z-μz)^2/σz^2)
μx = xMax/2.0
σx=xMax
μν=0.0
σν=0.1
uMat = zeros(Float64, Nxplus, Nνplus)
for νν=1:Nνplus
    for xx=1:Nxplus
        uMat[xx,νν] = fInitial(xs[xx], μx, σx, νs[νν], μν, σν)            
    end
end
u0 = reshape(uMat, nVerts)
u0[ghostMask.!=true] .= 0.0
integ = sum((W*u0)[ghostMask])
u0 ./= integ

#%%
# function model!(du, u, p, t)
#     du .= p*u
# end

# L = DiffEqArrayOperator(dropzeros(E*L1.+L2), update_func = update_func!)
L = MatrixOperator(Matrix(dropzeros(E*L1.+L2)), update_func! = update_func!)
# update_func!(L, u0, p, 0.0)
prob = ODEProblem(L, u0, (0.0, tMax), p)
# @time sol = solve(prob, LieRK4(), dt = 1/10, saveat=tMax/100.0, progress=true)
# @time sol = solve(prob, MagnusGL6(), dt = 0.001, saveat=tMax/100.0, progress=true)
# @time sol = solve(prob, RKMK2(), dt = 0.001, saveat=tMax/100.0, progress=true)
@time sol = solve(prob, MagnusAdapt4(), saveat=tMax/100.0, progress=true)



#%%

fig = Figure(size=(1000,1000))
ax = Axis3(fig[1, 1], aspect=:equal, azimuth=2.275π)
ax.xlabel = "x"
ax.ylabel = "ν"
ax.zlabel = "c"
uInternal = Observable(reshape(sol.u[1][ghostMask], (Nx, Nν)))
globalmin = minimum([minimum(u[ghostMask]) for u in sol.u])
globalmax = maximum([maximum(u[ghostMask]) for u in sol.u])
zlims!(ax, (globalmin, globalmax))
clims = (globalmin,globalmax)
surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
record(fig, datadir("surface2D.mp4"), 1:length(sol.t); framerate=10) do i
    uInternal[] .= reshape(sol.u[i][ghostMask], (Nx, Nν))
    uInternal[] = uInternal[]
end

fig2 = Figure(size=(1000,1000))
ax2 = Axis(fig2[1, 1])
ax2.xlabel = "x"
ax2.ylabel = "Fₑ"
lines!(ax2, xs[Nghost+1:end-Nghost], matFₑ[Nghost+1:end-Nghost])
display(fig2)
save(datadir("matF_e.png"), fig2)