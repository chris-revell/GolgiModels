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

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix

# General advection diffusion form:
# ċ = ∇⋅(D∇c - vc)
# where ċ is the time derivative of concentration
# ∇⋅ is the divergence
# D∇c is the gradient field of concentration multiplied by a (scalar) field of diffusivity, evaluated on edges
# vc is the concentration velocity vector field, also evaluated on edges. 
# Suppose vc is the velocity field enforced on edges multiplied by the mean concentration on vertices at either end of the edge.
# In other words, (vc)ⱼ = Āc.*v where velocity field v is expressed as a vector of size 2N
# Then we take the divergence simply by multiplying by Aᵀ
# so ċ = ∇⋅(D∇c - vc) -> Aᵀ*(D.*A*c .- v.*Ā*c) -> (AᵀDA - V*Ā)c where V is a diagonal vector of velocities and D is a diagonal vector of diffusivities 

# (1+αh)Ċ + β∇C = K₂K₄∇²C Gradients with respect to ν, dot with respect to time
# Ċ = (K₂K₄∇²C - β∇C)/(1+αh)

# ν ∈ (0.0,1.0)
# ∱ν over space = 1.0
# ν space discretised into N points  
Nx = 21
Ny = 21
Nν = 21
Nghost = 1 # number of ghost points on each side of the domain 
Nxplus = Nx+2*Nghost # number of discretised points including ghost points 
Nyplus = Ny+2*Nghost # number of discretised points including ghost points 
Nνplus = Nν+2*Nghost # number of discretised points including ghost points 
nVerts = Nxplus*Nyplus*Nνplus       # Total number of vertices 
nEdgesi = (Nxplus-1)*Nyplus*Nνplus  # Number of i-directed edges (x, in this case)
nEdgesj = Nxplus*(Nyplus-1)*Nνplus  # Number of j-directed edges (y, in this case)
nEdgesk = Nxplus*Nyplus*(Nνplus-1)  # Number of k-directed edges (ν, in this case)
nEdges = nEdgesi+nEdgesj+nEdgesk    # Total number of edges over all dimensions 
xMax = 1.0
yMax = 1.0
νMax = 1.0
xs = collect(range(0.0..xMax, Nxplus))
dx = xs[2]-xs[1]
ys = collect(range(0.0..yMax, Nyplus))
dy = ys[2]-ys[1]
νs = collect(range(0.0..νMax, Nνplus))
dν = νs[2]-νs[1]

α = 10.0
h = 1.0
β = 0.1
K₂ = 0.05
K₄ = 0.05  
tMax = 0.02

# Dx = K₂*K₄/(1+α*h)  # Diffusivity in x
# Dy = K₂*K₄/(1+α*h)  # Diffusivity in x
# Dν = K₂*K₄/(1+α*h)  # Diffusivity in ν


#%%

hFun(x, μx, σx, y, μy, σy) = exp(-(x-μx)^2/σx^2 -(y-μy)^2/σy)

#%%

A = makeIncidenceMatrix3D(Nxplus, Nyplus, Nνplus)
Ā = abs.(A)
Aᵀ = transpose(A)
Aᵤₚ = dropzeros((Ā-A)./2)




Gvec = vcat((dx^2).*collect(1:nEdgesi), (dy^2).*collect(1:nEdgesj), (dν^2).*collect(1:nEdgesk))
G = spdiagm(1.0./Gvec)



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
∇ₑ = l⁻¹*A      # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 
# L1 = a*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ) 
# L2 = D.*a*∇cdot*(hₑ*Pxy*∇ₑ)



# GDM⁻¹vec = vcat(dx.*reshape(D_i, edgesPerDirection)./dν, dν.*reshape(D_j, edgesPerDirection)./h)
# GDM⁻¹ = spdiagm(GDM⁻¹vec)
Aᵤₚ = dropzeros((Ā-A).÷2)
# G = dropzeros(M⁻¹*F)
# L = -W⁻¹*Aᵀ*GDM⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 
L = -W⁻¹*Aᵀ*G*D*F⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 

#%%

# Initial conditions using Gaussian
fInitial(x, μx, σx, y, μy, σy, z, μz, σz) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2 - (z-μz)^2/σz^2)
μx = 0.5
σx=10.0
μy = 0.5
σy=10.0
μν=0.0
σν=0.1
# Integrate Gaussian over domain 
# integ = hcubature(x -> fInitial(x[1], μx, σx, x[2], μy, σy, x[3], μν, σν), [0.0, 0.0, 0.0], [xMax, yMax, νMax])
uMat = zeros(Float64, Nxplus, Nyplus, Nνplus)
for νν=1:Nνplus
    for yy=1:Nyplus
        for xx=1:Nxplus
            # uMat[xx+Nghost,yy+Nghost,νν+Nghost] = fInitial(xs[xx+Nghost], μx, σx, ys[νν+Nghost], μy, σy, νs[νν+Nghost], μν, σν)            
            uMat[xx,yy,νν] = fInitial(xs[xx], μx, σx, ys[νν], μy, σy, νs[νν], μν, σν)            
        end
    end
end
u0 = reshape(uMat, nVerts)
u0[ghostMask.!=true] .= 0.0
integ = sum((W*u0)[ghostMask])
u0 ./= integ

#%%

p = L
function model!(du, u, p, t)
    du .= p*u
end

prob = ODEProblem(model!, u0, (0.0,tMax), p)

# prob = ODEProblem(MatrixOperator(L), u0, (0.0,tMax), p=())
t1 = @elapsed sol1 = solve(prob, Vern9(), saveat=tMax/100.0) ; nothing
# println("Finished 1")
# t2 = @elapsed sol2 = solve(prob, Trapezoid(), saveat=tMax/100.0) ; nothing
# println("Finished 2")
# sol = solve(prob, LinearExponential(), saveat=tMax/100.0)

#%%

fig = Figure(size=(1000,1000))
ax = Axis3(fig[1, 1], aspect=:equal, azimuth=2.275π)
ax.xlabel = "x"
ax.ylabel = "ν"
ax.zlabel = "c"
uInternal3D = reshape(sol1.u[1][ghostMask], (Nx, Ny, Nν))
uInternal = Observable(zeros(Nx, Nν))
globalmin = minimum([minimum(u[ghostMask]) for u in sol1.u])
globalmax = maximum([maximum(u[ghostMask]) for u in sol1.u])
zlims!(ax, (globalmin, globalmax))
clims = (globalmin,globalmax)
surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
record(fig, datadir("c_against_xν_at_halfy.mp4"), 1:length(sol1.t); framerate=10) do i
    uInternal3D .= reshape(sol1.u[i][ghostMask], (Nx, Ny, Nν))
    uInternal[] .= uInternal3D[:,Ny÷2,:]
    uInternal[] = uInternal[]
end

#%%

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "M, ∱cdxdy"
uInternal3D = reshape((W*sol1.u[end])[ghostMask], (Nx, Ny, Nν))
Mtmp = sum(uInternal3D, dims=1)
minima = Float64[]
maxima = Float64[]
for i=1:length(sol1.t)
    uInternal3D .= reshape(sol1.u[i][ghostMask], (Nx, Ny, Nν))
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
record(fig, datadir("Mvsν.mp4"), 1:length(sol1.t); framerate=10) do i
    uInternal3D .= reshape(sol1.u[i][ghostMask], (Nx, Ny, Nν))
    uInternal[] .= uInternal3D[:,Ny÷2,:]
    uInternal[] = uInternal[]
    Mtmp .= sum(uInternal3D, dims=1)
    M[] .= sum(Mtmp, dims=2)[1,1,:]
    M[] = M[]
end

#%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "x"
# ax.ylabel = "y"
# uInternal3D = reshape(sol.u[end][ghostMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nx,Ny))
# clims = (minimum(uInternal3D), maximum(uInternal3D))
# heatmap!(ax, xs[Nghost+1:end-Nghost], ys[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
# record(fig,"finalStateNuScan.mp4", 1:Nν; framerate=10) do i
#     uInternal[] .= uInternal3D[:,:,i]
#     uInternal[] = uInternal[]
# end

#%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "x"
# ax.ylabel = "y"
# uInternal3D = reshape(sol.u[1][ghostMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nx,Ny))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# clims = (globalmin,globalmax)
# clims = (minimum(uInternal3D), maximum(uInternal3D))
# heatmap!(ax, xs[Nghost+1:end-Nghost], ys[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
# record(fig, datadir("LargeNuTimeScan.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal3D .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))
#     uInternal[] .= uInternal3D[:,:,end]
#     uInternal[] = uInternal[]
#     ax.title = "xy profile of ν=1.0 at t=$(sol.t[i])"
# end

#%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "c"
# uInternal3D1 = reshape(sol.u[1][ghostMask], (Nx, Ny, Nν))
# uInternal3D2 = reshape(sol2.u[1][ghostMask], (Nx, Ny, Nν))
# uInternal1 = Observable(zeros(Nν))
# uInternal2 = Observable(zeros(Nν))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# ylims!(ax, (globalmin, globalmax))
# scatter!(ax, νs[Nghost+1:end-Nghost], uInternal1, marker=:circle, color=(:red,0.75), markersize=10)
# scatter!(ax, νs[Nghost+1:end-Nghost], uInternal2, marker=:cross, color=(:green,0.75), markersize=10)
# ax.title = "c against ν at x=0.5, y=0.5 at t=0.0"
# record(fig, datadir("NuProfileAtxyOverTime.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal3D1 .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))
#     uInternal3D2 .= reshape(sol2.u[i][ghostMask], (Nx, Ny, Nν))
#     uInternal1[] .= uInternal3D1[Nx÷2,Nx÷2,:]
#     uInternal2[] .= uInternal3D2[Nx÷2,Nx÷2,:]
#     uInternal1[] = uInternal1[]
#     uInternal2[] = uInternal2[]
#     ax.title = "c against ν at x=0.5, y=0.5 at t=$(sol.t[i])"
# end

#%%

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
# record(fig, datadir("xνOverTimeAty.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal[] .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))[:,Nyplus÷2,:]
#     uInternal[] = uInternal[]
#     ax.title = "c against x and ν at y=0.5 at t=$(sol.t[i])"
# end

#%%

# D_i2 = fill(0.0, (Nxplus-1, Nyplus, Nνplus))
# D_j2 = fill(0.0, (Nxplus, Nyplus-1, Nνplus))
# D_k2 = fill(0.0, (Nxplus, Nyplus, Nνplus-1))
# for k=2:Nνplus-1
#     for j=2:Nyplus-1
#         for i=2:Nxplus-1
#             D_i2[i,j,k] = dy*dν*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#             D_j2[i,j,k] = dx*dν*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#             D_k2[i,j,k] = dx*dy*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#         end
#     end
# end
# D_i2[:, :, 1] .= 0.0
# D_i2[:, :, end] .= 0.0
# D_i2[:, 1, :] .= 0.0
# D_i2[:, end, :] .= 0.0
# D_i2[1, :, :] .= 0.0
# D_i2[end, :, :] .= 0.0
# D_j2[:, :, 1] .= 0.0
# D_j2[:, :, end] .= 0.0
# D_j2[:, 1, :] .= 0.0
# D_j2[:, end, :] .= 0.0
# D_j2[1, :, :] .= 0.0
# D_j2[end, :, :] .= 0.0
# D_k2[:, :, 1] .= 0.0
# D_k2[:, :, end] .= 0.0
# D_k2[:, 1, :] .= 0.0
# D_k2[:, end, :] .= 0.0
# D_k2[1, :, :] .= 0.0
# D_k2[end, :, :] .= 0.0