using DifferentialEquations
using SparseArrays
using UnPack
using CairoMakie 
using HCubature
using FromFile
using DrWatson
using Printf
using LaTeXStrings
# using InvertedIndices
using GaussianRandomFields
using TerminalLoggers

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix

Nx = 11
Ny = 11
Nν = 11
Nghost = 1           # number of ghost points on each side of the domain 
Nxplus = Nx+2*Nghost # number of discretised points including ghost points 
Nyplus = Ny+2*Nghost # number of discretised points including ghost points 
Nνplus = Nν+2*Nghost # number of discretised points including ghost points 
xMax = 1.0; yMax = 1.0; νMax = 1.0
xs = collect(range(0.0..xMax, Nxplus))
dx = xs[2]-xs[1]
ys = collect(range(0.0..yMax, Nyplus))
dy = ys[2]-ys[1]
νs = collect(range(0.0..νMax, Nνplus))
dν = νs[2]-νs[1]

α = 1.0
δ = 1.0
β = 0.1
σ = 1.0
N = 100
K₂ = 0.5
K₃ = 1.0
K₄ = 0.5  
tMax = 0.01

# vx = 0.0            # Scalar advection speed in x
# vy = 0.0            # Scalar advection speed in y
# vν = β/(1+α*h)      # Scalar advection speed in ν
# Dx = K₂*K₄/(1+α*h)  # Diffusivity in x
# Dy = K₂*K₄/(1+α*h)  # Diffusivity in x
# Dν = K₂*K₄/(1+α*h)  # Diffusivity in ν


#%%

hFun(x, μx, σx, y, μy, σy) = exp(-(x-μx)^2/σx^2 -(y-μy)^2/σy)
# hFun(x, xhz, y, yhz) = 1+0.5*sin(2π*(x*xhz+y*yhz))

#%%

A = makeIncidenceMatrix3D(Nxplus, Nyplus, Nνplus)
Ā = abs.(A)
Aᵀ = transpose(A)
Aᵤₚ = dropzeros((Ā-A)./2)

nVerts = Nxplus*Nyplus*Nνplus       # Total number of vertices 
nEdgesi = (Nxplus-1)*Nyplus*Nνplus  # Number of i-directed edges (x, in this case)
nEdgesj = Nxplus*(Nyplus-1)*Nνplus  # Number of j-directed edges (y, in this case)
nEdgesk = Nxplus*Nyplus*(Nνplus-1)  # Number of k-directed edges (ν, in this case)
nEdges = nEdgesi+nEdgesj+nEdgesk    # Total number of edges over all dimensions 

# Ghost point mask is a 1D vector in which the value at component i 
# is true if vertex i in the flattened vector of vertices is an internal vertex
# but false if vertex i is a ghost vertex in the flattened vector of vertices 
# This can be used to exclude ghost points from calculations over the whole state vector
ghostMaskArray = fill(true, (Nxplus, Nyplus, Nνplus))
ghostMaskArray[:, :, 1] .= false
ghostMaskArray[:, :, end] .= false
ghostMaskArray[:, 1, :] .= false
ghostMaskArray[:, end, :] .= false
ghostMaskArray[1, :, :] .= false
ghostMaskArray[end, :, :] .= false
ghostMask = reshape(ghostMaskArray, nVerts)

# Vertex weights
# Forming a diagonal matrix of volumes around each vertex, divided by 2 at the periphery
matW = fill(dx*dy*dν, (Nxplus, Nyplus, Nνplus))
matW[:, :, 1] ./= 2.0
matW[:, :, end] ./= 2.0
matW[:, 1, :] ./= 2.0
matW[:, end, :] ./= 2.0
matW[1, :, :] ./= 2.0
matW[end, :, :] ./= 2.0
vecW = reshape(matW, nVerts)
W = spdiagm(vecW)
W⁻¹ = spdiagm(1.0./vecW)

# Diagonal matrix of volumes around each edge, divided by 2 at the periphery
# Matrix of i-directed edge weights  
F_i = fill(dx*dy*dν, (Nxplus-1, Nyplus, Nνplus))
F_i[:, :, 1] ./= 2.0; F_i[:, :, end] ./= 2.0; F_i[:, 1, :] ./= 2.0; F_i[:, end, :] ./= 2.0; F_i[1, :, :] ./= 2.0; F_i[end, :, :] ./= 2.0
# Matrix of j-directed edge weights  
F_j = fill(dx*dy*dν, (Nxplus, Nyplus-1, Nνplus))
F_j[:, :, 1] ./= 2.0; F_j[:, :, end] ./= 2.0; F_j[:, 1, :] ./= 2.0; F_j[:, end, :] ./= 2.0; F_j[1, :, :] ./= 2.0; F_j[end, :, :] ./= 2.0
# Matrix of k-directed edge weights  
F_k = fill(dx*dy*dν, (Nxplus, Nyplus, Nνplus-1))
F_k[:, :, 1] ./= 2.0; F_k[:, :, end] ./= 2.0; F_k[:, 1, :] ./= 2.0; F_k[:, end, :] ./= 2.0; F_k[1, :, :] ./= 2.0; F_k[end, :, :] ./= 2.0
Fvec = vcat(reshape(F_i, nEdgesi), reshape(F_j, nEdgesj), reshape(F_k, nEdgesk))
F = spdiagm(Fvec)
F⁻¹ = spdiagm(1.0./Fvec)

# Diagonal matrix of edge lengths
l_i = fill(dx, (Nxplus-1, Nyplus, Nνplus))
l_j = fill(dy, (Nxplus, Nyplus-1, Nνplus))
l_k = fill(dν, (Nxplus, Nyplus, Nνplus-1))
lvec = vcat(reshape(l_i, nEdgesi), reshape(l_j, nEdgesj), reshape(l_k, nEdgesk))
l = spdiagm(lvec)
l⁻¹ = spdiagm(1.0./lvec)

# Matrices for picking out ν and xy directions in derivatives 
Pν = spdiagm(vcat(zeros(nEdgesi), zeros(nEdgesj), ones(nEdgesk)))
Pxy = spdiagm(vcat(ones(nEdgesi), ones(nEdgesj), zeros(nEdgesk)))

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
mat_h = zeros(Nxplus, Nyplus, Nνplus)
for j=1:Nyplus
    for i=1:Nxplus
        mat_h[i, j, :] .= hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1)
    end
end
h_vec = reshape(mat_h, nVerts)
hᵥ = spdiagm(h_vec)
# hₑ = spdiagm(0.5.*Ā*h_vec)
hₑ = spdiagm(Aᵤₚ*h_vec)

# Prefactor diagonal matrix a to store 1/(1+αh) term 
a = spdiagm(1.0./(1.0 .+ α.*h_vec))

# Fₑ field over all vertices, for enxyme concentration calculation, derived from Gaussian random field 
cov = CovarianceFunction(2, Exponential(.5))
pts = range(0, stop=1, length=1001)
grf = GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=2001)
matFₑ = sample(grf)[1:Nxplus, 1:Nyplus]
matFₑ .= matFₑ.-minimum(matFₑ)

# function to find enzyme concentration at each vertex by integrating over ν
function enzymeE!(Niplus, Njplus, c, K₂, matE, vecE, E, Fₑ, stride)
    # Fₑ(x)*K₂/(K₂ + ∱c(x,ν,t)dν)
    for i=1:Niplus
        for j=1:Njplus
            cs = @view c[(1+(i-1)*stride[1] + (j-1)*stride[2]):stride[3]:end]
            matE[i,j,:] .= Fₑ[i,j]*K₂/(K₂ + 0.5*sum(cs[1:end-1].+cs[2:end]))
        end
    end
    vecE .= reshape(matE, nVerts)
    E .= spdiagm(vecE) 
    return nothing 
end

matE = zeros(Nxplus, Nyplus, Nνplus)
vecE = reshape(matE, nVerts)
E = spdiagm(vecE)

# (1+α)hċ + EβC_ν = EK₂K₄C_νν + D∇⟂⋅(h∇⟂C)
# D = αδN²(K₂+σK₃)

# (1+α)hċ = EK₂K₄C_νν - Eβ*C_ν + D∇⟂⋅(h∇⟂C)
# ċ = (EK₂K₄C_νν - Eβ*C_ν + D∇⟂⋅(h∇⟂C))/((1+α)h)
# ċ = (EK₂K₄C_νν - Eβ*C_ν)/((1+α)h) + D∇⟂⋅(h∇⟂C)/((1+α)h)
# ċ = E(K₂K₄C_νν - β*C_ν)/((1+α)h) + D∇⟂⋅(h∇⟂C)/((1+α)h)

# ċ = (E(C)*L1 + L2)*C
# L1 = ((1+α)h)⁻¹(K₂K₄C_νν - β*C_ν)
# L2 = ((1+α)h)⁻¹(D∇⟂⋅(h∇⟂))
# a = ((1+α)h)⁻¹ ⟹ L1 = a(K₂K₄C_νν - β*C_ν), L2 = aD∇⟂⋅(h∇⟂)


# E(x,t) = Fₑ(x)*K₂/(K₂ + ∱c(x,ν,t)dν)
# D = α*δ*N^2*(K₂+σ*K₃)

Dᵥ = spdiagm(fill(α*δ*N^2*(K₂+σ*K₃), nVerts))
Dₑ = spdiagm(fill(α*δ*N^2*(K₂+σ*K₃), nEdges))


H = Dₑ*Pxy*hₑ + K₂*K₄*Pν
L = a*(W⁻¹*Aᵀ*H*l⁻¹*A - Aᵀ*l⁻¹*Pν*Aᵤₚ)


Cνν = W⁻¹*Aᵀ*Pν*l⁻¹*A
Cν = Aᵀ*l⁻¹*Pν*Aᵤₚ
L1 = a*(K₂*K₄.*Cνν - β*Cν)
L2 = a*D*W⁻¹*Aᵀ*Pxy*hₑ*l⁻¹*A


#%%

# Initial conditions using Gaussian
fInitial(x, μx, σx, y, μy, σy, z, μz, σz) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2 - (z-μz)^2/σz^2)
μx = 0.5
σx=10.0
μy = 0.5
σy=10.0
μν=0.0
σν=0.1
uMat = zeros(Float64, Nxplus, Nyplus, Nνplus)
for νν=1:Nνplus
    for yy=1:Nyplus
        for xx=1:Nxplus
            uMat[xx,yy,νν] = fInitial(xs[xx], μx, σx, ys[νν], μy, σy, νs[νν], μν, σν)            
        end
    end
end
u0 = reshape(uMat, nVerts)
u0[ghostMask.!=true] .= 0.0
integ = sum((W*u0)[ghostMask])
u0 ./= integ

#%%

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

function model!(du, u, p, t)
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

    enzymeE!(Nxplus, Nyplus, u, K₂, matE, vecE, E, matFₑ, strides)

    du .= (E*L1 + L2)*u
end

prob = ODEProblem(model!, u0, (0.0,tMax), p)
# prob = ODEProblem(MatrixOperator(L), u0, (0.0,tMax), p)
sol = solve(prob, Rosenbrock23(autodiff=false), saveat=tMax/100.0, progress=true)
# sol = solve(prob, LinearExponential(), saveat=tMax/100.0)

#%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "M, ∱cdxdy"
# uInternal3D = reshape((W*sol.u[end])[ghostMask], (Nx, Ny, Nν))
# M = sum(uInternal3D, dims=1)
# M = sum(M, dims=2)
# lines!(ax, νs[1:Nghost:end-2*Nghost], M[1,1,:])
# ax.title = "Integral of c over x and y against ν at final time"
# save(datadir("finalνVsM.png"), fig)

# #%%

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

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "c"
# uInternal3D = reshape(sol.u[1][ghostMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nν))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# ylims!(ax, (globalmin, globalmax))
# lines!(ax, νs[Nghost+1:end-Nghost], uInternal)
# ax.title = "c against ν at x=0.5, y=0.5 at t=0.0"
# record(fig, datadir("NuProfileAtxyOverTime.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal3D .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))
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
# record(fig, datadir("xνOverTimeAty.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal[] .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))[:,Nyplus÷2,:]
#     uInternal[] = uInternal[]
#     ax.title = "c against x and ν at y=0.5 at t=$(sol.t[i])"
# end

#%%
# fig2 = Figure(size=(3250,2000), fontsize = 32)
# ax1 = Axis(fig2[1,1], aspect=1)
# ax1.xlabel = "h"
# ax1.ylabel = "x"
# heatmap!(ax1, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], zeros(Nplus,Nplus),colormap=:bwr)
# lines!(ax1, hFun.(xs[2:end-1], 0.5, 0.1), xs[2:end-1], linewidth=2)

# ax2 = Axis(fig2[1,2], aspect=1)
# ax2.xlabel = "ν"
# ax2.ylabel = "x"
# ax2.title = @sprintf "t = %.2f" 0.0
# # heatmap!(ax2, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
# heatmap!(ax2, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)
# ax3 = Axis(fig2[1,3], aspect=1)
# ax3.xlabel = "ν"
# ax3.ylabel = "x"
# ax3.title = @sprintf "t = %.2f" tMax/2.0
# # heatmap!(ax3, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
# heatmap!(ax3, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)
# ax4 = Axis(fig2[1,4], aspect=1)
# ax4.xlabel = "ν"
# ax4.ylabel = "x"
# ax4.title = @sprintf "t = %.2f" tMax
# # heatmap!(ax4, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
# heatmap!(ax4, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)

# linkyaxes!(ax1, ax2)

# ax5 = Axis(fig2[2,2], aspect=1)
# ax5.xlabel = "ν"
# ax5.ylabel = L"\int_x c(x,\nu)\partial x"
# ax5.title = @sprintf "t = %.2f" 0.0
# integ = zeros(Float64, N)
# u = reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
# for j=1:N
#     integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
# end
# lines!(ax5, νs[Nghost+1:end-Nghost], integ)
# ax6 = Axis(fig2[2,3], aspect=1)
# ax6.xlabel = "ν"
# ax6.ylabel = L"\int_x c(x,\nu)\partial x"
# ax6.title = @sprintf "t = %.2f" tMax/2.0
# u = reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
# for j=1:N
#     integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
# end
# lines!(ax6, νs[Nghost+1:end-Nghost], integ)
# ax7 = Axis(fig2[2,4], aspect=1)
# ax7.xlabel = "ν"
# ax7.ylabel = L"\int_x c(x,\nu)\partial x"
# ax7.title = @sprintf "t = %.2f" tMax
# u = reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
# for j=1:N
#     integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
# end
# lines!(ax7, νs[Nghost+1:end-Nghost], integ)


# display(fig2)
# #%%
# save("glycosylation.png", fig2)

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






# Gvec = vcat((dx^2).*collect(1:nEdgesi), (dy^2).*collect(1:nEdgesj), (dν^2).*collect(1:nEdgesk))
# G = spdiagm(1.0./Gvec)


# # Diffusivity field over edges 
# # Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
# D_i = fill(0.0, (Nxplus-1, Nyplus, Nνplus))
# D_j = fill(0.0, (Nxplus, Nyplus-1, Nνplus))
# D_k = fill(0.0, (Nxplus, Nyplus, Nνplus-1))
# for k=2:Nνplus-2
#     for j=2:Nyplus-1
#         for i=2:Nxplus-1
#             D_k[i,j,k] = dx*dy*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#         end
#     end
# end
# for k=2:Nνplus-1
#     for j=2:Nyplus-2
#         for i=2:Nxplus-1
#             D_j[i,j,k] = dx*dν*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#         end
#     end
# end
# for k=2:Nνplus-1
#     for j=2:Nyplus-1
#         for i=2:Nxplus-2
#             D_i[i,j,k] = dy*dν*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#         end
#     end
# end
# Dvec = vcat(reshape(D_i, nEdgesi), reshape(D_j, nEdgesj), reshape(D_k, nEdgesk))
# D = spdiagm(Dvec) # Diagonal matrix of advection velocities at each edge
