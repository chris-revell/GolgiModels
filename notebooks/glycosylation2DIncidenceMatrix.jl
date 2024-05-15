using DifferentialEquations
using SparseArrays
using UnPack
using CairoMakie 
using HCubature
using FromFile
using DrWatson
using Printf
using LaTeXStrings

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
N = 51
Nghost = 1 # number of ghost points on each side of the domain 
Nplus = N+2*Nghost # number of discretised points including ghost points 
xMax = 1.0
νMax = 1.0
xs = collect(range(0.0..xMax, Nplus))
dx = xs[2]-xs[1]
νs = collect(range(0.0..νMax, Nplus))
dν = νs[2]-νs[1]

α = 10.0
h = 1.0
β = 0.5
K₂ = 0.5
K₄ = 0.5  
tMax = 1.0

vx = 0.0            # Scalar advection velocity in x
vν = β/(1+α*h)      # Scalar advection velocity in ν
Dx = K₂*K₄/(1+α*h)  # Diffusivity in x
Dν = K₂*K₄/(1+α*h)  # Diffusivity in ν


#%%

hFun(x, μx, σx) = 2.0*exp(-(x-μx)^2/σx^2)+0.1
# hFun(x, μx, σx) = 1.0

#%%

A = makeIncidenceMatrix2D(Nplus)
Ā = abs.(A)
Aᵀ = transpose(A)

nEdges = size(A, 1)
nVerts = size(A, 2)
edgesPerDirection = nEdges÷2 #(Nplus-1)*Nplus

# Vertex weights
matE = fill(dx*dν, (Nplus, Nplus))
# Left edge
matE[:, 1] ./= 2.0
# Right edge
matE[:, end] ./= 2.0
# Top edge
matE[1, :] ./= 2.0
# Bottom edge
matE[end, :] ./= 2.0
vecE = reshape(matE, nVerts)
E⁻¹ = spdiagm(1.0./vecE)

# Edge weights 
F_i = fill(dx*dν, (Nplus-1, Nplus))
# Left edge, i-directed 
F_i[:, 1] ./= 2.0
# Right edge, i-directed 
F_i[:, end] ./= 2.0
# Top edge, i-directed 
F_i[1, :] ./= 2.0
# Bottom edge, i-directed 
F_i[end, :] ./= 2.0
# Matrix of i-directed edges  
F_j = fill(dx*dν, (Nplus, Nplus-1))
# Left edge, j-directed 
F_j[:, 1] ./= 2.0
# Right edge, j-directed 
F_j[:, end] ./= 2.0
# Top edge, j-directed 
F_j[1, :] ./= 2.0
# Bottom edge, j-directed 
F_j[end, :] ./= 2.0
Fvec = vcat(reshape(F_i, edgesPerDirection), reshape(F_j, edgesPerDirection))
F = spdiagm(Fvec)
F⁻¹ = spdiagm(1.0./Fvec)


# Edge weights 
vecM = fill(dx*dν, nEdges)
M⁻¹ = spdiagm(1.0./vecM)

# Velocity field 
V_i = vν.*ones(Float64, Nplus-1, Nplus)
for i=1:Nplus-1
    for j=1:Nplus
        V_i[i,j] *= β/(1+α*hFun(xs[j], 0.5, 0.1))
    end
end
# Left edge, i-directed 
V_i[:,1] .= 0.0
# Right edge, i-directed 
V_i[:,end] .= 0.0
# Top edge, i-directed 
V_i[1,:] .= 0.0
# Bottom edge, i-directed 
V_i[end,:] .= 0.0
# Matrix of j-directed edges  
V_j = vx.*ones(Float64, Nplus, Nplus-1)
for i=1:Nplus
    for j=1:Nplus-1
        V_j[i,j] *= β/(1+α*hFun(xs[j], 0.5, 0.1))
    end
end
# Left edge, j-directed 
V_j[:,1] .= 0.0
# Right edge, j-directed 
V_j[:,end] .= 0.0
# Top edge, j-directed 
V_j[1,:] .= 0.0
# Bottom edge, j-directed 
V_j[end,:] .= 0.0
Vvec = vcat(reshape(V_i, edgesPerDirection), reshape(V_j, edgesPerDirection))
V = spdiagm(Vvec)


# Diffusivity field over edges 
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
# Matrix of i-directed edges  
D_i = dx.*ones(Float64, Nplus-1, Nplus)   # Diffusive flux in ν is a function of grid edge length in x
for i=1:Nplus-1
    for j=1:Nplus
        D_i[i,j] *= K₂*K₄/(1+α*hFun(xs[j], 0.5, 0.1))  # Diffusivity in ν
    end
end
# Left edge, i-directed 
D_i[:,1] .= 0.0
# Right edge, i-directed 
D_i[:,end] .= 0.0
# Top edge, i-directed 
D_i[1,:] .= 0.0
# Bottom edge, i-directed 
D_i[end,:] .= 0.0
# Matrix of i-directed edges  
D_j = Dx*dν.*ones(Float64, Nplus, Nplus-1)   # Diffusive flux in ν is a function of grid edge length in x
for i=1:Nplus
    for j=1:Nplus-1
        D_j[i,j] *= K₂*K₄/(1+α*hFun(xs[j], 0.5, 0.1))  # Diffusivity in ν
    end
end
# Left edge, j-directed 
D_j[:,1] .= 0.0
# Right edge, j-directed 
D_j[:,end] .= 0.0
# Top edge, j-directed 
D_j[1,:] .= 0.0
# Bottom edge, j-directed 
D_j[end,:] .= 0.0
Dvec = vcat(reshape(D_i, edgesPerDirection), reshape(D_j, edgesPerDirection))
D = spdiagm(Dvec)



GDM⁻¹vec = vcat(dx.*reshape(D_i, edgesPerDirection)./dν, dν.*reshape(D_j, edgesPerDirection)./h)
GDM⁻¹ = spdiagm(GDM⁻¹vec)

Aᵤₚ = dropzeros((Ā-A).÷2)

G = dropzeros(M⁻¹*F)

# L = E⁻¹*Aᵀ*G*D*M⁻¹*A .- E⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 
L = -E⁻¹*Aᵀ*GDM⁻¹*A .+ E⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 

#%%

# Initial conditions using Gaussian
fInitial(x, μx, σx, y, μy, σy) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2)
μx = 0.5
σx=10.0
μν=0.0
σν=0.1
# Integrate Gaussian over domain 
integ = hcubature(x -> fInitial(x[1], μx, σx, x[2], μν, σν), [0.0, 0.0], [1.0, 1.0])
uMat = zeros(Float64, Nplus, Nplus)
for νν=1:N
    for xx=1:N
        uMat[νν+Nghost,xx+Nghost] = exp(-(xs[xx+Nghost]-μx)^2/σx^2 - (νs[νν+Nghost]-μν)^2/σν^2)
    end
end
uMat .= uMat./integ[1]
u0 = reshape(uMat, Nplus^2)

p = L

function model!(du, u, p, t)
    du .= p*u
end

prob = ODEProblem(model!, u0, (0.0,tMax), p)
sol = solve(prob, Trapezoid(), saveat=tMax/100.0)

#%%

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "x"
mins = [minimum(u) for u in sol.u[2:end]]
maxs = [maximum(u) for u in sol.u[2:end]]
uInternal = Observable(zeros(N,N))
clims = Observable((0.0,1.0))
heatmap!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:navia)
record(fig,"test2D2.mp4", 1:length(sol.t); framerate=10) do i
    uInternal[] .= reshape(sol.u[i],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
    uInternal[] = uInternal[]
    clims[] = (minimum(uInternal[]), maximum(uInternal[]))
    clims[] = clims[]
end

#%%

fig2 = Figure(size=(3250,2000), fontsize = 32)
ax1 = Axis(fig2[1,1], aspect=1)
ax1.xlabel = "h"
ax1.ylabel = "x"
heatmap!(ax1, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], zeros(Nplus,Nplus),colormap=:bwr)
lines!(ax1, hFun.(xs[2:end-1], 0.5, 0.1), xs[2:end-1], linewidth=2)

ax2 = Axis(fig2[1,2], aspect=1)
ax2.xlabel = "ν"
ax2.ylabel = "x"
ax2.title = @sprintf "t = %.2f" 0.0
# heatmap!(ax2, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
heatmap!(ax2, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)
ax3 = Axis(fig2[1,3], aspect=1)
ax3.xlabel = "ν"
ax3.ylabel = "x"
ax3.title = @sprintf "t = %.2f" tMax/2.0
# heatmap!(ax3, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
heatmap!(ax3, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)
ax4 = Axis(fig2[1,4], aspect=1)
ax4.xlabel = "ν"
ax4.ylabel = "x"
ax4.title = @sprintf "t = %.2f" tMax
# heatmap!(ax4, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
heatmap!(ax4, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)

linkyaxes!(ax1, ax2)

ax5 = Axis(fig2[2,2], aspect=1)
ax5.xlabel = "ν"
ax5.ylabel = "M"
ax5.title = @sprintf "t = %.2f" 0.0
integ = zeros(Float64, N)
u = reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
for j=1:N
    integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
end
lines!(ax5, νs[Nghost+1:end-Nghost], integ)
ax6 = Axis(fig2[2,3], aspect=1)
ax6.xlabel = "ν"
ax6.ylabel = "M"
ax6.title = @sprintf "t = %.2f" tMax/2.0
u = reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
for j=1:N
    integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
end
lines!(ax6, νs[Nghost+1:end-Nghost], integ)
ax7 = Axis(fig2[2,4], aspect=1)
ax7.xlabel = "ν"
ax7.ylabel = "M"
ax7.title = @sprintf "t = %.2f" tMax
u = reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
for j=1:N
    integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
end
lines!(ax7, νs[Nghost+1:end-Nghost], integ)


display(fig2)
#%%
save("glycosylation.png", fig2)