using DifferentialEquations
using SparseArrays
using UnPack
using CairoMakie 

# (1+αh)Ċ + β∇C = K₂K₄∇²C Gradients with respect to ν, dot with respect to time
# Ċ = (K₂K₄∇²C - β∇C)/(1+αh)

# ν ∈ (0.0,1.0)
# ∱ν = 1.0
# ν space discretised into N points 
N = 21
Nghost = 1 # number of ghost points on each side of the domain 
Nplus = N+2*Nghost # number of discretised points including ghost points 
x = collect(range(0.0..1.0, N))
dx = x[2]-x[1]
α = 1.0
h = 1.0
β = 3.0
K₂ = 3.0
K₄ = 1.0  
tMax = 3.0

#%%
# https://en.wikipedia.org/wiki/Laplacian_matrix Laplacian matrix for simple graph
# A adjacency matrix 
A = spzeros(Int64, Nplus, Nplus)
A[1,2] = 1
A[Nplus,Nplus-1] = 1
for i=2:Nplus-1
    A[i, i-1] = 1
    A[i, i+1] = 1
end
# deg degree matrix 
deg = spdiagm(sum(eachrow(A)))
∇² = (A.-deg)./dx^2

# https://en.wikipedia.org/wiki/Finite_difference_coefficient
∇ = spzeros(Float64, Nplus, Nplus)
∇[1,2] = 1
∇[Nplus,Nplus-1] = -1
for i=2:Nplus-1
    ∇[i, i-1] = -1
    ∇[i, i+1] = 1
end
∇ .= ∇./(2*dx)

#%%

# Express model as a matrix operator 
M = ∇².*K₂*K₄/(1+α*h) .- ∇.*β/(1+α*h)

# Initial conditions using Gaussian
u0 = zeros(Float64, Nplus)
y = exp.((-1.0.*x.^2)./0.1^2)
# Trapezium rule integration
integ = [0.0]
for i=1:length(y)-1
    integ[1] += dx*0.5*(y[i]+y[i+1])
end
u0[Nghost+1:end-Nghost] .= y./integ[1]

# Parameters named tuple
p = (D=K₂*K₄/(1+α*h),
    β=β/(1+α*h),
    M=M)

function boundaryConditions!(u, p)
    @unpack D, β, M = p 
    # D∇C - βC = 0 at periphery
    # D(u[x+1]-u[x-1])/(2dx) - βu[x] = 0 
    # u[x+1]-u[x-1] = 2dxβu[x]/D
    u[1] = u[3] - 2*dx*β*u[2]/D   
    u[end] = u[end-2] + 2*dx*β*u[end-1]/D    
end

function model!(du, u, p, t)
    @unpack M = p 
    boundaryConditions!(u, p)
    du .= M*u
end

prob = ODEProblem(model!, u0, (0.0,tMax), p)
sol = solve(prob, Trapezoid(), saveat=tMax/100.0)

#%%

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1])
ylims!(ax, (0.0,15.0))
uInternal = Observable(zeros(N))
lines!(ax, x, uInternal, color=:blue)
record(fig,"test.mp4", 1:length(sol.t); framerate=10) do i
    uInternal[] = sol.u[i][Nghost+1:end-Nghost]
    uInternal[] = uInternal[]
end


#%%

masses = zeros(length(sol.u))
for t=1:length(sol.u)
    for i=1+Nghost:Nplus-Nghost-1
        masses[t] += dx*0.5*(sol.u[t][i]+sol.u[t][i+1])
    end
end
fig2 = Figure(size=(1000,1000))
ax2 = CairoMakie.Axis(fig2[1, 1])
lines!(ax2, sol.t, masses)
ylims!(ax2, (0.95,1.05))
ax2.xlabel = "Time"
ax2.ylabel = "Mass"
display(fig2)



#%%

# fig = Figure(size=(1000, 1000))
# ax = CairoMakie.Axis(fig[1, 1])#, aspect=DataAspect())
# toPlot = hcat([u[Nghost+1:end-Nghost] for u in sol.u]...)
# heatmap!(ax, collect(range(0.0..1.0, N)), sol.t, toPlot, colormap=:batlow)
# resize_to_layout!(fig)
# display(fig)
# # save("$subFolder/$(fileName)_finalState.png", fig)


#%%
# # https://en.wikipedia.org/wiki/Laplacian_matrix Symmetric Laplacian via the incidence matrix
# # A incidence matrix 
# A = spzeros(Float64, νMax+3, νMax+4)
# for j=1:νMax+3
#     A[j,j] = -1
#     A[j,j+1] = 1
# end
# # Edge weights matrix 
# W = spdiagm(ones(νMax+3))
# ∇² = transpose(A)*W*A

# ∇ = A./2.0


#%%

# fig = Figure(size=(1000, 1000))
# ax = CairoMakie.Axis(fig[1, 1])#, aspect=DataAspect())
# mov = 
# lines!(ax, x, sol.u[1][Nghost+1:end-Nghost], color=:blue)
# lines!(ax, x, sol.u[end][Nghost+1:end-Nghost], color=:blue)
# resize_to_layout!(fig)
# display(fig)

#%%

# ∇f = spzeros(Float64, Nplus, Nplus)
# for i=1:N
#     ∇f[i, i] = -3
#     ∇f[i, i+1] = 2
#     ∇f[i, i+2] = 1
# end
# ∇f[N+1, N+1] = -3
# ∇f[N+1, N+2] = 1
# ∇f[N+2, N+2] = -3
# ∇f .= ∇f./(2*dx)
