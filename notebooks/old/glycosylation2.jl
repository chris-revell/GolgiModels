
using DifferentialEquations
using SparseArrays
using UnPack
using CairoMakie 
using LinearAlgebra

# (1+αh)Ċ + β∇C = K₂K₄∇²C Gradients with respect to ν, dot with respect to time
# Ċ = (K₂K₄∇²C - β∇C)/(1+αh)

# ν ∈ (0.0,1.0)
# ∱ν over space = 1.0
# ν space discretised into N points 
N = 101
Nghost = 1 # number of ghost points on each side of the domain 
Nplus = N+2*Nghost # number of discretised points including ghost points 
x = collect(range(0.0..1.0, N))
dx = x[2]-x[1]
α = 1.0
h = 1.0
β = 0.1
K₂ = 1.0
K₄ = 1.0

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
Ā = abs.(A)
# D degree matrix 
D = spdiagm(sum(eachrow(A)))
∇² = (D.-Ā)./dx^2

# https://en.wikipedia.org/wiki/Finite_difference_coefficient
∇ = spzeros(Float64, Nplus, Nplus)
∇[1,2] = 1
∇[Nplus,Nplus-1] = -1
for i=2:Nplus-1
    ∇[i, i-1] = -1
    ∇[i, i+1] = 1
end
∇ .= ∇./dx

#%%

# Express model as a matrix operator 
M = (K₂*K₄.*∇² .- β.*∇)./(1+α*h)

# Initial conditions
u0 = zeros(Float64, Nplus)
y = exp.((-1.0.*x.^2)./0.1^2)
# Trapezium rule integration
integ = [0.0]
for i=1:N-1
    integ[1] += dx*0.5*(y[i]+y[i+1])
end
u0[Nghost+1:end-Nghost] .= y./integ[1]

# Parameters named tuple
p = (α=α,
    h=h,
    β=β,
    K₂K₄=K₂*K₄,
    M=M)

function boundaryConditions!(u, p)
    @unpack β, K₂K₄ = p 
    # βC - K₂K₄∇C = 0 at periphery
    # (β - K₂K₄∇)*C = 0
    # u[1] = (2+β/K₂K₄)*u[2] - u[3]
    u[1] = - u[3]
    # u[end] = (2+β/K₂K₄)*u[end-1] - u[end-2]
    u[end] = - u[end-2]
end

function model!(du, u, p, t)
    @unpack M = p 
    boundaryConditions!(u, p)
    du .= M*u
end

prob = ODEProblem(model!, u0, (0.0,30.0), p)
sol = solve(prob)

#%%

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1])
ylims!(ax, (0.0,15.0))
uInternal = Observable(zeros(N))
lines!(ax, x, uInternal, color=:blue)
record(fig,"test.mp4", 1:length(sol.t); framerate=1) do i
    uInternal[] = sol.u[i][Nghost+1:end-Nghost]
    uInternal[] = uInternal[]
end

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