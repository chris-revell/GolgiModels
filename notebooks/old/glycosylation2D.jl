using DifferentialEquations
using SparseArrays
using UnPack
using CairoMakie 
using HCubature

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
N = 21
Nghost = 1 # number of ghost points on each side of the domain 
Nplus = N+2*Nghost # number of discretised points including ghost points 
xs = collect(range(0.0..1.0, N))
dx = xs[2]-xs[1]
νs = collect(range(0.0..1.0, N))
dν = νs[2]-νs[1]

α = 1.0
h = 1.0
β = 5.0
K₂ = 3.0
K₄ = 1.0  
tMax = 1.0

#%%
# State matrix uMat has dimensions Nplus x Nplus. 
# Column major order, so when flattened to a state vector u, u[1:Nplus] == uMat[:,1], u[1+Nplus:2*Nplus]==uMat[2,:] etc
# State vector has size Nplus*Nplus
# Let us say that in the state matrix, the first dimension is ν and the second is x.
# --> x -->
# \
# ν
# \
# v   
# So in the state vector, the first Nplus components correspond to 0<=ν<=1.0 for x=0.0; second Nplus components correspond to 0<=ν<=1.0 for x=dx, and so on 

# Function to convert index (i,j) to access uMat[i,j] where uMat is of size N1xN2 to index k, 
# which is the corresponding index of the same data point when uMat is flattened to a vector 
# matrixIndexToVectorIndex(i, j, N1, N2) = (j-1)*N1+i 
# ij_To_k(i, j, N1) = (j-1)*N1+i 
# # Function to convert index k to access u[k] where u is a vector of size N1*N2 to index (i,j) to access the same data 
# # in uMat, where uMat is u reshaped to be a 2D matrix of size N1xN2 
# # function vectorIndexToMatrixIndex(k, N1, N2)
# function k_To_ij(k, N1)
#     j = floor(Int64, (k-1)/N1)
#     i = k-j*N1
#     return (i,j)
# end

# https://en.wikipedia.org/wiki/Laplacian_matrix Laplacian matrix for simple graph
# A adjacency matrix 
function makeLaplacianMatrix(Nplus)
    Ā = spzeros(Int64, Nplus*Nplus, Nplus*Nplus)
    # Top left corner: k=1
    Ā[1, 2] = 1
    Ā[1, ij_To_k(1,2,Nplus)] = 1
    # Bottom left corner: k=Nplus
    Ā[Nplus, Nplus-1] = 1
    Ā[Nplus, ij_To_k(Nplus,2,Nplus)] = 1
    # Top right corner: k=ij_To_k(1,Nplus,Nplus)
    Ā[ij_To_k(1,Nplus,Nplus), ij_To_k(1,Nplus-1,Nplus)] = 1
    Ā[ij_To_k(1,Nplus,Nplus), ij_To_k(2,Nplus,Nplus)] = 1
    # Bottom right corner: k=ij_To_k(Nplus,Nplus,Nplus)
    Ā[ij_To_k(Nplus,Nplus,Nplus), ij_To_k(Nplus-1,Nplus,Nplus)] = 1
    Ā[ij_To_k(Nplus,Nplus,Nplus), ij_To_k(Nplus,Nplus-1,Nplus)] = 1
    # Left edge 
    for i=2:Nplus-1
        Ā[i,ij_To_k(i-1,1,Nplus)] = 1
        Ā[i,ij_To_k(i,2,Nplus)] = 1
        Ā[i,ij_To_k(i+1,1,Nplus)] = 1
    end
    # Right edge 
    for i=2:Nplus-1
        Ā[ij_To_k(i,Nplus,Nplus),ij_To_k(i-1,Nplus,Nplus)] = 1
        Ā[ij_To_k(i,Nplus,Nplus),ij_To_k(i,Nplus-1,Nplus)] = 1
        Ā[ij_To_k(i,Nplus,Nplus),ij_To_k(i+1,Nplus,Nplus)] = 1
    end
    # Top edge 
    for j=2:Nplus-1
        Ā[ij_To_k(1,j,Nplus),ij_To_k(1,j-1,Nplus)] = 1
        Ā[ij_To_k(1,j,Nplus),ij_To_k(2,j,Nplus)] = 1
        Ā[ij_To_k(1,j,Nplus),ij_To_k(1,j+1,Nplus)] = 1
    end
    # Bottom edge 
    for j=2:Nplus-1
        Ā[ij_To_k(Nplus,j,Nplus),ij_To_k(Nplus,j-1,Nplus)] = 1
        Ā[ij_To_k(Nplus,j,Nplus),ij_To_k(Nplus-1,j,Nplus)] = 1
        Ā[ij_To_k(Nplus,j,Nplus),ij_To_k(Nplus,j+1,Nplus)] = 1
    end
    # Internal adjacency
    for i=2:Nplus-1
        for j=2:Nplus-1
            k = ij_To_k(i,j,Nplus)
            Ā[k, ij_To_k(i-1,j,Nplus)] = 1
            Ā[k, ij_To_k(i+1,j,Nplus)] = 1
            Ā[k, ij_To_k(i,j-1,Nplus)] = 1
            Ā[k, ij_To_k(i,j+1,Nplus)] = 1
        end
    end
    # deg degree matrix 
    deg = spdiagm([nnz(Ā[i,:]) for i=1:Nplus*Nplus])
    ∇² = (Ā.-deg)./dx^2
    return ∇²
end 

# https://en.wikipedia.org/wiki/Finite_difference_coefficient
function make∇(Nplus)
    ∇ = spzeros(Float64, Nplus*Nplus, Nplus*Nplus)
    # Top left corner: k=1
    ∇[1, 2] = 1
    # ∇[1, ij_To_k(1,2,Nplus)] = 1
    # Bottom left corner: k=Nplus
    ∇[Nplus, Nplus-1] = -1
    # ∇[Nplus, ij_To_k(Nplus,2,Nplus)] = 1
    # Top right corner: k=ij_To_k(1,Nplus,Nplus)
    # ∇[ij_To_k(1,Nplus,Nplus), ij_To_k(1,Nplus-1,Nplus)] = -1
    ∇[ij_To_k(1,Nplus,Nplus), ij_To_k(2,Nplus,Nplus)] = 1
    # Bottom right corner: k=ij_To_k(Nplus,Nplus,Nplus)
    ∇[ij_To_k(Nplus,Nplus,Nplus), ij_To_k(Nplus-1,Nplus,Nplus)] = -1
    # ∇[ij_To_k(Nplus,Nplus,Nplus), ij_To_k(Nplus,Nplus-1,Nplus)] = -1
    # Left edge 
    for i=2:Nplus-1
        ∇[i,ij_To_k(i-1,1,Nplus)] = -1
        # ∇[i,ij_To_k(i,2,Nplus)] = 1
        ∇[i,ij_To_k(i+1,1,Nplus)] = 1
    end
    # Right edge 
    for i=2:Nplus-1
        ∇[ij_To_k(i,Nplus,Nplus),ij_To_k(i-1,Nplus,Nplus)] = -1
        # ∇[ij_To_k(i,Nplus,Nplus),ij_To_k(i,Nplus-1,Nplus)] = -1
        ∇[ij_To_k(i,Nplus,Nplus),ij_To_k(i+1,Nplus,Nplus)] = 1
    end
    # Top edge 
    for j=2:Nplus-1
        # ∇[ij_To_k(1,j,Nplus),ij_To_k(1,j-1,Nplus)] = -1
        ∇[ij_To_k(1,j,Nplus),ij_To_k(2,j,Nplus)] = 1
        # ∇[ij_To_k(1,j,Nplus),ij_To_k(1,j+1,Nplus)] = 1
    end
    # Bottom edge 
    for j=2:Nplus-1
        # ∇[ij_To_k(Nplus,j,Nplus),ij_To_k(Nplus,j-1,Nplus)] = -1
        ∇[ij_To_k(Nplus,j,Nplus),ij_To_k(Nplus-1,j,Nplus)] = -1
        # ∇[ij_To_k(Nplus,j,Nplus),ij_To_k(Nplus,j+1,Nplus)] = 1
    end
    # Internal adjacency
    for i=2:Nplus-1
        for j=2:Nplus-1
            k = ij_To_k(i,j,Nplus)
            ∇[k, ij_To_k(i-1,j,Nplus)] = -1
            ∇[k, ij_To_k(i+1,j,Nplus)] = 1
            # ∇[k, ij_To_k(i,j-1,Nplus)] = -1
            # ∇[k, ij_To_k(i,j+1,Nplus)] = 1
        end
    end
    ∇ .= ∇./(2*dx)
    return ∇
end

# V_i = ones(Float64, Nplus-1, Nplus)
# V_j = zeros(Float64, Nplus, Nplus-1)
# Vvec = vcat(reshape(V_i, Nplus*(Nplus-1)), reshape(V_j, Nplus*(Nplus-1)))
# V = spdiagm(Vvec)

#%%

# Express model as a matrix operator 
∇ = make∇(Nplus)
∇² = makeLaplacianMatrix(Nplus)
M = ∇².*K₂*K₄/(1+α*h) .- ∇.*β/(1+α*h)

# Initial conditions using Gaussian
fInitial(x, μx, σx, y, μy, σy) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2)
μx = 0.5
σx=0.1
μν=0.0
σν=0.1
# Integrate Gaussian over domain 
integ = hcubature(x -> fInitial(x[1], μx, σx, x[2], μν, σν), [0.0, 0.0], [1.0, 1.0])
uMat = zeros(Float64, Nplus, Nplus)
for νν=1:N
    for xx=1:N
        uMat[νν+Nghost,xx+Nghost] = exp(-(xs[xx]-μx)^2/σx^2 - (νs[νν]-μν)^2/σν^2)
    end
end
uMat .= uMat./integ[1]
u0 = reshape(uMat, Nplus^2)

# Parameters named tuple
p = (D=K₂*K₄/(1+α*h),
    β=β/(1+α*h),
    M=M)

function boundaryConditions!(u, p)
    @unpack D, β, M = p 
    # D∇C - βC = 0 at periphery
    # D(u[x+1]-u[x-1])/(2dx) - βu[x] = 0 
    # u[x+1]-u[x-1] = 2dxβu[x]/D
    # Left edge
    u[1:1:Nplus] = u[1+2*Nplus:1:3*Nplus] - 2*dx*β*u[1+Nplus:1:2*Nplus]/D   
    # Top edge
    u[1:Nplus:end] = u[3:Nplus:end] - 2*dx*β*u[2:Nplus:end]/D   
    # Right edge
    u[end-Nplus+1:1:end] = u[end-3*Nplus+1:1:end-2*Nplus] + 2*dx*β*u[end-2*Nplus+1:1:end-Nplus]/D   
    # Bottom edge
    u[Nplus:Nplus:end] = u[Nplus-2:Nplus:end] + 2*dx*β*u[Nplus-1:Nplus:end]/D
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
ax.xlabel = "ν"
ax.ylabel = "x"
mins = [minimum(u) for u in sol.u[2:end]]
maxs = [maximum(u) for u in sol.u[2:end]]
uInternal = Observable(zeros(N,N))
heatmap!(ax, xs, νs, uInternal, colorrange=(minimum(mins), maximum(maxs)), colormap=:navia)
record(fig,"test2D.mp4", 1:length(sol.t); framerate=10) do i
    uInternal[] = reshape(sol.u[i],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
    uInternal[] = uInternal[]
end



#%%

# masses = zeros(length(sol.u))
# for t=1:length(sol.u)
#     for i=1+Nghost:Nplus-Nghost-1
#         masses[t] += dx*0.5*(sol.u[t][i]+sol.u[t][i+1])
#     end
# end
# fig2 = Figure(size=(1000,1000))
# ax2 = CairoMakie.Axis(fig2[1, 1])
# lines!(ax2, sol.t, masses)
# ylims!(ax2, (0.95,1.05))
# ax2.xlabel = "Time"
# ax2.ylabel = "Mass"
# display(fig2)



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
