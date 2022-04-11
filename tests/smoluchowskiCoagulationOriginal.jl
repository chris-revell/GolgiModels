using ModelingToolkit, Catalyst, LinearAlgebra
using DiffEqBase, DiffEqJump
using Plots, SpecialFunctions

## Parameter
N = 10                       # maximum cluster size
Vₒ = (4π/3)*(10e-06*100)^3   # volume of a monomers in cm³
Nₒ = 1e-06/Vₒ                # initial conc. = (No. of init. monomers) / bulk volume
uₒ = 10000                   # No. of monomers initially
V = uₒ/Nₒ                    # Bulk volume of system in cm³

integ(x) = Int(floor(x))
n        = integ(N/2)
nr       = N%2 == 0 ? (n*(n + 1) - n) : (n*(n + 1)) # No. of forward reactions


# possible pairs of reactant multimers
pair = []
for i = 2:N
    push!(pair,[1:integ(i/2)  i .- (1:integ(i/2))])
end
pair = vcat(pair...)
vᵢ = @view pair[:,1]   # Reactant 1 indices
vⱼ = @view pair[:,2]   # Reactant 2 indices
volᵢ = Vₒ*vᵢ           # cm⁻³
volⱼ = Vₒ*vⱼ           # cm⁻³
sum_vᵢvⱼ = @. vᵢ + vⱼ  # Product index


# set i to  1 for additive kernel, 2  for constant
i = 1
if i==1
    B = 1.53e03                # s⁻¹
    kv = @. B*(volᵢ + volⱼ)/V  # dividing by volume as its a bi-molecular reaction chain
elseif i==2
    C = 1.84e-04               # cm³ s⁻¹
    kv = fill(C/V, nr)
end


# state variables are X, pars stores rate parameters for each rx
@parameters t
@variables k[1:nr]  X[1:N](t)
pars = Pair.(collect(k), kv)

# time-span
if i == 1
    tspan = (0. ,2000.)
elseif i == 2
    tspan = (0. ,350.)
end

 # initial condition of monomers
u₀    = zeros(Int64, N)
u₀[1] = uₒ
u₀map = Pair.(collect(X), u₀)   # map variable to its initial value


# vector to store the Reactions in
rx = []
for n = 1:nr
    # for clusters of the same size, double the rate
    if (vᵢ[n] == vⱼ[n])
        push!(rx, Reaction(k[n], [X[vᵢ[n]]], [X[sum_vᵢvⱼ[n]]], [2], [1]))
    else
        push!(rx, Reaction(k[n], [X[vᵢ[n]], X[vⱼ[n]]], [X[sum_vᵢvⱼ[n]]],
                           [1, 1], [1]))
    end
end
@named rs = ReactionSystem(rx, t, collect(X), collect(k))


# solving the system
jumpsys = convert(JumpSystem, rs)
dprob   = DiscreteProblem(jumpsys, u₀map, tspan, pars)
jprob   = JumpProblem(jumpsys, dprob, Direct(), save_positions=(false,false))
jsol    = solve(jprob, SSAStepper(), saveat = tspan[2]/30)


# Results for first three polymers...i.e. monomers, dimers and trimers
v_res = [1;2;3]

# comparison with analytical solution
# we only plot the stochastic solution at a small number of points
# to ease distinguishing it from the exact solution
t   = jsol.t
sol = zeros(length(v_res), length(t))
if i == 1
    ϕ = @. 1 - exp(-B*Nₒ*Vₒ*t)
    for j in v_res
        sol[j,:] = @. Nₒ*(1 - ϕ)*(((j*ϕ)^(j-1))/gamma(j+1))*exp(-j*ϕ)
    end
elseif i == 2
    ϕ = @. (C*Nₒ*t)
    for j in v_res
        sol[j,:] = @. 4Nₒ*((ϕ^(j-1))/((ϕ + 2)^(j+1)))
    end
end

# plotting normalised concentration vs analytical solution
default(lw=2, xlabel="Time (sec)")
scatter(ϕ, jsol(t)[1,:]/uₒ, label="X1 (monomers)", markercolor=:blue)
plot!(ϕ, sol[1,:]/Nₒ, line = (:dot,4,:blue), label="Analytical sol--X1")

scatter!(ϕ, jsol(t)[2,:]/uₒ, label="X2 (dimers)", markercolor=:orange)
plot!(ϕ, sol[2,:]/Nₒ, line = (:dot,4,:orange), label="Analytical sol--X2")

scatter!(ϕ, jsol(t)[3,:]/uₒ, label="X3 (trimers)", markercolor=:purple)
plot!(ϕ, sol[3,:]/Nₒ, line = (:dot,4,:purple), label="Analytical sol--X3",
      ylabel = "Normalized Concentration")

solMatrix = reduce(hcat,jsol.u)

for i=1:size(solMatrix)[2]
    display(sum(solMatrix[:,i].*collect(1:N)))
end
