using DrWatson
@quickactivate
using ModelingToolkit
using Catalyst
using LinearAlgebra
using DiffEqBase
using DiffEqJump
using Plots
using SpecialFunctions

## Parameter
N = 10                       # maximum cluster size
Vₒ = (4π/3)*(10e-06*100)^3   # volume of a monomers in cm³
Nₒ = 1e-06/Vₒ                # initial conc. = (No. of init. monomers) / bulk volume
uₒ = 10000                   # No. of monomers initially
V = uₒ/Nₒ                    # Bulk volume of system in cm³
B = 1.53e03                  # Reaction rate parameter, s⁻¹

# Find all possible monomer and oligomer pairs that result in oligomers with size <= N
pairs = Tuple[]
for i=1:N-2
   for j=i:N-1
       i+j<=N ? push!(pairs,(i,j)) : nothing
   end
end
nReactions = length(pairs)

# Additive kernel
# Reaction rate constants proportional to sum of reactant molecular volumes
rates = Float64[]
for i=1:nReactions
    push!(rates,B*(pairs[i][1]+pairs[i][2])/V)
end
# # Constant kernel
# C = 1.84e-04               # cm³ s⁻¹
# rates = fill(C/V, nReactions)


## state variables are X, pars stores rate parameters for each rx
@parameters t
@variables k[1:nReactions]  X[1:N](t)
pars = Pair.(collect(k), rates)

# time-span
tspan = (0. ,2000.)

 # initial condition of monomers
u₀    = zeros(Int64, N)
u₀[1] = uₒ
u₀map = Pair.(collect(X), u₀)    # map variable to its initial value

## vector to store the Reactions in
rx = []
for n = 1:nReactions
    # for clusters of the same size, double the rate
    if (pairs[n][1] == pairs[n][2])
        # Reaction(rate, [reactants], [products], [stoichiometry in], [stoichiometry out])
        push!(rx, Reaction(pars[n], [X[pairs[n][1]]], [X[pairs[n][1]+pairs[n][2]]], [2], [1]))
    else
        push!(rx, Reaction(pars[n], [X[pairs[n][1]], X[pairs[n][2]]], [X[pairs[n][1]+pairs[n][2]]], [1,1], [1]))
    end
end

# rx = []
# for n = 1:nReactions
#     # for clusters of the same size, double the rate
#     if (pairs[n][1] == pairs[n][2])
#         # Reaction(rate, reactant1, reactant2, [products], stoichiometry in, stoichiometry out )
#         push!(rx, Reaction(k[n], [X[vᵢ[n]]], [X[sum_vᵢvⱼ[n]]], [2], [1]))
#     else
#         push!(rx, Reaction(k[n], [X[vᵢ[n]], X[vⱼ[n]]], [X[sum_vᵢvⱼ[n]]],[1, 1], [1]))
#     end
# end
@named rs = ReactionSystem(rx, t, collect(X), collect(k))

## solving the system
jumpsys = convert(JumpSystem, rs)
dprob   = DiscreteProblem(jumpsys, u₀map, tspan, pars)
jprob   = JumpProblem(jumpsys, dprob, Direct(), save_positions=(false,false))
jsol    = solve(jprob, SSAStepper(), saveat = tspan[2]/30)

## Results for first three polymers...i.e. monomers, dimers and trimers
# Results for first three polymers...i.e. monomers, dimers and trimers
v_res = [1;2;3]

# comparison with analytical solution
# we only plot the stochastic solution at a small number of points
# to ease distinguishing it from the exact solution
t   = jsol.t
sol = zeros(length(v_res), length(t))
# if i == 1
    ϕ = @. 1 - exp(-B*Nₒ*Vₒ*t)
    for j in v_res
        sol[j,:] = @. Nₒ*(1 - ϕ)*(((j*ϕ)^(j-1))/gamma(j+1))*exp(-j*ϕ)
    end
# elseif i == 2
#     ϕ = @. (C*Nₒ*t)
#     for j in v_res
#         sol[j,:] = @. 4Nₒ*((ϕ^(j-1))/((ϕ + 2)^(j+1)))
#     end
# end

# plotting normalised concentration vs analytical solution
default(lw=2, xlabel="Time (sec)")
scatter(ϕ, jsol(t)[1,:]/uₒ, label="X1 (monomers)", markercolor=:blue)
plot!(ϕ, sol[1,:]/Nₒ, line = (:dot,4,:blue), label="Analytical sol--X1")

scatter!(ϕ, jsol(t)[2,:]/uₒ, label="X2 (dimers)", markercolor=:orange)
plot!(ϕ, sol[2,:]/Nₒ, line = (:dot,4,:orange), label="Analytical sol--X2")

scatter!(ϕ, jsol(t)[3,:]/uₒ, label="X3 (trimers)", markercolor=:purple)
plot!(ϕ, sol[3,:]/Nₒ, line = (:dot,4,:purple), label="Analytical sol--X3",
      ylabel = "Normalized Concentration")
