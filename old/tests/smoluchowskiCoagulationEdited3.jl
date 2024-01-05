using ModelingToolkit
using Catalyst
using LinearAlgebra
using DiffEqBase
using DiffEqJump
using CairoMakie
using SpecialFunctions

## Parameter
nMax = 10                       # maximum cluster size
Vₒ = (4π/3)*(10e-06*100)^3   # volume of a monomers in cm³
Nₒ = 1e-06/Vₒ                # initial conc. = (No. of init. monomers) / bulk volume
uₒ = 10000                   # No. of monomers initially
V = uₒ/Nₒ                    # Bulk volume of system in cm³
B = 1.53e03                # s⁻¹

# Find all possible monomer and oligomer pairs that result in oligomers with size <= nMax
pairs = Tuple[]
for i=1:nMax-2
    for j=i:nMax-1
        if i+j<=nMax
            push!(pairs,(i,j))
        end
    end
end
nReactions = length(pairs)

# Reaction rates proportional to sum of reactant volumes
rates = B*Vₒ.*[sum(a) for a in pairs]./V  # dividing by system volume as its a bi-molecular reaction chain
rates.*=rand()
# state variables are X, pars stores rate parameters for each rx
@parameters t
@variables k[1:nReactions]  X[1:nMax](t)

pars = Pair.(collect(k), rates)

# time-span
tspan = (0. ,2000.)

 # initial condition of monomers
u₀    = zeros(Int64, nMax)
u₀[1] = uₒ
u₀map = Pair.(collect(X), u₀)   # map variable to its initial value


# vector to store the Reactions in
rx = []
for (n,pair) in enumerate(pairs)
    push!(rx, Reaction(k[n], [X[pair[1]], X[pair[2]]], [X[sum(pair)]], [1, 1], [1]))
end
@named rs = ReactionSystem(rx, t, collect(X), collect(k))


# solving the system
jumpsys = convert(JumpSystem, rs)
dprob   = DiscreteProblem(jumpsys, u₀map, tspan, pars)
jumpProblem   = JumpProblem(jumpsys, dprob, Direct(), save_positions=(false,false))
jsol    = solve(jumpProblem, SSAStepper(), saveat = tspan[2]/30)

t   = jsol.t

solMatrix = reduce(hcat,jsol.u)

fig = Figure()#size=(00,1000))
ax = Axis(fig[1,1])
ax.xlabel="Time"
ax.ylabel="Oligomer count"
for i=1:nMax
    lines!(t, solMatrix[i,:], label="$i")
end
axislegend(ax,title="Oligomers")#, merge = merge, unique = unique)
# fig[1, 2] = Legend(fig, ax, "Oligomers", framevisible = false)
display(fig)

for i=1:size(solMatrix)[2]
    display(sum(solMatrix[:,i].*collect(1:nMax)))
end
