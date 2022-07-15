using ModelingToolkit
using Catalyst
using LinearAlgebra
using DiffEqBase
using DiffEqJump
using CairoMakie
using SpecialFunctions

## Parameter
nMax      = 10                      # maximum cluster size
nMonomers = 10000                   # No. of monomers initially


#
# Kb # Budding rate per patches of membrane in a compartment
# Km # Biochemical conversion rate for each patch of membrane (cis->medial->trans)
# Kf # Maximum fusion rate between two compartments
# J # Injection rate of cis-vesicles from the ER
kb = 1.0 # Budding rate normalized by the fusion rate (Kb/Kf )
km = 1.0 # Biochemical conversion rate normalized by the fusion rate (Km/Kf )
j = 1.0 # Injection rate normalized by the fusion rate (J/Kf )
# α # Fraction of active species (driving homotypic fusion) in the ER (cis) or the TGN (trans)
# αᵢ # When different between boundaries, a in the boundary i (i equals ER or TGN)


# Find all possible monomer and oligomer pairs that result in oligomers with size <= nMax
reactantsAggregating = Tuple[]
for i=1:nMax-2
    for j=i:nMax-1
        if i+j<=nMax
            push!(reactantsAggregating,(i,j))
        end
    end
end
nReactionsAggregating = length(reactantsAggregating)
# Reaction rates proportional to sum of reactant volumes
ratesAggregating = ones(nReactionsAggregating) # [sum(a) for a in reactantsAggregating]./V  # dividing by system volume as its a bi-molecular reaction chain

reactantsSplitting = collect(2:nMax)
nReactionsSplitting = length(reactantsSplitting)
ratesSplitting = ones(nReactionsSplitting) # (reactantsSplitting.^2)./V

rateInjection = [1.0]
rateRemoval = [1.0]

allRates = [rateInjection; ratesAggregating; ratesSplitting; rateRemoval]
nReactionsTotal = length(allRates)

# time-span
tspan = (0.0 ,100000.0)

# initial condition of monomers
u₀    = zeros(Int64, nMax)
# u₀[1] = nMonomers
# let
#     i=0
#     while i<nMonomers
#         index = ceil(Int64,nMax*rand())
#         u₀[index] += 1
#         i+=index
#     end
# end

# state variables are X, pars stores rate parameters for each rx
@parameters t
@variables k[1:nReactionsTotal]  X[1:nMax](t)
pars = Pair.(collect(k), allRates)
u₀map = Pair.(collect(X), u₀)   # map variable to its initial value

# vector to store the Reactions
rx = []
# Add injection reaction
push!(rx, Reaction(k[1], nothing, [X[1]], nothing, [1]))
# Add aggregation reactions
for (n,pair) in enumerate(reactantsAggregating)
    push!(rx, Reaction(k[n+1], [ X[pair[1]], X[pair[2]] ], [ X[sum(pair)] ], [1, 1], [1]))
end
# Add splitting reactions
for (n,reactant) in enumerate(reactantsSplitting)
    push!(rx, Reaction(k[n+nReactionsAggregating+1], [ X[reactant] ], [ X[reactant-1], X[1] ], [1], [1, 1]))
end
# Add removal reaction
push!(rx, Reaction(k[end], [X[1]], nothing, [1], nothing))


@named rs = ReactionSystem(rx, t, collect(X), collect(k))

# solving the system
jumpsys = convert(JumpSystem, rs)
dprob   = DiscreteProblem(jumpsys, u₀map, tspan, pars)
jumpProblem   = JumpProblem(jumpsys, dprob, Direct(), save_positions=(false,false))
jsol    = solve(jumpProblem, SSAStepper(), saveat = tspan[2]/30)
t = jsol.t

solMatrix = reduce(hcat,jsol.u)

fig = Figure()
ax = Axis(fig[1,1])
ax.xlabel="Time"
ax.ylabel="Oligomer count"
for i=1:nMax
    lines!(t, solMatrix[i,:], linewidth=4, label="$i")
end
lines!(t, solMatrix[nMax,:], linewidth=4, label="$nMax")
# allLines = [lines!(t, l) for l in eachrow(solMatrix)]
axislegend(ax,title="Oligomers")
display(fig)

for i=1:size(solMatrix)[2]
    display(sum(solMatrix[:,i].*collect(1:nMax)))
end
