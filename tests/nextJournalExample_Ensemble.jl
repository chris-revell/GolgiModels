#https://nextjournal.com/bebi5009/gillespie-julia

using Random # randexp()
using StatsBase # Weights() and sample()
using Plots
Plots.gr(fmt=:png)

#=
Stochastic chemical reaction: Gillespie Algorithm (direct method)
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3
=#
function ssa_direct(model, u0, tend, p, stoich; tstart=zero(tend))
    t = tstart   # Current time
    ts = [t]     # Time points
    u = copy(u0) # Current state
    us = copy(u) # Record of states
    while t < tend
        a = model(u, p, t)               # propensities
        dt = randexp() / sum(a)          # Time step
        du = sample(stoich, Weights(a))  # Choose the stoichiometry for the next reaction
        u .+= du  # Update state
        t += dt   # Update time

        us = [us u]  # Append state variable to record
        push!(ts, t) # Append time point to record
    end
    # Make column as variables, rows as observations
    us = collect(us')
    return (t = ts, u = us)
end

#=
Reaction of A <-> B with rate constants k1 & k2
=#
"Propensity model for this reaction"
model(u, p, t) = [p.k1 * u[1],  p.k2 * u[2]]

parameters = (k1=1.0, k2=0.5, stoich=[[-1, 1], [1, -1]])
u0 = [200, 0]
tend = 10.0

# Repeated simulations
Plots.plot()
for i in 1:50
    sol = ssa_direct(model, u0, tend, parameters, parameters.stoich)
    Plots.plot!(sol.t, sol.u, linecolor=[:blue :red], linealpha=0.1, label=false)
end

Plots.plot!(xlabel="time", ylabel="# of molecules", title = "SSA (1st reaction method) ensemble")
