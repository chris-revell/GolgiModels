#https://nextjournal.com/bebi5009/gillespie-julia

using Random # randexp()
using StatsBase # Weights() and sample()
using Plots
Plots.gr(fmt=:png)

#=
Stochastic chemical reaction: Gillespie Algorithm (first reaction method)
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3
=#
function ssa_first(model, u0, tend, p, stoich; tstart=zero(tend))
    t = tstart   # Current time
    ts = [t]     # Time points
    u = copy(u0) # Current state
    us = copy(u) # Record of states
    while t < tend
        a = model(u, p, t)  # propensities of reactions
        # dts from all reactions
        dts = randexp(length(a)) ./ a
        # Choose the reaction
        i = argmin(dts)
        dt = dts[i]
        du = stoich[i]
        # Update state and time
        u .+= du
        t += dt
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

sol2 = ssa_first(model, u0, tend, parameters, parameters.stoich)

Plots.plot(sol2.t, sol2.u, xlabel="time", ylabel="# of molecules", title = "SSA (1st reaction method)", label=["A" "B"])
