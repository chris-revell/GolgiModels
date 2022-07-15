module GolgiModel

using FromFile
using DrWatson
using ModelingToolkit
using Catalyst
using LinearAlgebra
using DifferentialEquations
using JumpProcesses
using CairoMakie
using Dates

@from "AllReactions.jl" using AllReactions
@from "Visualise.jl" using Visualise

function golgiModel(nMax,tMax,volume,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅)

    # nMax = maximum compartment size
    # nReactionsTotal =  Injection +
    #                   cis aggregation + cis splitting + 
    #                   cis to medial + medial to cis + 
    #                   medial aggregation + medial splitting + 
    #                   medial to trans + trans to medial +
    #                   trans aggregation + trans splitting
    #                   + removal 
    nReactionsTotal = 1+2*(nMax-1)+2+2*(nMax-1)+2+2*(nMax-1)+1

    # state variables are X, pars stores rate parameters for each reaction
    # species labels 1:nMax => cis vesicle size counts 
    #                nMax+1:2*nMax => medial vesicle size counts 
    #                2*nMax+1:3*nMax => trans vesicle size counts 
    @parameters t
    @variables k[1:nReactionsTotal]  X[1:3*nMax](t)
    
    reactants,products,rates,stoichiometryIn,stoichiometryOut = allReactions(nMax,X,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅)
    # vector to store the Reactions
    reactions = []
    for i=1:nReactionsTotal
        push!(reactions, Reaction(k[i], reactants[i], products[i], stoichiometryIn[i], stoichiometryOut[i]))
    end

    # initial condition of monomers
    u₀    = zeros(Int64, 3*nMax)
    # pair initial condition to state vector 
    u₀map = Pair.(collect(X), u₀)
    # pair rates to parameters vector 
    pars = Pair.(collect(k), rates)
    # time-span
    tspan = (0.0,tMax)

    # Set up reaction system object 
    @named system = ReactionSystem(reactions, t, collect(X), collect(k))

    # solving the system
    # jumpsys = convert(JumpSystem, system)
    discreteprob   = DiscreteProblem(system, u₀map, tspan, pars)
    jumpProblem   = JumpProblem(system, discreteprob, Direct(),save_positions=(false,false)) # Converts system to a set of MassActionJumps
    jsol    = solve(jumpProblem, SSAStepper(), saveat=tMax/100)

    params = @strdict nMax tMax
    fileName = savename(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"),params,"jld2",connector="")    
    safesave(datadir("sims",fileName),@strdict jsol params)
    visualise(nMax,jsol,params)

    return nothing

end

export golgiModel

end
    