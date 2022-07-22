module StochasticModel

using FromFile
using ModelingToolkit
using Catalyst
using LinearAlgebra
using DifferentialEquations
using JumpProcesses

@from "AllReactions.jl" using AllReactions

function stochasticModel(nMax,tMax,volume,p)
    
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
    @variables k[1:nReactionsTotal]  X[1:nMax*3](t)
    
    reactions, rates = allReactions(nReactionsTotal,k,X,volume,p)

    # initial condition of monomers
    u₀    = zeros(Int64, nMax*3)
    # pair initial condition to state vector 
    u₀map = Pair.(collect(X), u₀)
    # pair rates to parameters vector 
    pars = Pair.(collect(k), rates)
    # time-span
    tspan = (0.0,tMax)

    # Set up reaction system object 
    @named system = ReactionSystem(reactions, t, collect(X), collect(k))

    # solving the system    
    discreteprob  = DiscreteProblem(system, u₀map, tspan, pars)
    jumpProblem   = JumpProblem(system, discreteprob, Direct(),save_positions=(false,false)) # Converts system to a set of MassActionJumps
    stochasticSol = solve(jumpProblem, SSAStepper(), saveat=tMax/100)

    return stochasticSol

end

export stochasticModel

end