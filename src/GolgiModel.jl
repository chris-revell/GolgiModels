#
#  GolgiModel.jl
#  GolgiModels
#
#  Created by Christopher Revell on 31/01/2021.

module GolgiModel

using FromFile
using DrWatson
using CairoMakie
using Dates
using Statistics
using Catalyst
using ModelingToolkit
using DifferentialEquations
using JumpProcesses

@from "Visualise.jl" using Visualise
@from "DeterministicModel.jl" using DeterministicModel
@from "StochasticModel.jl" using StochasticModel
@from "AllReactions.jl" using AllReactions

function golgiModel(nMax,tMax,volume,k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀,k₁₁)

    # k₀  = ∅ToCis   
    # k₁  = cisAgg   
    # k₂  = cisSplit 
    # k₃  = cisToMed 
    # k₄  = medToCis 
    # k₅  = medAgg   
    # k₆  = medSplit 
    # k₇  = medToTran
    # k₈  = tranToMed
    # k₉  = tranAgg  
    # k₁₀ = tranSplit
    # k₁₁ = tranTo∅  
    # Package parameters into p array
    p = nMax,k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀,k₁₁

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
    u₀Stochastic    = zeros(Int64, nMax*3)
    u₀Deterministic = zeros(Float64, nMax*3)
    # pair initial condition to state vector 
    u₀MapS = Pair.(collect(X), u₀Stochastic)
    # pair rates to parameters vector 
    pars = Pair.(collect(k), rates)
    # time-span
    tspan = (0.0,tMax)

    # Set up reaction system object 
    @named system = ReactionSystem(reactions, t, collect(X), collect(k))

    # solving the system    
    @info "Solving stochastic model"
    discreteprob  = DiscreteProblem(system, u₀MapS, tspan, pars)
    jumpProblem   = JumpProblem(system, discreteprob, Direct(),save_positions=(false,false)) # Converts system to a set of MassActionJumps
    stochasticSol = solve(jumpProblem, SSAStepper(), saveat=tMax/100000)
    
    @info "Solving deterministic model"
    odeProblem = ODEProblem(system,u₀Deterministic,tspan,p)
    deterministicSol = solve(odeProblem, saveat=tMax/100000)
    
    windowLength = 1000
    stochasticTimeAverages = fill(zeros(nMax*3),101)
    stochasticTimeAverages[2:end] = [mean(stochasticSol.u[i-windowLength:i]) for i=windowLength+1:windowLength:length(stochasticSol.u)]

    params = @strdict nMax tMax volume k₀ k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁
    fileName = savename(Dates.format(Dates.now(),"mm-dd-HH-MM"),params,connector="")    
    safesave(datadir("sims","$fileName.jld2"),@strdict deterministicSol stochasticSol params)
    
    @info "Visualising results"
    visualise(fileName,nMax,volume,stochasticSol,stochasticTimeAverages,deterministicSol)

    return nothing

end

export golgiModel

end
    