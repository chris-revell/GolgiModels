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

@from "Visualise.jl" using Visualise
@from "DeterministicModel.jl" using DeterministicModel
@from "StochasticModel.jl" using StochasticModel

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
    
    @info "Solving deterministic model"
    deterministicSol = solveDeterministic(nMax,tMax,p)
    
    @info "Solving stochastic model"
    stochasticSol = stochasticModel(nMax,tMax,volume,p)
    windowLength = 1000
    stochasticTimeAverages = fill(zeros(nMax*3),101)
    stochasticTimeAverages[2:end] = [mean(stochasticSol.u[i-windowLength:i]) for i=windowLength+1:windowLength:length(stochasticSol.u)]

    params = @strdict nMax tMax volume k₀ k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁
    fileName = savename(Dates.format(Dates.now(),"mm-dd-HH-MM"),params,connector="")    
    safesave(datadir("sims","$fileName.jld2"),@strdict deterministicSol stochasticSol params)
    
    @info "Visualising results"
    visualise(fileName,nMax,stochasticSol,params,deterministicSol,stochasticTimeAverages,volume)

    return nothing

end

export golgiModel

end
    