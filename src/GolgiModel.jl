module GolgiModel

using FromFile
using DrWatson
using CairoMakie
using Dates

@from "Visualise.jl" using Visualise
@from "DeterministicModel.jl" using DeterministicModel
@from "StochasticModel.jl" using StochasticModel

function golgiModel(nMax,tMax,volume,∅ToCis,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅)

    # Package parameters into p array
    p = nMax,∅ToCis,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅
    
    @info "Solving deterministic model"
    deterministicSol = solveDeterministic(nMax,tMax,p)

    @info "Solving stochastic model"
    stochasticSol = stochasticModel(nMax,tMax,volume,p)

    params = @strdict nMax tMax
    fileName = savename(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"),params,"jld2",connector="")    
    safesave(datadir("sims",fileName),@strdict deterministicSol stochasticSol params)
    
    @info "Visualising results"
    visualise(nMax,stochasticSol,params,deterministicSol,volume)

    return nothing

end

export golgiModel

end
    