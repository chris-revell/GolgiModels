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

    params = @strdict nMax tMax volume ∅ToCis cisAgg cisSplit cisToMed medToCis medAgg medSplit medToTran tranToMed tranAgg tranSplit tranTo∅
    fileName = savename(Dates.format(Dates.now(),"mm-dd-HH-MM"),params,connector="")    
    safesave(datadir("sims","$fileName.jld2"),@strdict deterministicSol stochasticSol params)
    
    @info "Visualising results"
    visualise(fileName,nMax,stochasticSol,params,deterministicSol,volume)

    return nothing

end

export golgiModel

end
    