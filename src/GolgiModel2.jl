#
#  GolgiModel.jl
#  GolgiModels
#
#  Created by Christopher Revell on dd/mm/yyyy.

module GolgiModel2

using FromFile
using DrWatson
using Dates
using Statistics
using Catalyst
using ModelingToolkit
using DifferentialEquations
using JumpProcesses

@from "$(projectdir("src","Visualise.jl"))" using Visualise
@from "$(projectdir("src","AllReactions.jl"))" using AllReactions

# nMax= maximum compartment size
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

function golgiModel2(nMax,tMax,volume,k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀,k₁₁,nOutput)

    # Symbolic system parameters: time and rate constants 
    @parameters t k[1:12]
    # Symbolic system variables: Vector of number/concentration for cis, medial, and trans
    @variables C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax]
    
    # reactions = allReactions(nMax,C,M,T,K₀,K₁,K₂,K₃,K₄,K₅,K₆,K₇,K₈,K₉,K₁₀,K₁₁)
    
    # Set up reaction system object 
    # @named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], [K₀,K₁,K₂,K₃,K₄,K₅,K₆,K₇,K₇,K₈,K₉,K₁₀,K₁₁])
    
    system = allReactions(nMax,C,M,T,k,t)

    # Solving stochastic model
    @info "Solving stochastic model"
    # Map symbolic rate constants to values for stochastic model 
    p = Pair.(collect(k),[k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀,k₁₁])
    u₀Map = Pair.([collect(C); collect(M); collect(T)], zeros(Int32,3*nMax))
    # Convert to jump problem to solve 
    discreteprob  = DiscreteProblem(system, u₀Map, (0.0,tMax), p)
    jumpProblem   = JumpProblem(system, discreteprob, Direct(), save_positions=(false,false))
    integ = init(odeProblem,KenCarp3())
    stochasticSol = solve(jumpProblem, SSAStepper(), saveat=tMax/nOutput)
    

    @info "Solving deterministic model"
    # Map symbolic rate constants to values for stochastic model 
    p2 = Pair.(collect(k),[k₀/volume, k₁*volume, k₂, k₃, k₄, k₅*volume, k₆, k₇, k₈, k₉*volume, k₁₀, k₁₁])
    # Map symbolic state vectors to float vector for stochastic model 
    u₀Map = Pair.([collect(C); collect(M); collect(T)], zeros(Float32,3*nMax))
    odeProblem = ODEProblem(system,u₀Map,(0.0,tMax),p2)
    integ = init(odeProblem,KenCarp3())
    deterministicSol = solve!(integ) #solve(odeProblem, KenCarp3(), saveat=tMax/nOutput)
    
    # Calculate time average for 101 time points correspoding to 101 frames in the visualisation 
    windowLength = nOutput÷100
    stochasticTimeAverages = fill(zeros(nMax*3),100+1)
    stochasticTimeAverages[2:end] = [mean(stochasticSol.u[i-windowLength:i]) for i=windowLength+1:windowLength:length(stochasticSol.u)]

    # Save data to file 
    params = @strdict nMax tMax volume k₀ k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁ nOutput
    fileName = savename(Dates.format(Dates.now(),"mm-dd-HH-MM"),params,connector="")
    # @info "Saving data as $fileName.jld2"
    # safesave(datadir("sims","$fileName.jld2"),@strdict deterministicSol stochasticSol params)
    
    # Visualise results 
    @info "Visualising results; saving as $fileName.mp4"
    visualise(fileName,nMax,volume,stochasticSol,stochasticTimeAverages,deterministicSol,nOutput,windowLength)

    return nothing

end

export golgiModel2

end
    