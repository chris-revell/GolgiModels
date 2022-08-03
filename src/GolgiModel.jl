#
#  GolgiModel.jl
#  GolgiModels
#
#  Created by Christopher Revell on dd/mm/yyyy.

module GolgiModel

using FromFile
using DrWatson
using Dates
using Statistics
using Catalyst
using ModelingToolkit
using DifferentialEquations
using JumpProcesses

@from "Visualise.jl" using Visualise
@from "AllReactions.jl" using AllReactions

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

function golgiModel(nMax,tMax,volume,k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀,k₁₁,nOutput)

    # Symbolic system parameters: time and rate constants 
    @parameters t K₀ K₁ K₂ K₃ K₄ K₅ K₆ K₇ K₈ K₉ K₁₀ K₁₁
    # Symbolic system variables: Vector of number/concentration for cis, medial, and trans
    @variables C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax]
    
    reactions = allReactions(nMax,C,M,T,K₀,K₁,K₂,K₃,K₄,K₅,K₆,K₇,K₈,K₉,K₁₀,K₁₁)
    
    # Set up reaction system object 
    @named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], [K₀,K₁,K₂,K₃,K₄,K₅,K₆,K₇,K₇,K₈,K₉,K₁₀,K₁₁])
    
    
    # Solving stochastic model
    @info "Solving stochastic model"
    # Map symbolic rate constants to values for stochastic model 
    p = [:K₀=>k₀, :K₁=>k₁, :K₂=>k₂, :K₃=>k₃, :K₄=>k₄, :K₅=>k₅, :K₆=>k₆, :K₇=>k₇, :K₈=>k₈, :K₉=>k₉, :K₁₀=>k₁₀, :K₁₁=>k₁₁]
    # p = [:K₀=>k₀/volume, :K₁=>k₁*volume, :K₂=>k₂, :K₃=>k₃, :K₄=>k₄, :K₅=>k₅*volume, :K₆=>k₆, :K₇=>k₇, :K₈=>k₈, :K₉=>k₉*volume, :K₁₀=>k₁₀, :K₁₁=>k₁₁]
    # Map symbolic state vectors to integer vector for stochastic model 
    u₀ = zeros(Int64,3*nMax)
    u₀Map = Pair.([collect(C); collect(M); collect(T)],u₀)
    # Convert to jump problem to solve 
    discreteprob  = DiscreteProblem(system, u₀Map, (0.0,tMax), p)
    jumpProblem   = JumpProblem(system, discreteprob, Direct(),save_positions=(false,false)) # Converts system to a set of MassActionJumps
    stochasticSol = solve(jumpProblem, SSAStepper(), saveat=tMax/nOutput)
    

    @info "Solving deterministic model"
    # Map symbolic rate constants to values for stochastic model 
    p2 = [:K₀=>k₀/volume, :K₁=>k₁*volume, :K₂=>k₂, :K₃=>k₃, :K₄=>k₄, :K₅=>k₅*volume, :K₆=>k₆, :K₇=>k₇, :K₈=>k₈, :K₉=>k₉*volume, :K₁₀=>k₁₀, :K₁₁=>k₁₁]
    # Map symbolic state vectors to float vector for stochastic model 
    u₀ = zeros(Float64,3*nMax)
    u₀Map = Pair.([collect(C); collect(M); collect(T)],u₀)    
    odeProblem = ODEProblem(system,u₀Map,(0.0,tMax),p2)
    deterministicSol = solve(odeProblem, saveat=tMax/nOutput)
    
    # Calculate time average for 101 time points correspoding to 101 frames in the visualisation 
    windowLength = nOutput÷100
    stochasticTimeAverages = fill(zeros(nMax*3),100+1)
    stochasticTimeAverages[2:end] = [mean(stochasticSol.u[i-windowLength:i]) for i=windowLength+1:windowLength:length(stochasticSol.u)]

    # Save data to file 
    params = @strdict nMax tMax volume k₀ k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁ nOutput
    fileName = savename(Dates.format(Dates.now(),"mm-dd-HH-MM"),params,connector="")
    @info "Saving data as $fileName.jld2"
    safesave(datadir("sims","$fileName.jld2"),@strdict deterministicSol stochasticSol params)
    
    # Visualise results 
    @info "Visualising results; saving as $fileName.mp4"
    visualise(fileName,nMax,volume,stochasticSol,stochasticTimeAverages,deterministicSol,nOutput,windowLength)

    return nothing

end

export golgiModel

end
    