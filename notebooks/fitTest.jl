
using CairoMakie
using DifferentialEquations
using DrWatson
using UnPack
using Catalyst
using FromFile
using DiffEqCallbacks
using Distributions
using Turing
using DataStructures

nMax = 20           # Max compartment size /vesicles

# Catalyst system setup
# Symbolic system parameters: rate constants 
@parameters k[1:4]
@variables t
# Symbolic system variables: cis, medial, and trans compartment size counts 
@species C(t)[1:nMax]
# Use these parameters and variables to define a reaction system 
# Initialise reaction system (within array so it's a mutable object for ease of later updating)
reactions = []
push!(reactions, Reaction(k[1], nothing, [C[1]]))            # ∅->c₁
push!(reactions, Reaction(k[2], [C[1]], [C[2]], [2], [1]))   # 2c₁->c₂
push!(reactions, Reaction(k[3], [C[2]], [C[1]], [1], [2]))   # c₂->2c₁
for i=2:nMax-1
    push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]])) # c₁+cₙ->cₙ₊₁ for 2<=n<nMax
end
for i=3:nMax
    push!(reactions, Reaction(k[3], [C[i]], [C[i-1],C[1]]))  # cₙ->c₁+cₙ₋₁ for 3<=n<=nMax
end
push!(reactions, Reaction(k[4], [C[1]], nothing))             # c₁->∅
# Set up reaction system object. Collect symbolic state variables into a single vector.
@named system = ReactionSystem(reactions, t, collect(C), k, combinatoric_ratelaws=true)

pODE = Pair.(collect(k),zeros(Int64,4))
u₀MapODE = Pair.(collect(C), zeros(Float64,nMax))
# Create problem object
odeProblem = ODEProblem(system,u₀MapODE,(0.0,Inf),pODE)

@model function fitmodel(data, prob)
    # Prior
    σ ~ Uniform() #InverseGamma(2, 3) 
    k₂ ~ truncated(Normal(0.25,1.0),0,2)
    k₃ ~ truncated(Normal(0.5,1.0),0,2)
    k₄ ~ truncated(Normal(1.0,1.0),0,2)
    V ~ truncated(Normal(1000.0,500.0),500,1500)
    # Likelihood
    p = [1.0,k₂,k₃,k₄]
    prob = remake(prob, p=p)
    predicted = (solve(odeProblem,callback=TerminateSteadyState())).u[end].*V

    for i = 1:length(predicted)
        data[i] ~ Normal(predicted[i], σ)
    end
end

dummyData = last.(sort(collect(counter(ceil.(Int64,rand(Exponential(3),1000)))),by=first))
# dummyData = last.(sort(collect(counter(ceil.(Int64,1.0./rand(Uniform(0.01,1),1000)))),by=first))
if length(dummyData)<nMax
    append!(dummyData, zeros(Int64,20-length(dummyData)))
else
    dummyData = dummyData[1:nMax]
end

model = fitmodel(dummyData, odeProblem)
# chain = sample(model, NUTS(), MCMCThreads(), 1000, 4)
chain = sample(model, SMC(), 1000)
# @time chain = mapreduce(c -> sample(model, SMC(), 1000), chainscat, 1:4)
# @time chain = mapreduce(c -> sample(model, NUTS(.7), 1000), chainscat, 1:4)




pPredicted = Pair.(collect(k),[1.0,mean(chain[:,2,:]),mean(chain[:,3,:]),mean(chain[:,4,:])])
# Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
u₀MapPredicted = Pair.(collect(C), zeros(Float64,nMax))
# Create problem object
odePredicted = ODEProblem(system,u₀MapPredicted,(0.0,Inf),pPredicted)
solPredicted = solve(odePredicted,callback=TerminateSteadyState())
uPredicted = (solPredicted.u[end]).*mean(chain[:,5,:])


fig = Figure(); ax = Axis(fig[1,1])
CairoMakie.scatter!(ax, collect(1:nMax), dummyData, label="Dummy data")
CairoMakie.lines!(ax, collect(1:nMax), uPredicted, label="Model with fitted parameters")
display(fig)
save("tst.png",fig)


using StatsPlots
plot(chain)