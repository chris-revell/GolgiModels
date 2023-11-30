
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
using DataFrames
using StatsBase
using LinearAlgebra

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

pODE = Pair.(collect(k),zeros(Float64,4))
u₀MapODE = Pair.(collect(C), zeros(Float64,nMax))
# Create problem object
odeProblem = ODEProblem(system,u₀MapODE,(0.0,Inf),pODE)

##

dummyData = last.(sort(collect(counter(ceil.(Int64,rand(Exponential(3),1000)))),by=first))
# dummyData = last.(sort(collect(counter(ceil.(Int64,1.0./rand(Uniform(0.01,1),1000)))),by=first))
if length(dummyData)<nMax
    append!(dummyData, zeros(Int64,20-length(dummyData)))
else
    dummyData = dummyData[1:nMax]
end


##

@model function fitmodel(data, prob)
    # Prior
    σ ~ InverseGamma(2, 3)
    # k₂ ~ truncated(Normal(1.0,0.2),0.2,0.4)
    k₃ ~Uniform(1.0,2.0) #truncated(Normal(1.2,0.2),1.0,2.0)
    k₄ ~Uniform(0.0,1.2) #truncated(Normal(1.0,0.2),0.0,1.2)
    V ~ Uniform(200.0,400.0) #truncated(Normal(300,100),100.0,500.0)
    # Likelihood
    p = [1.0,1.0,k₃,k₄]
    prob = remake(prob, p=p)
    sol = solve(odeProblem,callback=TerminateSteadyState(min_t=100.0),save_idxs = [1])
    predicted = sol.u[end].*V

    if !(SciMLBase.successful_retcode(sol.retcode))
        Turing.@addlogprob! -Inf
        return
    elseif predicted[end]>predicted[1]
        Turing.@addlogprob! -Inf
        return
    end
    for i = 1:length(predicted)
        data[i] ~ Normal(predicted[i], σ^2)
    end
    return data
end

##

model = fitmodel(dummyData, odeProblem)

iterations = 10000
# iterator = NUTS(0.65)#HMCDA(0.15, 0.65) #HMCDA(0.15, 0.65) # SMC() # PG(10) # HMC(0.1, 5) # Gibbs(PG(10, :m) # HMC(0.1, 5, :s²)) # HMCDA(0.15, 0.65) # NUTS(0.65)
# chain = sample(model, iterator, iterations)
chain = sample(model, NUTS(), MCMCSerial(), iterations, 5)
chainDF = DataFrame(chain)


##

σHistogram = fit(Histogram, chainDF[!,:σ], nbins=iterations÷100)
σHistogramNormalized = normalize(h)#, mode=:density)
dx = σHistogram.edges[1].offset + σHistogram.edges[1].step/2
σMode = σHistogram.edges[1].offset + σHistogram.edges[1].step/2 + σHistogram.edges[1][findmax(σHistogram.weights)[2]]

k₃Histogram = fit(Histogram, chainDF[!,:k₃], nbins=iterations÷100)
k₃HistogramNormalized = normalize(h)#, mode=:density)
dx = k₃Histogram.edges[1].offset + k₃Histogram.edges[1].step/2
k₃Mode = k₃Histogram.edges[1].offset + k₃Histogram.edges[1].step/2 + k₃Histogram.edges[1][findmax(k₃Histogram.weights)[2]]

k₄Histogram = fit(Histogram, chainDF[!,:k₄], nbins=iterations÷100)
k₄HistogramNormalized = normalize(h)#, mode=:density)
dx = k₃Histogram.edges[1].offset + k₃Histogram.edges[1].step/2
k₄Mode = k₄Histogram.edges[1].offset + k₄Histogram.edges[1].step/2 + k₄Histogram.edges[1][findmax(k₄Histogram.weights)[2]]

VHistogram = fit(Histogram, chainDF[!,:V], nbins=iterations÷100)
VHistogramNormalized = normalize(h)#, mode=:density)
VMode = VHistogram.edges[1].offset + VHistogram.edges[1].step/2 + VHistogram.edges[1][findmax(VHistogram.weights)[2]]

Modes = (σ=0.0, k₃=k₃Mode, k₄=k₄Mode, V=VMode)
hists = (
    σ  σHistogram
    k₃HistogramNormalized
    k₄HistogramNormalized
    VHistogramNormalized

σ=0.0, k₃=k₃Mode, k₄=k₄Mode, V=VMode)

##

# @time chain = mapreduce(c -> sample(model, NUTS(.7), 1000), chainscat, 1:4)
# chainDF = DataFrame(chain)

pPredicted = Pair.(collect(k),[1.0,1.0,k₃Mode,k₄Mode])
# Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
u₀MapPredicted = Pair.(collect(C), zeros(Float64,nMax))
# Create problem object
odePredicted = ODEProblem(system,u₀MapPredicted,(0.0,Inf),pPredicted)
solPredicted = solve(odePredicted,callback=TerminateSteadyState())
uPredicted = (solPredicted.u[end]).*VMode

##

fig = Figure(resolution=(2000, 2000))
predictionFig = GridLayout(fig[1,1])
axPrediction = Axis(predictionFig[1,1],aspect=AxisAspect(2.5))
axPrediction.xlabel = "Compartment size"
axPrediction.ylabel = "Number of compartments"
CairoMakie.scatter!(axPrediction, collect(1:nMax), dummyData, label="Dummy data")
CairoMakie.lines!(axPrediction, collect(1:nMax), uPredicted, label="Model with fitted parameters")

params = names(chain, :parameters)
n_chains = length(chains(chain))
n_samples = length(chain)

chainFig = GridLayout(fig[2,1]) 



for (i, param) in enumerate(params)
    ax = Axis(chainFig[i, 1]; ylabel=string(param))
    for c in 1:n_chains
        values = chain[:, param, c]
        CairoMakie.scatter!(ax, 1:n_samples, values; label=string(i),markersize=4,markerspace=:pixel)#,color=(:black,0.01))
    end
    hlines!(ax,[Modes[param]])
    if i < length(params)
        hidexdecorations!(ax; grid=false)
    else
        ax.xlabel = "Iteration"
    end
end

for (i, param) in enumerate(params)
    ax = Axis(chainFig[i, 2]; ylabel=string(param))
    for c in 1:n_chains
        values = chain[:, param, c]
        # CairoMakie.density!(ax, values; label=string(i))
        barplot!(ax,hist;label=string(i))
    end
    vlines!(ax,[Modes[param]])
    if i == length(params)
        ax.xlabel = "Parameter estimate"
    end
end

rowsize!(fig.layout, 2, Relative(3/4))

display(fig)

##
# save("bayesianIterations.png",fig)

# N = 100_000
# x = rand(Gamma(6, 7), N)
# fig = Figure(); ax = Axis(fig[1,1])
# density!(ax, x)
# xlims!(ax,(0,100))
# display(fig)


# using StatsPlots
# StatsPlots.plot(chain)

