
using DifferentialEquations
using Distributions
using Turing
using DataFrames
using CairoMakie
using StatsBase
using LinearAlgebra
using Statistics
##
# Function describing rate of change of position of accelerating particle 
function fTest(du, u, p, t) 
    du .= p[1].*u
end
##
# Dummy data created by adding noise to a previous solution 
odeProblem = ODEProblem(fTest,[1.0],(0.0,1.5),[2.0])
sol = solve(odeProblem, saveat=0.01)
dummyData = first.(sol.u).+ rand(Normal(0.0,5.0),151)

##
# Model for inference 
@model function fitmodel(data, prob)
    # Prior
    σ ~ InverseGamma(1, 50)
    a ~ Normal(1.5,2.0)
    # Likelihood
    prob = remake(prob, p=[a])
    sol = solve(prob, saveat=0.01)
    predicted = first.(sol.u)

    for i = 1:length(predicted)
        data[i] ~ Normal(predicted[i], σ)
    end
    return data
end
##
# Sampling 
model = fitmodel(dummyData, odeProblem)
iterations = 1000
nChains = 4
chain = sample(model, NUTS(0.99), MCMCSerial(), iterations, nChains)
chainDF = DataFrame(chain)

##
# Solve with acceleration = mean of prior distribution 
odeProblem = ODEProblem(fTest,[1.0],(0.0,1.5),[mean(chainDF[!,:a])])
sol = solve(odeProblem, saveat=0.01)
##


hist = fit(Histogram, chainDF[!,:a], nbins=100)
histNormalised = normalize(hist)
dx = Float64(hist.edges[1][2] - hist.edges[1][1])
mode = Float64(hist.edges[1][findmax(hist.weights)[2]]+dx/2.0)
# mean = Statistics.mean(chainDF[!,:a])
maxWeight = maximum(hist.weights)
##

# Plot solution 
fig = Figure(resolution=(1000, 500), fontsize=18)
predictionFig = GridLayout(fig[1,1])
axPrediction = Axis(predictionFig[1,1],aspect=AxisAspect(2.5))
sca = CairoMakie.scatter!(axPrediction, sol.t, dummyData, label="Dummy data")
lin = CairoMakie.lines!(axPrediction, sol.t, first.(sol.u), label="Prediction")
# for i=1:length(hist.weights)
#     odeProblem = ODEProblem(fTest,[1.0],(0.0,1.5),[Float64(hist.edges[1][i]+dx/2.0)])
#     sol = solve(odeProblem, saveat=0.01)
#     CairoMakie.lines!(axPrediction, sol.t, first.(sol.u), color=(:black,0.5*hist.weights[i]/maxWeight))
# end

# ax2 = Axis(fig[1, 2])
# barplot!(ax2,hist,color=(Makie.wong_colors()[1],0.8))

Legend(predictionFig[1, 2],[lin, sca],["Prediction", "Data"])
# Plot samples and prior distributions
chainFig = GridLayout(fig[2,1]) 
for (i, param) in enumerate(names(chain, :parameters))
    ax1 = Axis(chainFig[i, 1]; ylabel=string(param))
    ax2 = Axis(chainFig[i, 2]; ylabel=string(param))
    for c in 1:length(chains(chain))
        values = chain[:, param, c]
        CairoMakie.scatter!(ax1, 1:length(chain), values; label=string(i),markersize=5,markerspace=:pixel,color=(:black,0.1))
        values = chain[:, param, c]
        CairoMakie.density!(ax2, values; label=string(i))        
    end
end
rowsize!(fig.layout, 2, Relative(3/4))
display(fig)



##
# save("bayesianIterations.png",fig)
