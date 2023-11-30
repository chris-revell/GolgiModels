
using CairoMakie
using DifferentialEquations
using DrWatson
using UnPack
using FromFile
using DiffEqCallbacks
using Distributions
using Turing
using DataStructures
using DataFrames

function fAcceleration(du, u, p, t) 
    du .= 0.5*p[1]*t^2
end

pODE = [2.0]
u₀ODE = [0.0]
# Create problem object
odeProblem = ODEProblem(fAcceleration,u₀ODE,(0.0,10.0),pODE)
sol = solve(odeProblem, saveat=0.1)

dummyData = first.(sol.u).+ rand(Uniform(-20,20),101)

##

@model function fitmodel(data, prob)
    # Prior
    σ ~ InverseGamma(2, 3)
    a ~ Uniform(1.0,5.0)#,0.1,5.0)
    
    p = [a]
    prob = remake(prob, p=p)
    sol = solve(odeProblem)
    predicted = first.(sol.u)

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
chain = sample(model, NUTS(), MCMCSerial(), iterations, 3)
chainDF = DataFrame(chain)


##

fig = Figure(resolution=(2000, 2000), fontsize=36)
predictionFig = GridLayout(fig[1,1])
axPrediction = Axis(predictionFig[1,1],aspect=AxisAspect(2.5))

pODE = [mean(chainDF[!,:a])]
u₀ODE = [0.0]
# Create problem object
odeProblem = ODEProblem(fAcceleration,u₀ODE,(0.0,10.0),pODE)
sol = solve(odeProblem, saveat=0.1)
sca = CairoMakie.scatter!(axPrediction, sol.t, dummyData, label="Dummy data")
lin = CairoMakie.lines!(axPrediction, sol.t, first.(sol.u), label="Model with fitted parameters")
# CairoMakie.lines!(axPrediction, collect(1:nMax), mean(chainDF[!,:k]).*[exp(-x/mean(chainDF[!,:θ])) for x=1:20], label="Model with fitted parameters")
Legend(predictionFig[1, 2],
    [lin, sca],
    ["Prediction", "Data"])

params = names(chain, :parameters)
n_chains = length(chains(chain))
n_samples = length(chain)

chainFig = GridLayout(fig[2,1]) 

for (i, param) in enumerate(params)
    ax = Axis(chainFig[i, 1]; ylabel=string(param))
    for c in 1:n_chains
        values = chain[:, param, c]
        CairoMakie.scatter!(ax, 1:n_samples, values; label=string(i),markersize=10,markerspace=:pixel,color=(:black,0.2))
    end

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
        CairoMakie.density!(ax, values; label=string(i))
    end

    if i == length(params)
        ax.xlabel = "Parameter estimate"
    end
end

# axes = [only(contents(chainFig[i, 2])) for i in 1:length(params)]
# linkxaxes!(axes...)

rowsize!(fig.layout, 2, Relative(3/4))

display(fig)

##
# save("bayesianIterations.png",fig)


# using StatsPlots
# StatsPlots.plot(chain)

