
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
using XLSX


##

measurements = DataFrame(XLSX.readtable(datadir("exp_raw","Nikki golgi data","AllNumbers.xlsx"),1))
measurements[!, [:Maturation,:CellType]] = convert.(String, measurements[!, [:Maturation,:CellType]])
headers = Symbol.(names(measurements))
floatHeaders = [h for h in headers if h∉[:Maturation,:CellType]]
measurements[!, floatHeaders] = convert.(Float64, measurements[!, 3:end])

##

# Delete row with outlier values
maxAreaInd = findmax(measurements[!,:Area])[2]
delete!(measurements, maxAreaInd)

##

medData = filter([:Maturation, :CellType] => (m, c) -> m == "Medial" && c == "WT", measurements, view=true)
traData = filter([:Maturation, :CellType] => (m, c) -> m == "Trans" && c == "WT", measurements, view=true)


(maxArea, maxAreaInd) = findmax(measurements[!,:Area])
step = maxArea*1.1/100
histbins = 0.0:step:maxArea*1.1
colors = (Cis=:red,Medial=:green,Trans=:blue)

medHist = fit(Histogram, medData[!,:Area], histbins)
traHist = fit(Histogram, traData[!,:Area], histbins)

##
nMax = length(medHist.weights)
@parameters k[1:8]
@variables t
@species C(t)[1:nMax] M(t)[1:nMax]
reactions = []
push!(reactions, Reaction(k[1], nothing, [C[1]]))            # ∅->c₁
push!(reactions, Reaction(k[2], [C[1]], [C[2]], [2], [1]))   # 2c₁->c₂
push!(reactions, Reaction(k[3], [C[2]], [C[1]], [1], [2]))   # c₂->2c₁
for i=2:nMax-1
    push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]])) # c₁+cₙ->cₙ₊₁ for 2<=n<nMax
end
for i=3:nMax
    push!(reactions, Reaction(k[3]*i^(2/3), [C[i]], [C[i-1],C[1]]))  # cₙ->c₁+cₙ₋₁ for 3<=n<=nMax
end
push!(reactions, Reaction(k[4], [C[1]], [M[1]]))             # c₁->m₁
push!(reactions, Reaction(k[5], [M[1]], [C[1]]))             # m₁->c₁
push!(reactions, Reaction(k[6]*2^(2/3), [M[1]], [M[2]], [2], [1]))   # 2m₁->m₂
push!(reactions, Reaction(k[7], [M[2]], [M[1]], [1], [2]))   # m₂->2m₁
for i=2:nMax-1
    push!(reactions, Reaction(k[6], [M[i],M[1]], [M[i+1]]))  # m₁+mₙ->mₙ₊₁ for 2<=n<nMax
end
for i=3:nMax
    push!(reactions, Reaction(k[7]*i^(2/3), [M[i]], [M[i-1],M[1]]))  # mₙ->m₁+mₙ₋₁ for 3<=n<=2nMax
end
push!(reactions, Reaction(k[8], [M[1]], nothing))             # m₁->∅
@named system = ReactionSystem(reactions, t, [collect(C); collect(M)], k, combinatoric_ratelaws=false)

pODE = Pair.(collect(k),ones(Float64,8))
u₀MapODE = Pair.([collect(C); collect(M)], zeros(Float64,2*nMax))
odeProblem = ODEProblem(system,u₀MapODE,(0.0,Inf),pODE)

##

@model function fitmodel(data, prob)
    # Prior
    σ ~ InverseGamma(1, 50)
    k₂Overk₃ ~ truncated(Normal(0.5,0.2),0.0,1.0)
    k₄ ~ truncated(Normal(1.0,0.2),0.0,1.2)
    k₅ ~ truncated(Normal(0.0,0.2),0.0,1.0)
    k₆Overk₇ ~ truncated(Normal(0.5,0.2),0.0,1.0)
    k₈ ~ truncated(Normal(1.0,0.2),0.0,1.2)
    V ~  Normal(150,50)
    # Likelihood
    prob = remake(prob, p=[1.0,k₂Overk₃,1.0, k₄, k₅, k₆Overk₇,1.0, k₈])
    sol = solve(prob,callback=TerminateSteadyState(min_t=100.0))
    predicted = sol.u[end].*V

    if !(SciMLBase.successful_retcode(sol.retcode))
        Turing.@addlogprob! -Inf
        return nothing
    elseif predicted[end]>predicted[1]
        Turing.@addlogprob! -Inf
        return nothing
    end
    for i = 1:length(predicted)
        data[i] ~ Normal(predicted[i], σ)
    end
    return data
end

##
model = fitmodel([medHist.weights...,traHist.weights...], odeProblem)
iterator = NUTS(0.65) 
iterations = 100
nChains = 1
chain = sample(model, iterator, MCMCSerial(), iterations, nChains)
chainDF = DataFrame(chain)

##

pPredicted = Pair.(collect(k), [1.0, mean(chainDF[!,:k₂Overk₃]), 1.0, mean(chainDF[!,:k₄]), mean(chainDF[!,:k₅]), mean(chainDF[!,:k₆Overk₇]),1.0, mean(chainDF[!,:k₈])])
# Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
u₀MapPredicted = Pair.([collect(C); collect(M)], zeros(Float64,2*nMax))
# Create problem object
odePredicted = ODEProblem(system,u₀MapPredicted,(0.0,Inf),pPredicted)
solPredicted = solve(odePredicted,callback=TerminateSteadyState())
uPredicted = (solPredicted.u[end]).*mean(chainDF[!,:V])

##
# Plot solution 
fig = Figure(size=(2000,2000))
axMed = Axis(fig[1,1],yscale=Makie.pseudolog10)
axTra = Axis(fig[1,2],yscale=Makie.pseudolog10)
barplot!(axMed, medHist,color=(colors[:Cis],0.8))#,direction=:x)
Label(fig[1,1,Top()],"Med",fontsize=32)
barplot!(axTra, traHist,color=(colors[:Medial],0.8))#,direction=:x)
Label(fig[1,2,Top()],"Tra",fontsize=32)
# display(fig)
linMed = CairoMakie.lines!(axMed, collect(1:nMax), uPredicted[1:nMax], label="Model with fitted parameters")
linTra = CairoMakie.lines!(axTra, collect(1:nMax), uPredicted[1+nMax:end], label="Model with fitted parameters")
# Legend(predictionFig[1, 2],[lin, sca],["Prediction", "Data"])
# Plot samples and prior distributions
chainFig = GridLayout(fig[2,1:2]) 
for (i, param) in enumerate(names(chain, :parameters))
    ax1 = Axis(chainFig[i, 1]; ylabel=string(param))
    ax2 = Axis(chainFig[i, 2]; ylabel=string(param))
    for c in 1:length(chains(chain))
        values = chain[:, param, c]
        CairoMakie.lines!(ax1, 1:length(chain), values; label=string(i),markersize=5,markerspace=:pixel,color=(Makie.wong_colors()[c],0.8))
        values = chain[:, param, c]
        CairoMakie.density!(ax2, values; label=string(i),color=(Makie.wong_colors()[c],0.8))        
    end
end
rowsize!(fig.layout, 2, Relative(3/4))
display(fig)
save("test2.png",fig)

##

hists = Dict()
for param names(chain, :parameters)
    hist = fit(Histogram, chainDF[!,param], nbins=iterations÷10)
    histNormalised = normalize(hist)
    dx = hist.edges[1].offset + hist.edges[1].step/2
    mode = Float64(hist.edges[1][findmax(hist.weights)[2]])
    mean = mean(chainDF[!,param])
    hists[param] = hist
    means[param] = mean
    modes[param] = mode
end

# ##

# fig = Figure(size=(2000, 2000))
# predictionFig = GridLayout(fig[1,1])
# axPrediction = Axis(predictionFig[1,1],aspect=AxisAspect(2.5))
# axPrediction.xlabel = "Compartment size"
# axPrediction.ylabel = "Number of compartments"

# params = names(chain, :parameters)
# n_chains = length(chains(chain))
# n_samples = length(chain)
# chainFig = GridLayout(fig[2,1]) 
# for (i, param) in enumerate(params)
#     ax = Axis(chainFig[i, 1]; ylabel=string(param))
#     for c in 1:n_chains
#         values = chain[:, param, c]
#         CairoMakie.scatter!(ax, 1:n_samples, values; label=string(i),markersize=10,markerspace=:pixel,color=(Makie.wong_colors()[i],0.2))
#     end
#     hlines!(ax,[modes[param]])
#     if i < length(params)
#         hidexdecorations!(ax; grid=false)
#     else
#         ax.xlabel = "Iteration"
#     end
# end

# for (i, param) in params
#     ax = Axis(chainFig[i, 2]; ylabel=string(param))
#     for c in 1:n_chains
#         values = chain[:, param, c]
#         # CairoMakie.density!(ax, values; label=string(i))
#         barplot!(ax,hists[Symbol(param)],label=string(i),color=(Makie.wong_colors()[i],0.8))
#     end
#     vlines!(ax,[modes[param]])
#     if i == length(params)
#         ax.xlabel = "Parameter estimate"
#     end
# end
# rowsize!(fig.layout, 2, Relative(3/4))
# display(fig)
