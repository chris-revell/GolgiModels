#
#  InteractiveGui.jl
#  GolgiModels
#
#  Created by Christopher Revell on dd/mm/yy.
#
#
# With components adapted from https://gist.github.com/Datseris/4b9d25a3ddb3936d3b83d3037f8188dd


using DynamicalSystems
using DifferentialEquations
using GLMakie
# using DataStructures: CircularBuffer
# using FileIO
using FromFile
using DrWatson
using ModelingToolkit
using Catalyst
using LinearAlgebra
using DifferentialEquations
using JumpProcesses
using UnPack

# @from "$(projectdir("src","AllReactions.jl"))" using AllReactions

# include(projectdir("scripts","testParameters.jl"))

nMax      = 100    #nMax
tMax      = 1000.0 #tMax
volume    = 1.0    #volume
k₀        = 1.0    #∅ToCis
k₁        = 0.5    #cisAgg
k₂        = 1.0    #cisSplit
k₃        = 1.0    #cisToMed
k₄        = 1.0    #medToCis
k₅        = 0.5    #medAgg
k₆        = 1.0    #medSplit
k₇        = 1.0    #medToTran
k₈        = 1.0    #tranToMed
k₉        = 0.5    #tranAgg
k₁₀       = 1.0    #tranSplit
k₁₁       = 1.0    #tranTo∅
nOutput   = 100000

function deterministicModel!(du,u,p,t)

    nMax,ks = p

    # Cis monomers: ∅->cis₁ + med₁->cis₁ + cis₂->2cis₁ - cis₁->med₁ - 2cis₁->cis₂, cisₙ->cis₁+cisₙ₋₁ for n>=3, cis₁+cisₙ->cisₙ₊₁ for n>=2 (n=1 in first line)
    du[1] = ks[1][] - ks[4][]*u[1] + ks[5][]*u[1+nMax] - ks[2][]*u[1]^2 - ks[2][]*u[1]*sum(u[2:nMax-1]) + 2*ks[3][]*u[2] + ks[3][]*sum(u[3:nMax])
    # Cis dimers: cis₁+cis₁->cis₂ - cis₂->cis₁ - cis₂+cis₁->cis₃ + cis₃->cis₂+cis₁
    du[2] = (ks[2][]*u[1]^2)/2 − ks[3][]*u[2] − ks[2][]*u[2]*u[1,1] + ks[3][]*u[3]
    # Cis n-omers: (cis₁+cisₙ₋₁->cisₙ) + (cisₙ₊₁->cisₙ+cis₁) - (cis₁+cisₙ->cisₙ₊₁) - (cisₙ->cisₙ₋₁+cis₁) for n=2:nMax-1
    du[3:nMax-1] .= [ks[2][]*u[n-1]*u[1] + ks[3][]*u[n+1] - ks[2][]*u[n]*u[1] - ks[3][]*u[n] for n=3:nMax-1]
    # Cis nMax-omers: (cis₁+cisₙ₋₁->cisₙ) - (cisₙ->cis₁+cisₙ₋₁)
    du[nMax] = ks[2][]*u[nMax-1]*u[1] - ks[3][]*u[nMax]

    # Med monomers: cis₁->med₁ + tran₁->med₁ - med₁->cis₁ - med₁->tran₁ + med₂->2med₁ - med₁->tran₁ - 2med₁->med₂, medₙ->med₁+medₙ₋₁ for n>=3 (n=1,2 in first line), med₁+medₙ->medₙ₊₁ for n>=2 (n=1 in first line)
    du[1+nMax] = ks[4][]*u[1] - ks[5][]*u[1+nMax] + ks[9][]*u[1+2*nMax] - ks[8][]*u[1+nMax] - ks[6][]*u[1+nMax]^2 - ks[6][]*u[1+nMax]*sum(u[2+nMax:2*nMax-1]) + 2*ks[6][]*u[2+nMax] + ks[7][]*sum(u[3+nMax:2*nMax])
    # Med dimers: med₁+med₁->med₂ - med₂->med₁ - med₂+med₁->med₃ + med₃->med₂+med₁
    du[2+nMax] = (ks[2][]*u[1+nMax]^2)/2 − ks[3][]*u[2+nMax] − ks[2][]*u[2+nMax]*u[1+nMax] + ks[3][]*u[3+nMax]
    # Med n-omers: (med₁+medₙ₋₁->medₙ) + (medₙ₊₁->medₙ+med₁) - (med₁+medₙ->medₙ₊₁) - (medₙ->medₙ₋₁+med₁) for n=2:nMax-1
    du[3+nMax:2*nMax-1] .= [ks[6][]*u[n-1+nMax]*u[1+nMax] + ks[7][]*u[n+1+nMax] - ks[6][]*u[n+nMax]*u[1+nMax] - ks[7][]*u[n+nMax] for n=3:nMax-1]
    # Med nMax-omers: (med₁+medₙ₋₁->medₙ) - (medₙ->med₁+medₙ₋₁)
    du[2*nMax] = ks[6][]*u[2*nMax-1]*u[1+nMax] - ks[7][]*u[2*nMax]

    # Tran monomers: med₁->tran₁ - tran₁->∅ - tran₁->med₁ + tran₂->2tran₁ - 2tran₁->tran₂, tranₙ->tran₁+tranₙ₋₁ for n>=3 (n=1,2 in first line), tran₁+tranₙ->tranₙ₊₁ for n>=2 (n=1 in first line)
    du[1+2*nMax] = ks[8][]*u[1+nMax] - ks[9][]*u[1+2*nMax] - ks[12][]*u[1+2*nMax] - ks[10][]*u[1+2*nMax]^2 - ks[10][]*u[1+2*nMax]*sum(u[2+2*nMax:3*nMax-1]) + 2*ks[11][]*u[2+2*nMax] + ks[11][]*sum(u[3+2*nMax:3*nMax])
    # Tran dimers: tran₁+tran₁->tran₂ - tran₂->tran₁ - tran₂+tran₁->tran₃ + tran₃->tran₂+tran₁
    du[2+2*nMax] = (ks[2][]*u[1+2*nMax]^2)/2 − ks[3][]*u[2+2*nMax] − ks[2][]*u[2+2*nMax]*u[1+2*nMax] + ks[3][]*u[3+2*nMax]
    # Tran n-omers: (tran₁+tranₙ₋₁->tranₙ) + (tranₙ₊₁->tranₙ+tran₁) - (tran₁+tranₙ->tranₙ₊₁) - (tranₙ->tranₙ₋₁+tran₁) for n=2:nMax-1
    du[3+2*nMax:3*nMax-1] .= [ks[10][]*u[n-1+2*nMax]*u[1+2*nMax] + ks[11][]*u[n+1+2*nMax] - ks[10][]*u[n+2*nMax]*u[1+2*nMax] - ks[11][]*u[n+2*nMax] for n=3:nMax-1]
    # Tran nMax-omers: (tran₁+tranₙ₋₁->tranₙ) - (tranₙ->tran₁+tranₙ₋₁)
    du[3*nMax] = ks[10][]*u[3*nMax-1]*u[1+2*nMax] - ks[11][]*u[3*nMax]

    return du

end

# Function to update figure based on system iteration
function animStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTranObservable,nMax)
    step!(integ)
    uDeterministic = integ.u
    xlims!(axCis,(0.0,max(5.0,maximum(uDeterministic))))
    xlims!(axMed,(0.0,max(5.0,maximum(uDeterministic))))
    xlims!(axTra,(0.0,max(5.0,maximum(uDeterministic))))
	# xlims!(axCis,(0.0,max(10,maximum(u[i][1:nMax]))))
	# xlims!(axMed,(0.0,max(10,maximum(u[i][1+nMax:2*nMax]))))
	# xlims!(axTra,(0.0,max(10,maximum(u[i][1+2*nMax:3*nMax]))))
	# stochasticCisObservable[]     .= uStochastic[cisRange]
	deterministicCisObservable[]  .= uDeterministic[1:nMax]
    deterministicCisObservable[]  = deterministicCisObservable[]
	# stochasticMedObservable[]     .= uStochastic[medRange]
	deterministicMedObservable[]  .= uDeterministic[1+nMax:2*nMax]
    deterministicMedObservable[]  = deterministicMedObservable[]
	# stochasticTranObservable[]    .= uStochastic[tranRange]
	deterministicTranObservable[] .= uDeterministic[1+2*nMax:3*nMax]
    deterministicTranObservable[] = deterministicTranObservable[]
end

# Initial conditions
u0   = zeros(nMax*3)

# Set up figure canvas
fig = Figure()
axCis = Axis(fig[1,1],aspect=0.55,yticksvisible=true,yticklabelsvisible=true)#,xminorgridvisible= true,yminorgridvisible=true)
axMed = Axis(fig[1,2],aspect=0.55,yticksvisible=false,yticklabelsvisible=false)#,xminorgridvisible=true,yminorgridvisible=true,xmajorgridvisible=true,ymajorgridvisible=true)
axTra = Axis(fig[1,3],aspect=0.55,yticksvisible=false,yticklabelsvisible=false)#,xminorgridvisible=true,yminorgridvisible=true,xmajorgridvisible=true,ymajorgridvisible=true)
xlims!(axCis,(0.0,5.0))
xlims!(axMed,(0.0,5.0))
xlims!(axTra,(0.0,5.0))
Label(fig[1,1,Bottom()],"Cis",fontsize=32)
Label(fig[1,2,Bottom()],"Medial",fontsize=32)
Label(fig[1,3,Bottom()],"Trans",fontsize=32)    
axCis.yticks = 0:10:nMax
axCis.ylabel = "Compartment size"


# Set up observable objects for cis results
# stochasticCisObservable             = Observable(u0[cisRange])
deterministicCisObservable          = Observable(u0[1:nMax].*volume)
# Set up observable objects for med results
# stochasticMedObservable             = Observable(u0[medRange])
deterministicMedObservable          = Observable(u0[1+nMax:2*nMax].*volume)
# Set up observable objects for tran results
# stochasticTranObservable            = Observable(u0[tranRange])
deterministicTranObservable         = Observable(u0[1+2*nMax:3*nMax].*volume)

# Set up bar and line plots for cis results
# barplot!(axCis, collect(1:nMax), stochasticCisObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:red,0.5))
lines!(axCis, deterministicCisObservable, collect(1:nMax), color=(:red,1.0), linewidth=4)
# Set up bar and line plots for med results
# barplot!(axMed, collect(1:nMax), stochasticMedObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:green,0.5))
lines!(axMed, deterministicMedObservable, collect(1:nMax), color=(:green,1.0), linewidth=4)
# Set up bar and line plots for tran results
# barplot!(axTra, collect(1:nMax), stochasticTranObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:blue,0.5))
lines!(axTra, deterministicTranObservable, collect(1:nMax), color=(:blue,1.0), linewidth=4)

# Set up parameter sliders
lsgrid = labelslidergrid!(
    fig,
    ["k₀","k₁","k₂","k₃","k₄","k₅","k₆","k₇","k₈","k₉","k₁₀","k₁₁"],
    [0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0];
    formats = [x -> "$(round(x, digits = 1))"],
)
fig[1, 4] = lsgrid.layout
# Set default slider values
defaults = ones(12)
set_close_to!.(lsgrid.sliders, defaults)
# Pull parameters from slider positions
ks = [s.value for s in lsgrid.sliders]

p = Dict{String,Any}()
@pack! p = nMax,ks

# odeProblem = ODEProblem(deterministicModel!,u0,(0.0,100000.0),p)

# Set up differential equation integrator
# dp = ContinuousDynamicalSystem(odeProblem)

# Set up integrator for each iteration
prob = ODEProblem(deterministicModel!,u0,(0.0,Inf),[nMax,ks])
integ = init(prob,Tsit5())
# integ = init(odeProblem, u0; alg = Tsit5(), adaptive = false, dt = 0.1)

# Add stop/start button to the top of the canvas
run = Button(fig[2,1:4]; label = "Start/Stop", tellwidth = false)

colsize!(fig.layout, 1, Relative(0.25))
colsize!(fig.layout, 2, Relative(0.25))
colsize!(fig.layout, 3, Relative(0.25))
colsize!(fig.layout, 4, Relative(0.25))
rowsize!(fig.layout, 1, Aspect(1, 2.0))
rowsize!(fig.layout, 2, Aspect(1, 0.1))
resize_to_layout!(fig)

isrunning = Observable(false)
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break # ensures computations stop if closed window
        animStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTranObservable,nMax)        
        # sleep(0.000001) # or `yield()` instead
    end
end

display(fig)
