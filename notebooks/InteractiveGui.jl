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

@from "$(projectdir("src","AllReactions.jl"))" using AllReactions

include(projectdir("scripts","testParameters.jl"))

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

# Function to update figure based on system iteration
function animstep!(integ,stochasticCisObservable,deterministicCisObservable,stochasticMedObservable,deterministicMedObservable,stochasticTranObservable,deterministicTranObservable,cisRange,medRange,tranRange)
    step!(integ)
    uDeterministic = integ.u
	# xlims!(axCis,(0.0,max(10,maximum(u[i][1:nMax]))))
	# xlims!(axMed,(0.0,max(10,maximum(u[i][1+nMax:2*nMax]))))
	# xlims!(axTra,(0.0,max(10,maximum(u[i][1+2*nMax:3*nMax]))))
	stochasticCisObservable[]     .= uStochastic[cisRange]
	deterministicCisObservable[]  .= uDeterministic[cisRange]
	stochasticMedObservable[]     .= uStochastic[medRange]
	deterministicMedObservable[]  .= uDeterministic[medRange]
	stochasticTranObservable[]    .= uStochastic[tranRange]
	deterministicTranObservable[] .= uDeterministic[tranRange]
end

# Initial conditions
u0   = zeros(nMax*3)

cisRange = 1:nMax
medRange = 1+nMax:2*nMax
tranRange = 1+2*nMax:3*nMax

# Set up figure canvas
fig = Figure()
axCis = Axis(fig[1,1],aspect=0.5)
axMed = Axis(fig[1,2],aspect=0.5,yticksvisible=false,yticklabelsvisible=false)
axTra = Axis(fig[1,3],aspect=0.5,yticksvisible=false,yticklabelsvisible=false)
Label(fig[1,1,Bottom()],"Cis",textsize=32)
Label(fig[1,2,Bottom()],"Medial",textsize=32)
Label(fig[1,3,Bottom()],"Trans",textsize=32)    
axCis.yticks = 0:10:nMax
axCis.ylabel = "Compartment size"
# #upperLim=10
# upperLim = maximum([max(maximum(deterministicSol.u[i]),maximum(stochasticSol.u[i])) for i=1:nOutput+1])
# xlims!(axCis,(0,upperLim))
# ylims!(axCis,(0,nMax))
# xlims!(axMed,(0,upperLim))
# ylims!(axMed,(0,nMax))
# xlims!(axTra,(0,upperLim))
# ylims!(axTra,(0,nMax))

# Set up observable objects for cis results
stochasticCisObservable             = Observable(u0[cisRange])
deterministicCisObservable          = Observable(u0[cisRange].*volume)
# Set up observable objects for med results
stochasticMedObservable             = Observable(u0[medRange])
deterministicMedObservable          = Observable(u0[medRange].*volume)
# Set up observable objects for tran results
stochasticTranObservable            = Observable(u0[tranRange])
deterministicTranObservable         = Observable(u0[tranRange].*volume)

# Set up bar and line plots for cis results
barplot!(axCis, collect(1:nMax), stochasticCisObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:red,0.5))
lines!(axCis, deterministicCisObservable, collect(1:nMax), color=(:blue,0.5), linewidth=2)
# Set up bar and line plots for med results
barplot!(axMed, collect(1:nMax), stochasticMedObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:green,0.5))
lines!(axMed, deterministicMedObservable, collect(1:nMax), color=(:blue,0.5), linewidth=2)
# Set up bar and line plots for tran results
barplot!(axTra, collect(1:nMax), stochasticTranObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:blue,0.5))
lines!(axTra, deterministicTranObservable, collect(1:nMax), color=(:blue,0.5), linewidth=2)


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
sliderobservables = [s.value for s in lsgrid.sliders]


# Symbolic system parameters: time and rate constants 
@parameters t K₀ K₁ K₂ K₃ K₄ K₅ K₆ K₇ K₈ K₉ K₁₀ K₁₁
# Symbolic system variables: Vector of number/concentration for cis, medial, and trans
@variables C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax]

reactions = allReactions(nMax,C,M,T,K₀,K₁,K₂,K₃,K₄,K₅,K₆,K₇,K₈,K₉,K₁₀,K₁₁)

# Set up reaction system object 
@named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], [K₀,K₁,K₂,K₃,K₄,K₅,K₆,K₇,K₇,K₈,K₉,K₁₀,K₁₁])

Ks = [:K₀, :K₁, :K₂, :K₃, :K₄, :K₅, :K₆, :K₇, :K₈, :K₉, :K₁₀, :K₁₁]
# Map symbolic rate constants to values for stochastic model 
p2 = [Ks[i] => sliderobservables[i][] for i=1:12]
# Map symbolic state vectors to float vector for stochastic model 
u₀ = zeros(Float64,3*nMax)
u₀Map = Pair.([collect(C); collect(M); collect(T)],u₀)    
odeProblem = ODEProblem(system,u₀Map,(0.0,10000.0),p2)

# Set up differential equation integrator
dp = ContinuousDynamicalSystem(odeProblem)

# Set up integrator for each iteration
integ = integrator(dp, u0; alg = Tsit5(), adaptive = false, dt = 0.1)


# Add stop/start button to the top of the canvas
run = Button(fig[1,1:3]; label = "Start/Stop", tellwidth = false)
isrunning = Observable(false)
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break # ensures computations stop if closed window
        animstep!(integ, dots, hands, sunPhaseLine)
        sleep(0.02) # or `yield()` instead
    end
end

display(fig)
