#
#  ClocksGui.jl
#  ClockOscillators
#
#  Created by Christopher Revell on 27/10/2021.
#
#
# With components adapted from https://gist.github.com/Datseris/4b9d25a3ddb3936d3b83d3037f8188dd

# ω0  Sunlight phase rate of change /day
# ωF  Fibroblast phase intrinsic rate of change /day
# ωM  Macrophage phase intrinsic rate of change /day
# ϕF  Fibroblast preferred phase offset in coupling to macrophage
# ϕM  Macrophage preferred phase offset in coupling to fibroblast
# μ   Amplitude of fibroblast phase coupling to macrophage (effect of macrophage on fibroblast)
# ν   Amplitude of macrophage phase coupling to fibroblast (effect of fibroblast on macrophage)
# αF  Amplitude of fibroblast phase coupling to sunlight (effect of sunlight on fibroblast)
# αM  Amplitude of macrophage phase coupling to sunlight (effect of sunlight on fibroblast)
# ψF  Fibroblast preferred phase offset from sunlight
# ψM  Macrophage preferred phase offset from sunlight


using DynamicalSystems
using DifferentialEquations
using GLMakie
using DataStructures: CircularBuffer
# using FileIO
using FromFile
using DrWatson
using ModelingToolkit
using Catalyst
using LinearAlgebra
using DifferentialEquations
using JumpProcesses

@from "$(projectdir())/src/AllReactions.jl" using AllReactions

nMax = 100

# Function to update figure based on system iteration
function animstep!(integ,nMax,uObservableCis,uObservableMed,uObservableTra)
    step!(integ)
    u = integ.u
	# xlims!(axCis,(0.0,max(10,maximum(u[i][1:nMax]))))
	# xlims!(axMed,(0.0,max(10,maximum(u[i][1+nMax:2*nMax]))))
	# xlims!(axTra,(0.0,max(10,maximum(u[i][1+2*nMax:3*nMax]))))
	uObservableCis[] .= u[i][1:nMax]
	uObservableMed[] .= u[i][1+nMax:2*nMax]
	uObservableTra[] .= u[i][1+2*nMax:3*nMax]
	uObservableCis[] = uObservableCis[]
	uObservableMed[] = uObservableMed[]
	uObservableTra[] = uObservableTra[]
end

# Initial conditions
u0   = zeros(nMax*3)

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
uObservableCis = Observable(u0[1:nMax])
uObservableMed = Observable(u0[1+nMax:2*nMax])
uObservableTra = Observable(u0[1+2*nMax:3*nMax])
barplot!(axCis, collect(1:nMax), uObservableCis, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:red)
barplot!(axMed, collect(1:nMax), uObservableMed, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:green)
barplot!(axTra, collect(1:nMax), uObservableTra, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:blue)

# Set up parameter sliders
lsgrid = labelslidergrid!(
    fig,
    ["volume","cisAgg","cisSplit","cisToMed","medToCis","medAgg","medSplit","medToTran","tranToMed","tranAgg","tranSplit","tranTo∅"],
    [0.0:0.01:1.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0, 0.0:0.1:2.0,0.0:0.1:2.0];
    formats = [x -> "$(round(x, digits = 1))"],
)
fig[1, 4] = lsgrid.layout
# Set default slider values
defaults = ones(12)
set_close_to!.(lsgrid.sliders, defaults)
# Pull parameters from slider positions
sliderobservables = [s.value for s in lsgrid.sliders]

nReactionsTotal = 1+2*(nMax-1)+2+2*(nMax-1)+2+2*(nMax-1)+1
@parameters t
@variables k[1:nReactionsTotal]  X[1:3*nMax](t)
reactants,products,rates,stoichiometryIn,stoichiometryOut = allReactions(nMax,X,sliderobservables...)
# vector to store the Reactions
reactions = []
for i=1:nReactionsTotal
	push!(reactions, Reaction(k[i], reactants[i], products[i], stoichiometryIn[i], stoichiometryOut[i]))
end
discreteprob   = DiscreteProblem(system, u₀map, tspan, pars)
jumpProblem   = JumpProblem(system, discreteprob, Direct(),save_positions=(false,false)) # Converts system to a set of MassActionJumps

# Set up differential equation integrator
prob = ODEProblem(model!,u0,(0.0,10000.0),sliderobservables)
dp = ContinuousDynamicalSystem(prob)
# Solve diffeq with constant step for smoother curves
diffeq = (alg = Tsit5(), adaptive = false, dt = δt)
# Set up integrator for each iteration
integ = integrator(dp, u0, SSAStepper(), saveat=tMax/100000)


# Add stop/start button to the top of the canvas
run = Button(ga[1,1:2]; label = "Start/Stop", tellwidth = false)
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
