#
#  InteractiveGui.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/23.
#
#
# With components adapted from https://gist.github.com/Datseris/4b9d25a3ddb3936d3b83d3037f8188dd

# Interactive parameters:
# k₁ : ∅->c₁
# k₂ : c₁+cₙ->cₙ₊₁
# k₃ : cₙ->c₁+cₙ₋₁
# k₄ : c₁->m₁
# k₅ : m₁->c₁
# k₆ : m₁+mₙ->mₙ₊₁
# k₇ : mₙ->m₁+mₙ₋₁
# k₈ : m₁->t₁
# k₉ : t₁->m₁
# k₁₀: t₁+tₙ->tₙ₊₁
# k₁₁: mₙ->m₁+mₙ₋₁
# k₁₂: t₁->∅

# Reaction rate is in number of reactions per time => mol s⁻¹
# For ODEs:
# 0th order reaction: Rate expression is of the form kₒ => kₒ has units of mol s⁻¹    #rate constant => kLv where k has units of /s
# 1st order reaction is of the form kₒ[x] where [x] is reactant concentration in mol*m⁻³ => kₒ has units of m³s⁻¹
# 2nd order reaction is of the form kₒ[x][x] where [x] is reactant concentration in mol*m⁻³ => kₒ has units of m⁶s⁻¹mol⁻¹

# For stochastic discrete problem:
# kₛ in units of mol s⁻¹ for all reaction orders
# 0th order reaction: kₛ=kₒv
# 1st order reaction: kₛ=kₒ/L 
# 2nd order reaction: kₛ=kₒ/L²v²

# using DynamicalSystems
using DifferentialEquations
using GLMakie
using DrWatson
using LinearAlgebra
using UnPack
using GeometryBasics
using FileIO
using FastBroadcast
using OrdinaryDiffEq
using Catalyst
using FromFile

@from "$(projectdir("src","AllReactions.jl"))" using AllReactions
# @from "$(projectdir("src","AnimStep.jl"))" using AnimStep
# @from "$(projectdir("src","ResetStep.jl"))" using ResetStep
# @from "$(projectdir("src","GuiFigureSetup.jl"))" using GuiFigureSetup

# Function to update figure based on system iteration
function animStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax,xLimTimeAv)
    step!(integ, 10.0)
    # Find time averaged maximum value to set xlim
    xLimTimeAv[1] = (xLimTimeAv[1]*19+maximum(integ.u))/20
    xlims!(axCis,(0.0,1.1*xLimTimeAv[1]))
    xlims!(axMed,(0.0,1.1*xLimTimeAv[1]))
    xlims!(axTra,(0.0,1.1*xLimTimeAv[1]))
	deterministicCisObservable[] .= integ.u[1:nMax]
    deterministicCisObservable[] = deterministicCisObservable[]
	deterministicMedObservable[] .= integ.u[1+nMax:2*nMax]
    deterministicMedObservable[] = deterministicMedObservable[]
	deterministicTraObservable[] .= integ.u[1+2*nMax:3*nMax]
    deterministicTraObservable[] = deterministicTraObservable[]
end

# Function to reset figure
function resetStep!(integ,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,nMax)
    reinit!(integ,erase_sol=true)
    stochasticCisObservable[] .= integ.u[1:nMax]
    stochasticCisObservable[] = stochasticCisObservable[]
	stochasticMedObservable[] .= integ.u[1+nMax:2*nMax]
    stochasticMedObservable[] = stochasticMedObservable[]
	stochasticTraObservable[] .= integ.u[1+2*nMax:3*nMax]
    stochasticTraObservable[] = stochasticTraObservable[]
end

# Set up figure canvas
fig = Figure(resolution=(1700,1500),fontsize=32)
axDiagram = Axis(fig[3,1:4],title="Model diagram",aspect=DataAspect())
image!(axDiagram,rotr90(load(joinpath("_research","model.png"))))
hidedecorations!(axDiagram)
hidespines!(axDiagram)
axCis = Axis(fig[1,1], aspect=0.55, ylabel = "Compartment size")
xlims!(axCis,(0,3))
axMed = Axis(fig[1,2], aspect=0.55, yticksvisible=false)
xlims!(axMed,(0,3))
axTra = Axis(fig[1,3], aspect=0.55, yticksvisible=false)
xlims!(axTra,(0,3))
Label(fig[1,1,Bottom()],"Cis",fontsize=32)
Label(fig[1,2,Bottom()],"Medial",fontsize=32)
Label(fig[1,3,Bottom()],"Trans",fontsize=32)

# Set up parameter sliders
parameterSliders = SliderGrid(
    fig[1,4],
    (label="k₁,  ∅ → c₁      " , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₂,  c₁+cₙ → cₙ₊₁" , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₃,  cₙ → c₁+cₙ₋₁" , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₄,  c₁ → m₁     " , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₅,  m₁ → c₁     " , range=0.0:1.0:120.0, startvalue=000.0, format="{:.2f}"),
    (label="k₆,  m₁+mₙ → mₙ₊₁" , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₇,  mₙ → m₁+mₙ₋₁" , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₈,  m₁ → t₁     " , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₉,  t₁ → m₁     " , range=0.0:1.0:120.0, startvalue=000.0, format="{:.2f}"),
    (label="k₁₀, t₁+tₙ → tₙ₊₁" , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₁₁, tₙ → t₁+tₙ₋₁" , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}"),
    (label="k₁₂, t₁ → ∅      " , range=0.0:1.0:120.0, startvalue=100.0, format="{:.2f}");
)

# Add stop/start button
run = Button(fig[2,1]; label = "Start/Stop", tellwidth = false)
reset = Button(fig[2,2]; label = "Reset", tellwidth = false)

colsize!(fig.layout, 1, Relative(0.25))
colsize!(fig.layout, 2, Relative(0.25))
colsize!(fig.layout, 3, Relative(0.25))
colsize!(fig.layout, 4, Relative(0.25))
rowsize!(fig.layout, 1, Aspect(1, 2.0))
rowsize!(fig.layout, 2, Aspect(1, 0.1))
resize_to_layout!(fig)

nMax    = 20             # Max compartment size
tMax    = Inf

# Catalyst system setup
# Symbolic system parameters: rate constants 
@parameters k[1:12] 
# Symbolic system variables: cis, medial, and trans compartment size counts 
@variables t C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax] 
# Use these parameters and variables to define a reaction system 
# vector to store the Reactions
system = allReactions(nMax,C,M,T,k,t)

# Map symbolic paramters to values. Collect symbolic parameters into a vector.
p = Pair.(collect(k),100.0.*ones(Float32,12))
# Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
u₀Map = Pair.([collect(C); collect(M); collect(T)], zeros(Float32,3*nMax))

# Create problem object
odeProblem = ODEProblem(system,u₀Map,(0.0,tMax),p)
# odeProblem = convert(ODESystem(),system,u₀Map,(0.0,tMax),p)
# Create integrator object
integ = init(odeProblem,KenCarp3())

xLimTimeAv = [5.0]

# Set up observable objects for cis, med, and trans results
deterministicCisObservable = Observable(zeros(Float32, nMax))
deterministicMedObservable = Observable(zeros(Float32, nMax))
deterministicTraObservable = Observable(zeros(Float32, nMax))
# Initialise plots
lines!(axCis, deterministicCisObservable, collect(1:nMax), color=(:red,1.0),   linewidth=6)
lines!(axMed, deterministicMedObservable, collect(1:nMax), color=(:green,1.0), linewidth=6)
lines!(axTra, deterministicTraObservable, collect(1:nMax), color=(:blue,1.0),  linewidth=6)

# Pull parameters from slider positions
kObservables = [s.value for s in parameterSliders.sliders]

# Set up button actions 
isrunning = Observable(false)
on(run.clicks) do clicks
    isrunning[] = !isrunning[]
end
on(reset.clicks) do clicks    
    resetStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
    isrunning[] = false
end

on(run.clicks) do clicks
    @async while isrunning[]       
        isopen(fig.scene) || break # ensures computations stop if closed window
        for i=1:12
            integ.p[i] = kObservables[i][]
        end        
        animStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax,xLimTimeAv)        
        sleep(0.1)
    end
end

display(fig)
