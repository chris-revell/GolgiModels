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
# using LinearAlgebra
using UnPack
using GeometryBasics
using FileIO
# using FastBroadcast
# using OrdinaryDiffEq
using Catalyst
using FromFile
using Format

@from "$(projectdir("src","AllReactions.jl"))" using AllReactions
@from "$(projectdir("src","GuiFigureSetup.jl"))" using GuiFigureSetup
# @from "$(projectdir("src","AnimStep.jl"))" using AnimStep
# @from "$(projectdir("src","ResetStep.jl"))" using ResetStep

# Function to update figure based on system iteration
function animStep!(integ,axCis,axMed,axTra,cisObservable,medObservable,traObservable,nMax,xLimTimeAv)
    step!(integ, 10.0)
	cisObservable[] .= integ.u[1:nMax]
    cisObservable[] = cisObservable[]
	medObservable[] .= integ.u[1+nMax:2*nMax]
    medObservable[] = medObservable[]
	traObservable[] .= integ.u[1+2*nMax:3*nMax]
    traObservable[] = traObservable[]
    # Find time averaged maximum value to set xlim
    # if integ.t>100.0
        xLimTimeAv[1] = (xLimTimeAv[1]*19+maximum(integ.u))/20
        xlims!(axCis,(0.0,1.1*xLimTimeAv[1]))
        xlims!(axMed,(0.0,1.1*xLimTimeAv[1]))
        xlims!(axTra,(0.0,1.1*xLimTimeAv[1]))
    # else
    #     xLimTimeAv[1] = (xLimTimeAv[1]*19+maximum(integ.u))/20
    # end
end

# Function to reset figure
function resetStep!(integ,axCis,axMed,axTra,cisObservable,medObservable,traObservable,nMax)
    reinit!(integ,erase_sol=true)
    cisObservable[] .= integ.u[1:nMax]
    cisObservable[] = cisObservable[]
	medObservable[] .= integ.u[1+nMax:2*nMax]
    medObservable[] = medObservable[]
	traObservable[] .= integ.u[1+2*nMax:3*nMax]
    traObservable[] = traObservable[]
    xlims!(axCis,(0.0,5.0))
    xlims!(axMed,(0.0,5.0))
    xlims!(axTra,(0.0,5.0))
end

nMax    = 20             # Max compartment size
tMax    = Inf
ksInit = [1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]

# Catalyst system setup
# Symbolic system parameters: rate constants 
@parameters k[1:12]
@variables t
# Symbolic system variables: cis, medial, and trans compartment size counts 
@species C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax] 
# Use these parameters and variables to define a reaction system 
# vector to store the Reactions
system = allReactions(nMax,C,M,T,k,t)

# Map symbolic paramters to values. Collect symbolic parameters into a vector.
p = Pair.(collect(k),ksInit)
# Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
u₀Map = Pair.([collect(C); collect(M); collect(T)], zeros(Float32,3*nMax))

# Create problem object
odeProblem = ODEProblem(system,u₀Map,(0.0,tMax),p)
# Create integrator object
integ = init(odeProblem,KenCarp3())

fig, axCis, axMed, axTra, parameterSliders, run, reset = guiFigureSetup(ksInit)

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
    resetStep!(integ,axCis,axMed,axTra,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
    isrunning[] = false
end

on(run.clicks) do clicks
    @async while isrunning[]       
        isopen(fig.scene) || break # ensures computations stop if closed window
        for i=1:12
            integ.p[i] = kObservables[i][]
        end        
        animStep!(integ,axCis,axMed,axTra,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax,xLimTimeAv)
        # axCis.title="t=$(format(integ.t, precision=1))"
        sleep(0.1)
    end
end

display(fig)