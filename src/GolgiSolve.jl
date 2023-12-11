#
#  GolgiSolve.jl
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

module GolgiSolve

# using PrecompileTools
using DifferentialEquations
using DrWatson
using UnPack
using Catalyst
using FromFile

include(srcdir("AllReactions.jl"))
include(srcdir("GuiFigureSetup.jl"))
include(srcdir("DwellTimes.jl"))
include(srcdir("AnimStep.jl"))
include(srcdir("RefreshObjects.jl"))
include(srcdir("HattedConstants.jl"))

function golgiSolve(;
    nMax=20,
    tMax=Inf,
    V=10,
    ks = [1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0],
    linearity = true
    )
    
    k̂ = [1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0]
    kStochFactors = [V, 1 / V, 1.0, 1.0, 1.0, 1 / V, 1.0, 1.0, 1.0, 1 / V, 1.0, 1.0]

    # Catalyst system setup
    # Symbolic system parameters: rate constants 
    @parameters k[1:12]
    @variables t
    # Symbolic system variables: cis, medial, and trans compartment size counts 
    @species C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax]
    # Use these parameters and variables to define a reaction system 
    # Initialise reaction system (within array so it's a mutable object for ease of later updating)
    system = allReactions(nMax,C,M,T,k,t)

    # Initialise ODE integrator (within array so it's a mutable object for ease of later updating)
    integODE = refreshODEs(nMax, C, M, T, k, ks, system)
    pODE = Pair.(collect(k),ks)
    # Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
    u₀MapODE = Pair.([collect(C); collect(M); collect(T)], zeros(Float64,3*nMax))
    # Create problem object
    odeProblem = ODEProblem(system,u₀MapODE,(0.0,tMax),pODE)
    solODE = solve(odeProblem)
    uODE = (solODE.u[end]).*V
 
    # Initialise stochastic integrator and accompanying objects (within arrays so they are mutable objects for ease of later updating)
    pStoch = Pair.(collect(k), kStochFactors.*ks)
    u₀MapStoch = Pair.([collect(C); collect(M); collect(T)], zeros(Int32, 3 * nMax))
    discreteProblem = DiscreteProblem(system, u₀MapStoch, (0.0, tMax), pStoch)
    jumpProblem = JumpProblem(system, discreteProblem, Direct(), save_positions=(false, false)) # Converts system to a set of MassActionJumps
    # integStoch = init(jumpProblem, SSAStepper())
    solStoch = solve(jumpProblem, SSAStepper())
    uStoch = solStoch.u[end]


    if linearity
        hattedConstantsLinear!(ks, k̂, solODE[end], nMax)
    else 
        hattedConstantsNonLinear!(ks, k̂, solStoch[end], nMax)
    end
    #   Improve this block; remove need for dwellTimesValues array
    dwellTimes = Float64[]
    for i = 1:7
        push!(dwellTimes,-dwellTimesFuns[i](k̂))
    end

    return uODE, uStoch, dwellTimes

end

# @compile_workload begin
#     golgiApp(displayFlag=false)
# end

export golgiSolve

end