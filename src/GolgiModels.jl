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

module GolgiModels

using PrecompileTools
using DifferentialEquations
using Makie
using WGLMakie; WGLMakie.activate!()
using DrWatson
using UnPack
using GeometryBasics
using FileIO
using Catalyst
using FromFile
using Format
using JSServe

include(srcdir("AllReactions.jl"))
include(srcdir("GuiFigureSetup.jl"))
include(srcdir("DwellTimes.jl"))
include(srcdir("AnimStep.jl"))
include(srcdir("RefreshIntegrators.jl"))
include(srcdir("HattedConstants.jl"))

function golgiApp(; displayFlag=true)

    JSServe.configure_server!(listen_port=9384, listen_url="0.0.0.0")

    nMax = 20           # Max compartment size /vesicles
    dt = 100.0          # Integration time interval /seconds
    V = 10              # μm³
    ksInit = [1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0]
    k̂ = [1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0]
    kStochFactors = [V, 1 / V, 1.0, 1.0, 1.0, 1 / V, 1.0, 1.0, 1.0, 1 / V, 1.0, 1.0]

    # Set up figure canvas
    fig, axCis, axMed, axTra, axDwell, parameterSliders, run, reset, linearityToggle, xLimTimeAv, linearityToggle = guiFigureSetup(ksInit)

    # Catalyst system setup
    # Symbolic system parameters: rate constants 
    @parameters k[1:12]
    @variables t
    # Symbolic system variables: cis, medial, and trans compartment size counts 
    @species C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax]
    # Use these parameters and variables to define a reaction system 

    # Initialise reaction system (within array so it's a mutable object for ease of later updating)
    system = [refreshSystem(nMax, C, M, T, k, t, true)]

    # Initialise ODE integrator (within array so it's a mutable object for ease of later updating)
    integODE = [refreshODEs(nMax, C, M, T, k, t, ksInit, system[1])]

    # Initialise stochastic integrator and accompanying objects (within arrays so they are mutable objects for ease of later updating)
    pStoch = Pair.(collect(k), kStochFactors.*ksInit)
    u₀MapStoch = Pair.([collect(C); collect(M); collect(T)], zeros(Int32, 3 * nMax))
    discreteProblem = [DiscreteProblem(system[1], u₀MapStoch, (0.0, Inf), pStoch)]
    jumpProblem = [JumpProblem(system[1], discreteProblem[1], Direct(), save_positions=(false, false))] # Converts system to a set of MassActionJumps
    integStoch = [init(jumpProblem[1], SSAStepper())]

    # Setup important observables
    # Observable value for whether or not system is currently running
    isrunning = Observable(false)
    # Obvservable vector of slider values 
    ksObservable = lift([s.value for s in parameterSliders.sliders]...) do values...
        [values...]
    end
    # Set up observable objects for cis, med, and trans results
    deterministicCisObservable = Observable(zeros(Float32, nMax))
    deterministicMedObservable = Observable(zeros(Float32, nMax))
    deterministicTraObservable = Observable(zeros(Float32, nMax))
    dwellTimeObservable = Observable(zeros(Float32, 7))
    # Set up observable objects for cis, med, and trans results
    stochasticCisObservable = Observable(zeros(Int32, nMax))
    stochasticMedObservable = Observable(zeros(Int32, nMax))
    stochasticTraObservable = Observable(zeros(Int32, nMax))
    # Observables for stochastic time averages
    stochTimeAvCisObservable = Observable(zeros(Float32, nMax))
    stochTimeAvMedObservable = Observable(zeros(Float32, nMax))
    stochTimeAvTraObservable = Observable(zeros(Float32, nMax))

    # 44444444444444444444444

    # Initialise plots
    lines!(axCis, deterministicCisObservable, collect(1:nMax), color=(:red, 1.0), linewidth=6)
    lines!(axMed, deterministicMedObservable, collect(1:nMax), color=(:green, 1.0), linewidth=6)
    lines!(axTra, deterministicTraObservable, collect(1:nMax), color=(:blue, 1.0), linewidth=6)
    barplot!(axCis, collect(1:nMax), stochasticCisObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:red)
    barplot!(axMed, collect(1:nMax), stochasticMedObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:green)
    barplot!(axTra, collect(1:nMax), stochasticTraObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:blue)
    lines!(axCis, stochTimeAvCisObservable, collect(1:nMax), color=(:black, 0.5), linestyle="--", linewidth=6)
    lines!(axMed, stochTimeAvMedObservable, collect(1:nMax), color=(:black, 0.5), linestyle="--", linewidth=6)
    lines!(axTra, stochTimeAvTraObservable, collect(1:nMax), color=(:black, 0.5), linestyle="--", linewidth=6)
    barplot!(axDwell, collect(1:7), dwellTimeObservable, direction=:y, bins=collect(0.5:1.0:nMax+0.5), color=:blue)
    dwellTimesValues = zeros(Float32, 7)


    # Set up button actions
    # Start/stop system when run button is clicked
    on(run.clicks) do clicks
        isrunning[] = !isrunning[]
    end
    # Reset system when reset button is clicked 
    on(reset.clicks) do clicks
        isrunning[] = false
        # Reset ODE integrator with new reaction system
        integODE[1] = refreshODEs(nMax, C, M, T, k, t, ksInit, system[1])
        # Reset stochastic integrator with new reaction system
        integStoch[1] = refreshStoch!(pStoch, u₀MapStoch, discreteProblem, jumpProblem, zeros(Int32, 3 * nMax), C, M, T, k, ksInit .* kStochFactors, system[1])
        # Return sliders to original values
        for i_slider in parameterSliders.sliders
            set_close_to!(i_slider, i_slider.startvalue[])
        end
        # Refresh all plots
        resetObservables(nMax, deterministicCisObservable, deterministicMedObservable, deterministicTraObservable, stochasticCisObservable, stochasticMedObservable, stochasticTraObservable, stochTimeAvCisObservable, stochTimeAvMedObservable, stochTimeAvTraObservable, dwellTimeObservable, xLimTimeAv)
    end
    # Switch linearity when linearityToggle is changed
    on(linearityToggle.active) do linearityValue
        # Reset reaction system dependent on value of linearity toggle
        system[1] = refreshSystem(nMax, C, M, T, k, t, linearityValue)
        # Reset ODE integrator with new reaction system
        integODE[1] = refreshODEs(nMax, C, M, T, k, t, ksInit, system[1])
        # Reset stochastic integrator with new reaction system
        integStoch[1] = refreshStoch!(pStoch, u₀MapStoch, discreteProblem, jumpProblem, zeros(Int32, 3 * nMax), C, M, T, k, ksInit .* kStochFactors, system[1])       
        resetObservables(nMax, deterministicCisObservable, deterministicMedObservable, deterministicTraObservable, stochasticCisObservable, stochasticMedObservable, stochasticTraObservable, stochTimeAvCisObservable, stochTimeAvMedObservable, stochTimeAvTraObservable, dwellTimeObservable, xLimTimeAv)
    end

    # Update integrators on changes to slider grid
    on(ksObservable) do ksVals
        integODE[1].p .= ksVals
        integStoch[1] = refreshStoch!(pStoch, u₀MapStoch, discreteProblem, jumpProblem, integStoch[1].u, C, M, T, k, kStochFactors .* ksVals, system[1])
    end

    on(xLimTimeAv) do x
        xlims!(axCis, (0.0, 1.1 * x))
        xlims!(axMed, (0.0, 1.1 * x))
        xlims!(axTra, (0.0, 1.1 * x))
    end

    

    # Start main loop when run button is clicked 
    on(run.clicks) do clicks
        @sync while isrunning[]
            isopen(fig.scene) || break # ensures computations stop if closed window

            animStepODE!(integODE[1], dt, axCis, axMed, axTra, deterministicCisObservable, deterministicMedObservable, deterministicTraObservable, nMax, V)
            animStepStoch!(integStoch[1], dt, axCis, axMed, axTra, stochasticCisObservable, stochasticMedObservable, stochasticTraObservable, stochTimeAvCisObservable, stochTimeAvMedObservable, stochTimeAvTraObservable, nMax)

            if linearityToggle.active[]
                hattedConstantsLinear!(integODE[1].p, k̂, integODE[1].u, nMax)
            else 
                hattedConstantsNonLinear!(integODE[1].p, k̂, integODE[1].u, nMax)
            end
            #   Improve this block; remove need for dwellTimesValues array
            for i = 1:7
                dwellTimesValues[i] = -dwellTimesFuns[i](k̂)
            end
            dwellTimeObservable[] .= dwellTimesValues
            dwellTimeObservable[] = dwellTimeObservable[]
            ylims!(axDwell, (0, maximum(dwellTimesValues)))
            ###########

            # Find time averaged maximum value to set xlim
            # xLimTimeAv[] = (xLimTimeAv[] * 19 + max(maximum(integStoch[1].u), maximum(integODE[1].u))) / 20
            xLimTimeAv[] = (xLimTimeAv[] * 19 + maximum(integODE[1].u)*V) / 20
            xLimTimeAv[] = xLimTimeAv[]            

            sleep(0.1)
        end
    end

    displayFlag == true ? display(fig) : nothing
end

@compile_workload begin
    golgiApp(displayFlag=false)
end

export golgiApp

end