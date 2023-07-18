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
using DrWatson
using UnPack
using GeometryBasics
using FileIO
using Catalyst
using FromFile
using Format
using GLMakie
using JSServe
using WGLMakie; WGLMakie.activate!()
import JSServe.TailwindDashboard as D

@from "$(projectdir("src","AllReactionsNonLinear.jl"))" using AllReactionsNonLinear
@from "$(projectdir("src","GuiFigureSetup.jl"))" using GuiFigureSetup

# Function to update figure based on system iteration
function animStep!(integ,dt,axCis,axMed,axTra,cisObservable,medObservable,traObservable,nMax)
    step!(integ, dt, true)
	cisObservable[] .= integ.u[1:nMax]
    cisObservable[] = cisObservable[]
	medObservable[] .= integ.u[1+nMax:2*nMax]
    medObservable[] = medObservable[]
	traObservable[] .= integ.u[1+2*nMax:3*nMax]
    traObservable[] = traObservable[]
end

# Function to reset figure
function resetStepODE!(integ,axCis,axMed,axTra,cisObservable,medObservable,traObservable,nMax)
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

function resetStepStoch!(pStoch,u₀MapStoch,nMax,discreteProblem,tMax,jumpProblem,integStoch,cisObservable,medObservable,traObservable)
    pStoch .= Pair.(collect(k),ksInit)
    # Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
    u₀MapStoch .= Pair.([collect(C); collect(M); collect(T)], zeros(Int32,3*nMax)) 
    # Create problem object
    discreteProblem  .= [DiscreteProblem(system, u₀MapStoch, (0.0,tMax), pStoch)]
    jumpProblem   .= [JumpProblem(system, discreteProblem[1], Direct(), save_positions=(false,false))] # Converts system to a set of MassActionJumps
    # Create integrator object
    integStoch .= [init(jumpProblem[1], SSAStepper())]
    cisObservable[] .= integStoch[1].u[1:nMax]
    cisObservable[] = cisObservable[]
	medObservable[] .= integStoch[1].u[1+nMax:2*nMax]
    medObservable[] = medObservable[]
	traObservable[] .= integStoch[1].u[1+2*nMax:3*nMax]
    traObservable[] = traObservable[]
end

app = App() do session
    nMax    = 20             # Max compartment size
    dt      = 100.0
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
    system = allReactionsNonLinear(nMax,C,M,T,k,t)


    # Map symbolic paramters to values. Collect symbolic parameters into a vector.
    pODE = Pair.(collect(k),ksInit)
    # Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
    u₀MapODE = Pair.([collect(C); collect(M); collect(T)], zeros(Float32,3*nMax))
    # Create problem object
    odeProblem = ODEProblem(system,u₀MapODE,(0.0,tMax),pODE)
    # Create integrator object
    integODE = init(odeProblem,KenCarp3())


    # Map symbolic paramters to values. Collect symbolic parameters into a vector.
    pStoch = Pair.(collect(k),ksInit)
    # Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
    u₀MapStoch = Pair.([collect(C); collect(M); collect(T)], zeros(Int32,3*nMax)) 
    # Create problem object
    discreteProblem  = [DiscreteProblem(system, u₀MapStoch, (0.0,tMax), pStoch)]
    jumpProblem   = [JumpProblem(system, discreteProblem[1], Direct(), save_positions=(false,false))] # Converts system to a set of MassActionJumps
    # Create integrator object
    integStoch = [init(jumpProblem[1], SSAStepper())]#, saveat=tMax/nOutput)

    sl_x = D.Slider("name 1", 0:0.01:10)
    sl_y = D.Slider("name 2", 0:0.01:10)

    # Set up figure canvas
    fig = Figure(resolution=(1700,1500),fontsize=32)
    axDiagram = Axis(fig[3,1:4],title="Model diagram",aspect=DataAspect())
    image!(axDiagram,rotr90(load(projectdir("supplementary","GolgiCompartmentModel.png"))))
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
    # parameterSliders = SliderGrid(
    #     fig[1,4],
    #     (label="k₁,  ∅ → c₁      " , range=0.0:0.01:1.2, startvalue=ksInit[1], format="{:.2f}"),
    #     (label="k₂,  c₁+cₙ → cₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[2], format="{:.2f}"),
    #     (label="k₃,  cₙ → c₁+cₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[3], format="{:.2f}"),
    #     (label="k₄,  c₁ → m₁     " , range=0.0:0.01:1.2, startvalue=ksInit[4], format="{:.2f}"),
    #     (label="k₅,  m₁ → c₁     " , range=0.0:0.01:1.2, startvalue=ksInit[5], format="{:.2f}"),
    #     (label="k₆,  m₁+mₙ → mₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[6], format="{:.2f}"),
    #     (label="k₇,  mₙ → m₁+mₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[7], format="{:.2f}"),
    #     (label="k₈,  m₁ → t₁     " , range=0.0:0.01:1.2, startvalue=ksInit[8], format="{:.2f}"),
    #     (label="k₉,  t₁ → m₁     " , range=0.0:0.01:1.2, startvalue=ksInit[9], format="{:.2f}"),
    #     (label="k₁₀, t₁+tₙ → tₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[10], format="{:.2f}"),
    #     (label="k₁₁, tₙ → t₁+tₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[11], format="{:.2f}"),
    #     (label="k₁₂, t₁ → ∅      " , range=0.0:0.01:1.2, startvalue=ksInit[12], format="{:.2f}");
    # )

    sl_x = D.Slider("k1", 0:0.01:1.2, value=1.0)
    sl_y = D.Slider("k2", 0:0.01:1.2, value=1.0)

    lift(sl_x.value, sl_y.value) do x, y
        integODE.p[1] = x
        integODE.p[2] = y
    end

    # Add stop/start button
    run = GLMakie.Button(fig[2,1]; label = "Start/Stop", tellwidth = false)
    reset = GLMakie.Button(fig[2,2]; label = "Reset", tellwidth = false)

    colsize!(fig.layout, 1, Relative(0.25))
    colsize!(fig.layout, 2, Relative(0.25))
    colsize!(fig.layout, 3, Relative(0.25))
    colsize!(fig.layout, 4, Relative(0.25))
    rowsize!(fig.layout, 1, Aspect(1, 2.0))
    rowsize!(fig.layout, 2, Aspect(1, 0.1))
    resize_to_layout!(fig)

    xLimTimeAv = [5.0]

    # Set up observable objects for cis, med, and trans results
    deterministicCisObservable = Observable(zeros(Float32, nMax))
    deterministicMedObservable = Observable(zeros(Float32, nMax))
    deterministicTraObservable = Observable(zeros(Float32, nMax))
    # Initialise plots
    lines!(axCis, deterministicCisObservable, collect(1:nMax), color=(:red,1.0),   linewidth=6)
    lines!(axMed, deterministicMedObservable, collect(1:nMax), color=(:green,1.0), linewidth=6)
    lines!(axTra, deterministicTraObservable, collect(1:nMax), color=(:blue,1.0),  linewidth=6)
    # Set up observable objects for cis, med, and trans results
    stochasticCisObservable = Observable(zeros(Int32, nMax))
    stochasticMedObservable = Observable(zeros(Int32, nMax))
    stochasticTraObservable = Observable(zeros(Int32, nMax))
    # Initialise plots
    barplot!(axCis, collect(1:nMax), stochasticCisObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:red)
    barplot!(axMed, collect(1:nMax), stochasticMedObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:green)
    barplot!(axTra, collect(1:nMax), stochasticTraObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:blue)

    # Observables for stochastic time averages
    stochTimeAvCisObservable = Observable(zeros(Float32, nMax))
    stochTimeAvMedObservable = Observable(zeros(Float32, nMax))
    stochTimeAvTraObservable = Observable(zeros(Float32, nMax))
    # Initialise plots
    lines!(axCis, stochTimeAvCisObservable, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=6)
    lines!(axMed, stochTimeAvMedObservable, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=6)
    lines!(axTra, stochTimeAvTraObservable, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=6)

    # Pull parameters from slider positions
    kObservables = [s.value for s in parameterSliders.sliders]

    # Set up button actions 
    isrunning = Observable(false)
    on(run.clicks) do clicks
        isrunning[] = !isrunning[]
    end
    on(reset.clicks) do clicks    
        resetStepODE!(integODE,axCis,axMed,axTra,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
        resetStepStoch!(pStoch,u₀MapStoch,nMax,discreteProblem,tMax,jumpProblem,integStoch,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable)
        isrunning[] = false
    end

    on(run.clicks) do clicks
        @async while isrunning[]       
            isopen(fig.scene) || break # ensures computations stop if closed window
            
            for i=1:12
                integODE.p[i] = kObservables[i][]
            end        
            animStep!(integODE,dt,axCis,axMed,axTra,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
            axCis.title="t=$(format(integODE.t, precision=1))"
            
            for i=1:12
                pStoch[i] = Pair(k[i],kObservables[i][])
            end 
            u₀MapStoch .= Pair.([collect(C); collect(M); collect(T)], integStoch[1].u) 
            discreteProblem[1] = DiscreteProblem(system, u₀MapStoch, (integStoch[1].t,Inf), pStoch)
            jumpProblem[1] = remake(jumpProblem[1],prob=discreteProblem[1])
            integStoch[1] = init(jumpProblem[end], SSAStepper())
            animStep!(integStoch[1],dt,axCis,axMed,axTra,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,nMax)        

            # Find time averaged maximum value to set xlim
            xLimTimeAv[1] = (xLimTimeAv[1]*19+max(maximum(integStoch[1].u),maximum(integODE.u)))/20
            xlims!(axCis,(0.0,1.1*xLimTimeAv[1]))
            xlims!(axMed,(0.0,1.1*xLimTimeAv[1]))
            xlims!(axTra,(0.0,1.1*xLimTimeAv[1]))

            stochTimeAvCisObservable[] .= (stochTimeAvCisObservable[].*19.0.+integStoch[1].u[1:nMax])./20.0
            stochTimeAvCisObservable[] = stochTimeAvCisObservable[]
            stochTimeAvMedObservable[] .= (stochTimeAvMedObservable[].*19.0.+integStoch[1].u[1+nMax:2*nMax])./20.0
            stochTimeAvMedObservable[] = stochTimeAvMedObservable[]
            stochTimeAvTraObservable[] .= (stochTimeAvTraObservable[].*19.0.+integStoch[1].u[1+2*nMax:3*nMax])./20.0
            stochTimeAvTraObservable[] = stochTimeAvTraObservable[]

            axCis.title="t=$(format(integODE.t, precision=1))"
            # axMed.title="t=$(format(integStoch[1].t, precision=1))"

            sleep(0.1)
        end
    end
    display(fig)
end
    