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

@from "$(projectdir("src","AllReactions.jl"))" using AllReactions
# @from "$(projectdir("src","GuiFigureSetup.jl"))" using GuiFigureSetup
@from "$(projectdir("src","DwellTimes.jl"))" using DwellTimes

# JSServe.browser_display()


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

function hattedConstants!(k,k̂,u,nMax)
    k̂.=k
    k̂[2] = k̂[2]*sum(u[1:nMax-1])
    k̂[6] = k̂[6]*sum(u[1+nMax:2*nMax-1])
    k̂[10] = k̂[10]*sum(u[1+2*nMax:3*nMax-1])
end

app = App() do 
    nMax    = 20             # Max compartment size
    dt      = 100.0
    tMax    = Inf
    ksInit = [1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]

    k̂ = [1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]

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


    # Set up figure canvas
    fig = Figure(resolution=(1500,1500),fontsize=32)

    grd1 = GridLayout(fig[1,1])
    grd2 = GridLayout(fig[2,1])
    grd3 = GridLayout(fig[3,1])

    axDiagram = Axis(grd1[1,1:3],title="Model diagram",aspect=DataAspect())
    image!(axDiagram,rotr90(load(projectdir("supplementary","GolgiCompartmentModel.png"))))
    hidedecorations!(axDiagram)
    hidespines!(axDiagram)

    axCis = Axis(grd2[1,1], ylabel = "Compartment size")#, aspect=0.3)
    xlims!(axCis,(0,3))
    axMed = Axis(grd2[1,2], yticksvisible=false)#, aspect=0.3)
    xlims!(axMed,(0,3))
    axTra = Axis(grd2[1,3], yticksvisible=false)#, aspect=0.3)
    xlims!(axTra,(0,3))

    axDwell = Axis(grd3[1:2,4],title = "Dwell Times")
    xlims!(axDwell,(1.5,7.5))
    ylims!(axDwell,(0,1))
    axDwell.xticks = (1:7, ["∅", "C₁", "C₊", "M₁", "M₊", "T₁", "T₊"])
    axDwell.ylabel = "Relative dwell time"

    axReducedDiagram = Axis(grd3[1,1:3],title="Reduced model",aspect=DataAspect())
    hidedecorations!(axReducedDiagram); hidespines!(axReducedDiagram)
    image!(axReducedDiagram,rotr90(load(projectdir("supplementary","GolgiCompartmentModel_reduced.png"))))

    Label(grd2[1,1,Bottom()],"Cis",fontsize=32)
    Label(grd2[1,2,Bottom()],"Medial",fontsize=32)
    Label(grd2[1,3,Bottom()],"Trans",fontsize=32)

    sliderstyle = JSServe.Asset(joinpath(@__DIR__,"scripts" "sliderstyle.css"))

    sliderk1 = D.Slider( "k1, " , 0.0:0.01:1.2, class="slider", value=ksInit[1])       #"k₁,  ∅ → c₁      " , 0.0:0.01:1.2, value=ksInit[1])
    sliderk2 = D.Slider( "k2, " , 0.0:0.01:1.2, class="slider", value=ksInit[2])       #"k₂,  c₁+cₙ → cₙ₊₁" , 0.0:0.01:1.2, value=ksInit[2])
    sliderk3 = D.Slider( "k3, " , 0.0:0.01:1.2, class="slider", value=ksInit[3])       #"k₃,  cₙ → c₁+cₙ₋₁" , 0.0:0.01:1.2, value=ksInit[3])
    sliderk4 = D.Slider( "k4, " , 0.0:0.01:1.2, class="slider", value=ksInit[4])       #"k₄,  c₁ → m₁     " , 0.0:0.01:1.2, value=ksInit[4])
    sliderk5 = D.Slider( "k5, " , 0.0:0.01:1.2, class="slider", value=ksInit[5])       #"k₅,  m₁ → c₁     " , 0.0:0.01:1.2, value=ksInit[5])
    sliderk6 = D.Slider( "k6, " , 0.0:0.01:1.2, class="slider", value=ksInit[6])       #"k₆,  m₁+mₙ → mₙ₊₁" , 0.0:0.01:1.2, value=ksInit[6])
    sliderk7 = D.Slider( "k7, " , 0.0:0.01:1.2, class="slider", value=ksInit[7])       #"k₇,  mₙ → m₁+mₙ₋₁" , 0.0:0.01:1.2, value=ksInit[7])
    sliderk8 = D.Slider( "k8, " , 0.0:0.01:1.2, class="slider", value=ksInit[8])       #"k₈,  m₁ → t₁     " , 0.0:0.01:1.2, value=ksInit[8])
    sliderk9 = D.Slider( "k9, " , 0.0:0.01:1.2, class="slider", value=ksInit[9])       #"k₉,  t₁ → m₁     " , 0.0:0.01:1.2, value=ksInit[9])
    sliderk10 = D.Slider("k10," , 0.0:0.01:1.2, class="slider", value=ksInit[10])       #"k₁₀, t₁+tₙ → tₙ₊₁" , 0.0:0.01:1.2, value=ksInit[10])
    sliderk11 = D.Slider("k11," , 0.0:0.01:1.2, class="slider", value=ksInit[11])       #"k₁₁, tₙ → t₁+tₙ₋₁" , 0.0:0.01:1.2, value=ksInit[11])
    sliderk12 = D.Slider("k12," , 0.0:0.01:1.2, class="slider", value=ksInit[12])       #"k₁₂, t₁ → ∅      " , 0.0:0.01:1.2, value=ksInit[12])

    lift(sliderk1.value,
    sliderk1.value,
    sliderk3.value,
    sliderk4.value, 
    sliderk5.value, 
    sliderk6.value, 
    sliderk7.value, 
    sliderk8.value, 
    sliderk9.value, 
    sliderk10.value,
    sliderk11.value,
    sliderk12.value) do a, b, c, d, e, f, g, h, i, j, k, l
        integODE.p[1] = a
        integODE.p[2] = b
        integODE.p[3] = c
        integODE.p[4] = d
        integODE.p[5] = e
        integODE.p[6] = f
        integODE.p[7] = g
        integODE.p[8] = h
        integODE.p[9] = i
        integODE.p[10] = j
        integODE.p[11] = k
        integODE.p[12] = l
    end

    # Add stop/start button
    run = GLMakie.Button(grd3[2,1]; label = "Start/Stop", tellwidth = false)
    reset = GLMakie.Button(grd3[2,2]; label = "Reset", tellwidth = false)

    xLimTimeAv = [5.0]

    # Set up observable objects for cis, med, and trans results
    deterministicCisObservable = Observable(zeros(Float32, nMax))
    deterministicMedObservable = Observable(zeros(Float32, nMax))
    deterministicTraObservable = Observable(zeros(Float32, nMax))
    dwellTimeObservable = Observable(zeros(Float32,7))
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

    barplot!(axDwell, collect(1:7), dwellTimeObservable, direction=:y, bins=collect(0.5:1.0:nMax+0.5), color=:blue)

    # Observables for stochastic time averages
    stochTimeAvCisObservable = Observable(zeros(Float32, nMax))
    stochTimeAvMedObservable = Observable(zeros(Float32, nMax))
    stochTimeAvTraObservable = Observable(zeros(Float32, nMax))
    # Initialise plots
    lines!(axCis, stochTimeAvCisObservable, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=6)
    lines!(axMed, stochTimeAvMedObservable, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=6)
    lines!(axTra, stochTimeAvTraObservable, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=6)

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

    dwellTimesValues = zeros(7)

    on(run.clicks) do clicks
        @async while isrunning[]
        # while isrunning[]
            isopen(fig.scene) || break # ensures computations stop if closed window

            animStep!(integODE,dt,axCis,axMed,axTra,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
            
            hattedConstants!(integODE.p,k̂,integODE.u,nMax)
            
            for i=1:7
                dwellTimesValues[i] = -dwellTimesFuns[i](k̂)
            end
            dwellTimeObservable[] .= dwellTimesValues
            dwellTimeObservable[] = dwellTimeObservable[]
            ylims!(axDwell,(0,maximum(dwellTimesValues)))

            axCis.title="t=$(format(integODE.t, precision=1))"
            
            for i=1:12
                pStoch[i] = Pair(k[i],integODE.p[i])
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

    # display(fig)
    return DOM.div(JSServe.TailwindCSS, sliderstyle,sliderk1,sliderk2,sliderk3,sliderk4,sliderk5,sliderk6,sliderk7,sliderk8,sliderk9,sliderk10,sliderk11,sliderk12, fig)
end
display(app)