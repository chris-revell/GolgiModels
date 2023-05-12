#
#  interactiveStochasticCatalyst.jl
#
#  Created by Christopher Revell on 28/04/23.
#
#
# With components adapted from https://gist.github.com/Datseris/4b9d25a3ddb3936d3b83d3037f8188dd

# Interactive parameters:
# k₁ : ∅->a₁
# k₂ : a₁+aₙ->aₙ₊₁
# k₃ : aₙ->a₁+aₙ₋₁
# k₄ : a₁->b₁
# k₅ : b₁->a₁
# k₆ : b₁+bₙ->bₙ₊₁
# k₇ : bₙ->b₁+bₙ₋₁
# k₈ : b₁->c₁
# k₉ : c₁->b₁
# k₁₀: c₁+cₙ->cₙ₊₁
# k₁₁: cₙ->c₁+cₙ₋₁
# k₁₂: c₁->∅

using DifferentialEquations
using GLMakie
using GeometryBasics
using Catalyst

# Find all possible reactions pairs that result in oligomers with size <= nMax
# Pass symbolic state vectors A, B, and C, and symbolic parameters k and t
function allReactions(nMax,A,B,C,k,t)
    # vector to store the Reactions
    reactions = []
    push!(reactions, Reaction(k[1], nothing, [A[1]]))            # ∅->c₁
    push!(reactions, Reaction(k[2], [A[1]], [A[2]], [2], [1]))   # 2c₁->c₂
    push!(reactions, Reaction(k[3], [A[2]], [A[1]], [1], [2]))   # c₂->2c₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[2], [A[i], A[1]], [A[i+1]])) # c₁+cₙ->cₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[3], [A[i]], [A[i-1],A[1]]))  # cₙ->c₁+cₙ₋₁ for 3<=n<=nMax
    end
    push!(reactions, Reaction(k[4], [A[1]], [B[1]]))             # c₁->m₁
    push!(reactions, Reaction(k[5], [B[1]], [A[1]]))             # m₁->c₁

    push!(reactions, Reaction(k[6], [B[1]], [B[2]], [2], [1]))   # 2m₁->m₂
    push!(reactions, Reaction(k[7], [B[2]], [B[1]], [1], [2]))   # m₂->2m₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[6], [B[i],B[1]], [B[i+1]]))  # m₁+mₙ->mₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[7], [B[i]], [B[i-1],B[1]]))  # mₙ->m₁+mₙ₋₁ for 3<=n<=2nMax
    end
    push!(reactions, Reaction(k[8], [B[1]], [C[1]]))             # m₁->t₁
    push!(reactions, Reaction(k[9], [C[1]], [B[1]]))             # t₁->m₁

    push!(reactions, Reaction(k[10], [C[1]], [C[2]], [2], [1]))  # 2t₁->t₂
    push!(reactions, Reaction(k[11], [C[2]], [C[1]], [1], [2]))  # t₂->2t₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[10], [C[i],C[1]], [C[i+1]])) # t₁+tₙ->tₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[11], [C[i]], [C[i-1],C[1]])) # tₙ->t₁+tₙ₋₁ for 3<=n<=2nMax
    end
    push!(reactions, Reaction(k[12], [C[1]], nothing))           # t₁->∅
    # Set up reaction system object. Collect symbolic state variables into a single vector.
    @named system = ReactionSystem(reactions, t, [collect(A); collect(B); collect(C)], k, combinatoric_ratelaws=false)
    return system
end

# Function to setup figure
function guiFigureSetup(ksInit)
    # Set up figure canvas
    fig = Figure(resolution=(1700,1500),fontsize=32)
    axA = Axis(fig[1,1], aspect=0.55, ylabel = "Aggregate size")
    xlims!(axA,(0,3))
    axB = Axis(fig[1,2], aspect=0.55, yticksvisible=false)
    xlims!(axB,(0,3))
    axC = Axis(fig[1,3], aspect=0.55, yticksvisible=false)
    xlims!(axC,(0,3))
    Label(fig[1,1,Bottom()],"A",fontsize=32)
    Label(fig[1,2,Bottom()],"B",fontsize=32)
    Label(fig[1,3,Bottom()],"A",fontsize=32)
    # Set up parameter sliders
    parameterSliders = SliderGrid(
        fig[1,4],
        (label="k₁,  ∅ → a₁      " , range=0.0:0.01:1.2, startvalue=ksInit[1], format="{:.2f}"),
        (label="k₂,  a₁+aₙ → aₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[2], format="{:.2f}"),
        (label="k₃,  aₙ → a₁+aₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[3], format="{:.2f}"),
        (label="k₄,  a₁ → b₁     " , range=0.0:0.01:1.2, startvalue=ksInit[4], format="{:.2f}"),
        (label="k₅,  b₁ → a₁     " , range=0.0:0.01:1.2, startvalue=ksInit[5], format="{:.2f}"),
        (label="k₆,  b₁+bₙ → bₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[6], format="{:.2f}"),
        (label="k₇,  bₙ → b₁+bₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[7], format="{:.2f}"),
        (label="k₈,  b₁ → c₁     " , range=0.0:0.01:1.2, startvalue=ksInit[8], format="{:.2f}"),
        (label="k₉,  c₁ → b₁     " , range=0.0:0.01:1.2, startvalue=ksInit[9], format="{:.2f}"),
        (label="k₁₀, c₁+cₙ → cₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[10], format="{:.2f}"),
        (label="k₁₁, cₙ → c₁+cₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[11], format="{:.2f}"),
        (label="k₁₂, c₁ → ∅      " , range=0.0:0.01:1.2, startvalue=ksInit[12], format="{:.2f}");
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

    return fig, axA, axB, axC, parameterSliders, run, reset
end

# Function to step system forwards in time and update figure 
function animStep!(integ,dt,axA,axB,axC,aObservable,bObservable,cObservable,nMax)
    # Step integrator forwards in time
    step!(integ, dt, true)
    # Update observables for each plot
	aObservable[] .= integ.u[1:nMax]
    aObservable[] = aObservable[]
	bObservable[] .= integ.u[1+nMax:2*nMax]
    bObservable[] = bObservable[]
	cObservable[] .= integ.u[1+2*nMax:3*nMax]
    cObservable[] = cObservable[]
end

# Function to reset ODE
function resetODE!(integ,axA,axB,axC,aObservable,bObservable,cObservable,nMax)
    reinit!(integ,erase_sol=true)
    aObservable[] .= integ.u[1:nMax]
    aObservable[] = aObservable[]
	bObservable[] .= integ.u[1+nMax:2*nMax]
    bObservable[] = bObservable[]
	cObservable[] .= integ.u[1+2*nMax:3*nMax]
    cObservable[] = cObservable[]
    xlims!(axA,(0.0,3.0))
    xlims!(axB,(0.0,3.0))
    xlims!(axC,(0.0,3.0))
end

# Function to reset stochastic integrator 
function resetStoch!(pStoch,u₀MapStoch,nMax,discreteProblem,jumpProblem,integStoch,aObservable,bObservable,cObservable)
    pStoch .= Pair.(collect(k),ksInit)
    # Map symbolic state vector zeros
    u₀MapStoch .= Pair.([collect(A); collect(B); collect(C)], zeros(Int32,3*nMax)) 
    # Reset problem object
    discreteProblem[1]  .= DiscreteProblem(system, u₀MapStoch, (0.0,Inf), pStoch)
    jumpProblem[1]   .= JumpProblem(system, discreteProblem[1], Direct(), save_positions=(false,false)) # Converts system to a set of MassActionJumps
    # Reset integrator object
    integStoch[1] = init(jumpProblem[1], SSAStepper())
    aObservable[] .= integStoch[1].u[1:nMax]
    aObservable[] = aObservable[]
	bObservable[] .= integStoch[1].u[1+nMax:2*nMax]
    bObservable[] = bObservable[]
	cObservable[] .= integStoch[1].u[1+2*nMax:3*nMax]
    cObservable[] = cObservable[]
end

nMax    = 20             # Max compartment size
dt      = 100.0          # Time step between GUI visualisation updates
ksInit = [1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]

# Catalyst system setup
@parameters k[1:12] # Rate constants 
@variables t        # Time 
@species A(t)[1:nMax] B(t)[1:nMax] C(t)[1:nMax] # Symbolic system variables: A, B, and C aggregate size counts 
# Create reaction system
system = allReactions(nMax,A,B,C,k,t)

# For ODEs:
# Map symbolic paramaters to values.
pODE = Pair.(collect(k),ksInit)
# Map symbolic state vectors to vector of initial values. 
u₀MapODE = Pair.([collect(A); collect(B); collect(C)], zeros(Float32,3*nMax))
# Create problem object
odeProblem = ODEProblem(system,u₀MapODE,(0.0,Inf),pODE)
# Create integrator object
integODE = init(odeProblem,KenCarp3())

# For stochastic problem:
# Map symbolic paramaters to values.
pStoch = Pair.(collect(k),ksInit)
# Map symbolic state vectors to vector of initial values. 
u₀MapStoch = Pair.([collect(A); collect(B); collect(C)], zeros(Int32,3*nMax)) 
# Create problem object
discreteProblem  = [DiscreteProblem(system, u₀MapStoch, (0.0,Inf), pStoch)]
jumpProblem   = [JumpProblem(system, discreteProblem[1], Direct(), save_positions=(false,false))]
# Create integrator object as the only component of a vector of length 1
integStoch = [init(jumpProblem[1], SSAStepper())]

# Create GUI figure
fig, axA, axB, axC, parameterSliders, run, reset = guiFigureSetup(ksInit)
xLimTimeAv = [5.0] # Stores a time average of max values for updating xlims

# Set up observable objects for a, b, and c ODE results
deterministicAObservable = Observable(zeros(Float32, nMax))
deterministicBObservable = Observable(zeros(Float32, nMax))
deterministicCObservable = Observable(zeros(Float32, nMax))
# Initialise line plots for ODEs
lines!(axA, deterministicAObservable, collect(1:nMax), color=(:red,1.0),   linewidth=6)
lines!(axB, deterministicBObservable, collect(1:nMax), color=(:green,1.0), linewidth=6)
lines!(axC, deterministicCObservable, collect(1:nMax), color=(:blue,1.0),  linewidth=6)
# Set up observable objects for a, b, and c stochastic results
stochasticAObservable = Observable(zeros(Int32, nMax))
stochasticBObservable = Observable(zeros(Int32, nMax))
stochasticCObservable = Observable(zeros(Int32, nMax))
# Initialise bar plots for stochastic results 
barplot!(axA, collect(1:nMax), stochasticAObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:red)
barplot!(axB, collect(1:nMax), stochasticBObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:green)
barplot!(axC, collect(1:nMax), stochasticCObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:blue)

# Pull parameters from slider positions
kObservables = [s.value for s in parameterSliders.sliders]

# Set up button actions 
isrunning = Observable(false)
on(run.clicks) do clicks
    # Start or stop when "run" button is clicked
    isrunning[] = !isrunning[]
end
on(reset.clicks) do clicks
    # Reset integrators when "reset" button is clicked 
    resetODE!(integODE,axA,axB,axC,deterministicAObservable,deterministicBObservable,deterministicCObservable,nMax)
    resetStoch!(pStoch,u₀MapStoch,nMax,discreteProblem,jumpProblem,integStoch,stochasticAObservable,stochasticBObservable,stochasticCObservable)
    isrunning[] = false
end

on(run.clicks) do clicks
    @async while isrunning[]       
        isopen(fig.scene) || break
        
        # Update ODE parameters according to sliders in GUI
        for i=1:12
            integODE.p[i] = kObservables[i][]
        end
        # Update ODE results and plots
        animStep!(integODE,dt,axA,axB,axC,deterministicAObservable,deterministicBObservable,deterministicCObservable,nMax)
        
        # Update stochastic parameters according to sliders in GUI
        # Call remake to update stochastic integrator.
        # NB discrete problem, jump problem, and stochastic integrator stored in vectors so they can be mutated, not reassigned
        for i=1:12
            pStoch[i] = Pair(k[i],kObservables[i][])
        end 
        u₀MapStoch .= Pair.([collect(A); collect(B); collect(C)], integStoch[1].u) 
        discreteProblem[1] = DiscreteProblem(system, u₀MapStoch, (integStoch[1].t,Inf), pStoch)
        jumpProblem[1] = remake(jumpProblem[1],prob=discreteProblem[1])
        integStoch[1] = init(jumpProblem[end], SSAStepper())
        # Update stochastic results and plots
        animStep!(integStoch[1],dt,axA,axB,axC,stochasticAObservable,stochasticBObservable,stochasticCObservable,nMax)        

        # Find time averaged maximum value to set xlim
        xLimTimeAv[1] = (xLimTimeAv[1]*19+max(maximum(integStoch[1].u),maximum(integODE.u)))/20
        xlims!(axA,(0.0,1.1*xLimTimeAv[1]))
        xlims!(axB,(0.0,1.1*xLimTimeAv[1]))
        xlims!(axC,(0.0,1.1*xLimTimeAv[1]))

        sleep(0.1)
    end
end

display(fig)