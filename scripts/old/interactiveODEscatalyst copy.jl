#
#  InteractiveGui.jl
#  GolgiModels
#
#  Created by Christopher Revell on dd/mm/yy.
#
#
# With components adapted from https://gist.github.com/Datseris/4b9d25a3ddb3936d3b83d3037f8188dd

# Interactive parameters:
# k₀ : ∅->c₁
# k₁ : c₁+cₙ->cₙ₊₁
# k₂ : cₙ->c₁+cₙ₋₁
# k₃ : c₁->m₁
# k₄ : m₁->c₁
# k₅ : m₁+mₙ->mₙ₊₁
# k₆ : mₙ->m₁+mₙ₋₁
# k₇ : m₁->t₁
# k₈ : t₁->m₁
# k₉ : t₁+tₙ->tₙ₊₁
# k₁₀: mₙ->m₁+mₙ₋₁
# k₁₁: t₁->∅

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

nMax      = 100    #nMax
volume    = 1.0    #volume

# Catalyst system setup

# Symbolic system parameters: time and rate constants 
@parameters t k₀ k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁
# Symbolic system variables: Vector of number/concentration for cis, medial, and trans
@variables C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax]

# Map symbolic rate constants to values for deterministic model 
kSymbols = [:k₀, :k₁, :k₂, :k₃, :k₄, :k₅, :k₆, :k₇, :k₈, :k₉, :k₁₀, :k₁₁]
p = kSymbols.=>ones(Float32,12)
# Map symbolic state vectors to float vector for stochastic model 
u₀ = zeros(Float32,3*nMax)
u₀Map = Pair.([collect(C); collect(M); collect(T)],u₀)    

# vector to store the Reactions
reactions = []
push!(reactions, Reaction(k₀, nothing, [C[1]]))                      #, nothing, [1]))    # Insertion into cis. Reactant ∅; product X[1]; rate: zeroth order kinetics
# push!(reactions, Reaction(k₁, [C[1]], [C[2]], [2], [1]))            # forming dimers
for i=1:nMax-1
    push!(reactions, Reaction(k₁, [C[i], C[1]], [C[i+1]]))           #, [1,1], [1])) # cis aggregation: second order kinetics
end
# push!(reactions, Reaction(k₃, [C[2]], [C[1]], [1], [2]))                # splitting dimers
for i=2:nMax
    push!(reactions, Reaction(k₂, [C[i]], [C[i-1],C[1]]))                #, [1], [1,1])) # cis splitting: first order kinetics        
end
push!(reactions, Reaction(k₃, [C[1]], [M[1]]))                           #, [1], [1])) # cis to medial: first order kinetics 
push!(reactions, Reaction(k₄, [M[1]], [C[1]]))                       #, [1], [1])) # medial to cis: first order kinetics 

# push!(reactions, Reaction(k₅, [M[1]], [M[2]], [2], [1]))        # forming dimers
for i=1:nMax-1
    push!(reactions, Reaction(k₅, [M[i],M[1]], [M[i+1]]))        #, [1,1], [1])) # med aggregation: second order kinetics
end
# push!(reactions, Reaction(k₆, [M[2]], [M[1]], [1], [2]))           # splitting dimers 
for i=2:nMax
    push!(reactions, Reaction(k₆, [M[i]], [M[i-1],M[1]]))           #, [1], [1,1])) # med splitting: first order kinetics                
end
push!(reactions, Reaction(k₇, [M[1]], [T[1]]))                          #, [1], [1])) # med to tran: first order kinetics        
push!(reactions, Reaction(k₈, [T[1]], [M[1]]))                           #, [1], [1])) # tran to med: first order kinetics            

# push!(reactions, Reaction(k₉, [T[1]], [T[2]], [2], [1]))        # forming dimers
for i=1:nMax-1
    push!(reactions, Reaction(k₉, [T[i],T[1]], [T[i+1]]))        #, [1,1], [1])) # tran aggregation: second order kinetics                
end
# push!(reactions, Reaction(k₁₀, [T[2]], [T[1]], [1], [2]))          # splitting dimers
for i=2:nMax
    push!(reactions, Reaction(k₁₀, [T[i]], [T[i-1],T[1]]))          #, [1], [1,1])) # tran splitting: first order kinetics                
end
push!(reactions, Reaction(k₁₁ , [T[1]], nothing))               #, [1], nothing)) # tran to ∅: first order kinetics            

# Set up reaction system object 
@named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], [k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₇,k₈,k₉,k₁₀,k₁₁])
# Create problem object
odeProblem = ODEProblem(system,u₀Map,(0.0,Inf),p)
# Create integrator object
integ = init(odeProblem,KenCarp3())

# Function to update figure based on system iteration
function animStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
    step!(integ, 10.0)
	deterministicCisObservable[] .= integ.u[1:nMax]
    deterministicCisObservable[] = deterministicCisObservable[]
	deterministicMedObservable[] .= integ.u[1+nMax:2*nMax]
    deterministicMedObservable[] = deterministicMedObservable[]
	deterministicTraObservable[] .= integ.u[1+2*nMax:3*nMax]
    deterministicTraObservable[] = deterministicTraObservable[]
end

function resetStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
    reinit!(integ,erase_sol=true)
    deterministicCisObservable[] .= integ.u[1:nMax]
    deterministicCisObservable[] = deterministicCisObservable[]
	deterministicMedObservable[] .= integ.u[1+nMax:2*nMax]
    deterministicMedObservable[] = deterministicMedObservable[]
	deterministicTraObservable[] .= integ.u[1+2*nMax:3*nMax]
    deterministicTraObservable[] = deterministicTraObservable[]
end


# Set up figure canvas
fig = Figure(resolution=(1700,1500),fontsize=32)
axDiagram = Axis(fig[3,1:4],title="Model diagram",aspect=DataAspect())
image!(axDiagram,rotr90(load(joinpath("supplementary","model.png"))))
hidedecorations!(axDiagram)
hidespines!(axDiagram)
axCis = Axis(fig[1,1], aspect=0.55, xminorgridvisible=true, yminorgridvisible=true)
axCis.yticks = 1:9:100
axCis.xticklabelsvisible = false
axMed = Axis(fig[1,2], aspect=0.55, xminorgridvisible=true, yminorgridvisible=true)
axMed.yticks = 1:9:100
axMed.yticklabelsvisible = false
axMed.xticklabelsvisible = false
axTra = Axis(fig[1,3], aspect=0.55, xminorgridvisible=true, yminorgridvisible=true)
axTra.yticks = 1:9:100
axTra.yticklabelsvisible = false
axTra.xticklabelsvisible = false
xlims!(axCis,(0.0,2.0))
xlims!(axMed,(0.0,2.0))
xlims!(axTra,(0.0,2.0))
Label(fig[1,1,Bottom()],"Cis concentration",fontsize=32)
Label(fig[1,2,Bottom()],"Medial concentration",fontsize=32)
Label(fig[1,3,Bottom()],"Trans concentration",fontsize=32)    
axCis.yticks = 0:10:nMax
axCis.ylabel = "Compartment size"


# Set up observable objects for cis results
deterministicCisObservable = Observable(u₀[1:nMax].*volume)
# Set up observable objects for med results
deterministicMedObservable = Observable(u₀[1+nMax:2*nMax].*volume)
# Set up observable objects for tran results
deterministicTraObservable = Observable(u₀[1+2*nMax:3*nMax].*volume)

yVals = collect(1:nMax)

lines!(axCis, deterministicCisObservable, yVals, color=(:red,1.0),   linewidth=6)
lines!(axMed, deterministicMedObservable, yVals, color=(:green,1.0), linewidth=6)
lines!(axTra, deterministicTraObservable, yVals, color=(:blue,1.0),  linewidth=6)

# Set up parameter sliders
parameterSliders = SliderGrid(
    fig[1,4],
    (label="k₀,  ∅ → c₁      " , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₁,  c₁+cₙ → cₙ₊₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₂,  cₙ → c₁+cₙ₋₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₃,  c₁ → m₁     " , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₄,  m₁ → c₁     " , range=0.0:0.1:2.0, startvalue=0.0, format="{:.2f}"),
    (label="k₅,  m₁+mₙ → mₙ₊₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₆,  mₙ → m₁+mₙ₋₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₇,  m₁ → t₁     " , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₈,  t₁ → m₁     " , range=0.0:0.1:2.0, startvalue=0.0, format="{:.2f}"),
    (label="k₉,  t₁+tₙ → tₙ₊₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₁₀, tₙ → t₁+tₙ₋₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₁₁, t₁ → ∅      " , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}");
)
# Pull parameters from slider positions
kObservables = [s.value for s in parameterSliders.sliders]

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

isrunning = Observable(false)
# isreset   = Observable(false)
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
        animStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)        
        sleep(0.1) # yield()
    end
end

display(fig)
