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

nMax      = 100    #nMax
volume    = 1.0    #volume

function deterministicModel!(du,u,p,t)

    @unpack nMax,ks = p

    # Cis monomers: ∅->cis₁ + med₁->cis₁ + cis₂->2cis₁ - cis₁->med₁ - 2cis₁->cis₂, cisₙ->cis₁+cisₙ₋₁ for n>=3, cis₁+cisₙ->cisₙ₊₁ for n>=2 (n=1 in first line)
    du[1] = ks[1] - ks[4]*u[1] + ks[5]*u[1+nMax] - 2*ks[2]*u[1]^2 - ks[2]*u[1]*sum(@view u[2:nMax-1]) + 2*ks[3]*u[2]^2 + ks[3]*sum(@view u[3:nMax])
    # Cis dimers: cis₁+cis₁->cis₂ - cis₂->cis₁ - cis₂+cis₁->cis₃ + cis₃->cis₂+cis₁
    du[2] = (ks[2]*u[1]^2)/2 − ks[3]*u[2]^2 − ks[2]*u[2]*u[1,1] + ks[3]*u[3]^2
    
    # Cis n-omers: (cis₁+cisₙ₋₁->cisₙ) + (cisₙ₊₁->cisₙ+cis₁) - (cis₁+cisₙ->cisₙ₊₁) - (cisₙ->cisₙ₋₁+cis₁) for n=2:nMax-1
    for n=3:nMax-1
        du[n] = ks[2]*u[n-1]*u[1] + ks[3]*u[n+1]^2 - ks[2]*u[n]*u[1] - ks[3]*u[n]^2
    end
    # du[3:nMax-1] .= [ks[2]*u[n-1]*u[1] + ks[3]*u[n+1]^2 - ks[2]*u[n]*u[1] - ks[3]*u[n]^2 for n=3:nMax-1]
    
    # Cis nMax-omers: (cis₁+cisₙ₋₁->cisₙ) - (cisₙ->cis₁+cisₙ₋₁)
    du[nMax] = ks[2]*u[nMax-1]*u[1] - ks[3]*u[nMax]^2

    # Med monomers: cis₁->med₁ + tran₁->med₁ - med₁->cis₁ - med₁->tran₁ + med₂->2med₁ - med₁->tran₁ - 2med₁->med₂, medₙ->med₁+medₙ₋₁ for n>=3 (n=1,2 in first line), med₁+medₙ->medₙ₊₁ for n>=2 (n=1 in first line)
    du[1+nMax] = ks[4]*u[1] - ks[5]*u[1+nMax] + ks[9]*u[1+2*nMax] - ks[8]*u[1+nMax] - 2*ks[6]*u[1+nMax]^2 - ks[6]*u[1+nMax]*sum(@view u[2+nMax:2*nMax-1]) + 2*ks[6]*u[2+nMax] + ks[7]*sum(@view u[3+nMax:2*nMax])
    # Med dimers: med₁+med₁->med₂ - med₂->med₁ - med₂+med₁->med₃ + med₃->med₂+med₁
    du[2+nMax] = (ks[2]*u[1+nMax]^2)/2 − ks[3]*u[2+nMax] − ks[2]*u[2+nMax]*u[1+nMax] + ks[3]*u[3+nMax]

    # Med n-omers: (med₁+medₙ₋₁->medₙ) + (medₙ₊₁->medₙ+med₁) - (med₁+medₙ->medₙ₊₁) - (medₙ->medₙ₋₁+med₁) for n=2:nMax-1
    for n=3+nMax:2*nMax-1
        du[n] = ks[6]*u[n-1]*u[1+nMax] + ks[7]*u[n+1] - ks[6]*u[n]*u[1+nMax] - ks[7]*u[n]
    end
    # du[3+nMax:2*nMax-1] .= [ks[6]*u[n-1+nMax]*u[1+nMax] + ks[7]*u[n+1+nMax] - ks[6]*u[n+nMax]*u[1+nMax] - ks[7]*u[n+nMax] for n=3:nMax-1]

    # Med nMax-omers: (med₁+medₙ₋₁->medₙ) - (medₙ->med₁+medₙ₋₁)
    du[2*nMax] = ks[6]*u[2*nMax-1]*u[1+nMax] - ks[7]*u[2*nMax]

    # Tran monomers: med₁->tran₁ - tran₁->∅ - tran₁->med₁ + tran₂->2tran₁ - 2tran₁->tran₂, tranₙ->tran₁+tranₙ₋₁ for n>=3 (n=1,2 in first line), tran₁+tranₙ->tranₙ₊₁ for n>=2 (n=1 in first line)
    du[1+2*nMax] = ks[8]*u[1+nMax] - ks[9]*u[1+2*nMax] - ks[12]*u[1+2*nMax] - 2*ks[10]*u[1+2*nMax]^2 - ks[10]*u[1+2*nMax]*sum(@view u[2+2*nMax:3*nMax-1]) + 2*ks[11]*u[2+2*nMax] + ks[11]*sum(@view u[3+2*nMax:3*nMax])
    # Tran dimers: tran₁+tran₁->tran₂ - tran₂->tran₁ - tran₂+tran₁->tran₃ + tran₃->tran₂+tran₁
    du[2+2*nMax] = (ks[2]*u[1+2*nMax]^2)/2 − ks[3]*u[2+2*nMax] − ks[2]*u[2+2*nMax]*u[1+2*nMax] + ks[3]*u[3+2*nMax]
    
    # Tran n-omers: (tran₁+tranₙ₋₁->tranₙ) + (tranₙ₊₁->tranₙ+tran₁) - (tran₁+tranₙ->tranₙ₊₁) - (tranₙ->tranₙ₋₁+tran₁) for n=2:nMax-1
    for n=3+2*nMax:3*nMax-1
        du[n] = ks[10]*u[n-1]*u[1+2*nMax] + ks[11]*u[n+1] - ks[10]*u[n]*u[1+2*nMax] - ks[11]*u[n]
    end
    # du[3+2*nMax:3*nMax-1] .= [ks[10]*u[n-1+2*nMax]*u[1+2*nMax] + ks[11]*u[n+1+2*nMax] - ks[10]*u[n+2*nMax]*u[1+2*nMax] - ks[11]*u[n+2*nMax] for n=3:nMax-1]
    
    # Tran nMax-omers: (tran₁+tranₙ₋₁->tranₙ) - (tranₙ->tran₁+tranₙ₋₁)
    du[3*nMax] = ks[10]*u[3*nMax-1]*u[1+2*nMax] - ks[11]*u[3*nMax]

    return du
end

# Function to update figure based on system iteration
function animstep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
    step!(integ, 5.0)
	deterministicCisObservable[]  .= integ.u[1:nMax]
    deterministicCisObservable[]  = deterministicCisObservable[]
	deterministicMedObservable[]  .= integ.u[1+nMax:2*nMax]
    deterministicMedObservable[]  = deterministicMedObservable[]
	deterministicTraObservable[] .= integ.u[1+2*nMax:3*nMax]
    deterministicTraObservable[] = deterministicTraObservable[]
end

function resetstep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
    reinit!(integ,erase_sol=true)
    deterministicCisObservable[] .= integ.u[1:nMax]
    deterministicCisObservable[] = deterministicCisObservable[]
	deterministicMedObservable[] .= integ.u[1+nMax:2*nMax]
    deterministicMedObservable[] = deterministicMedObservable[]
	deterministicTraObservable[] .= integ.u[1+2*nMax:3*nMax]
    deterministicTraObservable[] = deterministicTraObservable[]
end

# Initial conditions
u0   = zeros(Float32,nMax*3)

# Set up figure canvas
fig = Figure(resolution=(1700,1500),fontsize=32)
axDiagram = Axis(fig[3,1:4],title="Model diagram",aspect=DataAspect())
image!(axDiagram,rotr90(load(joinpath("_research","model.png"))))
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
deterministicCisObservable = Observable(u0[1:nMax].*volume)
# Set up observable objects for med results
deterministicMedObservable = Observable(u0[1+nMax:2*nMax].*volume)
# Set up observable objects for tran results
deterministicTraObservable = Observable(u0[1+2*nMax:3*nMax].*volume)

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
    (label="k₄,  m₁ → c₁     " , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₅,  m₁+mₙ → mₙ₊₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₆,  mₙ → m₁+mₙ₋₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₇,  m₁ → t₁     " , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₈,  t₁ → m₁     " , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₉,  t₁+tₙ → tₙ₊₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₁₀, tₙ → t₁+tₙ₋₁" , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}"),
    (label="k₁₁, t₁ → ∅      " , range=0.0:0.1:2.0, startvalue=1.0, format="{:.2f}");
)
# Pull parameters from slider positions
kObservables = [s.value for s in parameterSliders.sliders]
ks = ones(Float32,12)

# Setup parameters for ODE
p = @dict nMax ks
p = NamedTuple([pair for pair in p])

# Set up integrator for each iteration
prob = ODEProblem(deterministicModel!, u0, (0.0,1000), p)
integ = init(prob,Tsit5())

# solve!(integ)

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
    resetstep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)
    isrunning[] = false
end

on(run.clicks) do clicks
    @async while isrunning[]       
        isopen(fig.scene) || break # ensures computations stop if closed window
        for i=1:12
            integ.p.ks[i] = kObservables[i][]
        end        
        animstep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax)        
        sleep(0.01) # yield()
    end
end

display(fig)
