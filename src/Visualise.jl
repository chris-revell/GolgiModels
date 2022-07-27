module Visualise

using CairoMakie
using DrWatson
using Printf
using Dates
using Statistics

function visualise(fileName,nMax,volume,stochasticSol,stochasticTimeAverages,deterministicSol)

    cisRange = 1:nMax
    medRange = 1+nMax:2*nMax
    tranRange = 1+2*nMax:3*nMax

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
    upperLim=10# upperLim = maximum([max(maximum(deterministicSol.u[i]),maximum(stochasticSol.u[i])) for i=1:101])
    xlims!(axCis,(0,upperLim))
    ylims!(axCis,(0,nMax))
    xlims!(axMed,(0,upperLim))
    ylims!(axMed,(0,nMax))
    xlims!(axTra,(0,upperLim))
    ylims!(axTra,(0,nMax))

    # Set up observable objects for cis results
    stochasticCisObservable             = Observable(stochasticSol.u[1][cisRange])
    stochasticCisObservableTimeAverage  = Observable(stochasticTimeAverages[1][cisRange])
    deterministicCisObservable          = Observable(deterministicSol.u[1][cisRange].*volume)
    
    # Set up observable objects for med results
    stochasticMedObservable             = Observable(stochasticSol.u[1][medRange])
    stochasticMedObservableTimeAverage  = Observable(stochasticTimeAverages[1][tranRange])
    deterministicMedObservable          = Observable(deterministicSol.u[1][medRange].*volume)
    
    # Set up observable objects for tran results
    stochasticTranObservable            = Observable(stochasticSol.u[1][tranRange])
    stochasticTranObservableTimeAverage = Observable(stochasticTimeAverages[1][tranRange])
    deterministicTranObservable         = Observable(deterministicSol.u[1][tranRange].*volume)
    
    # Set up bar and line plots for cis results
    barplot!(axCis, collect(1:nMax), stochasticCisObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:red,0.5))
    lines!(axCis, stochasticCisObservableTimeAverage, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=2)
    lines!(axCis, deterministicCisObservable, collect(1:nMax), color=(:blue,0.5), linewidth=2)
    
    # Set up bar and line plots for med results
    barplot!(axMed, collect(1:nMax), stochasticMedObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:green,0.5))
    lines!(axMed, stochasticMedObservableTimeAverage, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=2)
    lines!(axMed, deterministicMedObservable, collect(1:nMax), color=(:blue,0.5), linewidth=2)
    
    # Set up bar and line plots for tran results
    barplot!(axTra, collect(1:nMax), stochasticTranObservable, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=(:blue,0.5))
    lines!(axTra, stochasticTranObservableTimeAverage, collect(1:nMax), color=(:black,0.5), linestyle="--", linewidth=2)
    lines!(axTra, deterministicTranObservable, collect(1:nMax), color=(:blue,0.5), linewidth=2)

    windowLength = 1000
    tSteps = range(1+windowLength,length(stochasticSol.t),step=windowLength)    

    record(fig,datadir("sims","$fileName.mp4"),tSteps; framerate=10) do i
            
        stochasticCisObservable[]             .= stochasticSol.u[i][cisRange]
        deterministicCisObservable[]          .= deterministicSol.u[i][cisRange].*volume
        stochasticCisObservableTimeAverage[]  .= mean(stochasticSol.u[1:i])[cisRange] #stochasticTimeAverages[(i-1)÷windowLength][cisRange] 
        
        stochasticMedObservable[]             .= stochasticSol.u[i][medRange]
        deterministicMedObservable[]          .= deterministicSol.u[i][medRange].*volume
        stochasticMedObservableTimeAverage[]  .= mean(stochasticSol.u[1:i])[medRange] #stochasticTimeAverages[(i-1)÷windowLength][medRange] 
        
        stochasticTranObservable[]            .= stochasticSol.u[i][tranRange]
        deterministicTranObservable[]         .= deterministicSol.u[i][tranRange].*volume
        stochasticTranObservableTimeAverage[] .= mean(stochasticSol.u[1:i])[tranRange] #stochasticTimeAverages[(i-1)÷windowLength][tranRange] 

        stochasticCisObservable[]             = stochasticCisObservable[]
        deterministicCisObservable[]          = deterministicCisObservable[]
        stochasticCisObservableTimeAverage[]  = stochasticCisObservableTimeAverage[]
        stochasticMedObservable[]             = stochasticMedObservable[]
        deterministicMedObservable[]          = deterministicMedObservable[]
        stochasticMedObservableTimeAverage[]  = stochasticMedObservableTimeAverage[]
        stochasticTranObservable[]            = stochasticTranObservable[]
        deterministicTranObservable[]         = deterministicTranObservable[]
        stochasticTranObservableTimeAverage[] = stochasticTranObservableTimeAverage[]

    end
    
    save(datadir("sims","$fileName.png"),fig)

end

export visualise

end