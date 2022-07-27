module Visualise

using CairoMakie
using DrWatson
using Printf
using Dates
using Statistics

function visualise(fileName,nMax,stochasticSol,params,deterministicSol,stochasticTimeAverages,volume)

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
    stochasticCisObservable             = Observable(stochasticSol.u[1][1:nMax])
    stochasticCisObservableTimeAverage  = Observable(stochasticTimeAverages[1][1:nMax])
    deterministicCisObservable          = Observable(deterministicSol.u[1][:,1].*volume)
    
    # Set up observable objects for med results
    stochasticMedObservable             = Observable(stochasticSol.u[1][1+nMax:2*nMax])
    stochasticMedObservableTimeAverage  = Observable(stochasticTimeAverages[1][1+nMax:2*nMax])
    deterministicMedObservable          = Observable(deterministicSol.u[1][:,2].*volume)
    
    # Set up observable objects for tran results
    stochasticTranObservable            = Observable(stochasticSol.u[1][1+2*nMax:3*nMax])
    stochasticTranObservableTimeAverage = Observable(stochasticTimeAverages[1][1+2*nMax:3*nMax])
    deterministicTranObservable         = Observable(deterministicSol.u[1][:,3].*volume)
    
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
    # animationFilename = savename(Dates.format(Dates.now(),"mm-dd-HH-MM"),params,"mp4",connector="")

    record(fig,datadir("sims","$fileName.mp4"),tSteps; framerate=10) do i
            
        stochasticCisObservable[]             .= stochasticSol.u[i][1:nMax]
        deterministicCisObservable[]          .= deterministicSol.u[i][:,1].*volume
        stochasticCisObservableTimeAverage[]  .= mean(stochasticSol.u[1:i])[1:nMax]
        
        stochasticMedObservable[]             .= stochasticSol.u[i][1+nMax:2*nMax]
        deterministicMedObservable[]          .= deterministicSol.u[i][:,2].*volume
        stochasticMedObservableTimeAverage[]  .= mean(stochasticSol.u[1:i])[1+nMax:2*nMax]
        
        stochasticTranObservable[]            .= stochasticSol.u[i][1+2*nMax:3*nMax]
        deterministicTranObservable[]         .= deterministicSol.u[i][:,3].*volume
        stochasticTranObservableTimeAverage[] .= mean(stochasticSol.u[1:i])[1+2*nMax:3*nMax]

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