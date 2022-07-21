module Visualise

using CairoMakie
using DrWatson
using Printf
using Dates

function visualise(nMax,stochasticSol,params,deterministicSol,volume)

    fig = Figure()
    axCis = Axis(fig[1,1],aspect=0.5)
    axMed = Axis(fig[1,2],aspect=0.5,yticksvisible=false,yticklabelsvisible=false)
    axTra = Axis(fig[1,3],aspect=0.5,yticksvisible=false,yticklabelsvisible=false)
    Label(fig[1,1,Bottom()],"Cis",textsize=32)
    Label(fig[1,2,Bottom()],"Medial",textsize=32)
    Label(fig[1,3,Bottom()],"Trans",textsize=32)    
    axCis.yticks = 0:10:nMax
    axCis.ylabel = "Compartment size"    

    uObservableCis = Observable(stochasticSol.u[1][1:nMax])
    detUObservableCis = Observable(deterministicSol.u[1][1:nMax].*volume)
    uObservableMed = Observable(stochasticSol.u[1][1+nMax:2*nMax])
    detUObservableMed = Observable(deterministicSol.u[1][1+nMax:2*nMax].*volume)
    uObservableTra = Observable(stochasticSol.u[1][1+2*nMax:3*nMax])
    detUObservableTra = Observable(deterministicSol.u[1][1+2*nMax:3*nMax].*volume)
    barplot!(axCis, collect(1:nMax), uObservableCis, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:red)
    lines!(axCis, detUObservableCis, collect(1:nMax), color=:black)
    barplot!(axMed, collect(1:nMax), uObservableMed, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:green)
    lines!(axMed, detUObservableMed, collect(1:nMax), color=:black)
    barplot!(axTra, collect(1:nMax), uObservableTra, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:blue)
    lines!(axTra, detUObservableTra, collect(1:nMax), color=:black)

    tSteps = range(1,length(stochasticSol.t),step=1)

    animationFilename = savename(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"),params,"mp4",connector="")   

    record(fig,datadir("sims",animationFilename),tSteps; framerate=10) do i
        
        CairoMakie.xlims!(axCis,(0.0,max(10,maximum(stochasticSol.u[i][1:nMax]))))
        CairoMakie.xlims!(axMed,(0.0,max(10,maximum(stochasticSol.u[i][1+nMax:2*nMax]))))
        CairoMakie.xlims!(axTra,(0.0,max(10,maximum(stochasticSol.u[i][1+2*nMax:3*nMax]))))
        uObservableCis[] .= stochasticSol.u[i][1:nMax]
        detUObservableCis[] .= deterministicSol.u[i][1:nMax].*volume
        uObservableMed[] .= stochasticSol.u[i][1+nMax:2*nMax]
        detUObservableMed[] .= deterministicSol.u[i][1+nMax:2*nMax].*volume
        uObservableTra[] .= stochasticSol.u[i][1+2*nMax:3*nMax]
        detUObservableTra[] .= deterministicSol.u[i][1+2*nMax:3*nMax].*volume
        
        uObservableCis[] = uObservableCis[]
        detUObservableCis[] = detUObservableCis[]
        uObservableMed[] = uObservableMed[]
        detUObservableMed[] = detUObservableMed[]
        uObservableTra[] = uObservableTra[]
        detUObservableTra[] = detUObservableTra[]
    end

end

export visualise

end