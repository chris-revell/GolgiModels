module Visualise

using CairoMakie
using DrWatson
using Printf
using Dates

function visualise(nMax,jsol,params)

    fig = Figure()
    axCis = Axis(fig[1,1],aspect=0.5)
    axMed = Axis(fig[1,2],aspect=0.5,yticksvisible=false,yticklabelsvisible=false)
    axTra = Axis(fig[1,3],aspect=0.5,yticksvisible=false,yticklabelsvisible=false)
    # hidedecorations!(axMed)
    # hidedecorations!(axTra)    
    Label(fig[1,1,Bottom()],"Cis",textsize=32)
    Label(fig[1,2,Bottom()],"Medial",textsize=32)
    Label(fig[1,3,Bottom()],"Trans",textsize=32)    
    axCis.yticks = 0:10:nMax
    axCis.ylabel = "Compartment size"
    # axMed.yticks = 0:10:nMax
    # axMed.ylabel = "Compartment size"
    # axTra.yticks = 1:10:nMax
    # axTra.ylabel = "Compartment size"
    # ax.title = "Size distribution at t=$(@sprintf("%.2f", jsol.t[end]))"
    # CairoMakie.ylims!(ax,(0.5,nMax+0.5))
    # CairoMakie.xlims!(ax,(0.0,4.0))
    # save(datadir("test.png"),fig)
    # display(fig)

    uObservableCis = Observable(jsol.u[1][1:nMax])
    uObservableMed = Observable(jsol.u[1][1+nMax:2*nMax])
    uObservableTra = Observable(jsol.u[1][1+2*nMax:3*nMax])
    barplot!(axCis, collect(1:nMax), uObservableCis, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:red)
    # barplot!(axCis, collect(1:nMax), -uObservableCis, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:red)
    barplot!(axMed, collect(1:nMax), uObservableMed, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:green)
    # barplot!(axMed, collect(1:nMax), -uObservableMed, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:green)
    barplot!(axTra, collect(1:nMax), uObservableTra, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:blue)
    # barplot!(axTra, collect(1:nMax), -uObservableTra, direction=:x, bins=collect(0.5:1.0:nMax+0.5), color=:blue)

    tSteps = range(1,length(jsol.t),step=1)

    animationFilename = savename(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"),params,"mp4",connector="")   

    record(fig,datadir("sims",animationFilename),tSteps; framerate=10) do i
        # ax.title = "t=$(@sprintf("%.2f", jsol.t[i]))"
        display(i)
        CairoMakie.xlims!(axCis,(0.0,max(10,maximum(jsol.u[i][1:nMax]))))
        CairoMakie.xlims!(axMed,(0.0,max(10,maximum(jsol.u[i][1+nMax:2*nMax]))))
        CairoMakie.xlims!(axTra,(0.0,max(10,maximum(jsol.u[i][1+2*nMax:3*nMax]))))
        uObservableCis[] .= jsol.u[i][1:nMax]
        uObservableMed[] .= jsol.u[i][1+nMax:2*nMax]
        uObservableTra[] .= jsol.u[i][1+2*nMax:3*nMax]
        uObservableCis[] = uObservableCis[]
        uObservableMed[] = uObservableMed[]
        uObservableTra[] = uObservableTra[]
    end

end

export visualise

end