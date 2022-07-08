module Visualise

using CairoMakie
using DrWatson
using Printf
using Dates

function visualise(nMax,jsol,params)

    fig = Figure()
    ax = Axis(fig[1,1])
    ax.xticks = (1:1:3, ["cis", "medial", "trans"])
    ax.yticks = 1:10:nMax
    ax.xlabel = "Maturation state"
    ax.ylabel = "Compartment size"
    ax.title = "Size distribution at t=$(@sprintf("%.2f", jsol.t[end]))"
    CairoMakie.ylims!(ax,(0.5,nMax+0.5))
    CairoMakie.xlims!(ax,(0.0,4.0))
    # save(datadir("test.png"),fig)
    # display(fig)

    uObservableCis = Observable(jsol.u[1][1:nMax])
    uObservableMed = Observable(jsol.u[1][1+nMax:2*nMax])
    uObservableTra = Observable(jsol.u[1][1+2*nMax:3*nMax])
    hist!(ax, uObservableCis, scale_to=-0.4, offset=1, direction=:x,bins=collect(0.5:1.0:7.5),color=:red)
    hist!(ax, uObservableCis, scale_to=0.4, offset=1, direction=:x,bins=collect(0.5:1.0:7.5),color=:red)
    hist!(ax, uObservableMed, scale_to=-0.4, offset=2, direction=:x,bins=collect(0.5:1.0:7.5),color=:green)
    hist!(ax, uObservableMed, scale_to=0.4, offset=2, direction=:x,bins=collect(0.5:1.0:7.5),color=:green)
    hist!(ax, uObservableTra, scale_to=-0.4, offset=3, direction=:x,bins=collect(0.5:1.0:7.5),color=:blue)
    hist!(ax, uObservableTra, scale_to=0.4, offset=3, direction=:x,bins=collect(0.5:1.0:7.5),color=:blue)

    tSteps = range(1,length(jsol.t),step=1)

    animationFilename = savename(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"),params,"mp4",connector="")   

    record(fig,datadir("sims",animationFilename),tSteps; framerate=10) do i
        ax.title = "t=$(@sprintf("%.2f", jsol.t[i]))"
        display(i)
        uObservableCis[] = jsol.u[i][1:nMax]
        uObservableMed[] = jsol.u[i][1+nMax:2*nMax]
        uObservableTra[] = jsol.u[i][1+2*nMax:3*nMax]
        uObservableCis[] = uObservableCis[]
        uObservableMed[] = uObservableMed[]
        uObservableTra[] = uObservableTra[]
    end

end

export visualise

end