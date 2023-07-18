using DifferentialEquations
using Catalyst
using JSServe
using WGLMakie; WGLMakie.activate!()
import JSServe.TailwindDashboard as D


App() do 

    nMax    = 20             # Max compartment size
    dt      = 100.0
    tMax    = Inf
    ksInit = [1.0,1.0,1.0,1.0]
    @parameters k[1:4]
    @variables t 
    @species C(t)[1:nMax]

    reactions = []
    push!(reactions, Reaction(k[1], nothing, [C[1]]))
    push!(reactions, Reaction(k[2]*2^(2/3), [C[1]], [C[2]], [2], [1])) 
    push!(reactions, Reaction(k[3], [C[2]], [C[1]], [1], [2]))
    for i=2:nMax-1
        push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]]))
    end
    for i=3:nMax
        push!(reactions, Reaction(k[3]*i^(2/3), [C[i]], [C[i-1],C[1]]))
    end
    push!(reactions, Reaction(k[4], [C[1]], nothing))
    @named system = ReactionSystem(reactions, t, collect(C), k, combinatoric_ratelaws=false)

    pODE = Pair.(collect(k),ksInit)
    u₀MapODE = Pair.(collect(C), zeros(Float32,nMax))
    odeProblem = ODEProblem(system,u₀MapODE,(0.0,tMax),pODE)
    integODE = init(odeProblem,KenCarp3())

    fig = Figure(resolution=(1700,1500),fontsize=32)
    axCis = Axis(fig[1,1], aspect=0.55, ylabel = "Compartment size")
    deterministicCisObservable = Observable(zeros(Float32, nMax))
    lines!(axCis, deterministicCisObservable, collect(1:nMax), color=(:red,1.0),   linewidth=6)

    run = Makie.Button(fig[2,1]; label = "Start/Stop", tellwidth = false)

    sl_x = D.Slider("k1", 0:0.01:1.2)
    sl_y = D.Slider("k2", 0:0.01:1.2)

    lift(sl_x.value, sl_y.value) do x, y
        integODE.p[1] = x
        integODE.p[2] = y
    end

    isrunning = Observable(false)
    on(run.clicks) do clicks
        isrunning[] = !isrunning[]
    end
    on(run.clicks) do clicks
        @async while isrunning[]       
            isopen(fig.scene) || break
            step!(integODE, dt, true)
            deterministicCisObservable[] .= integODE.u[1:nMax]
            deterministicCisObservable[] = deterministicCisObservable[]
            sleep(0.1)
        end
    end 

    DOM.div(JSServe.TailwindCSS, sl_x, sl_y, fig)
end