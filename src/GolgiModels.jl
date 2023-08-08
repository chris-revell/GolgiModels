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

module GolgiModels

    using PrecompileTools
    using DifferentialEquations
    using Makie
    using WGLMakie; WGLMakie.activate!()
    using DrWatson
    using UnPack
    using GeometryBasics
    using FileIO
    using Catalyst
    using FromFile
    using Format
    using JSServe  

    # function vecObsToFloat(a,b)
    #     for i=1:length(a)
    #         b[i] = a[i][]
    #     end
    # end

    # function pairSymbolToObservable(pairs,sym,obs)
    #     for i=1:length(pairs)
    #         pairs[i] = Pair(sym[i],kStochFactors[i]*kObservables[i][])
    #     end 

    # quickactivate("GolgiModels")
    include(srcdir("AllReactions.jl"))
    include(srcdir("GuiFigureSetup.jl"))
    include(srcdir("DwellTimes.jl"))
    include(srcdir("AnimStep.jl"))
    include(srcdir("RefreshIntegrators.jl"))
    include(srcdir("HattedConstants.jl"))

    function golgiApp(;displayFlag=true)

# 11111111111111111111111        
        # quickactivate("GolgiModels")
         # TODO: handle this by environment variables instead
        JSServe.configure_server!(listen_port=9384, listen_url="0.0.0.0")

        nMax          = 20             # Max compartment size /vesicles
        dt            = 100.0          # Integration time interval /seconds
        V             = 10 # μm³
        ksInit        = [1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]
        k̂             = [1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]
        kStochFactors = [V, 1/V, 1.0, 1.0, 1.0, 1/V, 1.0, 1.0, 1.0, 1/V, 1.0, 1.0]
# 11111111111111111111111


# 22222222222222222222222        
        # Set up figure canvas
        fig, axCis, axMed, axTra, axDwell, parameterSliders, run, reset, linearityToggle, xLimTimeAv, linearityToggle = guiFigureSetup(ksInit)

# 22222222222222222222222
    
# 33333333333333333333333
        # Catalyst system setup
        # Symbolic system parameters: rate constants 
        @parameters k[1:12]
        @variables t
        # Symbolic system variables: cis, medial, and trans compartment size counts 
        @species C(t)[1:nMax] M(t)[1:nMax] T(t)[1:nMax] 
        # Use these parameters and variables to define a reaction system 
        
        system = [refreshSystem(nMax,C,M,T,k,t,true)]
        integODE = [refreshODEs(nMax,C,M,T,k,t,ksInit,system[1])]

        pStoch           = Pair.(collect(k),ksInit)
        u₀MapStoch       = Pair.([collect(C); collect(M); collect(T)], zeros(Int32,3*nMax)) 
        discreteProblem  = [DiscreteProblem(system[1], u₀MapStoch, (0.0,Inf), pStoch)]
        jumpProblem      = [JumpProblem(system[1], discreteProblem[1], Direct(), save_positions=(false,false))] # Converts system to a set of MassActionJumps
        integStoch       = [init(jumpProblem[1], SSAStepper())]
        # integStoch       = [refreshStoch(nMax,C,M,T,k,ksInit,system[1])]
        
# 33333333333333333333333

       
# 44444444444444444444444
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

        dwellTimesValues = zeros(Float32,7)

# 44444444444444444444444        
        
# 55555555555555555555555
        # Set up button actions 
        isrunning = Observable(false)
        on(run.clicks) do clicks
            isrunning[] = !isrunning[]
        end
        on(reset.clicks) do clicks    
            isrunning[] = false            
            # resetStepODE!(integODE,nMax,linearityToggleVal)
            # resetStepStoch!(pStoch,u₀MapStoch,nMax,discreteProblem,tMax,jumpProblem,integStoch,C,M,T,k,ksInit,system,linearityToggleVal)
            # resetObservables(axCis,axMed,axTra,nMax,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,stochTimeAvCisObservable,stochTimeAvMedObservable,stochTimeAvTraObservable,dwellTimeObservable)
            system[1] = refreshSystem(nMax,C,M,T,k,t,linearityToggle.active[])
            integODE[1] = refreshODEs(nMax,C,M,T,k,t,ksInit,system[1])            
            integStoch[1] = refreshStoch!(pStoch,u₀MapStoch,discreteProblem,jumpProblem,zeros(Int32,3*nMax),C,M,T,k,ksInit.*kStochFactors,system[1])
            for i_slider in parameterSliders.sliders
                set_close_to!(i_slider, i_slider.startvalue[])
            end
            resetObservables(axCis,axMed,axTra,nMax,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,stochTimeAvCisObservable,stochTimeAvMedObservable,stochTimeAvTraObservable,dwellTimeObservable)
        end
# 55555555555555555555555
        

        kObsVec = lift([s.value for s in parameterSliders.sliders]...) do values...
            [values...]
        end


        on(run.clicks) do clicks
            while isrunning[]
                isopen(fig.scene) || break # ensures computations stop if closed window

                integODE[1].p .= kObsVec[] 
                animStepODE!(integODE[1],dt,axCis,axMed,axTra,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax,V)
 
                
                hattedConstants!(integODE[1].p,k̂,integODE[1].u,nMax)
                for i=1:7
                    dwellTimesValues[i] = -dwellTimesFuns[i](k̂)
                end
                dwellTimeObservable[] .= dwellTimesValues
                dwellTimeObservable[] = dwellTimeObservable[]
                ylims!(axDwell,(0,maximum(dwellTimesValues)))
                

                # integStoch[1] = refreshStoch(nMax,C,M,T,k,kStochFactors.*kObsVec[],system[1])
                # pStoch = Pair.(collect(kSyms),kStochFactors.*kObsVec[])
                # discreteProblem  = [DiscreteProblem(system, integStoch[1].u, (0.0,Inf), pStoch)]
                # jumpProblem   = [JumpProblem(system, discreteProblem[1], Direct(), save_positions=(false,false))] # Converts system to a set of MassActionJumps
                integStoch[1] = refreshStoch!(pStoch,u₀MapStoch,discreteProblem,jumpProblem,integStoch[1].u,C,M,T,k,kStochFactors.*kObsVec[],system[1])
                animStepStoch!(integStoch[1],dt,axCis,axMed,axTra,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,stochTimeAvCisObservable,stochTimeAvMedObservable,stochTimeAvTraObservable,nMax)
                
                # Find time averaged maximum value to set xlim
                xLimTimeAv[1] = (xLimTimeAv[1]*19+max(maximum(integStoch[1].u),maximum(integODE[1].u)))/20
                xlims!(axCis,(0.0,1.1*xLimTimeAv[1]))
                xlims!(axMed,(0.0,1.1*xLimTimeAv[1]))
                xlims!(axTra,(0.0,1.1*xLimTimeAv[1]))                

                sleep(0.1)
            end
        end

        displayFlag==true ? display(fig) : nothing
    end

    @compile_workload begin
        golgiApp(displayFlag=false)
    end

export golgiApp

end