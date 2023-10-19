#
#  RefreshIntegrators.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

using OrdinaryDiffEq
using ForwardDiff

function refreshSystem(nMax,C,M,T,k,t,linearityToggleVal)
    if linearityToggleVal
        system = allReactions(nMax,C,M,T,k,t)
    else
        system = allReactionsNonLinear(nMax,C,M,T,k,t)
    end
end

function refreshODEs(nMax,C,M,T,k,t,ks,system)
    # Map symbolic paramters to values. Collect symbolic parameters into a vector.
    pODE = Pair.(collect(k),ks)
    # Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
    u₀MapODE = Pair.([collect(C); collect(M); collect(T)], zeros(Float64,3*nMax))
    # Create problem object
    odeProblem = ODEProblem(system,u₀MapODE,(0.0,Inf),pODE)
    # Create integrator object
    integODE = init(odeProblem,Rodas4P(),maxiters=1e6) #KenCarp4())
    
    return integODE
end

function refreshStoch!(pStoch,u₀MapStoch,discreteProblem,jumpProblem,u,C,M,T,kSyms,kVals,system)
    pStoch .= Pair.(collect(kSyms),kVals)
    u₀MapStoch .= Pair.([collect(C); collect(M); collect(T)], u) 
    discreteProblem[1] = DiscreteProblem(system, u₀MapStoch, (0.0,Inf), pStoch)
    jumpProblem[1]     = JumpProblem(system, discreteProblem[1], Direct(), save_positions=(false,false))
    return init(jumpProblem[1], SSAStepper())
end

function resetObservables(nMax,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,stochTimeAvCisObservable,stochTimeAvMedObservable,stochTimeAvTraObservable,dwellTimeObservable,xLimTimeAv)

    stochasticCisObservable[] .= zeros(Int32,nMax)
    stochasticCisObservable[] = stochasticCisObservable[]
	stochasticMedObservable[] .= zeros(Int32,nMax)
    stochasticMedObservable[] = stochasticMedObservable[]
	stochasticTraObservable[] .= zeros(Int32,nMax)
    stochasticTraObservable[] = stochasticTraObservable[]

    stochTimeAvCisObservable[] .= zeros(Float64,nMax)
    stochTimeAvCisObservable[] = stochTimeAvCisObservable[]
    stochTimeAvMedObservable[] .= zeros(Float64,nMax)
    stochTimeAvMedObservable[] = stochTimeAvMedObservable[]
    stochTimeAvTraObservable[] .= zeros(Float64,nMax)
    stochTimeAvTraObservable[] = stochTimeAvTraObservable[]

    deterministicCisObservable[] .= zeros(Float64,nMax)
    deterministicCisObservable[] = deterministicCisObservable[]
	deterministicMedObservable[] .= zeros(Float64,nMax)
    deterministicMedObservable[] = deterministicMedObservable[]
	deterministicTraObservable[] .= zeros(Float64,nMax)
    deterministicTraObservable[] = deterministicTraObservable[]

    dwellTimeObservable[] .= zeros(Float64,7)
    dwellTimeObservable[] = dwellTimeObservable[]

    xLimTimeAv[] = 5.0
    xLimTimeAv[] = xLimTimeAv[]

end

# export resetStep!

# end