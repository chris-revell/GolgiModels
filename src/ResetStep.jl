#
#  ResetStep.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

# module ResetStep

# Function to reset figure
function resetStepODE!(integODE,nMax,ksInit, linearityToggleVal)
    if linearityToggleVal

    else
        system[1] = allReactionsNonLinear(nMax,C,M,T,k,t)
        # Map symbolic paramters to values. Collect symbolic parameters into a vector.
        pODE = Pair.(collect(k),ksInit)
        # Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
        u₀MapODE = Pair.([collect(C); collect(M); collect(T)], zeros(Float32,3*nMax))
        # Create problem object
        odeProblem = ODEProblem(system[1],u₀MapODE,(0.0,tMax),pODE)
        # Create integrator object
        integODE[] = init(odeProblem,KenCarp3())

    # reinit!(integ,erase_sol=true)
    
end

function resetStepStoch!(pStoch,u₀MapStoch,nMax,discreteProblem,tMax,jumpProblem,integStoch,C,M,T,k,ksInit,system,linearityToggleVal)
    pStoch .= Pair.(collect(k),ksInit)
    # Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
    u₀MapStoch .= Pair.([collect(C); collect(M); collect(T)], zeros(Int32,3*nMax)) 
    # Create problem object
    discreteProblem  .= [DiscreteProblem(system, u₀MapStoch, (0.0,tMax), pStoch)]
    jumpProblem   .= [JumpProblem(system, discreteProblem[1], Direct(), save_positions=(false,false))] # Converts system to a set of MassActionJumps
    # Create integrator object
    integStoch .= [init(jumpProblem[1], SSAStepper())]
end

function resetObservables(axCis,axMed,axTra,nMax,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,stochTimeAvCisObservable,stochTimeAvMedObservable,stochTimeAvTraObservable,dwellTimeObservable)

    stochasticCisObservable[] .= zeros(Int32,nMax)
    stochasticCisObservable[] = cisObservable[]
	stochasticMedObservable[] .= zeros(Int32,nMax)
    stochasticMedObservable[] = medObservable[]
	stochasticTraObservable[] .= zeros(Int32,nMax)
    stochasticTraObservable[] = traObservable[]

    stochTimeAvCisObservable[] .= zeros(Float32,nMax)
    stochTimeAvCisObservable[] = stochTimeAvCisObservable[]
    stochTimeAvMedObservable[] .= zeros(Float32,nMax)
    stochTimeAvMedObservable[] = stochTimeAvMedObservable[]
    stochTimeAvTraObservable[] .= zeros(Float32,nMax)
    stochTimeAvTraObservable[] = stochTimeAvTraObservable[]

    deterministicCisObservable[] .= zeros(Float32,nMax)
    deterministicCisObservable[] = cisObservable[]
	deterministicMedObservable[] .= zeros(Float32,nMax)
    deterministicMedObservable[] = medObservable[]
	deterministicTraObservable[] .= zeros(Float32,nMax)
    deterministicTraObservable[] = traObservable[]

    xlims!(axCis,(0.0,5.0))
    xlims!(axMed,(0.0,5.0))
    xlims!(axTra,(0.0,5.0))
    dwellTimeObservable[] .= zeros(Float32,7)
    dwellTimeObservable[] = dwellTimeObservable[]

end

# export resetStep!

# end