#
#  ResetStep.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

# module ResetStep

# Function to reset figure
function resetStepODE!(integ,axCis,axMed,axTra,cisObservable,medObservable,traObservable,nMax,dwellTimeObservable)
    reinit!(integ,erase_sol=true)
    cisObservable[] .= integ.u[1:nMax]
    cisObservable[] = cisObservable[]
	medObservable[] .= integ.u[1+nMax:2*nMax]
    medObservable[] = medObservable[]
	traObservable[] .= integ.u[1+2*nMax:3*nMax]
    traObservable[] = traObservable[]
    xlims!(axCis,(0.0,5.0))
    xlims!(axMed,(0.0,5.0))
    xlims!(axTra,(0.0,5.0))
    dwellTimeObservable[] .= zeros(Float32,7)
    dwellTimeObservable[] = dwellTimeObservable[]
end

function resetStepStoch!(pStoch,u₀MapStoch,nMax,discreteProblem,tMax,jumpProblem,integStoch,cisObservable,medObservable,traObservable,stochTimeAvCisObservable,stochTimeAvMedObservable,stochTimeAvTraObservable,C,M,T,k,ksInit,system)
    pStoch .= Pair.(collect(k),ksInit)
    # Map symbolic state vector to vector of values. Collect symbolic state variables into a single vector.
    u₀MapStoch .= Pair.([collect(C); collect(M); collect(T)], zeros(Int32,3*nMax)) 
    # Create problem object
    discreteProblem  .= [DiscreteProblem(system, u₀MapStoch, (0.0,tMax), pStoch)]
    jumpProblem   .= [JumpProblem(system, discreteProblem[1], Direct(), save_positions=(false,false))] # Converts system to a set of MassActionJumps
    # Create integrator object
    integStoch .= [init(jumpProblem[1], SSAStepper())]
    cisObservable[] .= integStoch[1].u[1:nMax]
    cisObservable[] = cisObservable[]
	medObservable[] .= integStoch[1].u[1+nMax:2*nMax]
    medObservable[] = medObservable[]
	traObservable[] .= integStoch[1].u[1+2*nMax:3*nMax]
    traObservable[] = traObservable[]

    stochTimeAvCisObservable[] .= zeros(Float32,nMax)
    stochTimeAvCisObservable[] = stochTimeAvCisObservable[]
    stochTimeAvMedObservable[] .= zeros(Float32,nMax)
    stochTimeAvMedObservable[] = stochTimeAvMedObservable[]
    stochTimeAvTraObservable[] .= zeros(Float32,nMax)
    stochTimeAvTraObservable[] = stochTimeAvTraObservable[]
end


# export resetStep!

# end