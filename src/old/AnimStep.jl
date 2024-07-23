#
#  AnimStep.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

# module AnimStep

# Function to update figure based on system iteration

function animStepODE!(integODE,dt,axCis,axMed,axTra,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax,V)
    step!(integODE, dt, true)
	deterministicCisObservable[] .= integODE.u[1:nMax].*V
    deterministicCisObservable[] = deterministicCisObservable[]
	deterministicMedObservable[] .= integODE.u[1+nMax:2*nMax].*V
    deterministicMedObservable[] = deterministicMedObservable[]
	deterministicTraObservable[] .= integODE.u[1+2*nMax:3*nMax].*V
    deterministicTraObservable[] = deterministicTraObservable[]
end


function animStepStoch!(integStoch,dt,axCis,axMed,axTra,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,stochTimeAvCisObservable,stochTimeAvMedObservable,stochTimeAvTraObservable,nMax)
    step!(integStoch, dt, true)
	stochasticCisObservable[] .= integStoch.u[1:nMax]
    stochasticCisObservable[] = stochasticCisObservable[]
	stochasticMedObservable[] .= integStoch.u[1+nMax:2*nMax]
    stochasticMedObservable[] = stochasticMedObservable[]
	stochasticTraObservable[] .= integStoch.u[1+2*nMax:3*nMax]
    stochasticTraObservable[] = stochasticTraObservable[]

    stochTimeAvCisObservable[] .= (stochTimeAvCisObservable[].*19.0.+integStoch.u[1:nMax])./20.0
    stochTimeAvCisObservable[] = stochTimeAvCisObservable[]
    stochTimeAvMedObservable[] .= (stochTimeAvMedObservable[].*19.0.+integStoch.u[1+nMax:2*nMax])./20.0
    stochTimeAvMedObservable[] = stochTimeAvMedObservable[]
    stochTimeAvTraObservable[] .= (stochTimeAvTraObservable[].*19.0.+integStoch.u[1+2*nMax:3*nMax])./20.0
    stochTimeAvTraObservable[] = stochTimeAvTraObservable[]
end

# export animStep!

# end