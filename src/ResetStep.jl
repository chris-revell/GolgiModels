#
#  ResetStep.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

module ResetStep

using GLMakie

# Function to reset figure
function resetStep!(integ,stochasticCisObservable,stochasticMedObservable,stochasticTraObservable,nMax)
    reinit!(integ,erase_sol=true)
    stochasticCisObservable[] .= integ.u[1:nMax]
    stochasticCisObservable[] = stochasticCisObservable[]
	stochasticMedObservable[] .= integ.u[1+nMax:2*nMax]
    stochasticMedObservable[] = stochasticMedObservable[]
	stochasticTraObservable[] .= integ.u[1+2*nMax:3*nMax]
    stochasticTraObservable[] = stochasticTraObservable[]
end


export resetStep!

end