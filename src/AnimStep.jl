#
#  AnimStep.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

module AnimStep

using GLMakie

# Function to update figure based on system iteration
function animStep!(integ,deterministicCisObservable,deterministicMedObservable,deterministicTraObservable,nMax,xLimTimeAv)
    step!(integ, 10.0)
    # Find time averaged maximum value to set xlim
    xLimTimeAv[1] = (xLimTimeAv[1]*19+maximum(integ.u))/20
    xlims!(axCis,(0.0,1.1*xLimTimeAv[1]))
    xlims!(axMed,(0.0,1.1*xLimTimeAv[1]))
    xlims!(axTra,(0.0,1.1*xLimTimeAv[1]))
	deterministicCisObservable[] .= integ.u[1:nMax]
    deterministicCisObservable[] = deterministicCisObservable[]
	deterministicMedObservable[] .= integ.u[1+nMax:2*nMax]
    deterministicMedObservable[] = deterministicMedObservable[]
	deterministicTraObservable[] .= integ.u[1+2*nMax:3*nMax]
    deterministicTraObservable[] = deterministicTraObservable[]
end

export animStep!

end