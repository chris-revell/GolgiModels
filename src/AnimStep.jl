#
#  AnimStep.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

module AnimStep

using GLMakie

# Function to update figure based on system iteration
function animStep!(integ,axCis,axMed,axTra,cisObservable,medObservable,traObservable,nMax,xLimTimeAv)
    step!(integ, 10.0)
	cisObservable[] .= integ.u[1:nMax]
    cisObservable[] = cisObservable[]
	medObservable[] .= integ.u[1+nMax:2*nMax]
    medObservable[] = medObservable[]
	traObservable[] .= integ.u[1+2*nMax:3*nMax]
    traObservable[] = traObservable[]
end

export animStep!

end