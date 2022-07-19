module Deterministic

using DifferentialEquations

function deterministicModel(du,u,p,t)

    ∅ToCis,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅ = p

    du[1] = ∅ToCis + medToCis*u[1+nMax] - cisAgg*u[1]^2 - cisToMed*u[1]
    du[2:nMax-1] = [cisAgg*u[i-1]*u[1] + cisSplit*u[i+1] - cisAgg*u[i]*u[1] - cisSplit*u[i] for i=2:nMax-1]
    du[nMax] = cisAgg*u[nMax-1]*u[1] - cisSplit*u[nMax]
    du[1+nMax] = cisToMed*u[1] + tranToMed*u[1+2*nMax] - medAgg*u[1+nMax] - medToTran*u[1+nMax]
    du[2+nMax:2*nMax-1] = [cisAgg*u[nMax+i-1]*u[1+nMax] + cisSplit*u[nMax+i+1] - cisAgg*u[nMax+i]*u[1+nMax] - cisSplit*u[nMax+i] for i=2:nMax-1]
    du[nMax] = cisAgg*u[nMax-1]*u[1] - cisSplit*u[nMax]
    
end

function determinstic()


end

export deterministic

end