module DeterministicModel

using DifferentialEquations

function deterministicModel!(du,u,p,t)

    nMax,∅ToCis,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅ = p

    # Cis monomers: ∅->cis₁ + med₁->cis₁ + cis₂->2cis₁ - cis₁->med₁ - 2cis₁->cis₂
    du[1] = ∅ToCis + medToCis*u[1+nMax] - cisToMed*u[1] + 2*cisSplit*u[2] - 2*cisAgg*u[1]^2 
    # Cis monomers: cisₙ->cis₁+cisₙ₋₁ for n>=3 (n=1,2 in first line)
    du[1] += sum(cisSplit.*u[3:nMax])
    # Cis monomers: cis₁+cisₙ->cisₙ₊₁ for n>=2 (n=1 in first line)
    du[1] -= sum(cisAgg.*u[2:nMax-1].*u[1])
    # Cis n-omers: (cis₁+cisₙ₋₁->cisₙ) + (cisₙ₊₁->cisₙ+cis₁) - (cis₁+cisₙ->cisₙ₊₁) - (cisₙ->cisₙ₋₁+cis₁) for n=2:nMax-1
    du[2:nMax-1] = [cisAgg*u[n-1]*u[1] + cisSplit*u[n+1] - cisAgg*u[n]*u[1] - cisSplit*u[n] for n=2:nMax-1]
    # Cis nMax-omers: (cis₁+cisₙ₋₁->cisₙ) - (cisₙ->cis₁+cisₙ₋₁)
    du[nMax] = cisAgg*u[nMax-1]*u[1] - cisSplit*u[nMax]

    # Med monomers: cis₁->med₁ + tran₁->med₁ - med₁->cis₁ - med₁->tran₁ + med₂->2med₁ - med₁->tran₁ - 2med₁->med₂
    du[1+nMax] = cisToMed*u[1] + tranToMed*u[1+2*nMax] - medToCis*u[1+nMax] - medToTran*u[1+nMax] + 2*medSplit*u[2+nMax] - 2*medAgg*u[1+nMax]^2 
    # Med monomers: medₙ->med₁+medₙ₋₁ for n>=3 (n=1,2 in first line)
    du[1+nMax] += sum(medSplit.*u[3+nMax:2*nMax])
    # Med monomers: med₁+medₙ->medₙ₊₁ for n>=2 (n=1 in first line)
    du[1+nMax] -= sum(medAgg.*u[2+nMax:2*nMax-1].*u[1+nMax])
    # Med n-omers: (med₁+medₙ₋₁->medₙ) + (medₙ₊₁->medₙ+med₁) - (med₁+medₙ->medₙ₊₁) - (medₙ->medₙ₋₁+med₁) for n=2:nMax-1
    du[2+nMax:2*nMax-1] = [medAgg*u[nMax+n-1]*u[1+nMax] + medSplit*u[nMax+n+1] - medAgg*u[nMax+n]*u[1+nMax] - medSplit*u[n+nMax] for n=2:nMax-1]
    # Med nMax-omers: (med₁+medₙ₋₁->medₙ) - (medₙ->med₁+medₙ₋₁)
    du[2*nMax] = medAgg*u[2*nMax-1]*u[1] - medSplit*u[2*nMax]

    # Tran monomers: med₁->tran₁ - tran₁->∅ - tran₁->med₁ + tran₂->2tran₁ - 2tran₁->tran₂
    du[1+2*nMax] = medToTran*u[1+nMax] - tranTo∅*u[1+2*nMax] - tranToMed*u[1+2*nMax] + 2*tranSplit*u[2+2*nMax] - 2*tranAgg*u[1+2*nMax]^2 
    # Tran monomers: tranₙ->tran₁+tranₙ₋₁ for n>=3 (n=1,2 in first line)
    du[1+2*nMax] += sum(tranSplit.*u[3+2*nMax:3*nMax])
    # Tran monomers: tran₁+tranₙ->tranₙ₊₁ for n>=2 (n=1 in first line)
    du[1+2*nMax] -= sum(tranAgg.*u[2+2*nMax:3*nMax-1].*u[1+2*nMax])
    # Tran n-omers: (tran₁+tranₙ₋₁->tranₙ) + (tranₙ₊₁->tranₙ+tran₁) - (tran₁+tranₙ->tranₙ₊₁) - (tranₙ->tranₙ₋₁+tran₁) for n=2:nMax-1
    du[2+2*nMax:3*nMax-1] = [tranAgg*u[2*nMax+n-1]*u[1+2*nMax] + tranSplit*u[2*nMax+n+1] - tranAgg*u[2*nMax+n]*u[1+2*nMax] - tranSplit*u[n+2*nMax] for n=2:nMax-1]
    # Tran nMax-omers: (tran₁+tranₙ₋₁->tranₙ) - (tranₙ->tran₁+tranₙ₋₁)
    du[3*nMax] = tranAgg*u[3*nMax-1]*u[1] - tranSplit*u[3*nMax]

    return du

end

function solveDeterministic(nMax,tMax,p)

    # Setup initial condition of deterministic model
    u0 = zeros(3*nMax)
    # Create ODEProblem object for deterministic model
    prob = ODEProblem(deterministicModel!,u0,(0.0,tMax),p)
    # Solve deterministic model 
    deterministicSol = solve(prob,saveat=(tMax/100))

    return deterministicSol

end

export solveDeterministic#deterministicModel!

end