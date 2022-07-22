module DeterministicModel

using DifferentialEquations

function deterministicModel!(du,u,p,t)

    nMax,∅ToCis,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅ = p

    # Cis monomers: ∅->cis₁ + med₁->cis₁ + cis₂->2cis₁ - cis₁->med₁ - 2cis₁->cis₂
    du[1,1] = ∅ToCis + medToCis*u[1,2] - cisToMed*u[1,1] + 2*cisSplit*u[2,1] - 2*cisAgg*u[1,1]^2 
    # Cis monomers: cisₙ->cis₁+cisₙ₋₁ for n>=3 (n=1,2 in first line)
    du[1,1] += sum(cisSplit.*u[3:nMax,1])
    # Cis monomers: cis₁+cisₙ->cisₙ₊₁ for n>=2 (n=1 in first line)
    du[1,1] -= sum(cisAgg.*u[2:nMax-1,1].*u[1,1])
    # Cis n-omers: (cis₁+cisₙ₋₁->cisₙ) + (cisₙ₊₁->cisₙ+cis₁) - (cis₁+cisₙ->cisₙ₊₁) - (cisₙ->cisₙ₋₁+cis₁) for n=2:nMax-1
    du[2:nMax-1,1] = [cisAgg*u[n-1,1]*u[1,1] + cisSplit*u[n+1,1] - cisAgg*u[n,1]*u[1,1] - cisSplit*u[n,1] for n=2:nMax-1]
    # Cis nMax-omers: (cis₁+cisₙ₋₁->cisₙ) - (cisₙ->cis₁+cisₙ₋₁)
    du[nMax,1] = cisAgg*u[nMax-1,1]*u[1,1] - cisSplit*u[nMax,1]

    # Med monomers: cis₁->med₁ + tran₁->med₁ - med₁->cis₁ - med₁->tran₁ + med₂->2med₁ - med₁->tran₁ - 2med₁->med₂
    du[1,2] = cisToMed*u[1,1] + tranToMed*u[1,3] - medToCis*u[1,2] - medToTran*u[1,2] + 2*medSplit*u[2,2] - 2*medAgg*u[1,2]^2 
    # Med monomers: medₙ->med₁+medₙ₋₁ for n>=3 (n=1,2 in first line)
    du[1,2] += sum(medSplit.*u[3:nMax,2])
    # Med monomers: med₁+medₙ->medₙ₊₁ for n>=2 (n=1 in first line)
    du[1,2] -= sum(medAgg.*u[2:nMax-1,2].*u[1,2])
    # Med n-omers: (med₁+medₙ₋₁->medₙ) + (medₙ₊₁->medₙ+med₁) - (med₁+medₙ->medₙ₊₁) - (medₙ->medₙ₋₁+med₁) for n=2:nMax-1
    du[2:nMax-1,2] = [medAgg*u[n-1,2]*u[1,2] + medSplit*u[n+1,2] - medAgg*u[n,2]*u[1,2] - medSplit*u[n,2] for n=2:nMax-1]
    # Med nMax-omers: (med₁+medₙ₋₁->medₙ) - (medₙ->med₁+medₙ₋₁)
    du[nMax,2] = medAgg*u[nMax-1,2]*u[1,2] - medSplit*u[nMax,2]

    # Tran monomers: med₁->tran₁ - tran₁->∅ - tran₁->med₁ + tran₂->2tran₁ - 2tran₁->tran₂
    du[1,3] = medToTran*u[1,2] - tranTo∅*u[1,3] - tranToMed*u[1,3] + 2*tranSplit*u[2,3] - 2*tranAgg*u[1,3]^2 
    # Tran monomers: tranₙ->tran₁+tranₙ₋₁ for n>=3 (n=1,2 in first line)
    du[1,3] += sum(tranSplit.*u[3:nMax,3])
    # Tran monomers: tran₁+tranₙ->tranₙ₊₁ for n>=2 (n=1 in first line)
    du[1,3] -= sum(tranAgg.*u[2:nMax-1,3].*u[1,3])
    # Tran n-omers: (tran₁+tranₙ₋₁->tranₙ) + (tranₙ₊₁->tranₙ+tran₁) - (tran₁+tranₙ->tranₙ₊₁) - (tranₙ->tranₙ₋₁+tran₁) for n=2:nMax-1
    du[2:nMax-1,3] = [tranAgg*u[n-1,3]*u[1,3] + tranSplit*u[n+1,3] - tranAgg*u[n,3]*u[1,3] - tranSplit*u[n,3] for n=2:nMax-1]
    # Tran nMax-omers: (tran₁+tranₙ₋₁->tranₙ) - (tranₙ->tran₁+tranₙ₋₁)
    du[nMax,3] = tranAgg*u[nMax-1,3]*u[1,3] - tranSplit*u[nMax,3]

    return du

end

function solveDeterministic(nMax,tMax,p)

    # Setup initial condition of deterministic model
    u0 = zeros(nMax,3)
    # Create ODEProblem object for deterministic model
    prob = ODEProblem(deterministicModel!,u0,(0.0,tMax),p)
    # Solve deterministic model 
    deterministicSol = solve(prob,saveat=(tMax/100))

    return deterministicSol

end

export solveDeterministic#deterministicModel!

end