module DeterministicModel

using DifferentialEquations

function deterministicModel!(du,u,p,t)

    nMax,k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀,k₁₁ = p

    # Cis monomers: ∅->cis₁ + med₁->cis₁ + cis₂->2cis₁ - cis₁->med₁ - 2cis₁->cis₂, cisₙ->cis₁+cisₙ₋₁ for n>=3, cis₁+cisₙ->cisₙ₊₁ for n>=2 (n=1 in first line)
    du[1,1] = k₀ - k₃*u[1,1] + k₄*u[1,2] - 2*k₁*u[1,1]^2 - k₁*u[1,1]*sum(u[2:nMax-1,1]) + 2*k₂*u[2,1] + k₂*sum(u[3:nMax,1])
    # Cis n-omers: (cis₁+cisₙ₋₁->cisₙ) + (cisₙ₊₁->cisₙ+cis₁) - (cis₁+cisₙ->cisₙ₊₁) - (cisₙ->cisₙ₋₁+cis₁) for n=2:nMax-1
    du[2:nMax-1,1] .= [k₁*u[n-1,1]*u[1,1] + k₂*u[n+1,1] - k₁*u[n,1]*u[1,1] - k₂*u[n,1] for n=2:nMax-1]
    # Cis nMax-omers: (cis₁+cisₙ₋₁->cisₙ) - (cisₙ->cis₁+cisₙ₋₁)
    du[nMax,1] = k₁*u[nMax-1,1]*u[1,1] - k₂*u[nMax,1]

    # Med monomers: cis₁->med₁ + tran₁->med₁ - med₁->cis₁ - med₁->tran₁ + med₂->2med₁ - med₁->tran₁ - 2med₁->med₂, medₙ->med₁+medₙ₋₁ for n>=3 (n=1,2 in first line), med₁+medₙ->medₙ₊₁ for n>=2 (n=1 in first line)
    du[1,2] = k₃*u[1,1] - k₄*u[1,2] + k₈*u[1,3] - k₇*u[1,2] - 2*k₅*u[1,2]^2 - k₅*u[1,2]*sum(u[2:nMax-1,2]) + 2*k₆*u[2,2] + k₆*sum(u[3:nMax,2]) 
    # Med n-omers: (med₁+medₙ₋₁->medₙ) + (medₙ₊₁->medₙ+med₁) - (med₁+medₙ->medₙ₊₁) - (medₙ->medₙ₋₁+med₁) for n=2:nMax-1
    du[2:nMax-1,2] .= [k₅*u[n-1,2]*u[1,2] + k₆*u[n+1,2] - k₅*u[n,2]*u[1,2] - k₆*u[n,2] for n=2:nMax-1]
    # Med nMax-omers: (med₁+medₙ₋₁->medₙ) - (medₙ->med₁+medₙ₋₁)
    du[nMax,2] = k₅*u[nMax-1,2]*u[1,2] - k₆*u[nMax,2]

    # Tran monomers: med₁->tran₁ - tran₁->∅ - tran₁->med₁ + tran₂->2tran₁ - 2tran₁->tran₂, tranₙ->tran₁+tranₙ₋₁ for n>=3 (n=1,2 in first line), tran₁+tranₙ->tranₙ₊₁ for n>=2 (n=1 in first line)
    du[1,3] = k₇*u[1,2] - k₈*u[1,3] - k₁₁*u[1,3] - 2*k₉*u[1,3]^2 - k₉*u[1,3]*sum(u[2:nMax-1,3]) + 2*k₁₀*u[2,3] + k₁₀*sum(u[3:nMax,3])
    # Tran n-omers: (tran₁+tranₙ₋₁->tranₙ) + (tranₙ₊₁->tranₙ+tran₁) - (tran₁+tranₙ->tranₙ₊₁) - (tranₙ->tranₙ₋₁+tran₁) for n=2:nMax-1
    du[2:nMax-1,3] .= [k₉*u[n-1,3]*u[1,3] + k₁₀*u[n+1,3] - k₉*u[n,3]*u[1,3] - k₁₀*u[n,3] for n=2:nMax-1]
    # Tran nMax-omers: (tran₁+tranₙ₋₁->tranₙ) - (tranₙ->tran₁+tranₙ₋₁)
    du[nMax,3] = k₉*u[nMax-1,3]*u[1,3] - k₁₀*u[nMax,3]

    return du

end

function solveDeterministic(nMax,tMax,p)

    # Setup initial condition of deterministic model
    u0 = zeros(nMax,3)
    # Create ODEProblem object for deterministic model
    prob = ODEProblem(deterministicModel!,u0,(0.0,tMax),p)
    # Solve deterministic model 
    deterministicSol = solve(prob,saveat=tMax/100000)

    return deterministicSol

end

export solveDeterministic#deterministicModel!

end