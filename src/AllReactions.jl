module AllReactions

using Catalyst

function allReactions(nMax,X,∅ToCis,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅)

    # Find all possible reactions pairs that result in oligomers with size <= nMax
    reactants = []
    products = []
    rates = Float64[]
    stoichiometryIn = []
    stoichiometryOut = []

    push!(reactants,nothing)          # Insertion into cis reactant: ∅
    push!(products,[X[1]])            # Insertion into cis product: single vesicle in cis bucket X[1]
    push!(rates,∅ToCis*volume)        # Insertion rate: zeroth order kinetics
    push!(stoichiometryIn,nothing)
    push!(stoichiometryOut,[1])
    for i=1:nMax-1
        push!(reactants,[X[i],X[1]])  # cis aggregation reactants 
        push!(products,[X[i+1]])      # cis aggregation products
        push!(rates,cisAgg)           # cis aggregation rate: second order kinetics
        push!(stoichiometryIn,[1,1])
        push!(stoichiometryOut,[1])
    end
    for i=2:nMax
        push!(reactants,[X[i]])       # cis splitting reactants
        push!(products,[X[i-1],X[1]]) # cis splitting products 
        push!(rates,cisSplit)         # cis splitting rates: first order kinetics
        push!(stoichiometryIn,[1])
        push!(stoichiometryOut,[1,1])
    end
    push!(reactants,[X[1]])           # cis to medial reactant
    push!(products,[X[nMax+1]])       # cis to medial product 
    push!(rates,cisToMed)             # cis to medial rate: first order kinetics 
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,[1])
    push!(reactants,[X[1+nMax]])      # medial to cis reactant
    push!(products,[X[1]])            # medial to cis product 
    push!(rates,medToCis)             # medial to cis rate: first order kinetics 
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,[1])

    for i=nMax+1:2*nMax-1
        push!(reactants,[X[i],X[1+nMax]])  # medial aggregation reactants 
        push!(products,[X[i+1]])           # medial aggregation products
        push!(rates,medAgg)                # medial aggregation rate: second order kinetics 
        push!(stoichiometryIn,[1,1])
        push!(stoichiometryOut,[1])
    end
    for i=nMax+2:2*nMax
        push!(reactants,[X[i]])            # medial splitting reactants
        push!(products,[X[i-1],X[1+nMax]]) # medial splitting products 
        push!(rates,medSplit)              # medial splitting rate: first order kinetics 
        push!(stoichiometryIn,[1])
        push!(stoichiometryOut,[1,1])
    end
    push!(reactants,[X[nMax+1]])           # medial to trans reactant
    push!(products,[X[2*nMax+1]])          # medial to trans product 
    push!(rates,medToTran)                 # medial to trans rate: first order kinetics 
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,[1])
    push!(reactants,[X[1+2*nMax]])         # trans to medial reactant
    push!(products,[X[1+nMax]])            # trans to medial product 
    push!(rates,tranToMed)                 # trans to medial rate: first order kinetics 
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,[1])

    for i=2*nMax+1:3*nMax-1
        push!(reactants,[X[i],X[1+2*nMax]]) # trans aggregation reactants 
        push!(products,[X[i+1]])            # trans aggregation products
        push!(rates,tranAgg)                # trans aggregation rate: second order kinetics 
        push!(stoichiometryIn,[1,1])
        push!(stoichiometryOut,[1]) 
    end
    for i=2*nMax+2:3*nMax
        push!(reactants,[X[i]])             # trans splitting reactants
        push!(products,[X[i-1],X[1+2*nMax]])# trans splitting products 
        push!(rates,tranSplit)              # trans splitting rate: first order kinetics
        push!(stoichiometryIn,[1])
        push!(stoichiometryOut,[1,1])
    end

    # Removal 
    push!(reactants,[X[2*nMax+1]])          # Removal from trans reactant: single vesicle in trans
    push!(products,nothing)                 # Removal from trans product: ∅
    push!(rates,tranTo∅)                    # Removal from trans rate: first order kinetics 
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,nothing)

    return reactants,products,rates,stoichiometryIn,stoichiometryOut

end

export allReactions

end