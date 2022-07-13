module AllReactions

using Catalyst

function allReactions(nMax,X)

    # Find all possible reactions pairs that result in oligomers with size <= nMax
    reactants = []
    products = []
    rates = Float64[]
    stoichiometryIn = []
    stoichiometryOut = []

    push!(reactants,nothing) # Insertion into cis reactant: ∅
    push!(products,[X[1]]) # Insertion into cis product: single vesicle in cis bucket 
    push!(rates,1.0)
    push!(stoichiometryIn,nothing)
    push!(stoichiometryOut,[1])
    for i=1:nMax-1
        push!(reactants,[X[i],X[1]]) # cis aggregation reactants 
        push!(products,[X[i+1]])  # cis aggregation products
        push!(rates,1.0)       # cis aggregation rates 
        push!(stoichiometryIn,[1,1])
        push!(stoichiometryOut,[1])
    end
    for i=2:nMax
        push!(reactants,[X[i]])    # cis splitting reactants
        push!(products,[X[i-1],X[1]]) # cis splitting products 
        push!(rates,1.0)        # cis splitting rates 
        push!(stoichiometryIn,[1])
        push!(stoichiometryOut,[1,1])
    end
    push!(reactants,[X[1]])      # cis to medial reactant
    push!(products,[X[nMax+1]])  # cis to medial product 
    push!(rates,1.0)          # cis to medial rate
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,[1])
    push!(reactants,[X[1+nMax]]) # medial to cis reactant
    push!(products,[X[1]])       # medial to cis product 
    push!(rates,1.0)          # medial to cis rate 
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,[1])

    for i=nMax+1:2*nMax-1
        push!(reactants,[X[i],X[1]]) # cis aggregation reactants 
        push!(products,[X[i+1]])  # cis aggregation products
        push!(rates,1.0)       # cis aggregation rates 
        push!(stoichiometryIn,[1,1])
        push!(stoichiometryOut,[1])
    end
    for i=nMax+2:2*nMax
        push!(reactants,[X[i]])    # cis splitting reactants
        push!(products,[X[i-1],X[1]]) # cis splitting products 
        push!(rates,1.0)        # cis splitting rates 
        push!(stoichiometryIn,[1])
        push!(stoichiometryOut,[1,1])
    end
    push!(reactants,[X[nMax+1]])      # cis to medial reactant
    push!(products,[X[2*nMax+1]])  # cis to medial product 
    push!(rates,1.0)          # cis to medial rate 
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,[1])
    push!(reactants,[X[1+2*nMax]]) # medial to cis reactant
    push!(products,[X[1+nMax]])       # medial to cis product 
    push!(rates,1.0)          # medial to cis rate 
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,[1])

    for i=2*nMax+1:3*nMax-1
        push!(reactants,[X[i],X[1]]) # cis aggregation reactants 
        push!(products,[X[i+1]])  # cis aggregation products
        push!(rates,1.0)       # cis aggregation rates 
        push!(stoichiometryIn,[1,1])
        push!(stoichiometryOut,[1]) 
    end
    for i=2*nMax+2:3*nMax
        push!(reactants,[X[i]])    # cis splitting reactants
        push!(products,[X[i-1],X[1]]) # cis splitting products 
        push!(rates,1.0)        # cis splitting rates 
        push!(stoichiometryIn,[1])
        push!(stoichiometryOut,[1,1])
    end

    # Removal 
    push!(reactants,[X[2*nMax+1]])# Removal from trans reactant: single vesicle in trans bucket
    push!(products,nothing)         # Removal from trans product: ∅
    push!(rates,1.0)
    push!(stoichiometryIn,[1])
    push!(stoichiometryOut,nothing)

    return reactants,products,rates,stoichiometryIn,stoichiometryOut

end

export allReactions

end