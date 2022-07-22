module AllReactions

using Catalyst

# Find all possible reactions pairs that result in oligomers with size <= nMax    
function allReactions(nReactionsTotal,k,X,volume,p)

    nMax,∅ToCis,cisAgg,cisSplit,cisToMed,medToCis,medAgg,medSplit,medToTran,tranToMed,tranAgg,tranSplit,tranTo∅ = p

    # Vector of reaction rates 
    rates = Float64[]
    # vector to store the Reactions
    reactions = []

    push!(reactions, Reaction(∅ToCis*volume, nothing, [X[1+0*nMax]], nothing, [1]))    # Insertion into cis. Reactant ∅; product X[1]; rate: zeroth order kinetics
    push!(rates,∅ToCis*volume)
    for i=1:nMax-1
        push!(reactions, Reaction(cisAgg/volume, [X[i+0*nMax],X[1+0*nMax]], [X[i+1+0*nMax]], [1,1], [1])) # cis aggregation: second order kinetics
        push!(rates,cisAgg/volume)
    end
    for i=2:nMax
        push!(reactions, Reaction(cisSplit, [X[i+0*nMax]], [X[i-1+0*nMax],X[1+0*nMax]], [1], [1,1])) # cis splitting: first order kinetics
        push!(rates,cisSplit)
    end
    push!(reactions, Reaction(cisToMed, [X[1+0*nMax]], [X[1+1*nMax]], [1], [1])) # cis to medial: first order kinetics 
    push!(rates,cisToMed)
    
    push!(reactions, Reaction(medToCis, [X[1+1*nMax]], [X[1+0*nMax]], [1], [1])) # medial to cis: first order kinetics 
    push!(rates,medToCis)
    for i=1:nMax-1
        push!(reactions, Reaction(medAgg/volume, [X[i+1*nMax],X[1+1*nMax]],[X[i+1+1*nMax]], [1,1], [1])) # med aggregation: second order kinetics
        push!(rates,medAgg/volume)
    end
    for i=2:nMax
        push!(reactions, Reaction(medSplit, [X[i+1*nMax]],[X[i-1+1*nMax],X[1+1*nMax]], [1], [1,1])) # med splitting: first order kinetics        
        push!(rates,medSplit)
    end
    push!(reactions, Reaction(medToTran, [X[1+1*nMax]], [X[1+2*nMax]], [1], [1])) # med to tran: first order kinetics        
    push!(rates,medToTran)

    push!(reactions, Reaction(tranToMed, [X[1+2*nMax]], [X[1+1*nMax]], [1], [1])) # tran to med: first order kinetics        
    push!(rates,tranToMed)
    
    for i=1:nMax-1
        push!(reactions, Reaction(tranAgg/volume, [X[i+2*nMax],X[1+2*nMax]], [X[i+1+2*nMax]], [1,1], [1])) # tran aggregation: second order kinetics        
        push!(rates,tranAgg/volume)
    end
    for i=2:nMax
        push!(reactions, Reaction(tranSplit, [X[i+2*nMax]], [X[i-1+2*nMax],X[1+2*nMax]], [1], [1,1])) # tran splitting: first order kinetics        
        push!(rates,tranSplit)
    end
    
    push!(reactions, Reaction(tranTo∅, [X[1+2*nMax]], nothing, [1], nothing)) # tran to ∅: first order kinetics        
    push!(rates,tranTo∅)
    
    return reactions, rates

end

export allReactions

end