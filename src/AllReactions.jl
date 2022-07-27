module AllReactions

using Catalyst

# Find all possible reactions pairs that result in oligomers with size <= nMax    
function allReactions(nReactionsTotal,k,X,volume,p)

    nMax,k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀,k₁₁ = p

    # Vector of reaction rates 
    rates = Float64[]
    # vector to store the Reactions
    reactions = []

    push!(reactions, Reaction(k₀*volume, nothing, [X[1+0*nMax]], nothing, [1]))    # Insertion into cis. Reactant ∅; product X[1]; rate: zeroth order kinetics
    push!(rates,k₀*volume)
    for i=1:nMax-1
        push!(reactions, Reaction(k₁/volume, [X[i+0*nMax],X[1+0*nMax]], [X[i+1+0*nMax]], [1,1], [1])) # cis aggregation: second order kinetics
        push!(rates,k₁/volume)
    end
    for i=2:nMax
        push!(reactions, Reaction(k₂, [X[i+0*nMax]], [X[i-1+0*nMax],X[1+0*nMax]], [1], [1,1])) # cis splitting: first order kinetics
        push!(rates,k₂)
    end
    push!(reactions, Reaction(k₃, [X[1+0*nMax]], [X[1+1*nMax]], [1], [1])) # cis to medial: first order kinetics 
    push!(rates,k₃)
    
    push!(reactions, Reaction(k₄, [X[1+1*nMax]], [X[1+0*nMax]], [1], [1])) # medial to cis: first order kinetics 
    push!(rates,k₄)
    for i=1:nMax-1
        push!(reactions, Reaction(k₅/volume, [X[i+1*nMax],X[1+1*nMax]],[X[i+1+1*nMax]], [1,1], [1])) # med aggregation: second order kinetics
        push!(rates,k₅/volume)
    end
    for i=2:nMax
        push!(reactions, Reaction(k₆, [X[i+1*nMax]],[X[i-1+1*nMax],X[1+1*nMax]], [1], [1,1])) # med splitting: first order kinetics        
        push!(rates,k₆)
    end
    push!(reactions, Reaction(k₇, [X[1+1*nMax]], [X[1+2*nMax]], [1], [1])) # med to tran: first order kinetics        
    push!(rates,k₇)

    push!(reactions, Reaction(k₈, [X[1+2*nMax]], [X[1+1*nMax]], [1], [1])) # tran to med: first order kinetics        
    push!(rates,k₈)
    
    for i=1:nMax-1
        push!(reactions, Reaction(k₉/volume, [X[i+2*nMax],X[1+2*nMax]], [X[i+1+2*nMax]], [1,1], [1])) # tran aggregation: second order kinetics        
        push!(rates,k₉/volume)
    end
    for i=2:nMax
        push!(reactions, Reaction(k₁₀, [X[i+2*nMax]], [X[i-1+2*nMax],X[1+2*nMax]], [1], [1,1])) # tran splitting: first order kinetics        
        push!(rates,k₁₀)
    end
    
    push!(reactions, Reaction(k₁₁, [X[1+2*nMax]], nothing, [1], nothing)) # tran to ∅: first order kinetics        
    push!(rates,k₁₁)
    
    return reactions, rates

end

export allReactions

end