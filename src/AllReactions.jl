#
#  AllReactions.jl
#  GolgiModels
#
#  Created by Christopher Revell on dd/mm/yyyy.

module AllReactions

using Catalyst
# using UnPack

# Find all possible reactions pairs that result in oligomers with size <= nMax    
function allReactions(nMax,C,M,T,K₀,K₁,K₂,K₃,K₄,K₅,K₆,K₇,K₈,K₉,K₁₀,K₁₁)

    # @unpack :K₀ :K₁ :K₂ :K₃ :K₄ :K₅ :K₆ :K₇ :K₈ :K₉ :K₁₀ :K₁₁ = p

    # vector to store the Reactions
    reactions = []

    push!(reactions, Reaction(K₀, nothing, [C[1]]))                      #, nothing, [1]))    # Insertion into cis. Reactant ∅; product X[1]; rate: zeroth order kinetics
    for i=1:nMax-1
        push!(reactions, Reaction(K₁, [C[i], C[1]], [C[i+1]]))           #, [1,1], [1])) # cis aggregation: second order kinetics
    end
    for i=2:nMax
        push!(reactions, Reaction(K₂, [C[i]], [C[i-1],C[1]]))                #, [1], [1,1])) # cis splitting: first order kinetics        
    end
    push!(reactions, Reaction(K₃, [C[1]], [M[1]]))                           #, [1], [1])) # cis to medial: first order kinetics 
    
    push!(reactions, Reaction(K₄, [M[1]], [C[1]]))                       #, [1], [1])) # medial to cis: first order kinetics 
    for i=1:nMax-1
        push!(reactions, Reaction(K₅, [M[i],M[1]], [M[i+1]]))        #, [1,1], [1])) # med aggregation: second order kinetics
    end
    for i=2:nMax
        push!(reactions, Reaction(K₆, [M[i]], [M[i-1],M[1]]))           #, [1], [1,1])) # med splitting: first order kinetics                
    end
    push!(reactions, Reaction(K₇, [M[1]], [T[1]]))                          #, [1], [1])) # med to tran: first order kinetics        
    
    push!(reactions, Reaction(K₈, [T[1]], [M[1]]))                           #, [1], [1])) # tran to med: first order kinetics            
    for i=1:nMax-1
        push!(reactions, Reaction(K₉, [T[i],T[1]], [T[i+1]]))        #, [1,1], [1])) # tran aggregation: second order kinetics                
    end
    for i=2:nMax
        push!(reactions, Reaction(K₁₀, [T[i]], [T[i-1],T[1]]))          #, [1], [1,1])) # tran splitting: first order kinetics                
    end
    
    push!(reactions, Reaction(K₁₁ , [T[1]], nothing))               #, [1], nothing)) # tran to ∅: first order kinetics            
    
    return reactions

end

export allReactions

end