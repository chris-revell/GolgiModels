#
#  AllReactions.jl
#  GolgiModels
#
#  Created by Christopher Revell on dd/mm/yyyy.

module AllReactions

using Catalyst
# using UnPack

# Find all possible reactions pairs that result in oligomers with size <= nMax    
function allReactions(nMax,C,M,T,k₀,k₁,k₂,k₃,k₄,k₅,k₆,k₇,k₈,k₉,k₁₀,k₁₁)

    # vector to store the Reactions
    reactions = []
    
    push!(reactions, Reaction(k₀, nothing, [C[1]]))            # ∅->c₁
    push!(reactions, Reaction(k₁, [C[1]], [C[2]], [2], [1]))   # 2c₁->c₂
    push!(reactions, Reaction(k₃, [C[2]], [C[1]], [1], [2]))   # c₂->2c₁
    for i=2:nMax-1
        push!(reactions, Reaction(k₁, [C[i], C[1]], [C[i+1]])) # c₁+cₙ->cₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k₂, [C[i]], [C[i-1],C[1]]))  # cₙ->c₁+cₙ₋₁ for 3<=n<=nMax
    end
    push!(reactions, Reaction(k₃, [C[1]], [M[1]]))             # c₁->m₁
    push!(reactions, Reaction(k₄, [M[1]], [C[1]]))             # m₁->c₁

    push!(reactions, Reaction(k₅, [M[1]], [M[2]], [2], [1]))   # 2m₁->m₂
    push!(reactions, Reaction(k₆, [M[2]], [M[1]], [1], [2]))   # m₂->2m₁
    for i=2:nMax-1
        push!(reactions, Reaction(k₅, [M[i],M[1]], [M[i+1]]))  # m₁+mₙ->mₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k₆, [M[i]], [M[i-1],M[1]]))  # mₙ->m₁+mₙ₋₁ for 3<=n<=2nMax
    end
    push!(reactions, Reaction(k₇, [M[1]], [T[1]]))             # m₁->t₁
    push!(reactions, Reaction(k₈, [T[1]], [M[1]]))             # t₁->m₁

    push!(reactions, Reaction(k₉, [T[1]], [T[2]], [2], [1]))   # 2t₁->t₂
    push!(reactions, Reaction(k₁₀, [T[2]], [T[1]], [1], [2]))  # t₂->2t₁
    for i=2:nMax-1
        push!(reactions, Reaction(k₉, [T[i],T[1]], [T[i+1]]))  # t₁+tₙ->tₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k₁₀, [T[i]], [T[i-1],T[1]])) # tₙ->t₁+tₙ₋₁ for 3<=n<=2nMax
    end
    push!(reactions, Reaction(k₁₁ , [T[1]], nothing))          # t₁->∅ 
    
    return reactions

end

export allReactions

end