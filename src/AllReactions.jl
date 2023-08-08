#
#  AllReactions.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

# module AllReactions

using Catalyst

# Find all possible reactions pairs that result in oligomers with size <= nMax
# Pass symbolic state vectors C, M, and T, and symbolic parameters k and t
function allReactions(nMax,C,M,T,k,t)

    # vector to store the Reactions
    reactions = []
    push!(reactions, Reaction(k[1], nothing, [C[1]]))            # ∅->c₁
    push!(reactions, Reaction(k[2], [C[1]], [C[2]], [2], [1]))   # 2c₁->c₂
    push!(reactions, Reaction(k[3], [C[2]], [C[1]], [1], [2]))   # c₂->2c₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]])) # c₁+cₙ->cₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[3], [C[i]], [C[i-1],C[1]]))  # cₙ->c₁+cₙ₋₁ for 3<=n<=nMax
    end
    push!(reactions, Reaction(k[4], [C[1]], [M[1]]))             # c₁->m₁
    push!(reactions, Reaction(k[5], [M[1]], [C[1]]))             # m₁->c₁

    push!(reactions, Reaction(k[6], [M[1]], [M[2]], [2], [1]))   # 2m₁->m₂
    push!(reactions, Reaction(k[7], [M[2]], [M[1]], [1], [2]))   # m₂->2m₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[6], [M[i],M[1]], [M[i+1]]))  # m₁+mₙ->mₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[7], [M[i]], [M[i-1],M[1]]))  # mₙ->m₁+mₙ₋₁ for 3<=n<=2nMax
    end
    push!(reactions, Reaction(k[8], [M[1]], [T[1]]))             # m₁->t₁
    push!(reactions, Reaction(k[9], [T[1]], [M[1]]))             # t₁->m₁

    push!(reactions, Reaction(k[10], [T[1]], [T[2]], [2], [1]))  # 2t₁->t₂
    push!(reactions, Reaction(k[11], [T[2]], [T[1]], [1], [2]))  # t₂->2t₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[10], [T[i],T[1]], [T[i+1]])) # t₁+tₙ->tₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[11], [T[i]], [T[i-1],T[1]])) # tₙ->t₁+tₙ₋₁ for 3<=n<=2nMax
    end
    push!(reactions, Reaction(k[12], [T[1]], nothing))           # t₁->∅
    # Set up reaction system object. Collect symbolic state variables into a single vector.
    @named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], k, combinatoric_ratelaws=true)
    
    return system

end

# Find all possible reactions pairs that result in oligomers with size <= nMax
# Pass symbolic state vectors C, M, and T, and symbolic parameters k and t
function allReactionsNonLinear(nMax,C,M,T,k,t)

    # vector to store the Reactions
    reactions = []
    push!(reactions, Reaction(k[1], nothing, [C[1]]))            # ∅->c₁
    push!(reactions, Reaction(k[2]*2^(2/3), [C[1]], [C[2]], [2], [1]))   # 2c₁->c₂
    push!(reactions, Reaction(k[3], [C[2]], [C[1]], [1], [2]))   # c₂->2c₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]])) # c₁+cₙ->cₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[3]*i^(2/3), [C[i]], [C[i-1],C[1]]))  # cₙ->c₁+cₙ₋₁ for 3<=n<=nMax
    end
    push!(reactions, Reaction(k[4], [C[1]], [M[1]]))             # c₁->m₁
    push!(reactions, Reaction(k[5], [M[1]], [C[1]]))             # m₁->c₁

    push!(reactions, Reaction(k[6]*2^(2/3), [M[1]], [M[2]], [2], [1]))   # 2m₁->m₂
    push!(reactions, Reaction(k[7], [M[2]], [M[1]], [1], [2]))   # m₂->2m₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[6], [M[i],M[1]], [M[i+1]]))  # m₁+mₙ->mₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[7]*i^(2/3), [M[i]], [M[i-1],M[1]]))  # mₙ->m₁+mₙ₋₁ for 3<=n<=2nMax
    end
    push!(reactions, Reaction(k[8], [M[1]], [T[1]]))             # m₁->t₁
    push!(reactions, Reaction(k[9], [T[1]], [M[1]]))             # t₁->m₁

    push!(reactions, Reaction(k[10]*2^(2/3), [T[1]], [T[2]], [2], [1]))  # 2t₁->t₂
    push!(reactions, Reaction(k[11], [T[2]], [T[1]], [1], [2]))  # t₂->2t₁
    for i=2:nMax-1
        push!(reactions, Reaction(k[10], [T[i],T[1]], [T[i+1]])) # t₁+tₙ->tₙ₊₁ for 2<=n<nMax
    end
    for i=3:nMax
        push!(reactions, Reaction(k[11]*i^(2/3), [T[i]], [T[i-1],T[1]])) # tₙ->t₁+tₙ₋₁ for 3<=n<=2nMax
    end
    push!(reactions, Reaction(k[12], [T[1]], nothing))           # t₁->∅
    # Set up reaction system object. Collect symbolic state variables into a single vector.
    @named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], k, combinatoric_ratelaws=false)
    
    return system

end

# export allReactions

# end