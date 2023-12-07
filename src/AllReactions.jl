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
    # Source 
    push!(reactions, Reaction(k[1], nothing, [C[1]]))            # ∅->C₁
    # Aggregation 
    push!(reactions, Reaction(k[2], [C[1]], [C[2]], [2], [1]))   # 2C₁->C₂
    push!(reactions, Reaction(k[6], [M[1]], [M[2]], [2], [1]))   # 2M₁->M₂
    push!(reactions, Reaction(k[10], [T[1]], [T[2]], [2], [1]))  # 2T₁->T₂
    for i=2:nMax-1
        push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]])) # C₁+Cₙ->Cₙ₊₁ for 2<=n<nMax
        push!(reactions, Reaction(k[6], [M[i],M[1]], [M[i+1]]))  # M₁+Mₙ->Mₙ₊₁ for 2<=n<nMax
        push!(reactions, Reaction(k[10], [T[i],T[1]], [T[i+1]])) # T₁+Tₙ->Tₙ₊₁ for 2<=n<nMax
    end
    # Splitting 
    push!(reactions, Reaction(k[3], [C[2]], [C[1]], [1], [2]))   # C₂->2C₁
    push!(reactions, Reaction(k[7], [M[2]], [M[1]], [1], [2]))   # M₂->2M₁
    push!(reactions, Reaction(k[11], [T[2]], [T[1]], [1], [2]))  # T₂->2T₁
    for i=3:nMax
        push!(reactions, Reaction(k[3], [C[i]], [C[i-1],C[1]]))  # Cₙ->C₁+Cₙ₋₁ for 3<=n<=nMax
        push!(reactions, Reaction(k[7], [M[i]], [M[i-1],M[1]]))  # Mₙ->M₁+Mₙ₋₁ for 3<=n<=2nMax
        push!(reactions, Reaction(k[11], [T[i]], [T[i-1],T[1]])) # Tₙ->T₁+Tₙ₋₁ for 3<=n<=2nMax
    end
    # Maturation 
    push!(reactions, Reaction(k[4], [C[1]], [M[1]]))             # C₁->M₁
    push!(reactions, Reaction(k[5], [M[1]], [C[1]]))             # M₁->C₁
    push!(reactions, Reaction(k[8], [M[1]], [T[1]]))             # M₁->T₁
    push!(reactions, Reaction(k[9], [T[1]], [M[1]]))             # T₁->M₁
    # Sink 
    push!(reactions, Reaction(k[12], [T[1]], nothing))           # T₁->∅
    # Set up reaction system object. Collect symbolic state variables into a single vector.
    @named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], k, combinatoric_ratelaws=true)
    return system
end

# Find all possible reactions pairs that result in oligomers with size <= nMax
# Pass symbolic state vectors C, M, and T, and symbolic parameters k and t
function allReactionsHeterogeneous(nMax,C,M,T,k,t)

    # vector to store the Reactions
    reactions = []
    # Source 
    push!(reactions, Reaction(k[1], nothing, [C[1]]))            # ∅->C₁
    # Aggregation 
    push!(reactions, Reaction(k[2], [C[1]], [C[2]], [2], [1]))   # 2C₁->C₂
    push!(reactions, Reaction(k[6], [M[1]], [M[2]], [2], [1]))   # 2M₁->M₂
    push!(reactions, Reaction(k[10], [T[1]], [T[2]], [2], [1]))  # 2T₁->T₂
    for i=2:nMax-1
        push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]])) # C₁+Cₙ->Cₙ₊₁ for 2<=n<nMax
        push!(reactions, Reaction(k[6], [M[i],M[1]], [M[i+1]]))  # M₁+Mₙ->Mₙ₊₁ for 2<=n<nMax
        push!(reactions, Reaction(k[10], [T[i],T[1]], [T[i+1]])) # T₁+Tₙ->Tₙ₊₁ for 2<=n<nMax
    end
    # Splitting 
    push!(reactions, Reaction(k[3]*2^(2/3), [C[2]], [C[1]], [1], [2]))   # C₂->2C₁
    push!(reactions, Reaction(k[7]*2^(2/3), [M[2]], [M[1]], [1], [2]))   # M₂->2M₁
    push!(reactions, Reaction(k[11]*2^(2/3), [T[2]], [T[1]], [1], [2]))  # T₂->2T₁
    for i=3:nMax
        push!(reactions, Reaction(k[3]*i^(2/3), [C[i]], [C[i-1],C[1]]))  # Cₙ->C₁+Cₙ₋₁ for 3<=n<=nMax
        push!(reactions, Reaction(k[7]*i^(2/3), [M[i]], [M[i-1],M[1]]))  # Mₙ->M₁+Mₙ₋₁ for 3<=n<=2nMax
        push!(reactions, Reaction(k[11]*i^(2/3), [T[i]], [T[i-1],T[1]])) # Tₙ->T₁+Tₙ₋₁ for 3<=n<=2nMax
    end
    # Maturation 
    push!(reactions, Reaction(k[4], [C[1]], [M[1]]))             # C₁->M₁
    push!(reactions, Reaction(k[5], [M[1]], [C[1]]))             # M₁->C₁
    push!(reactions, Reaction(k[8], [M[1]], [T[1]]))             # M₁->T₁
    push!(reactions, Reaction(k[9], [T[1]], [M[1]]))             # T₁->M₁
    # Sink 
    push!(reactions, Reaction(k[12], [T[1]], nothing))           # T₁->∅
    # Set up reaction system object. Collect symbolic state variables into a single vector.
    @named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], k, combinatoric_ratelaws=true)
    return system
end

# export allReactions

# end