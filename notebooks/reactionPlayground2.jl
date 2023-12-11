using Catalyst

nMax = 5
@parameters k[1:12]
@variables t
@species C₁(t) C₂(t) C₃(t) C₄(t) C₅(t) M₁(t) M₂(t) M₃(t) M₄(t) M₅(t) T₁(t) T₂(t) T₃(t) T₄(t) T₅(t)

C = [C₁, C₂, C₃, C₄, C₅]
M = [M₁, M₂, M₃, M₄, M₅]
T = [T₁, T₂, T₃, T₄, T₅]

reactions = []

# # Source 
# push!(reactions, Reaction(k[1], nothing, [C[1]]))            # ∅->C₁
# # Aggregation 
# push!(reactions, Reaction(k[2], [C[1]], [C[2]], [2], [1]))   # 2C₁->C₂
# push!(reactions, Reaction(k[6], [M[1]], [M[2]], [2], [1]))   # 2M₁->M₂
# push!(reactions, Reaction(k[10], [T[1]], [T[2]], [2], [1]))  # 2T₁->T₂
# for i=2:nMax-1
#     push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]])) # C₁+Cₙ->Cₙ₊₁ for 2<=n<nMax
#     push!(reactions, Reaction(k[6], [M[i],M[1]], [M[i+1]]))  # M₁+Mₙ->Mₙ₊₁ for 2<=n<nMax
#     push!(reactions, Reaction(k[10], [T[i],T[1]], [T[i+1]])) # T₁+Tₙ->Tₙ₊₁ for 2<=n<nMax
# end
# # Splitting 
# push!(reactions, Reaction(k[3], [C[2]], [C[1]], [1], [2]))   # C₂->2C₁
# push!(reactions, Reaction(k[7], [M[2]], [M[1]], [1], [2]))   # M₂->2M₁
# push!(reactions, Reaction(k[11], [T[2]], [T[1]], [1], [2]))  # T₂->2T₁
# for i=3:nMax
#     push!(reactions, Reaction(k[3], [C[i]], [C[i-1],C[1]]))  # Cₙ->C₁+Cₙ₋₁ for 3<=n<=nMax
#     push!(reactions, Reaction(k[7], [M[i]], [M[i-1],M[1]]))  # Mₙ->M₁+Mₙ₋₁ for 3<=n<=2nMax
#     push!(reactions, Reaction(k[11], [T[i]], [T[i-1],T[1]])) # Tₙ->T₁+Tₙ₋₁ for 3<=n<=2nMax
# end
# # Maturation 
# push!(reactions, Reaction(k[4], [C[1]], [M[1]]))             # C₁->M₁
# push!(reactions, Reaction(k[5], [M[1]], [C[1]]))             # M₁->C₁
# push!(reactions, Reaction(k[8], [M[1]], [T[1]]))             # M₁->T₁
# push!(reactions, Reaction(k[9], [T[1]], [M[1]]))             # T₁->M₁
# # Sink 
# push!(reactions, Reaction(k[12], [T[1]], nothing))           # T₁->∅


# Source 
push!(reactions, Reaction(k[1], nothing, [C[1]]))  # ∅->C₁
# Aggregation
push!(reactions, Reaction(k[2], [C[1]], [C[2]], [2], [1]))   # 2C₁->C₂
push!(reactions, Reaction(k[6], [M[1]], [M[2]], [2], [1]))   # 2M₁->M₂
push!(reactions, Reaction(k[10], [T[1]], [T[2]], [2], [1]))   # 2T₁->T₂
for i=2:nMax-1
    push!(reactions, Reaction(k[2], [C[i], C[1]], [C[i+1]])) # C₁+Cₙ->Cₙ₊₁ for 2<=n<nMax
    push!(reactions, Reaction(k[6], [M[i], M[1]], [M[i+1]])) # M₁+Mₙ->Mₙ₊₁ for 2<=n<nMax
    push!(reactions, Reaction(k[10], [T[i], T[1]], [T[i+1]])) # T₁+Tₙ->Tₙ₊₁ for 2<=n<nMax
end
# Splitting
push!(reactions, Reaction(k[3], [C[2]], [C[1]], [1], [2]))   # C₂->2C₁
push!(reactions, Reaction(k[7], [M[2]], [M[1]], [1], [2]))   # M₂->2M₁
push!(reactions, Reaction(k[11], [T[2]], [T[1]], [1], [2]))   # T₂->2T₁
for i=3:nMax
    push!(reactions, Reaction(k[3], [C[i]], [C[i-1],C[1]]))  # Cₙ->C₁+Cₙ₋₁ for 3<=n<=nMax
    push!(reactions, Reaction(k[7], [M[i]], [M[i-1],M[1]]))  # Mₙ->M₁+Mₙ₋₁ for 3<=n<=nMax
    push!(reactions, Reaction(k[11], [T[i]], [T[i-1],T[1]]))  # Tₙ->T₁+Tₙ₋₁ for 3<=n<=nMax
end
# Maturation
for i=2:nMax
    push!(reactions, Reaction(k[4], [C[i]], [M[i]]))  # Cₙ->Mₙ
    push!(reactions, Reaction(k[8], [M[i]], [T[i]]))  # Mₙ->Tₙ
end
# Retrograde transport
push!(reactions, Reaction(k[5], [M[1]], [C[1]]))      # M₁->C₁
push!(reactions, Reaction(k[9], [T[1]], [M[1]]))      # T₁->M₁
# Sink
push!(reactions, Reaction(k[12], [T[1]], nothing)) # T₁->∅


@named system = ReactionSystem(reactions, t, [collect(C); collect(M); collect(T)], k, combinatoric_ratelaws=true)

grph = Graph(system)
savegraph(grph,"martinSystemGraph.png")