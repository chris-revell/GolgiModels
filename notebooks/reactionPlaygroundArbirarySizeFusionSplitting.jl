using Catalyst

nMax = 5
@parameters k[1:12]
@variables t
@species C₁(t) C₂(t) C₃(t) C₄(t) C₅(t) M₁(t) M₂(t) M₃(t) M₄(t) M₅(t) T₁(t) T₂(t) T₃(t) T₄(t) T₅(t)

C = [C₁, C₂, C₃, C₄, C₅]
M = [M₁, M₂, M₃, M₄, M₅]
T = [T₁, T₂, T₃, T₄, T₅]

reactionsVec = []

# Source 
push!(reactionsVec, Reaction(k[1], nothing, [C[1]]))  # ∅->C₁
# Homotypic aggregation
for i=1:floor(Int64,nMax/2)
    push!(reactionsVec, Reaction(k[2], [C[i]], [C[2*i]], [2], [1]))   # 2Cᵢ->C₂ᵢ
    push!(reactionsVec, Reaction(k[6], [M[i]], [M[2*i]], [2], [1]))   # 2Mᵢ->M₂ᵢ
    push!(reactionsVec, Reaction(k[10], [T[i]], [T[2*i]], [2], [1]))   # 2Tᵢ->T₂ᵢ
end
# Heterotypic aggregation
for i=1:nMax-1
    for j=i+1:nMax-i
        push!(reactionsVec, Reaction(k[2], [C[i], C[j]], [C[i+j]])) # Cᵢ+Cⱼ->Cᵢ₊ⱼ for all i+j<=nMax
        push!(reactionsVec, Reaction(k[6], [M[i], M[j]], [M[i+j]])) # Mᵢ+Mⱼ->Mᵢ₊ⱼ for all i+j<=nMax
        push!(reactionsVec, Reaction(k[10], [T[i], T[j]], [T[i+j]])) # Tᵢ+Tⱼ->Tᵢ₊ⱼ for all i+j<=nMax
    end
end
# Homotypic splitting
for i in [i for i in 1:nMax if iseven(i)]
    push!(reactionsVec, Reaction(k[3], [C[i]], [C[i÷2]], [1], [2]))   # C₂ᵢ->2Cᵢ
    push!(reactionsVec, Reaction(k[7], [M[i]], [M[i÷2]], [1], [2]))   # M₂ᵢ->2Mᵢ
    push!(reactionsVec, Reaction(k[11], [T[i]], [T[i÷2]], [1], [2]))   # T₂ᵢ->2Tᵢ
end
# Heterotypic splitting
for i=1:nMax
    for j=1:floor(Int64,(i-1)/2)
        push!(reactionsVec, Reaction(k[3], [C[i]], [C[i-j],C[j]]))  # Cᵢ->Cᵢ₋ⱼ+Cⱼ for 3<=i<=nMax
        push!(reactionsVec, Reaction(k[7], [M[i]], [M[i-j],M[j]]))  # Mᵢ->Mᵢ₋ⱼ+Mⱼ for 3<=i<=nMax
        push!(reactionsVec, Reaction(k[11], [T[i]], [T[i-j],T[j]]))  # Tᵢ->Tᵢ₋ⱼ+Tⱼ for 3<=i<=nMax
    end
end
# Maturation
for i=1:nMax
    push!(reactionsVec, Reaction(k[4], [C[i]], [M[i]]))  # Cᵢ->Mᵢ
    push!(reactionsVec, Reaction(k[8], [M[i]], [T[i]]))  # Mᵢ->Tᵢ
end
# Retrograde transport
push!(reactionsVec, Reaction(k[5], [M[1]], [C[1]]))      # M₁->C₁
push!(reactionsVec, Reaction(k[9], [T[1]], [M[1]]))      # T₁->M₁
# Sink
push!(reactionsVec, Reaction(k[12], [T[1]], nothing)) # T₁->∅


@named system = ReactionSystem(reactionsVec, t, [collect(C); collect(M); collect(T)], k, combinatoric_ratelaws=true)

grph = Graph(system)
savegraph(grph,"ArbitrarySizeFusionSplitting.png")