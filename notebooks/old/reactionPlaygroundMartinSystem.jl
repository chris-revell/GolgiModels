using Catalyst

nMax = 5
@parameters k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁ k₁₂
@variables t
@species C₁(t) C₂(t) C₃(t) C₄(t) C₅(t) M₁(t) M₂(t) M₃(t) M₄(t) M₅(t) T₁(t) T₂(t) T₃(t) T₄(t) T₅(t)

C = [C₁, C₂, C₃, C₄, C₅]
M = [M₁, M₂, M₃, M₄, M₅]
T = [T₁, T₂, T₃, T₄, T₅]

reactionsVec = []

# Source 
push!(reactionsVec, Reaction(k₁, nothing, [C[1]]))  # ∅->C₁
# Aggregation
push!(reactionsVec, Reaction(k₂, [C[1]], [C[2]], [2], [1]))   # 2C₁->C₂
push!(reactionsVec, Reaction(k₆, [M[1]], [M[2]], [2], [1]))   # 2M₁->M₂
push!(reactionsVec, Reaction(k₁₀, [T[1]], [T[2]], [2], [1]))   # 2T₁->T₂
for i=2:nMax-1
    push!(reactionsVec, Reaction(k₂, [C[i], C[1]], [C[i+1]])) # C₁+Cₙ->Cₙ₊₁ for 2<=n<nMax
    push!(reactionsVec, Reaction(k₆, [M[i], M[1]], [M[i+1]])) # M₁+Mₙ->Mₙ₊₁ for 2<=n<nMax
    push!(reactionsVec, Reaction(k₁₀, [T[i], T[1]], [T[i+1]])) # T₁+Tₙ->Tₙ₊₁ for 2<=n<nMax
end
# Splitting
push!(reactionsVec, Reaction(k₃, [C[2]], [C[1]], [1], [2]))   # C₂->2C₁
push!(reactionsVec, Reaction(k₇, [M[2]], [M[1]], [1], [2]))   # M₂->2M₁
push!(reactionsVec, Reaction(k₁₁, [T[2]], [T[1]], [1], [2]))   # T₂->2T₁
for i=3:nMax
    push!(reactionsVec, Reaction(k₃, [C[i]], [C[i-1],C[1]]))  # Cₙ->C₁+Cₙ₋₁ for 3<=n<=nMax
    push!(reactionsVec, Reaction(k₇, [M[i]], [M[i-1],M[1]]))  # Mₙ->M₁+Mₙ₋₁ for 3<=n<=nMax
    push!(reactionsVec, Reaction(k₁₁, [T[i]], [T[i-1],T[1]]))  # Tₙ->T₁+Tₙ₋₁ for 3<=n<=nMax
end
# Maturation
for i=2:nMax
    push!(reactionsVec, Reaction(k₄, [C[i]], [M[i]]))  # Cₙ->Mₙ
    push!(reactionsVec, Reaction(k₈, [M[i]], [T[i]]))  # Mₙ->Tₙ
end
# Retrograde transport
push!(reactionsVec, Reaction(k₅, [M[1]], [C[1]]))      # M₁->C₁
push!(reactionsVec, Reaction(k₉, [T[1]], [M[1]]))      # T₁->M₁
# Sink
push!(reactionsVec, Reaction(k₁₂, [T[1]], nothing)) # T₁->∅


@named system = ReactionSystem(reactionsVec, t, [collect(C); collect(M); collect(T)], [k₁, k₂, k₃, k₄, k₅, k₆, k₇, k₈, k₉, k₁₀, k₁₁, k₁₂], combinatoric_ratelaws=true)

##

grph = Graph(system)

##

incidencemat(system)

complexstoichmat(system)

##

M = netstoichmat(system)
rxs = reactions(system)
ν = oderatelaw.(rxs)

odes = M*ν

##

import HomotopyContinuation

psA = Pair.(collect(k), ones(12))
psD = Pair.([k₁, k₂, k₃, k₄, k₅, k₆, k₇, k₈, k₉, k₁₀, k₁₁, k₁₂], ones(12))
psC = Pair(k,ones(12))

wilhelm_2009_model = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
end
psB = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5]

hc_steady_states(system, psD)
