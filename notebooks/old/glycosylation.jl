using Catalyst
using DifferentialEquations

nMax = 5           # Max compartment size /vesicles
V = 10              # μm³

ksInit = [1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0]
k̂ = [1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0]
kStochFactors = [V, 1 / V, 1.0, 1.0, 1.0, 1 / V, 1.0, 1.0, 1.0, 1 / V, 1.0, 1.0]

# Catalyst system setup
# @parameters k[0:5]
@parameters k₀ k₁ k₂ k₃ k₄ k₅
@variables t
@species C(t)[1:nMax] Q(t)[1:nMax] E(t) S(t) 
# Use these parameters and variables to define a reaction system 

reactions = []
# Source 
push!(reactions, Reaction(k₀, nothing, [C[1]]))            # ∅->C₁
# Enzyme complex formation 
for n=1:nMax
    push!(reactions, Reaction(k₁, [C[n], E], [Q[n]])) # Cₙ+E->Qₙ
end
# Enzyme complex splitting without polymerisation
for n=1:nMax
    push!(reactions, Reaction(k₂, [Q[n]], [C[n], E])) # Qₙ->Cₙ+E
end
# Complex polymerisation
for n=1:nMax-1
    push!(reactions, Reaction(k₃, [Q[n], S], [C[n+1]])) # Qₙ+S->Cₙ₊₁
end
# Polymer depolymerisation into complex
for n=1:nMax-1
    push!(reactions, Reaction(k₄, [C[n+1]], [Q[n], S])) # Cₙ₊₁->Qₙ+S
end
# Removal of polymer 
for n=1:nMax
    push!(reactions, Reaction(k₅, [C[n]], nothing))           # T₁->∅
end 

# Set up reaction system object. Collect symbolic state variables into a single vector.
@named system = ReactionSystem(reactions, t, [collect(C); collect(Q); [S, E]], [k₀, k₁, k₂, k₃, k₄, k₅], combinatoric_ratelaws=true)

Graph(system)


# @parameters k₀ k₁ k₂ k₃ k₄ k₅
# @variables t
# @species C₁(t) C₂(t) C₃(t) C₄(t) C₅(t) Q₁(t) Q₂(t) Q₃(t) Q₄(t) Q₅(t) S(t) E(t)
# reactions = []
# push!(reactions, Reaction(k₀, nothing, [C₁]))            # ∅->C₁
# # Enzyme complex formation 
# push!(reactions, Reaction(k₁, [C₁, E], [Q₁])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₁, [C₂, E], [Q₂])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₁, [C₃, E], [Q₃])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₁, [C₄, E], [Q₄])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₁, [C₅, E], [Q₅])) # Cₙ+E->Qₙ

# push!(reactions, Reaction(k₂, [Q₁], [C₁, E])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₂, [Q₂], [C₂, E])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₂, [Q₃], [C₃, E])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₂, [Q₄], [C₄, E])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₂, [Q₅], [C₅, E])) # Cₙ+E->Qₙ

# push!(reactions, Reaction(k₃, [Q₁, S], [C₂])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₃, [Q₂, S], [C₃])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₃, [Q₃, S], [C₄])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₃, [Q₄, S], [C₅])) # Cₙ+E->Qₙ

# push!(reactions, Reaction(k₄, [C₂], [Q₁, S])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₄, [C₃], [Q₂, S])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₄, [C₄], [Q₃, S])) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₄, [C₅], [Q₄, S])) # Cₙ+E->Qₙ

# push!(reactions, Reaction(k₅, [C₁], nothing)) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₅, [C₂], nothing)) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₅, [C₃], nothing)) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₅, [C₄], nothing)) # Cₙ+E->Qₙ
# push!(reactions, Reaction(k₅, [C₅], nothing)) # Cₙ+E->Qₙ

# # Set up reaction system object. Collect symbolic state variables into a single vector.
# @named system = ReactionSystem(reactions, t, [C₁,C₂,C₃,C₄,C₅,Q₁,Q₂,Q₃,Q₄,Q₅,S,E], [k₀, k₁, k₂, k₃, k₄, k₅], combinatoric_ratelaws=true)

# Graph(system)