### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ ece8a9a2-1427-11ee-1e57-c9fb83e28f88
begin
	using DrWatson
	DrWatson.quickactivate(".")
end

# ╔═╡ df623d7d-ce61-4500-a30c-d8a1db7ceac8
begin
	using Symbolics
	using Catalyst
	using LinearAlgebra
	using SparseArrays
	using Latexify
	using FileIO
	using Images
	using FromFile
end

# ╔═╡ dbb036e4-c94e-42c8-a308-a46e091649ba
load("$(homedir())/Postdoc/Code/GolgiModels/supplementary/GolgiCompartmentModel.png")

# ╔═╡ 6932679f-b5ee-414b-9775-08e2dac966fa
# Parameters
begin
	nMax    = 5             # Max compartment size
	dt      = 100.0
	tMax    = Inf
	ksInit = [1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]
end

# ╔═╡ fbe20878-7e02-48cf-847f-6afdd93e1910
begin
	@parameters k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁ k₁₂
	@variables t
	@species ∅₁ C₁(t) C₂(t) C₃(t) C₄(t) C₅(t) M₁(t) M₂(t) M₃(t) M₄(t) M₅(t) T₁(t) T₂(t) T₃(t) T₄(t) T₅(t) ∅₂
end

# ╔═╡ df094732-0238-404f-8d69-fc864acad6fe
begin
	k = [k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁ k₁₂]
	C = [C₁ C₂ C₃ C₄ C₅]
	M = [M₁ M₂ M₃ M₄ M₅]
	T = [T₁ T₂ T₃ T₄ T₅]
end

# ╔═╡ 5a3bb540-b5c0-46de-b324-7238786ba8f5
begin
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
    
end	

# ╔═╡ bb3f1646-79fc-4fd9-94ee-38496a4ffa88
@named system = ReactionSystem(reactions)#, t, [C; M; T], k, combinatoric_ratelaws=false)


# ╔═╡ f2cc1da8-7831-4665-8ca4-05f625a64c07
grph = Graph(system)

# ╔═╡ 0147de42-053b-4399-ab65-71942bfe0d54
savegraph(grph, "fullGraph.png")#, fmt="png")

# ╔═╡ Cell order:
# ╠═ece8a9a2-1427-11ee-1e57-c9fb83e28f88
# ╠═df623d7d-ce61-4500-a30c-d8a1db7ceac8
# ╠═dbb036e4-c94e-42c8-a308-a46e091649ba
# ╠═6932679f-b5ee-414b-9775-08e2dac966fa
# ╠═fbe20878-7e02-48cf-847f-6afdd93e1910
# ╠═df094732-0238-404f-8d69-fc864acad6fe
# ╠═5a3bb540-b5c0-46de-b324-7238786ba8f5
# ╠═bb3f1646-79fc-4fd9-94ee-38496a4ffa88
# ╠═f2cc1da8-7831-4665-8ca4-05f625a64c07
# ╠═0147de42-053b-4399-ab65-71942bfe0d54
