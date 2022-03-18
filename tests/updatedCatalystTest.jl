using DrWatson
@quickactivate
using Catalyst
using DiffEqJump
using Plots

golgiMaturationModel = @reaction_network begin
    β, cis --> trans
    ν, trans --> cis
end β ν

p     = (0.01,0.01)
u₀    = [100,0]
tspan = (0.0,500.0)
prob  = DiscreteProblem(golgiMaturationModel, u₀, tspan, p)

jump_prob = JumpProblem(golgiMaturationModel, prob, Direct())

sol = solve(jump_prob, SSAStepper())

plot(sol)

Graph(golgiMaturationModel)


# using Latexify
# display(latexify(golgiMaturationModel))
# eq = latexify(golgiMaturationModel)
# using CairoMakie
# using LaTeXStrings
# fig = Figure()
# ax = Axis(fig[1,1])
# hidedecorations!(ax)
# hidespines!(ax)
# ax.title = eq
# display(fig)
# using MathJaxRenderer
# display(write("tmp.svg",Math(eq)))
#
