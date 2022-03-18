using Catalyst
using DiffEqJump
using Plots

sir_model = @reaction_network begin
    β, S + I --> 2I
    ν, I --> R
    l, R --> S
end β ν l

p     = (0.001,0.01,0.01)
u₀    = [100,1,0]
tspan = (0.0,500.0)
prob  = DiscreteProblem(sir_model, u₀, tspan, p)

jump_prob = JumpProblem(sir_model, prob, Direct())

sol = solve(jump_prob, SSAStepper())

plot(sol)
