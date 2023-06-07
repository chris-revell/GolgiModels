using Catalyst, DifferentialEquations, Plots

rs = @reaction_network begin
  c1, S + E --> SE
  c2, SE --> S + E
  c3, SE --> P + E
end
p = (:c1 => 0.00166, :c2 => 0.0001, :c3 => 0.1)
tspan = (0., 100.)
u0 = [:S => 301., :E => 100., :SE => 0., :P => 0.]

# solve ODEs
oprob = ODEProblem(rs, u0, tspan, p)
osol  = solve(oprob, Tsit5())

# solve JumpProblem
u0 = [:S => 301, :E => 100, :SE => 0, :P => 0]
dprob = DiscreteProblem(rs, u0, tspan, p)
jprob = JumpProblem(rs, dprob, Direct())
jsol = solve(jprob, SSAStepper())

Plots.plot(Plots.plot(osol; title = "Reaction Rate Equation ODEs"),
     Plots.plot(jsol; title = "Stochastic Chemical Kinetics Jump Processes");
     layout = (2, 1))