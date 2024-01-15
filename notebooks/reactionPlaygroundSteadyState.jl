using Catalyst
include(srcdir("AllReactions.jl"))

nMax = 5
@parameters k₁ k₂ k₃ k₄ k₅ k₆ k₇ k₈ k₉ k₁₀ k₁₁ k₁₂
@variables t
@species C₁(t) C₂(t) C₃(t) C₄(t) C₅(t) M₁(t) M₂(t) M₃(t) M₄(t) M₅(t) T₁(t) T₂(t) T₃(t) T₄(t) T₅(t)

C = [C₁, C₂, C₃, C₄, C₅]
M = [M₁, M₂, M₃, M₄, M₅]
T = [T₁, T₂, T₃, T₄, T₅]

##

system = allReactions(nMax,C,M,T,[k₁, k₂, k₃, k₄, k₅, k₆, k₇, k₈, k₉, k₁₀, k₁₁, k₁₂],t)

##

grph = Graph(system)

##

incidence = incidencemat(system)
compstoichmat = complexstoichmat(system)

##

M = netstoichmat(system)
rxs = reactions(system)
ν = oderatelaw.(rxs)

odes = M*ν

##

# import HomotopyContinuation

# psA = Pair.(collect(k), ones(12))
# psD = Pair.([k₁, k₂, k₃, k₄, k₅, k₆, k₇, k₈, k₉, k₁₀, k₁₁, k₁₂], ones(12))
# psC = Pair(k,ones(12))

# wilhelm_2009_model = @reaction_network begin
#     k1, Y --> 2X
#     k2, 2X --> X + Y
#     k3, X + Y --> Y
#     k4, X --> 0
# end
# psB = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5]

# hc_steady_states(system, psD)

##
using NonlinearSolve
u_guess = [(C.=>zeros(nMax)); (M.=>zeros(nMax)); (T.=>zeros(nMax))]
p = [:k₁=>1.0, :k₂=>1.0, :k₃=>1.0, :k₄=>1.0, :k₅=>1.0, :k₆=>1.0, :k₇=>1.0, :k₈=>1.0, :k₉=>1.0, :k₁₀=>1.0, :k₁₁=>1.0, :k₁₂=>1.0]
nl_prob = NonlinearProblem(system, u_guess, p)
sol = solve(nl_prob)