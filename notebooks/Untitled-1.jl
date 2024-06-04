

using DifferentialEquations, SciMLOperators, SparseArrays

_A = [2 -1; -3 -5] / 5
A = MatrixOperator(_A)
prob = ODEProblem(A, [1.0, -1.0], (1.0, 6.0))
sol = solve(prob, LinearExponential())

# function update_func(A, u, p, t)
#     A[1, 1] = cos(t)
#     A[2, 1] = sin(t)
#     A[1, 2] = -sin(t)
#     A[2, 2] = cos(t)
# end
# A = DiffEqArrayOperator(ones(2, 2), update_func = update_func)
# prob = ODEProblem(A, ones(2), (1.0, 6.0))
# sol = solve(prob, MagnusGL6(), dt = 1 / 10)


function update_func(A, u, p, t)
    A[1, 1] = 0
    A[2, 1] = sin(u[1])
    A[1, 2] = -1
    A[2, 2] = 0
end
A = DiffEqArrayOperator(sparse(ones(2, 2)), update_func = update_func)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol = solve(prob, LieRK4(), dt = 1 / 4)



function update_func(A, u, p, t)
    A[1, 1] = 0
    A[2, 1] = 1
    A[1, 2] = -2 * (1 - cos(u[2]) - u[2] * sin(u[2]))
    A[2, 2] = 0
end
A = DiffEqArrayOperator(ones(2, 2), update_func = update_func)
prob = ODEProblem(A, ones(2), (30, 150.0))
sol = solve(prob, MagnusAdapt4())


# CayleyEuler - First order method using Cayley transformations.  (Doesn't work because it tries to take the inverse of a sparse matrix)
# LieEuler - First order Lie Euler method.                          (Seems incredibly slow)
# RKMK2 - Second order Runge–Kutta–Munthe-Kaas method.
# RKMK4 - Fourth order Runge–Kutta–Munthe-Kaas method.
# LieRK4 - Fourth order Lie Runge-Kutta method.
# CG2 - Second order Crouch–Grossman method.
# CG4a - Fourth order Crouch-Grossman method.
# MagnusAdapt4 



#%%

using SciMLOperators
using LinearAlgebra, FFTW

n = 256
L = 2π

dx = L / n
x = range(start = -L / 2, stop = L / 2 - dx, length = n) |> Array
u = @. sin(5x)cos(7x);
du = @. 5cos(5x)cos(7x) - 7sin(5x)sin(7x);

k = rfftfreq(n, 2π * n / L) |> Array
m = length(k)
P = plan_rfft(x)

fwd(u, p, t) = P * u
bwd(u, p, t) = P \ u

fwd(du, u, p, t) = mul!(du, P, u)
bwd(du, u, p, t) = ldiv!(du, P, u)

F = FunctionOperator(fwd, x, im * k;
    T = ComplexF64, op_adjoint = bwd,
    op_inverse = bwd,
    op_adjoint_inverse = fwd, islinear = true
)

ik = im * DiagonalOperator(k)
Dx = F \ ik * F

Dx = cache_operator(Dx, x)

@show ≈(Dx * u, du; atol = 1e-8)
@show ≈(mul!(copy(u), Dx, u), du; atol = 1e-8)




using OrdinaryDiffEq, Test, DiffEqDevTools
using LinearAlgebra, Random

# Linear exponential solvers
A = MatrixOperator([2.0 -1.0; -1.0 2.0])
u0 = ones(2)
prob = ODEProblem(A, u0, (0.0, 1.0))
solve(prob, LinearExponential(krylov = :off))

sol1 = solve(prob, LinearExponential(krylov = :off))(1.0)
sol2 = solve(prob, LinearExponential(krylov = :simple))(1.0)
sol3 = solve(prob, LinearExponential(krylov = :adaptive))(1.0)
sol4 = solve(prob, Rosenbrock23(), reltol = 1e-12, abstol = 1e-12)(1.0)
sol_analytic = exp(1.0 * Matrix(A)) * u0

@test isapprox(sol1, sol_analytic, rtol = 1e-10)
@test isapprox(sol2, sol_analytic, rtol = 1e-10)
@test isapprox(sol3, sol_analytic, rtol = 1e-10)
@test isapprox(sol4, sol_analytic, rtol = 1e-8)



# u' = A(t)u solvers
function update_func!(A, u, p, t)
    A[1, 1] = 0
    A[2, 1] = sin(u[1])
    A[1, 2] = -1
    A[2, 2] = 0
end
A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (10, 50.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern9(), dt = 1 / 4)
sol2 = solve(prob, OrdinaryDiffEq.RKMK2(), dt = 1 / 4)
dts = 1 ./ 2 .^ (10:-1:5)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, RKMK2(), test_setup)
@test sim.𝒪est[:l2]≈2 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern9(), dt = 1 / 4)
sol2 = solve(prob, OrdinaryDiffEq.RKMK4(), dt = 1 / 4)
dts = (0.38) .^ (6:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, RKMK4(), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.22

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern9(), dt = 1 / 4)
sol2 = solve(prob, OrdinaryDiffEq.LieRK4(), dt = 1 / 4)
dts = 1 ./ 2 .^ (7:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, LieRK4(), test_setup)
@test sim.𝒪est[:l2]≈5 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern9(), dt = 1 / 4)
sol2 = solve(prob, OrdinaryDiffEq.CG2(), dt = 1 / 4)
dts = 1 ./ 2 .^ (7:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, CG2(), test_setup)
@test sim.𝒪est[:l2]≈2 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 20.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern6(), dt = 1 / 8)
sol2 = solve(prob, OrdinaryDiffEq.CG3(), dt = 1 / 8)
dts = 1 ./ 2 .^ (10:-1:3)
test_setup = Dict(:alg => Vern6(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, CG3(), test_setup)
@test sim.𝒪est[:l2]≈3 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol1 = solve(prob, Vern9(), dt = 1 / 4)
sol2 = solve(prob, CG4a(), dt = 1 / 4)
dts = (0.38) .^ (6:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, CG4a(), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.28

function update_func!(A, u, p, t)
    A[1, 1] = 0
    A[2, 1] = 1
    A[1, 2] = -2 * (1 - cos(u[2]) - u[2] * sin(u[2]))
    A[2, 2] = 0
end
A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (30, 150.0))
dts = 1 ./ 2 .^ (7:-1:1)
test_setup = Dict(:alg => Tsit5(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusAdapt4(), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.2