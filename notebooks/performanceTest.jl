using LightGraphs
using DifferentialEquations
using Plots
using Distributions
using SparseArrays
using LinearAlgebra

@info "Initializing system..."
degree = 5
kdist = Poisson(degree)
N = Int(100)
ks = 1
while sum(ks) % 2 != 0
    global ks = rand(kdist, N)
end

g = random_configuration_model(N, ks)

i₀ = rand(1:N, round(Int, N*0.1)) 
x₀ = [ones(N); zeros(N); zeros(N)]
x₀[i₀] .= 0
x₀[i₀ .+ N] .= 1
tspan = (0.0, 100.0)
p = [0.25, 0.05, degree, adjacency_matrix(g)]

function agg_sol(sol, N)
    agg = reshape(sol[:, :], (N, 3, :))
    agg = mean(agg; dims=1)
    agg = dropdims(agg; dims=1)
    agg = permutedims(agg, (2, 1))
    return agg
end


function dNᵢ(i)

    function rate(u, p, t)
        u[i]*(p[1]/p[3])*dot(p[4][i, :], u[(N+1):2N])
    end

    function affect!(integrator)
        integrator.u[i] -= 1
        integrator.u[N+i] += 1
    end

    return rate, affect!

end

function dMᵢ(i)

    function rate(u, p, t)
        u[N+i]*p[2]
    end

    function affect!(integrator)
        integrator.u[N+i] -= 1

        integrator.u[2N+i] += 1
    end

    return rate, affect!

end

@info "Collecting jumps..."
jumps = ConstantRateJump[]
for i in 1:N
    push!(jumps, ConstantRateJump(dNᵢ(i)...))
end

for i in 1:N
    push!(jumps, ConstantRateJump(dMᵢ(i)...))
end

@info "Building dependency graph..."
vtoj = Vector{Vector{Int64}}(undef, 3N)
for i in 1:N
    vtoj[i] = [i]
    vtoj[N+i] = [findnz(p[4][i, :])[1]; N+i]
    vtoj[2N+i] = []
end
vtoj

jtov = Vector{Vector{Int64}}(undef, length(jumps))
for i in 1:N
    jtov[i] = [i, N+i]
    jtov[N+i] = [N+i, 2N+i]
end
jtov

@info "Building discrete problem..."
discrete_prob = DiscreteProblem(x₀, tspan, p)

@info "Building jump problem..."
# see https://github.com/SciML/ModelingToolkit.jl/blob/f12f472f630fd85a6fab4ca547b1c679217c33df/src/systems/jumps/jumpsystem.jl
# see https://diffeq.sciml.ai/stable/types/jump_types/
jump_prob = JumpProblem(discrete_prob, RSSA(), JumpSet(jumps...);
    vartojumps_map=vtoj, jumptovars_map=jtov)

@info "Solving problem..."
jump_sol = solve(jump_prob)

@info "Plotting..."
plot(agg_sol(jump_sol, N));
ylims!(0, 1);
plot!(legend=:right);
savefig("test3.png")

@info "Done."