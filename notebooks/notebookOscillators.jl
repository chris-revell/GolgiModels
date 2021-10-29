### A Pluto.jl notebook ###
# v0.16.3

using Markdown
using InteractiveUtils

# ╔═╡ 6446adac-a3e2-4961-bea5-620118149b50
begin
    import Pkg
    Pkg.add(PackageSpec(url="https://github.com/Pocket-titan/DarkMode"))
    import DarkMode
    DarkMode.enable()
    # OR DarkMode.Toolbox(theme="default")
end

# ╔═╡ 81bfcb8c-3345-11ec-21ec-bf37d08f624f
begin
	using DifferentialEquations
	using Plots	
	using PlutoUI
	#using DarkMode; DarkMode.Toolbox(theme="default")
end;

# ╔═╡ 369a25d1-6c20-4b22-9f66-a845364b3a0d
function model!(du, u, p, t)
	ω0, ωF, ωM, ϕF, ϕM, μ, ν, αF, αM, ψF, ψM = p
	du[1] = ω0
	du[2] = ωF + μ*sin(u[3]-u[2]-ϕF) + αF*sin(u[1]-u[2]-ψF)
	du[3] = ωM + ν*sin(u[2]-u[3]-ϕM) + αM*sin(u[1]-u[3]-ψM)
end;

# ╔═╡ ecf11f46-b4fa-46c5-8c10-b833adcbd72b
begin
	θ  = rand(3).%(2π) # Initial phases (Sunlight, fibroblast, macrophage)
	tMax = 10.0 	   # Total time
	ω0 = 2π            # Sunlight phase rate of change /day
	ωF = 2π            # Fibroblast phase intrinsic rate of change /day
	ωM = 2π            # Macrophage phase intrinsic rate of change /day
	ϕF = 3π/5          # Fibroblast phase offset in coupling to macrophage
	ϕM = 2π/5 		   # Macrophage phase offset in coupling to fibroblast
	μ  = 1.0           # Amplitude of fibroblast phase coupling to macrophage
	ν  = 0.0           # Amplitude of macrophage phase coupling to fibroblast
	αF = 10.0          # Amplitude of fibroblast phase coupling to sunlight
	αM = 10.0          # Amplitude of macrophage phase coupling to sunlight
	ψF = π/8           # Fibroblast phase offset from sunlight
	ψM = π/4          # Macrophage phase offset from sunlight
	p = [ω0, ωF, ωM, ϕF, ϕM, μ, ν, αF, αM, ψF, ψM]
	
	prob = ODEProblem(model!,θ,(0.0,tMax),p)
	sol = solve(prob, Tsit5(),saveat=0.01)

	anim = @animate for u in sol.u
		scatter(sin.(u),cos.(u),xlims=(-1.1,1.1),ylims=(-1.1,1.1),aspect_ratio=:equal,marker=[:star8,:circle,:circle],color=[:red,:green,:blue],ms=10,series_annotations=["Light","Fib","Mac"],legend=:none)
	end
	gif(anim,"test.gif",fps=3)
	
end

# ╔═╡ Cell order:
# ╠═6446adac-a3e2-4961-bea5-620118149b50
# ╠═81bfcb8c-3345-11ec-21ec-bf37d08f624f
# ╠═369a25d1-6c20-4b22-9f66-a845364b3a0d
# ╠═ecf11f46-b4fa-46c5-8c10-b833adcbd72b
