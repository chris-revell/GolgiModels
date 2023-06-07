using DrWatson
@quickactivate
using DynamicalSystems
using DifferentialEquations
using GLMakie
using DataStructures: CircularBuffer
using GeometryBasics

# Function defining ODEs for model
function model!(du, u, p, t)
	k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13 = p
	du[1] = k1[]*u[2] - k2[]*u[1] + 1.0
	du[2] = k2[]*u[1] + k4[]*u[4] - k1[]*u[2] - k6[]*u[2] - k5[]*u[2]
	du[3] = k6[]*u[2] + k7[]*u[4] - k11[]*u[3] - k10[]*u[3]
	du[4] = k3[]*u[1] + k5[]*u[2] + k8[]*u[5] - k4[]*u[4] - k7[]*u[4] -k9[]*u[4]
	du[5] = k9[]*u[4] + k10[]*u[3] - k8[]*u[5] - k12[]*u[5]
	du[6] = k11[]*u[3] - k13[]*u[6]
	du[7] = k12[]*u[5] - k13[]*u[7]
end

# Function to iterate system state
function progress_for_one_step!(integ)
    step!(integ)
    u = integ.u
    return u
end

# Function to update figure based on system iteration
function animStep!(integ, positions, circles)
	u = progress_for_one_step!(integ)
	circles[] = [
		Circle(positions[1],0.1*sqrt(u[1])),
		Circle(positions[2],0.1*sqrt(u[2])),
		Circle(positions[3],0.1*sqrt(u[3])),
		Circle(positions[4],0.1*sqrt(u[4])),
		Circle(positions[5],0.1*sqrt(u[5])),
		Circle(positions[6],0.1*sqrt(u[6])),
		Circle(positions[7],0.1*sqrt(u[7]))
	]
end

# Define positions of each organelle in diagram
nPoints = 7
positions = [
	Point2(0.0,0.0),
	Point2(0.0,0.5),
	Point2(0.0,1.0),
	Point2(0.5,0.5),
	Point2(0.5,1.0),
	Point2(0.0,1.5),
	Point2(0.5,1.5)
]
# Set up circles at defined positions with as observable list
circles = Observable([
	Circle(positions[1],0.01),
	Circle(positions[2],0.01),
	Circle(positions[3],0.01),
	Circle(positions[4],0.01),
	Circle(positions[5],0.01),
	Circle(positions[6],0.01),
	Circle(positions[7],0.01)
])
# Set up labels or circle
labels = [
	"ER",
	"Golgi",
	"Secretory vesicle",
	"Lysosome",
	"Endosome",
	"Monomer secretion",
	"Fibril formation"
]
# Set up lines between interacting organelles
lines = [
	[positions[1], positions[1].+(positions[2].-positions[1])./2.0 .+ Point2(0.05,0.0),positions[2]],
	[positions[1], positions[1].+(positions[2].-positions[1])./2.0 .- Point2(0.05,0.0),positions[2]],
	[positions[1],positions[4]],
	[positions[2],positions[3]],
	[positions[2],positions[2].+(positions[4].-positions[2])./2.0 .+ Point2(0.0,0.05),positions[4]],
	[positions[2],positions[2].+(positions[4].-positions[2])./2.0 .- Point2(0.0,0.05),positions[4]],
	[positions[4],positions[4].+(positions[5].-positions[4])./2.0 .+ Point2(0.05,0.0),positions[5]],
	[positions[4],positions[4].+(positions[5].-positions[4])./2.0 .- Point2(0.05,0.0),positions[5]],
	[positions[3],positions[4]],
	[positions[3],positions[5]],
	[positions[3],positions[6]],
	[positions[5],positions[7]]
]

# Set up figure canvas
fig = Figure(resolution=(1000,1100))
grid = fig[1,1] = GridLayout()
ax = Axis(grid[2,1],aspect=DataAspect())
xlims!(ax,(-0.3,0.9))
ylims!(ax,(-0.3,1.9))

hidedecorations!(ax)
for l in lines
	lines!(ax,l;color=:black)
end
text!(labels,position=positions,color=:black,align=(:left,:baseline),textsize=0.05,space=:data)
poly!(circles)

# Set up parameter sliders
lsgrid = labelslidergrid!(
    fig,
    ["k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", "k12", "k13"],
    [0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0];
    formats = x -> "$(round(x, digits = 1))", #[x -> "$(round(x, digits = 1))$s" for s in ["π", "π", "π", "π", "π", "", "", "", "", "π", "π"]],
    width = 350,
    tellheight = false
)
grid[2, 2] = lsgrid.layout
# Set default slider values
defaults = ones(13)
set_close_to!.(lsgrid.sliders, defaults)
# Pull parameters from slider positions [ω0, ωF, ωM, ϕF, ϕM, μ, ν, αF, αM, ψF, ψM]
kObservables = [s.value for s in lsgrid.sliders]


# Set up differential equation integrator to solve concentrations/circle radii
u0 = zeros(nPoints)
prob = ODEProblem(model!,u0,(0.0,10000.0),kObservables)
dp = ContinuousDynamicalSystem(prob)
# Solve diffeq with constant step for smoother curves
diffeq = (alg = Tsit5(), adaptive = false, dt = 0.005)
# Set up integrator for each iteration
integ = integrator(dp, u0; diffeq...)


run = Button(grid[1,1:2]; label = "Start/Stop")#, tellwidth = false)
isrunning = Observable(false)
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break # ensures computations stop if closed window
        animStep!(integ, positions, circles)
        sleep(0.02) # or `yield()` instead
    end
end

display(fig)
