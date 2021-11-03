using DifferentialEquations
using GLMakie

nPoints = 7;

fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
xlims!(ax,(-0.5,1.1))
ylims!(ax,(-0.5,2.1))
hidedecorations!(ax)


positions = [
	Point2f0(0.0,0.0),
	Point2f0(0.0,0.5),
	Point2f0(0.0,1.0),
	Point2f0(0.5,0.5),
	Point2f0(0.5,1.0),
	Point2f0(0.0,1.5),
	Point2f0(0.5,1.5)
]

markerSizes = 40.0.*ones(nPoints);

labels = [
	"ER",
	"Golgi",
	"Secretory vesicle",
	"Lysosome",
	"Endosome",
	"Monomer secretion",
	"Fibril formation"
]

lines = [
	[positions[1], positions[1].+(positions[2].-positions[1])./2.0 .+ Point2f0(0.05,0.0),positions[2]],
	[positions[1], positions[1].+(positions[2].-positions[1])./2.0 .- Point2f0(0.05,0.0),positions[2]],
	[positions[1],positions[4]],
	[positions[2],positions[3]],
	[positions[2],positions[2].+(positions[4].-positions[2])./2.0 .+ Point2f0(0.0,0.05),positions[4]],
	[positions[2],positions[2].+(positions[4].-positions[2])./2.0 .- Point2f0(0.0,0.05),positions[4]],
	[positions[4],positions[4].+(positions[5].-positions[4])./2.0 .+ Point2f0(0.05,0.0),positions[5]],
	[positions[4],positions[4].+(positions[5].-positions[4])./2.0 .- Point2f0(0.05,0.0),positions[5]],
	[positions[3],positions[4]],
	[positions[3],positions[5]],
	[positions[3],positions[6]],
	[positions[5],positions[7]]
]

for l in lines
	lines!(ax,l;color=:black)
end

scatter!(ax,positions;markersize=markerSizes);

text!(
    labels,
    position = positions,
    color = :black,
    align = (:left, :baseline),
    textsize = 0.05,
    space = :data
)

lsgrid = labelslidergrid!(
    fig,
    ["k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", "k12", "k13"],
    [0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0, 0:0.01:1.0];
    formats = x -> "$(round(x, digits = 1))", #[x -> "$(round(x, digits = 1))$s" for s in ["π", "π", "π", "π", "π", "", "", "", "", "π", "π"]],
    width = 350,
    tellheight = false
)
fig[1, 2] = lsgrid.layout
# Set default slider values
defaults = ones(13)
set_close_to!.(lsgrid.sliders, defaults)


display(fig)
