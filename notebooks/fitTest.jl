using Distributions
using DataStructures
using GeometryBasics
using CairoMakie

nMax = 20

probDist = Exponential(2)

allSizes = ceil.(Int64,rand(probDist,1000))

c = counter(allSizes)

pts = [Point2(key,value) for (key,value) in collect(c)]

# fig = Figure(resolution=(500,500))
# ax = Axis(fig[1,1])
# scatter!(ax,pts)
# display(fig)


