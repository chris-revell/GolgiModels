using XLSX
using DataFrames
using CairoMakie
using StatsBase
using GeometryBasics

##

# Volume of ellipsoid (4/3)πa²b
# Area of ellipse πab 
# a = (feret length)/2
# b = 2*area / (π*feret length)

##

measurements = DataFrame(XLSX.readtable(datadir("exp_raw","Nikki golgi data","AllNumbers.xlsx"),1))
measurements[!, [:Maturation,:CellType]] = convert.(String, measurements[!, [:Maturation,:CellType]])
headers = Symbol.(names(measurements))
floatHeaders = [h for h in headers if h∉[:Maturation,:CellType]]
measurements[!, floatHeaders] = convert.(Float64, measurements[!, 3:end])

##

# Delete row with outlier values
maxVolInd = findmax(measurements[!,:Area])[2]
delete!(measurements, maxVolInd)

##

fig = Figure(size=(1000,2000))
axPerim = Axis(fig[1,1], aspect=DataAspect())
axPerim.xlabel = "Measured compartment perimeters"
axPerim.ylabel = "Perimeters calculated from area and feret length,\nassuming compartment is an ellipse"
axVol = Axis(fig[2,1],yscale = Makie.pseudolog10)
axVol.yticks = ([1.0,10.0,100.0],string.([1.0,10.0,100.0]))
axVol.xlabel = "Volume"
axVol.ylabel = "Frequency"

maxPerim = maximum(measurements[!,:Perimeter])
xlims!(axPerim,(0.0,maxPerim)); ylims!(axPerim,(0.0,maxPerim))
lines!(axPerim,[0.0,maxPerim],[0.0,maxPerim],color=:black,linestyle=:dash)
maxVolInd = findmax(measurements[!,:Area])[2]
maxVol = (4.0/3.0)*measurements[maxVolInd,:Area]*measurements[maxVolInd,:Feret]/2.0
histbins = 0.0:10.0:maxVol*1.1


for cellType in ["WT"]#unique(measurements[!,:CellType])#["WT", "1G6"]
    for maturity in ["Trans","Medial"]
        filteredData = filter([:Maturation, :CellType] => (m, c) -> m == maturity && c == cellType, measurements, view=true)

        a = filteredData[!,:Feret]./2.0
        b = filteredData[!,:Area]./(π.*a)
        h = ((a.-b).^2)./((a.+b).^2)
        ellipseCircumferences = π.*(a.+b).*(1.0.+3.0.*h./(10.0.+sqrt.(4.0.-3.0.*h)))
        volumes = (4.0/3.0)*π.*(a.^2).*b

        scatter!(axPerim,filteredData[!,:Perimeter],ellipseCircumferences, label="$cellType, $maturity, $(length(a))")

        hist = fit(Histogram, volumes, histbins)
        barplot!(axVol, hist, label = "$cellType, $maturity, $(length(a))")        
    end
end

axislegend(axPerim)
axislegend(axVol)
display(fig)



##
