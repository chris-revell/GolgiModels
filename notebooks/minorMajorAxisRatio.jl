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

measurements = dropmissing(DataFrame(XLSX.readtable(datadir("exp_raw","Nikki golgi data","AllNumbers2.xlsx"),1)))
measurements[!, [:Maturation,:CellType]] = convert.(String, measurements[!, [:Maturation,:CellType]])
headers = Symbol.(names(measurements))
floatHeaders = [h for h in headers if h∉[:Maturation,:CellType]]
# floatHeaders = [:MinFeret,:Area]
measurements[!, floatHeaders] = convert.(Float64, measurements[!, floatHeaders])

##

fig = Figure(resolution=(1000,2000),fontsize=36)
Label(fig[1,1,Top()],"Major/minor axis ratios from data")
ax = Axis(fig[1,1],aspect=DataAspect())
lines!(ax,[0,20],[0,20],linestyle=:dash,color=(:black,0.5))
ax.ylabel = "b"
ax.xlabel = "a"
ax2 = Axis(fig[2,1], aspect=AxisAspect(1.0))
ax2.xlabel = "a/b ratio"
ax2.ylabel = "Frequency density"

for (ic,cellType) in enumerate(unique(measurements[!,:CellType]))
    for (im,maturity) in enumerate(unique(measurements[!,:Maturation]))
        filteredData = filter([:Maturation, :CellType] => (m, c) -> m == maturity && c == cellType, measurements, view=true)

        a = filteredData[!,:Feret]./2.0
        b = filteredData[!,:Area]./(π.*a)
        
        scatter!(ax,a,b,color=(Makie.wong_colors()[im],0.25))
        
        # hist = fit(Histogram,a./b)
        density!(ax2,a./b,color=(Makie.wong_colors()[im],0.25))
    end
end

display(fig)

##

save("minorMajorAxisRatio.png",fig)



##
filteredData = filter([:Maturation, :CellType] => (m, c) -> m == "Trans" && c == "WT", measurements, view=true)
a = filteredData[!,:Feret]./2.0
b = filteredData[!,:Area]./(π.*a)
rTransWT = mean(b./a)

##
filteredData = filter([:Maturation, :CellType] => (m, c) -> m == "Medial" && c == "WT", measurements, view=true)
a = filteredData[!,:Feret]./2.0
b = filteredData[!,:Area]./(π.*a)
rMedialWT = mean(b./a)


##
filteredData = filter([:Maturation, :CellType] => (m, c) -> m == "Cis" && c == "WT", measurements, view=true)
a = filteredData[!,:Feret]./2.0
b = filteredData[!,:Area]./(π.*a)
rCisWT = mean(b./a)


# all_a = measurements[!,:Feret]./2.0
# all_b = measurements[!,:Area]./(π.*all_a)

# mean(all_a./all_b)