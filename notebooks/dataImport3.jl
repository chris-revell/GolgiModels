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

fig_a = Figure(size=(1000,1000))
ax_a = Axis(fig_a[1,1])
ax_a.ylabel = "MinFeret"
ax_a.xlabel = "Frequency"
step = 0.06
histbins = 0.0:step:maximum(measurements[!,:MinFeret])/3.0
hist = fit(Histogram, measurements[!,:MinFeret], histbins)
barplot!(ax_a, hist)
vlines!(ax_a,[0.06],color=:black,linestyle=:dash)
display(fig_a)
save("MinFeretHist.png",fig_a)
##

fig_b = Figure(size=(1000,1000))
ax_b = Axis(fig_a[1,1])
ax_b.ylabel = "Area"
ax_b.xlabel = "Frequency"
step = π*0.06^2
histbins = 0.0:step:maximum(measurements[!,:Area])
hist = fit(Histogram, measurements[!,:Area], histbins)
barplot!(ax_b, hist)
vlines!(ax_b,[π*0.06^2],color=:black,linestyle=:dash)
display(fig_b)
save("AreaHist.png",fig_b)

##
# Delete row with outlier values
# maxVolInd = findmax(measurements[!,:Area])[2]
# delete!(measurements, maxVolInd)

##

fig = Figure(size=(2000,3000))

maxVolInd = findmax(measurements[!,:Area])[2]
maxVol = (4.0/3.0)*measurements[maxVolInd,:Area]*measurements[maxVolInd,:Feret]/2.0
step = maxVol*1.1/100
histbins = 0.0:step:maxVol*1.1
colors = (Cis=:red,Medial=:green,Trans=:blue)
axes = Axis[]

typicalVesicleVolume = (4π/3)*0.03^3

for (ic,cellType) in enumerate(unique(measurements[!,:CellType]))#["WT", "1G6"]
    for (im,maturity) in enumerate(["Cis","Medial","Trans"])
        filteredData = filter([:Maturation, :CellType] => (m, c) -> m == maturity && c == cellType, measurements, view=true)

        ax = Axis(fig[ic,im],xscale = Makie.pseudolog10)
        ax.xticks = ([1.0,10.0,100.0],string.([1.0,10.0,100.0]))
        ax.ylabel = "Volume"
        ax.xlabel = "Frequency"

        a = filteredData[!,:Feret]./2.0
        b = filteredData[!,:Area]./(π.*a)
        h = ((a.-b).^2)./((a.+b).^2)
        ellipseCircumferences = π.*(a.+b).*(1.0.+3.0.*h./(10.0.+sqrt.(4.0.-3.0.*h)))
        volumes = (4.0/3.0)*π.*(a.^2).*b

        @show median(volumes)/typicalVesicleVolume

        hist = fit(Histogram, volumes, histbins)
        barplot!(ax, hist, label = "$cellType, $maturity, $(length(a))",color=colors[Symbol(maturity)],direction=:x)
        barLocations = collect(hist.edges[1]).+Float64(hist.edges[1].step)/2
        barLabels = string.(collect(1:hist.edges[1].len))
        ax.yticks = (barLocations[1:10:end],barLabels[1:10:end])
        Label(fig[ic,im,Top()],"$cellType, $maturity, ($(length(a)) compartments)",fontsize=32)

        push!(axes,ax)
    
    end
end

linkyaxes!(axes...)
linkxaxes!(axes...)
display(fig)

##

save("allhistograms.png",fig)



##
