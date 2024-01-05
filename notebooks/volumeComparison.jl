using XLSX
using DataFrames
using CairoMakie
using StatsBase
using GeometryBasics

##

measuredVolumes = XLSX.readxlsx(datadir("exp_raw","Nikki golgi data","20231004 Analysis for model.xlsx"))[5]

transWTvolumes = measuredVolumes["I4:I1016"]
trans1G12volumes = measuredVolumes["J4:J999"]
trans1G6volumes = measuredVolumes["K4:K896"]
trans2H7volumes = measuredVolumes["L4:L993"]
trans2B8volumes = measuredVolumes["M4:M399"]

medialWTvolumes = measuredVolumes["A4:A2305"]
medial1G12volumes = measuredVolumes["B4:B176"]
medial1G6volumes = measuredVolumes["C4:C1543"]
medial2H7volumes = measuredVolumes["D4:D1691"]
medial2B8volumes = measuredVolumes["E4:E664"]

##

measurements = dropmissing(DataFrame(XLSX.readtable(datadir("exp_raw","Nikki golgi data","AllNumbers2.xlsx"),1)))
measurements[!, [:Maturation,:CellType]] = convert.(String, measurements[!, [:Maturation,:CellType]])
headers = Symbol.(names(measurements))
floatHeaders = [h for h in headers if h∉[:Maturation,:CellType]]
# floatHeaders = [:MinFeret,:Area]
measurements[!, floatHeaders] = convert.(Float64, measurements[!, floatHeaders])

##

# fig = Figure(size=(2000,3000))

maxVolInd = findmax(measurements[!,:Area])[2]
maxVol = (4.0/3.0)*measurements[maxVolInd,:Area]*measurements[maxVolInd,:Feret]/2.0
step = maxVol*1.1/1000
histbins = 0.0:step:maxVol*1.1
colors = (Cis=:red,Medial=:green,Trans=:blue)
axes = Axis[]

# typicalVesicleVolume = (4π/3)*0.03^3

##

fig = Figure(size=(1000,2000))
ax1 = Axis(fig[1,1],aspect=AxisAspect(1.0),yscale=Makie.pseudolog10)
ax2 = Axis(fig[2,1],aspect=AxisAspect(1.0),yscale=Makie.pseudolog10)
filteredData = filter([:Maturation, :CellType] => (m, c) -> m == "Trans" && c == "WT", measurements, view=true)
a = filteredData[!,:Feret]./2.0
b = filteredData[!,:Area]./(π.*a)
h = ((a.-b).^2)./((a.+b).^2)
ellipseCircumferences = π.*(a.+b).*(1.0.+3.0.*h./(10.0.+sqrt.(4.0.-3.0.*h)))
volumes = (4.0/3.0)*π.*(a.^2).*b
# hist = fit(Histogram, volumes, nbins = 100)#, histbins)
# barplot!(ax, hist, label = "Calculated")
density!(ax1,volumes,color=(:blue,0.5))
Label(fig[1,1,Top()], "Calculated")
density!(ax2,Float64.(transWTvolumes[:,1]),color=(:red,0.5))
Label(fig[2,1,Top()], "Measured")
# hist = fit(Histogram, Float64.(transWTvolumes[:,1]), nbins = 100)#, histbins)
# barplot!(ax, hist, label = "Calculated")

linkyaxes!(ax1,ax2)
linkxaxes!(ax1,ax2)

display(fig)

##


# for (ic,cellType) in enumerate(unique(measurements[!,:CellType]))#["WT", "1G6"]
#     for (im,maturity) in enumerate(["Cis","Medial","Trans"])
#         filteredData = filter([:Maturation, :CellType] => (m, c) -> m == maturity && c == cellType, measurements, view=true)

#         ax = Axis(fig[ic,im],xscale = Makie.pseudolog10)
#         ax.xticks = ([1.0,10.0,100.0],string.([1.0,10.0,100.0]))
#         ax.ylabel = "Volume"
#         ax.xlabel = "Frequency"

#         a = filteredData[!,:Feret]./2.0
#         b = filteredData[!,:Area]./(π.*a)
#         h = ((a.-b).^2)./((a.+b).^2)
#         ellipseCircumferences = π.*(a.+b).*(1.0.+3.0.*h./(10.0.+sqrt.(4.0.-3.0.*h)))
#         volumes = (4.0/3.0)*π.*(a.^2).*b

#         @show median(volumes)/typicalVesicleVolume

#         hist = fit(Histogram, volumes, histbins)
#         barplot!(ax, hist, label = "$cellType, $maturity, $(length(a))",color=colors[Symbol(maturity)],direction=:x)
#         barLocations = collect(hist.edges[1]).+Float64(hist.edges[1].step)/2
#         barLabels = string.(collect(1:hist.edges[1].len))
#         ax.yticks = (barLocations[1:10:end],barLabels[1:10:end])
#         Label(fig[ic,im,Top()],"$cellType, $maturity, ($(length(a)) compartments)",fontsize=32)

#         push!(axes,ax)
    
#     end
# end

# linkyaxes!(axes...)
# linkxaxes!(axes...)
# # display(fig)

# ##

# save("allhistograms.png",fig)



##
