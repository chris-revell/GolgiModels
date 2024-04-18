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

measurementsMed = XLSX.readxlsx(datadir("exp_raw","Nikki golgi data","20231004 Analysis for model.xlsx"))[4]
areasWTMedStr = Vector(measurementsMed["D4:D586"][:,1])
deleteat!(areasWTMedStr, findall(ismissing,areasWTMedStr))
areasWTMed = Float64.(areasWTMedStr)
perimetersWTMedStr = Vector(measurementsMed["E4:E586"][:,1])
deleteat!(perimetersWTMedStr, findall(ismissing,perimetersWTMedStr))
perimetersWTMed = Float64.(perimetersWTMedStr)
# Take feret length to be 2*a (width of ellipse)
feretsWTMedStr = Vector(measurementsMed["G4:G586"][:,1])
deleteat!(feretsWTMedStr, findall(ismissing,feretsWTMedStr))
feretsWTMed = Float64.(feretsWTMedStr)
circularitiesWTMed = 4π.*areasWTMed./perimetersWTMed.^2
aMed = feretsWTMed./2.0
bMed = areasWTMed./(π.*aMed)
hMed = ((aMed.-bMed).^2)./((aMed.+bMed).^2)
ellipseCircumferencesMed = π.*(aMed.+bMed).*(1.0.+3.0.*hMed./(10.0.+sqrt.(4.0.-3.0.*hMed)))
volumesMed = (4.0/3.0).*(aMed.^2).*bMed

##


measurementsMed = XLSX.readxlsx(datadir("exp_raw","Nikki golgi data","20231004 Analysis for model.xlsx"))[3]
areasWTMedStr = Vector(measurementsMed["D4:D586"][:,1])
deleteat!(areasWTMedStr, findall(ismissing,areasWTMedStr))
areasWTMed = Float64.(areasWTMedStr)
perimetersWTMedStr = Vector(measurementsMed["E4:E586"][:,1])
deleteat!(perimetersWTMedStr, findall(ismissing,perimetersWTMedStr))
perimetersWTMed = Float64.(perimetersWTMedStr)
# Take feret length to be 2*a (width of ellipse)
feretsWTMedStr = Vector(measurementsMed["G4:G586"][:,1])
deleteat!(feretsWTMedStr, findall(ismissing,feretsWTMedStr))
feretsWTMed = Float64.(feretsWTMedStr)
circularitiesWTMed = 4π.*areasWTMed./perimetersWTMed.^2
aMed = feretsWTMed./2.0
bMed = areasWTMed./(π.*aMed)
hMed = ((aMed.-bMed).^2)./((aMed.+bMed).^2)
ellipseCircumferencesMed = π.*(aMed.+bMed).*(1.0.+3.0.*hMed./(10.0.+sqrt.(4.0.-3.0.*hMed)))
volumesMed = (4.0/3.0).*(aMed.^2).*bMed


measurementsTrans = XLSX.readxlsx(datadir("exp_raw","Nikki golgi data","20231004 Analysis for model.xlsx"))[3]
areasWTTrans = Float64.(Vector(measurementsTrans["D4:D1328"][:,1]))
perimetersWTTrans = Float64.(Vector(measurementsTrans["E4:E1328"][:,1]))
# Take feret length to be 2*a (width of ellipse)
feretsWTTrans = Float64.(Vector(measurementsTrans["G4:G1328"][:,1]))
circularitiesWTTrans = 4π.*areasWTTrans./perimetersWTTrans.^2
aTrans = feretsWTTrans./2.0
bTrans = areasWTTrans./(π.*aTrans)
hTrans = ((aTrans.-bTrans).^2)./((aTrans.+bTrans).^2)
ellipseCircumferencesTrans = π.*(aTrans.+bTrans).*(1.0.+3.0.*hTrans./(10.0.+sqrt.(4.0.-3.0.*hTrans)))
volumesTrans = (4.0/3.0).*(aTrans.^2).*bTrans

##
figPerim = Figure()
axPerim = Axis(figPerim[1,1])
scatter!(axPerim,perimetersWTTrans,ellipseCircumferencesTrans, label="Trans")
scatter!(axPerim,perimetersWTMed,ellipseCircumferencesMed, label="Medial")
minPerim = minimum([ellipseCircumferencesTrans...,perimetersWTTrans...])
maxPerim = maximum([ellipseCircumferencesTrans...,perimetersWTTrans...])
lines!(axPerim,[Point2(minPerim,minPerim),Point2(maxPerim,maxPerim)],linestyle=:dash,color=:black)
axPerim.xlabel = "Measured compartment perimeters"
axPerim.ylabel = "Perimeters calculated from area and feret length,\nassuming compartment is an ellipse"
axislegend(axPerim)
display(figPerim)


##
figVol = Figure()
axVol = Axis(figVol[1,1],yscale = Makie.pseudolog10)
histTrans = fit(Histogram, volumesTrans, nbins=100)
histMed = fit(Histogram, volumesMed, histTrans.edges[1])
axVol.yticks = ([1.0,10.0,100.0],string.([1.0,10.0,100.0]))
axVol.xlabel = "Volume"
axVol.ylabel = "Frequency"
barplot!(axVol,histTrans, label="Trans ($(length(areasWTTrans)) compartments)", color=(:red,0.5))
barplot!(axVol,histMed, label="Medial ($(length(areasWTMed)) compartments)", color=(:green,0.5))
axislegend(axVol)
display(figVol)





##

# using Symbolics
# @variables A a b e h C r
# expr_eccentricity = e ~ sqrt(1-b^2/a^2)
# exprArea = A ~ π*a*b
# expr_h = h ~ (a-b)^2/(a+b)^2
# expr_perim = C ~ π*(a+b)*(1+3*h/(10+sqrt(4-3*h)))
# circularity = r ~ 4π*A/C^2
# Symbolics.solve_for([expr_eccentricity,exprArea,expr_h,expr_perim,circularity], [e])