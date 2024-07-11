using CairoMakie 
using DrWatson
using DataFrames
using XLSX

df = DataFrame(XLSX.readtable(datadir("exp_raw", "ModellingData.xlsx"), 2))

simpsonsRule(xs) = sum(xs[1:end-1].+xs[2:end])./2.0

wtVec = filter(x->x.KO=="WT", df)[!,"normalised intensity"]
wt = reverse(wtVec./simpsonsRule(wtVec))
GMAP210_1Vec = filter(x->x.KO=="GMAP210 KO1", df)[!,"normalised intensity"]
GMAP210_1 = reverse(GMAP210_1Vec./simpsonsRule(GMAP210_1Vec))
GMAP210_2Vec = filter(x->x.KO=="GMAP210 KO2", df)[!,"normalised intensity"]
GMAP210_2 = reverse(GMAP210_2Vec./simpsonsRule(GMAP210_2Vec))
Golgin160_2Vec = filter(x->x.KO=="Golgin160 KO2", df)[!,"normalised intensity"]
Golgin160_2 = reverse(Golgin160_2Vec./simpsonsRule(Golgin160_2Vec))

νs = collect(1:length(wt))

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, νs, wt, label="wt")
lines!(ax, νs, GMAP210_1, label="GMAP210_1")
lines!(ax, νs, GMAP210_2, label="GMAP210_2")
lines!(ax, νs, Golgin160_2, label="Golgin160_2")
ax.xlabel = "ν"
ax.ylabel = "Normalised intensity"
axislegend(ax)
save("intensityData.png", fig)
display(fig)