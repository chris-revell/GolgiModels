using DrWatson
using FromFile
@quickactivate "ClockOscillators"

include("scripts/testParameters.jl")
@from "src/GolgiModel.jl" using GolgiModel
