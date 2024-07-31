#
#  GolgiModels.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/23.
#
#


module GolgiModels

using DrWatson
using PrecompileTools
using FromFile

@from "$(srcdir("GolgiApp.jl"))" using GolgiApp
@from "$(srcdir("GolgiSolve.jl"))" using GolgiSolve

# @compile_workload begin
#     golgiApp(displayFlag=false)
# end

export golgiApp
export golgiSolve
# 
end