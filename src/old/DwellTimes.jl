#
#  DwellTimes.jl
#  GolgiModels
#
#  Created by Christopher Revell on 22/05/2023.

# module DwellTimes

using Symbolics
using LinearAlgebra
using SparseArrays

@variables k₁ k̂₂ k̂₃ k₄ k₅ k̂₆ k̂₇ k₈ k₉ k̂₁₀ k̂₁₁ k₁₂

K₇ = sparse(Diagonal([-k₁, -k̂₂-k₄, -k̂₃, -k₅-k̂₆-k₈, -k̂₇, -k₉-k̂₁₀-k₁₂, -k̂₁₁]))
K₇[2,1] = k₁
K₇[2,3] = k̂₃
K₇[2,4] = k₅
K₇[3,2] = k̂₂

K₇[4,2] = k₄
K₇[4,5] = k̂₇
K₇[4,6] = k₉
K₇[5,4] = k̂₆

K₇[6,4] = k₈
K₇[6,7] = k̂₁₁
K₇[7,6] = k̂₁₀

K₇⁻¹ = simplify.(inv(Matrix(K₇)))


dwellTimesExpr = build_function.(K₇⁻¹[:,1], fill([k₁, k̂₂, k̂₃, k₄, k₅, k̂₆, k̂₇, k₈, k₉, k̂₁₀, k̂₁₁, k₁₂],7))
# Base.remove_linenums!.(f_expr)
dwellTimesFuns = eval.(dwellTimesExpr)

# export dwellTimesFuns

# end