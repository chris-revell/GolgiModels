
using Symbolics

nMax = 100

@variables C[1:nMax] M[1:nMax] k[1:8]

exp1 = k[1] + k[5]M[1] ~ k[4]*C[1] 
exp2 = k[4]*C[1] ~ k[8]M[1] + k[5]M[1]

C₁Exp = simplify(Symbolics.solve_for( substitute(exp1,Dict(M[1]=> Symbolics.solve_for(exp2,M[1]))) ,C[1]))

M₁Exp = simplify(Symbolics.solve_for( substitute(exp2,Dict(C[1]=> Symbolics.solve_for(exp1,C[1]))) ,M[1]))



Symbolics.Valu
# exp1 = C[1] ~ (k[1]/k[12])
# exp2 = M[1] ~ (k[1]/k[12])

##

exp3 = k[2]C[]