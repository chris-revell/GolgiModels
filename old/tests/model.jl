using DifferentialEquations


# Function defining ODEs for model
function model!(du, u, p, t)
	k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13 = p
	du[1] = k1[]*u[2] - k2[]*u[1]
	du[2] = k2[]*u[1] + k4[]*u[4] - k1[]*u[2] - k6[]*u[2] - k5[]*u[2]
	du[3] = k6[]*u[2] + k7[]*u[4] - k11[]*u[3] - k10[]*u[3]
	du[4] = k3[]*u[1] + k5[]*u[2] + k8[]*u[5] - k4[]*u[4] - k7[]*u[4] -k9[]*u[4]
	du[5] = k9[]*u[4] + k10[]*u[3] - k8[]*u[5] - k12[]*u[5]
	du[6] = k11[]*u[3] - k13[]*u[6]
	du[7] = k12[]*u[5] - k13[]*u[7]
end


p = rand(13)

u0 = zeros(7)
