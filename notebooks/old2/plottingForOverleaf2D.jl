using JLD2
using CairoMakie

fig2 = Figure(size=(3250,2000), fontsize = 32)
ax1 = Axis(fig2[1,1], aspect=1, xlabel = "h", ylabel = "x")
xlims!(ax1, (0.0,2.0))
ylims!(ax1, (minimum(xs),maximum(xs)))
lines!(ax1, hFun.(xs[2:end-1], μxh, σxh), xs[2:end-1], linewidth=4)

ax2 = Axis(fig2[1,2], aspect=1, xlabel = "ν", ylabel = "x")
ax2.title = @sprintf "t = %.2f" 0.0
heatmap!(ax2, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], rotr90(reshape(sol.u[1],(Nxplus,Nνplus)))[1+Nghost:Nxplus-Nghost,1+Nghost:Nνplus-Nghost], colormap=:batlow)
ax3 = Axis(fig2[1,3], aspect=1, xlabel = "ν", ylabel = "x")
ax3.title = @sprintf "t = %.2f" tMax/2.0
heatmap!(ax3, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], rotr90(reshape(sol.u[30],(Nxplus,Nνplus)))[1+Nghost:Nxplus-Nghost,1+Nghost:Nνplus-Nghost], colormap=:batlow)
ax4 = Axis(fig2[1,4], aspect=1, xlabel = "ν", ylabel = "x")
ax4.title = @sprintf "t = %.2f" tMax
heatmap!(ax4, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], rotr90(reshape(sol.u[60],(Nxplus,Nνplus)))[1+Nghost:Nxplus-Nghost,1+Nghost:Nνplus-Nghost], colormap=:batlow)

linkyaxes!(ax1, ax2, ax3, ax4)

ax5 = Axis(fig2[2,2], aspect=1, xlabel = "ν", ylabel = "M(t;ν)")
integ = zeros(Float64, Nνplus)
u = reshape(sol.u[1],(Nxplus,Nνplus))
for j=1:Nνplus
    integ[j] = (sum(u[2:end, j])+sum(u[1:end-1, j]))*dx/2.0
end
lines!(ax5, νs[Nghost+1:end-Nghost], integ[Nghost+1:end-Nghost])
xlims!(ax5, (minimum(νs),maximum(νs)))

linkxaxes!(ax2, ax5)

ax6 = Axis(fig2[2,3], aspect=1, xlabel = "ν", ylabel = "M(t;ν)")
integ = zeros(Float64, Nνplus)
u = reshape(sol.u[30],(Nxplus,Nνplus))
for j=1:Nνplus
    integ[j] = (sum(u[2:end, j])+sum(u[1:end-1, j]))*dx/2.0
end
lines!(ax6, νs[Nghost+1:end-Nghost], integ[Nghost+1:end-Nghost])
xlims!(ax6, (minimum(νs),maximum(νs)))

linkxaxes!(ax3, ax6)

ax7 = Axis(fig2[2,4], aspect=1, xlabel = "ν", ylabel = "M(t;ν)")
integ = zeros(Float64, Nνplus)
u = reshape(sol.u[60],(Nxplus,Nνplus))
for j=1:Nνplus
    integ[j] = (sum(u[2:end, j])+sum(u[1:end-1, j]))*dx/2.0
end
lines!(ax7, νs[Nghost+1:end-Nghost], integ[Nghost+1:end-Nghost])
xlims!(ax7, (minimum(νs),maximum(νs)))

linkxaxes!(ax4, ax7)

linkxaxes!(ax5, ax6, ax7)

display(fig2)


save("glycosylation2D_hVariation.png", fig2)