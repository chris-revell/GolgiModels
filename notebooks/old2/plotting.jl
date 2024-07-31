using JLD2
using CairoMakie

# JLD2.@load "sol.jld2" sol


fig = Figure(size=(1000,1000))
ax = Axis3(fig[1, 1], aspect=:equal, azimuth=2.275π)
ax.xlabel = "x"
ax.ylabel = "ν"
ax.zlabel = "c"
uInternal3D = reshape(sol.u[1][ghostMask], (Nx, Ny, Nν))
globalmin = minimum([minimum(u[ghostMask]) for u in sol.u])
globalmax = maximum([maximum(u[ghostMask]) for u in sol.u])
uInternal = Observable(zeros(Nx,Nν))

zlims!(ax, (globalmin, globalmax))
clims = (globalmin,globalmax)

surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
record(fig, datadir("surface.mp4"), 1:length(sol.t); framerate=10) do i
    uInternal[] .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))[:,Nyplus÷2,:]
    uInternal[] = uInternal[]
end


# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "M, ∱cdxdy"
# uInternal3D = reshape((W*sol.u[end])[ghostMask], (Nx, Ny, Nν))
# M = sum(uInternal3D, dims=1)
# M = sum(M, dims=2)
# lines!(ax, νs[1:Nghost:end-2*Nghost], M[1,1,:])
# ax.title = "Integral of c over x and y against ν at final time"
# save(datadir("finalνVsM.png"), fig)

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "x"
# ax.ylabel = "y"
# uInternal3D = reshape(sol.u[1][ghostMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nx,Ny))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# clims = (globalmin,globalmax)
# clims = (minimum(uInternal3D), maximum(uInternal3D))
# heatmap!(ax, xs[Nghost+1:end-Nghost], ys[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
# record(fig, datadir("LargeNuTimeScan.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal3D .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))
#     uInternal[] .= uInternal3D[:,:,end]
#     uInternal[] = uInternal[]
#     ax.title = "xy profile of ν=1.0 at t=$(sol.t[i])"
# end

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "c"
# uInternal3D = reshape(sol.u[1][ghostMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nν))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# ylims!(ax, (globalmin, globalmax))
# lines!(ax, νs[Nghost+1:end-Nghost], uInternal)
# ax.title = "c against ν at x=0.5, y=0.5 at t=0.0"
# record(fig, datadir("NuProfileAtxyOverTime.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal3D .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))
#     uInternal[] .= uInternal3D[Nx÷2,Nx÷2,:]
#     uInternal[] = uInternal[]
#     ax.title = "c against ν at x=0.5, y=0.5 at t=$(sol.t[i])"
# end



# #%%


#%%
# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "x"
# ax.ylabel = "ν"
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# uInternal = Observable(zeros(Nx,Nν))
# clims = (globalmin,globalmax)
# ax.title = "c against x and ν at y=0.5 at t=0.0"
# heatmap!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
# record(fig, datadir("xνOverTimeAty.mp4"), 1:length(sol.t); framerate=10) do i
#     uInternal[] .= reshape(sol.u[i][ghostMask], (Nx, Ny, Nν))[:,Nyplus÷2,:]
#     uInternal[] = uInternal[]
#     ax.title = "c against x and ν at y=0.5 at t=$(sol.t[i])"
# end

#%%
# fig2 = Figure(size=(3250,2000), fontsize = 32)
# ax1 = Axis(fig2[1,1], aspect=1)
# ax1.xlabel = "h"
# ax1.ylabel = "x"
# heatmap!(ax1, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], zeros(Nplus,Nplus),colormap=:bwr)
# lines!(ax1, hFun.(xs[2:end-1], 0.5, 0.1), xs[2:end-1], linewidth=2)

# ax2 = Axis(fig2[1,2], aspect=1)
# ax2.xlabel = "ν"
# ax2.ylabel = "x"
# ax2.title = @sprintf "t = %.2f" 0.0
# # heatmap!(ax2, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
# heatmap!(ax2, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)
# ax3 = Axis(fig2[1,3], aspect=1)
# ax3.xlabel = "ν"
# ax3.ylabel = "x"
# ax3.title = @sprintf "t = %.2f" tMax/2.0
# # heatmap!(ax3, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
# heatmap!(ax3, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)
# ax4 = Axis(fig2[1,4], aspect=1)
# ax4.xlabel = "ν"
# ax4.ylabel = "x"
# ax4.title = @sprintf "t = %.2f" tMax
# # heatmap!(ax4, νs[Nghost+1:end-Nghost], xs[Nghost+1:end-Nghost], transpose(reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]), colormap=:batlow)
# heatmap!(ax4, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost], colormap=:batlow)

# linkyaxes!(ax1, ax2)

# ax5 = Axis(fig2[2,2], aspect=1)
# ax5.xlabel = "ν"
# ax5.ylabel = L"\int_x c(x,\nu)\partial x"
# ax5.title = @sprintf "t = %.2f" 0.0
# integ = zeros(Float64, N)
# u = reshape(sol.u[1],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
# for j=1:N
#     integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
# end
# lines!(ax5, νs[Nghost+1:end-Nghost], integ)
# ax6 = Axis(fig2[2,3], aspect=1)
# ax6.xlabel = "ν"
# ax6.ylabel = L"\int_x c(x,\nu)\partial x"
# ax6.title = @sprintf "t = %.2f" tMax/2.0
# u = reshape(sol.u[50],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
# for j=1:N
#     integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
# end
# lines!(ax6, νs[Nghost+1:end-Nghost], integ)
# ax7 = Axis(fig2[2,4], aspect=1)
# ax7.xlabel = "ν"
# ax7.ylabel = L"\int_x c(x,\nu)\partial x"
# ax7.title = @sprintf "t = %.2f" tMax
# u = reshape(sol.u[100],(Nplus,Nplus))[1+Nghost:Nplus-Nghost,1+Nghost:Nplus-Nghost]
# for j=1:N
#     integ[j] = (sum(u[j,2:end])+sum(u[j,1:end-1]))*dx/2.0
# end
# lines!(ax7, νs[Nghost+1:end-Nghost], integ)


# display(fig2)
# #%%
# save("glycosylation.png", fig2)

