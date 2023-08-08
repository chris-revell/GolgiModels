#
#  GuiFigureSetup.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

# module GuiFigureSetup

using WGLMakie
using FileIO
using DrWatson

# Function to setup figure
function guiFigureSetup(ksInit)
    
    fig = Figure(resolution=(2000,1500),fontsize=32)

    grd1 = GridLayout(fig[1,1])
    grd2 = GridLayout(fig[2,1])
    grd3 = GridLayout(fig[3,1])

    axDiagram = Axis(grd1[1,1],title="Model diagram",aspect=DataAspect())
    image!(axDiagram,rotr90(load(projectdir("supplementary","GolgiCompartmentModel.png"))))
    hidedecorations!(axDiagram)
    hidespines!(axDiagram)

    axCis = Axis(grd2[1,1], ylabel = "Compartment size")
    # xlims!(axCis,(0,3))
    axMed = Axis(grd2[1,2], yticksvisible=false)
    # xlims!(axMed,(0,3))
    axTra = Axis(grd2[1,3], yticksvisible=false)
    # xlims!(axTra,(0,3))

    axDwell = Axis(grd3[1:2,4],title = "Dwell Times")
    xlims!(axDwell,(1.5,7.5))
    ylims!(axDwell,(0,1))
    axDwell.xticks = (1:7, ["∅", "C₁", "C₊", "M₁", "M₊", "T₁", "T₊"])
    axDwell.ylabel = "Relative dwell time"

    axReducedDiagram = Axis(grd3[1,1:3],title="Reduced model",aspect=DataAspect())
    hidedecorations!(axReducedDiagram); hidespines!(axReducedDiagram)
    image!(axReducedDiagram,rotr90(load(projectdir("supplementary","GolgiCompartmentModel_reduced.png"))))

    Label(grd2[1,1,Bottom()],"Cis",fontsize=32)
    Label(grd2[1,2,Bottom()],"Medial",fontsize=32)
    Label(grd2[1,3,Bottom()],"Trans",fontsize=32)

    # Set up parameter sliders
    parameterSliders = SliderGrid(
        grd1[1,2],
        (label="k₁,  ∅ → c₁      " , range=0.0:0.01:1.2, startvalue=ksInit[1], format="{:.2f}"),
        (label="k₂,  c₁+cₙ → cₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[2], format="{:.2f}"),
        (label="k₃,  cₙ → c₁+cₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[3], format="{:.2f}"),
        (label="k₄,  c₁ → m₁     " , range=0.0:0.01:1.2, startvalue=ksInit[4], format="{:.2f}"),
        (label="k₅,  m₁ → c₁     " , range=0.0:0.01:1.2, startvalue=ksInit[5], format="{:.2f}"),
        (label="k₆,  m₁+mₙ → mₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[6], format="{:.2f}"),
        (label="k₇,  mₙ → m₁+mₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[7], format="{:.2f}"),
        (label="k₈,  m₁ → t₁     " , range=0.0:0.01:1.2, startvalue=ksInit[8], format="{:.2f}"),
        (label="k₉,  t₁ → m₁     " , range=0.0:0.01:1.2, startvalue=ksInit[9], format="{:.2f}"),
        (label="k₁₀, t₁+tₙ → tₙ₊₁" , range=0.0:0.01:1.2, startvalue=ksInit[10], format="{:.2f}"),
        (label="k₁₁, tₙ → t₁+tₙ₋₁" , range=0.0:0.01:1.2, startvalue=ksInit[11], format="{:.2f}"),
        (label="k₁₂, t₁ → ∅      " , range=0.0:0.01:1.2, startvalue=ksInit[12], format="{:.2f}"),
        width = 700,
    )

    # Add stop/start button
    run = Makie.Button(grd3[2,1]; label = "Start/Stop", tellwidth = false)
    reset = Makie.Button(grd3[2,2]; label = "Reset", tellwidth = false)

    linearityToggle = Toggle(fig, active = true)
    toggleLabel = Label(fig, lift(x -> x ? "Linear" : "Nonlinear", linearityToggle.active))
    grd3[2, 3] = grid!(hcat(linearityToggle, toggleLabel), tellheight = false)

    rowsize!(fig.layout,2,Relative(0.25))
    rowsize!(fig.layout,3,Relative(0.25))
    resize_to_layout!(fig)

    xLimTimeAv = [5.0]
    xlims!(axCis,(0.0,1.1*xLimTimeAv[1]))
    xlims!(axMed,(0.0,1.1*xLimTimeAv[1]))
    xlims!(axTra,(0.0,1.1*xLimTimeAv[1]))                

    # Pull parameters from slider positions
    # kObservables = [s.value for s in parameterSliders.sliders]

    return fig, axCis, axMed, axTra, axDwell, parameterSliders, run, reset, linearityToggle, xLimTimeAv, linearityToggle#, kObservables
end


# export guiFigureSetup

# end