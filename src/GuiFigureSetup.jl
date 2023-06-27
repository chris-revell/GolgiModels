#
#  GuiFigureSetup.jl
#  GolgiModels
#
#  Created by Christopher Revell on 28/04/2023.

module GuiFigureSetup

using GLMakie
using FileIO
using DrWatson

# Function to setup figure
function guiFigureSetup(ksInit)
    # Set up figure canvas
    fig = Figure(resolution=(1700,1500),fontsize=32)
    axDiagram = Axis(fig[3,1:4],title="Model diagram",aspect=DataAspect())
    image!(axDiagram,rotr90(load(projectdir("supplementary","model.png"))))
    hidedecorations!(axDiagram)
    hidespines!(axDiagram)
    axCis = Axis(fig[1,1], aspect=0.55, ylabel = "Compartment size")
    xlims!(axCis,(0,3))
    axMed = Axis(fig[1,2], aspect=0.55, yticksvisible=false)
    xlims!(axMed,(0,3))
    axTra = Axis(fig[1,3], aspect=0.55, yticksvisible=false)
    xlims!(axTra,(0,3))
    Label(fig[1,1,Bottom()],"Cis",fontsize=32)
    Label(fig[1,2,Bottom()],"Medial",fontsize=32)
    Label(fig[1,3,Bottom()],"Trans",fontsize=32)

    # Set up parameter sliders
    parameterSliders = SliderGrid(
        fig[1,4],
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
        (label="k₁₂, t₁ → ∅      " , range=0.0:0.01:1.2, startvalue=ksInit[12], format="{:.2f}");
    )

    # Add stop/start button
    run = Button(fig[2,1]; label = "Start/Stop", tellwidth = false)
    reset = Button(fig[2,2]; label = "Reset", tellwidth = false)

    colsize!(fig.layout, 1, Relative(0.25))
    colsize!(fig.layout, 2, Relative(0.25))
    colsize!(fig.layout, 3, Relative(0.25))
    colsize!(fig.layout, 4, Relative(0.25))
    rowsize!(fig.layout, 1, Aspect(1, 2.0))
    rowsize!(fig.layout, 2, Aspect(1, 0.1))
    resize_to_layout!(fig)

    return fig, axCis, axMed, axTra, parameterSliders, run, reset
end


export guiFigureSetup

end