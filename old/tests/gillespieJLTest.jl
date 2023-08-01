using Gillespie
using Gadfly
using Random

function F(x,parms)
  (S,I,R) = x
  (beta,gamma) = parms
  infection = beta*S*I
  recovery = gamma*I
  [infection,recovery]
end

x0 = [9999,1,0]

nu = [[-1 1 0];[0 -1 1]]

parms = [0.1/10000.0,0.05]
tf = 1000.0

Random.seed!(1236)

result = ssa(x0,F,nu,parms,tf);

data = ssa_data(result)

plot(data,
  layer(x="time",y="x1",Geom.step,Theme(default_color=colorant"red")),
  layer(x="time",y="x2",Geom.step,Theme(default_color=colorant"blue")),
  layer(x="time",y="x3",Geom.step,Theme(default_color=colorant"green")),
  Guide.xlabel("Time"),
  Guide.ylabel("Number"),
  Guide.manual_color_key("Population",
                            ["S", "I", "R"],
                            ["red", "blue", "green"]),
  Guide.title("SIR epidemiological model"))
