import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider,Button
from scipy.integrate import odeint
from math import cos,pi,exp
import networkx as nx

nNodes = 5
nEdges = 8

Graph = nx.DiGraph()

nodes = ["ER", "Golgi", "Endosome", "Extracellular trimer", "Extracellular fibril"]
Graph.add_nodes_from(nodes)
edges = [["ER","Golgi"],["Golgi","ER"],["Golgi","Extracellular trimer"],["Golgi","Extracellular fibril"],["Extracellular trimer","Endosome"],["Endosome","Extracellular fibril"],["Extracellular trimer","Extracellular fibril"]]
edgeLabels = {("ER","Golgi"):"$k_1$",("Golgi","ER"):"$k_2$",("Golgi","Extracellular trimer"):"$k_3$",("Golgi","Extracellular fibril"):"$k_4$",("Extracellular trimer","Endosome"):"$k_5$",("Endosome","Extracellular fibril"):"$k_6$",("Extracellular trimer","Extracellular fibril"):"$k_7$"}
Graph.add_edges_from(edges)
fixedpos = nx.shell_layout(Graph)

ks = np.zeros(nEdges)+0.5

def SecretoryPathwayODEs(y,t):
    C_ER = y[0]
    C_G  = y[1]
    C_E  = y[2]
    C_T  = y[3]
    C_F  = y[4]

    dC_ERdt = ks[0] - ks[1]*C_ER + ks[2]*C_G
    dC_Gdt  = ks[1]*C_ER - ks[2]*C_G - ks[3]*C_G - ks[4]*C_G
    dC_Edt  = ks[5]*C_T - ks[6]*C_E
    dC_Tdt  = ks[3]*C_G - ks[5]*C_T - ks[7]*C_T
    dC_Fdt  = ks[4]*C_G + ks[6]*C_E + ks[7]*C_T - ks[0]*C_F

    x = np.array([dC_ERdt,dC_Gdt ,dC_Edt ,dC_Tdt ,dC_Fdt])/100

    return x


# Set up plot canvas
fig,ax=plt.subplots()
plt.subplots_adjust(left=0.1,bottom=0.5)

sliders = []
axSliders = []
#Create sliders to vary k values
for i in range(nEdges):
    axSliders.append(plt.axes([0.1,0.05*i,0.8,0.03]))
    axSliders[-1].set_xticks([0.0,1.0])
    sliders.append(Slider(axSliders[-1],"$k_{}$".format(i),valinit=0.5,valmin=0,valmax=1))
    axSliders[-1].xaxis.set_visible(True)
    axSliders[-1].set_xticks([0,1])


onOff = False  # Global variable to switch dynamic plotting on and off

time = 0
dt = 1
ts = np.arange(0,96,dt)
output = np.zeros([20,nNodes])

nx.draw_networkx_nodes(Graph,ax=ax,pos=fixedpos,node_size=output[-1,:6]*5000)
nx.draw_networkx_labels(Graph,ax=ax,pos=fixedpos)
nx.draw_networkx_edges(Graph,ax=ax,pos=fixedpos,width=ks*10,arrowstyle="fancy")
nx.draw_networkx_edge_labels(Graph,ax=ax,pos=fixedpos,edge_labels=edgeLabels,font_color='red')


# Function to continually update the plot so long as onOff=True
def updateGraph(val):
    global onOff
    global time
    global output
    # Set onOff = True if currently False
    if onOff:
        pass
    else:
        onOff = True
    # Extract k values from sliders
    for i in range(nEdges):
        ks[i] = sliders[i].val
        widths = np.append(ks[1:5],[ks[6],ks[5],ks[7:]])
    while onOff:
        trange = np.linspace(time,time+20*dt,num=20,endpoint=False)
        output = odeint(SecretoryPathwayODEs,output[-1,:],trange)
        time = time+20*dt
        ax.cla()
        nx.draw_networkx_nodes(Graph,ax=ax,pos=fixedpos,node_size=output[-1,:6]*5000)
        nx.draw_networkx_labels(Graph,ax=ax,pos=fixedpos)
        nx.draw_networkx_edges(Graph,ax=ax,pos=fixedpos,width=widths*10,arrowstyle="fancy")
        nx.draw_networkx_edge_labels(Graph,ax=ax,pos=fixedpos,edge_labels=edgeLabels,font_color='red')
        # Stop briefly before repeating loop
        plt.pause(0.2)

# Call updatedLine in response to slider updates
for i in range(nEdges):
    sliders[i].on_changed(updateGraph)

# Create start and stop buttons
axButton1 = plt.axes([0.1,0.9,0.1,0.1])
btn1 = Button(axButton1,'Start')
axButton2 = plt.axes([0.25,0.9,0.1,0.1])
btn2 = Button(axButton2,'Stop')

# Define function to set onOff=False and stop continuous plot updating
def stopPlot(val):
    global onOff
    onOff = False

# Set responses to button clicks. Button 1 calls updateGraph with current slider values.
btn1.on_clicked(updateGraph)
# Button 2 calls function to stop continuous plotting
btn2.on_clicked(stopPlot)

# Show plot in GUI window
plt.show()
