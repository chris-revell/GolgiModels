from sys import argv
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import cv2
from skimage.util import invert
from skimage.morphology import medial_axis
from astropy import modeling

# Function to import image and generate bw image
def read_cell(file):
    image = cv2.imread(file)
    gray_image = invert(cv2.cvtColor(image, cv2.COLOR_BGR2GRAY))
    #blur = cv2.GaussianBlur(gray_image, (5, 5), 0)
    #ret3, bw_image = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    bw_image = gray_image
    return image,bw_image

# Function to generate medial axis
def med_axis(bw_image):
    skel, distance = medial_axis(bw_image, return_distance=True)
    dist_on_skel = distance * skel/2.0
    return dist_on_skel

# Create black and white image
originalimage,bw = read_cell(argv[1])

# Create medial axis from imported image
maxis = med_axis(bw)

# Visually inspect axis to avoid loops
#plt.imshow(bw+100*maxis,cmap="gray_r")
#plt.show()

# Array to store the indices of all neighbouring pixels within the medial axis for each pixel within the medial axis.
neighbours = np.zeros((np.shape(maxis)[0],np.shape(maxis)[1],18),dtype='int')
# Routing to find all such neighbours
for x in range(1,np.shape(maxis)[0]-1):
    for y in range(1,np.shape(maxis)[1]-1):
        if maxis[x,y] > 0:
            # if maxis[x,y] > 0 then pixel x,y is part of the medial axis.
            # Now test all pixels surrounding the current pixel under consideration.
            for i in range(-1,2):
                for j in range(-1,2):
                    if (i,j) != (0,0) and maxis[x+i,y+j]>0:
                        # 0th component of 3rd dimension stores the number of neighbours pixel x,y has. Remainder of components store the indices of those neighbours.
                        neighbours[x,y,0] = neighbours[x,y,0]+1
                        neighbours[x,y,int(2*neighbours[x,y,0]-1)] = x+i
                        neighbours[x,y,int(2*neighbours[x,y,0])]   = y+j

# Find all tips (1 neighbouring pixel) and branch points (3 neighbour pixels) within the medial axis
tips = []
branches = []
pseudobranches = []
for x in range(np.shape(maxis)[0]):
    for y in range(np.shape(maxis)[1]):
        if neighbours[x,y,0] == 1:
            tips.append((x,y))
        elif neighbours[x,y,0] == 3:
            notBranch = 0
            for i in range(3):
                for j in range(3):
                    if i==j:
                        pass
                    else:
                        xdif = neighbours[x,y,2*i+1]-neighbours[x,y,2*j+1]
                        ydif = neighbours[x,y,2*i+2]-neighbours[x,y,2*j+2]
                        if (xdif==0 and abs(ydif)==1) or (ydif==0 and abs(xdif)==1):
                            notBranch = 1
            if notBranch == 0:
                branches.append((x,y))
            else:
                pseudobranches.append((x,y))

# For each tip of the medial axis, find all possible paths to other tips in order to find the longes path through the axis
maxLength = 0
route = []
widths = []
for t in tips:

    lastPixel = t
    currentPixel = (neighbours[lastPixel[0],lastPixel[1],1],neighbours[lastPixel[0],lastPixel[1],2])
    length = sqrt((lastPixel[0]-currentPixel[0])**2+(lastPixel[1]-currentPixel[1])**2)
    pixelsSinceLastBranch = [t]
    blockedPixels = np.zeros(np.shape(bw))
    route = [t,currentPixel]
    widths = [maxis[t[0],t[1]],maxis[currentPixel[0],currentPixel[1]]]

    while blockedPixels[t[0],t[1]] == 0:

        if currentPixel in tips:

            #pathLengths.append(length)
            if length > maxLength:
                maxLength = length
                longestRoute = route
                longestRouteWidths = widths
            for pixel in pixelsSinceLastBranch:
                blockedPixels[pixel[0],pixel[1]] = 1
            pixelsSinceLastBranch.clear()
            lastPixel = t
            currentPixel = (neighbours[lastPixel[0],lastPixel[1],1],neighbours[lastPixel[0],lastPixel[1],2])
            length = sqrt((lastPixel[0]-currentPixel[0])**2+(lastPixel[1]-currentPixel[1])**2)
            pixelsSinceLastBranch = [t,currentPixel]
            route = [t,currentPixel]
            widths = [maxis[t[0],t[1]],maxis[currentPixel[0],currentPixel[1]]]

        elif currentPixel in branches:
            for i in range(4):
                if i == 3:
                    #set pixels back to last branch to blocked and return to tip
                    for pixel in pixelsSinceLastBranch:
                        blockedPixels[pixel[0],pixel[1]] = 1
                    currentPixel = t
                else:
                    pixelConsidered = (neighbours[currentPixel[0],currentPixel[1],2*i+1],neighbours[currentPixel[0],currentPixel[1],2*i+2])
                    if pixelConsidered == lastPixel or blockedPixels[pixelConsidered[0],pixelConsidered[1]] == 1:
                        pass
                    else:
                        newPixel = pixelConsidered
                        break
            if currentPixel == t:
                #End this loop iteration if we have returned to the start
                pixelsSinceLastBranch.clear()
                lastPixel = t
                currentPixel = (neighbours[lastPixel[0],lastPixel[1],1],neighbours[lastPixel[0],lastPixel[1],2])
                length = sqrt((lastPixel[0]-currentPixel[0])**2+(lastPixel[1]-currentPixel[1])**2)
                pixelsSinceLastBranch = [t,currentPixel]
                route = [t,currentPixel]
                widths = [maxis[t[0],t[1]],maxis[currentPixel[0],currentPixel[1]]]
                continue
            else:
                length = length + sqrt((newPixel[0]-currentPixel[0])**2+(newPixel[1]-currentPixel[1])**2)
                lastPixel = currentPixel
                currentPixel = newPixel
                pixelsSinceLastBranch.clear()
                pixelsSinceLastBranch.append(currentPixel)
                route.append(currentPixel)
                widths.append(maxis[currentPixel[0],currentPixel[1]])

        elif currentPixel in pseudobranches:
            for i in range(3):
                pixelConsidered = (neighbours[currentPixel[0],currentPixel[1],2*i+1],neighbours[currentPixel[0],currentPixel[1],2*i+2])
                if pixelConsidered == lastPixel or pixelConsidered in pseudobranches:
                    pass
                else:
                    newPixel = pixelConsidered
                    break
            length = length + sqrt((newPixel[0]-currentPixel[0])**2+(newPixel[1]-currentPixel[1])**2)
            lastPixel = currentPixel
            currentPixel = newPixel
            pixelsSinceLastBranch.append(currentPixel)
            route.append(currentPixel)
            widths.append(maxis[currentPixel[0],currentPixel[1]])

        else:
            for i in range(2):
                pixelConsidered = (neighbours[currentPixel[0],currentPixel[1],2*i+1],neighbours[currentPixel[0],currentPixel[1],2*i+2])
                if pixelConsidered == lastPixel:
                    pass
                else:
                    newPixel = pixelConsidered
            length = length + sqrt((newPixel[0]-currentPixel[0])**2+(newPixel[1]-currentPixel[1])**2)
            lastPixel = currentPixel
            currentPixel = newPixel
            pixelsSinceLastBranch.append(currentPixel)
            route.append(currentPixel)
            widths.append(maxis[currentPixel[0],currentPixel[1]])

# Place longest axis route on image and express distance along path and width as a function of distance as xs and ys arrays.
pathImage = np.zeros(np.shape(bw))
xs = np.zeros(len(longestRoute))
ys = np.zeros(len(longestRoute))
dist=0
xs[0] = dist
ys[0] = longestRouteWidths[0]
for i,pixel in enumerate(longestRoute[1:]):
    pathImage[pixel[0],pixel[1]] = 1
    dist = dist + sqrt((longestRoute[i+1][0]-longestRoute[i][0])**2+(longestRoute[i+1][1]-longestRoute[i][1])**2)
    xs[i+1] = dist
    ys[i+1] = longestRouteWidths[i+1]
# Find the location of the widest point in the axis
widestPixel = longestRoute[np.argmax(longestRouteWidths)]

# Normalise cell area to 1
area = np.trapz(ys,xs)
xs = xs/sqrt(area)
ys = ys/sqrt(area)

# Find normalised max values
widestValue = np.max(ys)
pathLength = np.max(xs)

# Fit Gaussian to normalised width data
fitter = modeling.fitting.LevMarLSQFitter()
model = modeling.models.Gaussian1D()   # depending on the data you need to give some initial values
fitted_model = fitter(model, xs, ys)


# Save cell image with spline and widest point 
img = (bw+1000*pathImage)
h, w = np.shape(img)
fig0 = plt.figure(figsize=(w/100.,h/100.), dpi=100)
ax0 = fig0.add_axes([0,0,1,1])
ax0.axis('off')
ax0.imshow(img, zorder=0,  extent=[0,w,  h, 0], interpolation="nearest",cmap="gray_r")
ax0.scatter(widestPixel[1],widestPixel[0],color="red",s=2)
plt.savefig(argv[1][:-4]+"WithAxis.png", bbox_inches=0, dpi=100)


# Plot width data with fitted Gaussian
fig1,ax1 = plt.subplots(1,1,figsize=(6,6))
ax1.plot(xs,ys,label="Raw data")
ax1.plot(xs,fitted_model(xs),label="Gaussian:\nA={:04.2f}\n$\sigma$={:04.2f}\n$\mu$={:04.2f}".format(fitted_model.amplitude.value,fitted_model.stddev.value,fitted_model.mean.value))
ax1.set_xlabel("Axis distance")
ax1.set_ylabel("Width")
ax1.legend(loc="best")
fig1.savefig(argv[1][:-4]+"GaussianFit.png",bbox_inches='tight',dpi=300)

# Save fit values to file
outfile = open(argv[1][:-4]+"Parameters.txt","w")
outfile.write("Amplitude = "+str(fitted_model.amplitude.value)+"\n")
outfile.write("Mean = "+str(fitted_model.mean.value)+"\n")
outfile.write("Standard deviation = "+str(fitted_model.stddev.value)+"\n")
outfile.write("Ratio std/amp = "+str(fitted_model.stddev.value/fitted_model.amplitude.value)+"\n")
outfile.write("Normalised max width = "+str(widestValue)+"\n")
outfile.write("Normalised axis length = "+str(pathLength)+"\n")
outfile.close()

# Find cell number
parts = (argv[1]).split("_")
for n,i in enumerate(parts):
    if "/Cell" in i:
        Nc = parts[n][-1]
        break
    else:
        pass
# Determine whether before or in metaphase
if "Before" in argv[1]:
    M = 0
else:
    M= 1


# Output details to stdout
print(Nc, ", ", M ,",",fitted_model.mean.value,",",fitted_model.amplitude.value,",",fitted_model.stddev.value,",",widestValue,",",pathLength)
#%%