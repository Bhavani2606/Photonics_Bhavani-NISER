import meep as mp
from meep import mpb
import matplotlib.pyplot as plt
import math
import numpy as np

phi = (math.sqrt(5)+1)/2 #diophantine number
N = 13  # number of bilayers
thickness_list = [] #to store the thickness of the various bilayers
A = 1
t0 = 0.125

def thickness_list(A, t0, N, phi): #to return the thickness of the various bilayers
    thickness_list = []
    for i in range (0, N):
        t = t0*(1+A*math.cos(2*math.pi*phi*i))
        thickness_list.append(t)
    return thickness_list

# to define a function that takes thickness of the bilayers and other constants as inputs and calculate all the frequencies of the given band.
def PhQc_freqlist(A, t0, N, phi, bands):

    t_list1 = np.array(thickness_list(A, t0, N, phi))#to store the thickness of the various bilayers fro a particular set of parameters
    T = np.sum(t_list1)# to calculate the total thickness of the heterostructure
    # print ("Total thickness = ", 2*T)
    geometry_lattice = mp.Lattice(size=mp.Vector3(2*T), basis1 = mp.Vector3(1))
    geometry = [] #to store the heterostructure geometry
    x = 0 #to store the x-coordinate of the starting point of the first block
    for i in range(N):
        geometry.append(mp.Block(size = mp.Vector3(t_list1[i], mp.inf, mp.inf), center = mp.Vector3(x+0.5*t_list1[i], 0, 0), material = mp.Medium(epsilon = 3.5)))
        geometry.append(mp.Block(size = mp.Vector3(t_list1[i], mp.inf, mp.inf), center = mp.Vector3(x+1.5*t_list1[i], 0, 0), material = mp.Medium(epsilon = 1.5)))
        x += 2*t_list1[i]

    k_points = [mp.Vector3(), mp.Vector3(0.5), mp.Vector3()]
    k_points = mp.interpolate(16, k_points)
    num_bands = bands
    resolution = 32
    ms = mpb.ModeSolver(geometry = geometry, geometry_lattice = geometry_lattice, k_points = k_points, resolution = resolution, num_bands = num_bands)
    ms.run_tm()
    tm_freqs = ms.all_freqs
    tm_freqs = np.array(tm_freqs)
    return tm_freqs[0]

# A = 0.1
plot_y_list = []
for A in np.arange(0.0, 1.0, 0.05):
    freq_list = PhQc_freqlist(A, t0, N, phi, 60)
    plot_y_list.append(freq_list)
#to plot multiple data points in a single graph
plot_y_list = np.array(plot_y_list)
plot_y_list = np.transpose(plot_y_list)
# plot_y_list = plot_y_list
print(plot_y_list)
fig, ax = plt.subplots()
A_list = np.arange(0.0, 1.0, 0.05)

# Plot bands
#to plot the data point of each list in plot_y_list in a single plot vs A
for y in plot_y_list:
    ax.plot(A_list, y, color='blue')
# ax.plot(tm_freqs, color='blue')
# ax.set_xlim([A_list[0], A_list[-1]])
# #set limit of y-axis so that the graph is more readable
# ax.set_ylim([0, ])
ax.axis('auto')
ax.set_autoscale_on(True)
ax.set_ylabel('frequency (c/a)', size=16)
ax.set_xlabel('A', size=16)
plt.show()
