import meep as mp
from meep import mpb
import matplotlib.pyplot as plt
import math
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

phi = 144/89 #diophantine number
N = 89  # number of bilayers
thickness_list = [] #to store the thickness of the various bilayers
A = 1
t0 = 0.125

def thickness_list(A, t0, N, phi): #to return thte thickness of the various bilayers
    thickness_list = []
    for i in range (0, N):
        t = t0*(1+A*math.cos(2*math.pi*phi*i))
        thickness_list.append(t)
    return thickness_list



# to define a function that takes thickness of the bilayers other constants as inputs and calcultes all the frequencies of hte given band.
def PhQc_freqlist(A, t0, N, phi, bands):

    t_list1 = np.array(thickness_list(A, t0, N, phi))#to store the thickness of the various bilayers fro a particular set of parameters
    T = np.sum(t_list1)# to calculate the total thickness of the heterostructure
    print ("Total thickness = ", 2*T)
    geometry_lattice = mp.Lattice(size=mp.Vector3(2*T))
    geometry = [] #to store the heterostructure geometry
    x = 0 #to store the x-coordinate of the starting point of the first block
    for i in range(N):
        geometry.append(mp.Block(size = mp.Vector3(t_list1[i], mp.inf, mp.inf), center = mp.Vector3(x+0.5*t_list1[i], 0, 0), material = mp.Medium(epsilon = 3.5)))
        geometry.append(mp.Block(size = mp.Vector3(t_list1[i], mp.inf, mp.inf), center = mp.Vector3(x+1.5*t_list1[i], 0, 0), material = mp.Medium(epsilon = 1.5)))
        x += 2*t_list1[i]

    k_points = [mp.Vector3(), mp.Vector3(0.5), mp.Vector3()]# defining the ends of brillouin zone
    k_points = mp.interpolate(16, k_points)
    num_bands = bands
    resolution = 36
    ms = mpb.ModeSolver(geometry = geometry, geometry_lattice = geometry_lattice, k_points = k_points, resolution = resolution, num_bands = num_bands)
    ms.run_tm()
    tm_freqs = ms.all_freqs
    tm_freqs = np.array(tm_freqs)
    IPR_array = []#to store the IPR values for different bands
    for i in range(1, bands+1):
        field = ms.get_hfield(i)
        field = field[:, 0, 0, 1]
        field_2 = (field*np.conj(field))
        IPR = (np.sum(field_2**2))/((np.sum(field_2))**2)
        IPR_array.append(IPR)
    return tm_freqs[0], IPR_array

# To plot the eigen-spectrum of the heterostructure for different values of A
plot_y_list = []
IPR_list = []
for A in np.arange(0.0, 1.0, 0.01):
    freq_list, IPR_arr = PhQc_freqlist(A, t0, N, phi, 200)
    IPR_list.append(np.real(IPR_arr))
    plot_y_list.append(freq_list)

# print(np.shape(IPR_list), np.shape(plot_y_list))
#to plot multiple data points in a single graph
plot_y_list = np.array(plot_y_list)
IPR_list = np.log(np.array(IPR_list))
# plot_y_list = plot_y_list/(2*math.pi)
plot_y_list = np.transpose(plot_y_list)
IPR_list = np.transpose(IPR_list)

IPR_max, IPR_min = np.max(IPR_list), np.min(IPR_list)
# print(plot_y_list)
A_list = np.arange(0.0, 1.0, 0.01)

fig, ax = plt.subplots()

# color_data = np.linspace(IPR_min, IPR_max, 1000)  # This will determine the color variation

# Create a colormap
colormap = cm.get_cmap('viridis')  # You can choose other colormaps such as 'plasma', 'inferno', etc.

# Normalize color data
# norm = mcolors.Normalize(vmin=min(color_data), vmax=max(color_data))

# Plot bands
# Scatter plot for multiple y values, see https://stackoverflow.com/a/34280815/2261298
#to plot the data point of each list in plot_y_list in a single plot vs A
for i in range (0, len(plot_y_list)):
    color_data = IPR_list[i]
    norm = mcolors.Normalize(vmin=min(color_data), vmax=max(color_data))
    
    for j in range(len(A_list) - 1):
        # Calculate color for this segment
        color = colormap(norm(IPR_list[i][j]))
        
        # Plot the segment of the line
        ax.scatter(A_list[j], plot_y_list[i][j], color=color, lw=2)

# Add a colorbar to indicate the color mapping
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation='vertical')
cbar.set_label('IPR')

# ax.plot(tm_freqs, color='blue')
ax.set_xlim([A_list[0], A_list[-1]])
ax.set_ylim([0, 3])
ax.set_ylabel('frequency (c/a)', size=16)
ax.set_xlabel('A', size=16)
plt.show()
