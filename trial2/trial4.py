import math
import meep as mp
from meep import mpb
import matplotlib.pyplot as plt
import numpy as np

geometry_lattice = mp.Lattice(size=mp.Vector3(13), basis1 = mp.Vector3(1, 0))

x = 0
t = 0.125
N = 13
phi = (math.sqrt(5)+1)/2

def thickness_list(A, t0, N, phi): #to return thte thickness of the various bilayers
    thickness_list = []
    for i in range (0, N):
        t = t0*(1+A*math.cos(2*math.pi*phi*i))
        thickness_list.append(t)
    return thickness_list

t_list1 = thickness_list(0, t, 13, phi)
geometry = []
for i in range(N):
    geometry.append(mp.Block(size = mp.Vector3(t_list1[i], mp.inf, mp.inf), center = mp.Vector3(x+0.5*t_list1[i], 0, 0), material = mp.Medium(epsilon = 3.5)))
    geometry.append(mp.Block(size = mp.Vector3(t_list1[i], mp.inf, mp.inf), center = mp.Vector3(x+1.5*t_list1[i], 0, 0), material = mp.Medium(epsilon = 1.5)))
    x += 2*t_list1[i]
k_points = [mp.Vector3(), mp.Vector3(0.5), mp.Vector3()]
k_points = mp.interpolate(8, k_points)
num_bands = 8
resolution = 8
ms = mpb.ModeSolver(geometry = geometry, geometry_lattice = geometry_lattice, k_points = k_points, resolution = resolution, num_bands = num_bands)
ms.run_tm()
tm_freqs = ms.all_freqs
tm_gaps = ms.gap_list
ms.run_te()
te_freqs = ms.all_freqs
te_gaps = ms.gap_list
# output = ms.output_epsilon()

# epsilon = ms.get_epsilon_point(mp.Vector3(1))
fig, ax = plt.subplots()
x = range(len(tm_freqs))
# Plot bands
# Scatter plot for multiple y values, see https://stackoverflow.com/a/34280815/2261298
for xz, tmz, tez in zip(x, tm_freqs, te_freqs):
    ax.scatter([xz]*len(tmz), tmz, color='blue')
    # ax.scatter([xz]*len(tez), tez, color='red', facecolors='none')
ax.plot(tm_freqs, color='blue')
# ax.plot(te_freqs, color='red')
ax.set_ylim([0, 0.5])
ax.set_xlim([x[0], x[-1]])
print("frequency values = ", tm_freqs)
# Plot gaps
# for gap in tm_gaps:
#     if gap[0] > 1:
#         ax.fill_between(x, gap[1], gap[2], color='blue', alpha=0.2)

# for gap in te_gaps:
#     if gap[0] > 1:
#         ax.fill_between(x, gap[1], gap[2], color='red', alpha=0.2)


# # Plot labels
ax.text(12, 0.75, 'TM bands', color='blue', size=15)
# ax.text(13.05, 0.8, 'TE bands', color='red', size=15)

# points_in_between = (len(tm_freqs) - 4) / 3
# tick_locs = [i*points_in_between+i for i in range(4)]
# tick_labs = ['Γ', 'X', 'M', 'Γ']
# ax.set_xticks(tick_locs)
# ax.set_xticklabels(tick_labs, size=16)
ax.set_ylabel('frequency (c/a)', size=16)
ax.grid(True)

plt.show()
# print("Epsilon values = ", epsilon)
