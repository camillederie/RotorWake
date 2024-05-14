"""
Plotting function
"""
from create_geometry import BladeGeometry
from lifting_line_solver import LiftingLineSolver
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from read_BEMdata_into_Python import read_matlab_data
import numpy as np


phase_shifts = np.linspace(0, 180, 4)
colors = ["lawngreen", "deepskyblue", "orangered", "darkviolet"]
# start plots
fig_rotor, ax_rotor = plt.subplots(1, 1, dpi=150)

nspan = 20
ntheta = 200
nblades = 3
spacing = 'equal'
nrotor = 2

first_blades_right = []
first_blade_left = []

for color, phase_shift in zip(colors, phase_shifts):


    prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, phase_diff=phase_shift,
                             double_rotor=True,
                             n_theta=ntheta, spacing=spacing, a=0.33, xshift=0, yshift=100, zshift=0)
    blade = prop_geo.bladepanels

    theta = np.linspace(0, 2 * np.pi)
    y = prop_geo.span_arr[0] * np.cos(theta) * prop_geo.radius
    z = prop_geo.span_arr[0] * np.sin(theta) * prop_geo.radius

    ax_rotor.plot(y, z, color='k')
    ax_rotor.plot(y + 100, z, color='k')
    # split blade panels
    blades_split = np.split(blade, 6)

    for i, blade_i in enumerate(blades_split):
        y = blade_i[:, [1, 4, 7, 10]].flatten()
        z = blade_i[:, [2, 5, 8, 11]].flatten()

        ax_rotor.plot(y, z, c='k')

    first_blades_right.append(blades_split[3])
    first_blade_left.append(blades_split[0])

ax_rotor.plot(first_blade_left[0][:, [1, 4, 7, 10]].flatten(),
              first_blade_left[0][:, [2, 5, 8, 11]].flatten(),
              c=colors[0])
for phase_shift, blade_i, color in zip(phase_shifts, first_blades_right, colors):
    ax_rotor.plot(blade_i[:, [1, 4, 7, 10]].flatten(),
                  blade_i[:, [2, 5, 8, 11]].flatten(),
                  c=color, label=r"$\gamma$: {}$^\circ$".format(phase_shift))

ax_rotor.set_xlabel('x (m)', fontsize=14)
ax_rotor.set_ylabel('y (m)', fontsize=14)
ax_rotor.legend(prop={"size": 14}, loc="lower center")
# ax_rotor.grid(True)
# ax_rotor.set_aspect("equal", adjustable="box")
fig_rotor.suptitle("Double Rotor Configuration")