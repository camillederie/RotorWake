"""
Plotting function
"""
from create_geometry import BladeGeometry
from lifting_line_solver import LiftingLineSolver
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from read_BEMdata_into_Python import read_matlab_data
import numpy as np

nspan = 20
ntheta = 200
nblades = 3
spacing = "equal"
nrotor = 1

prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta, spacing=spacing,
                         a=0, double_rotor=False, phase_diff=0)

blade = prop_geo.bladepanels
rings = prop_geo.filaments
solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6, n_iter=200)

data = solver.run_solver()
omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius
[CP_LLM, CT_LLM, CP_LLM2, CT_LLM2] = solver.CP_and_CT(np.resize(data[0], data[2].shape), np.resize(data[1], data[2].shape), data[2],
                                    np.resize(data[3], data[2].shape), np.resize(data[4], data[2].shape),
                                    solver.geo.v_inf, omega, solver.geo.radius, nblades)

c=["lawngreen", "deepskyblue", "orangered", "darkviolet"]
# =============================================================================
# Rotor performance plots
# =============================================================================
def plottingFunction(prop_geo, data):
    plt.plot(prop_geo.centerPoints, data[-1][:len(prop_geo.centerPoints)])
    plt.show()
    return


# plottingFunction(solver,prop_geo,data)

import random

# fig = plt.figure()
# ax = Axes3D(fig)


# =============================================================================
# BLADE VIZ
# =============================================================================
fig = plt.figure(dpi=150)
ax1 = fig.add_subplot(111, projection="3d")
# ax2 = fig.add_subplot(132, projection="3d")
# ax3 = fig.add_subplot(133, projection="3d")
# axes = [ax1, ax2, ax3]

# split blade panels
blades_split = np.split(blade, 3)

for blade_i in (blades_split):
    x = blade_i[:, [0, 9]]
    y = blade_i[:, [1, 10]]
    z = blade_i[:, [2, 11]]

    ax1.plot_surface(x, y, z, color='k')

# =============================================================================
# Draw rotor hub
# =============================================================================
theta = np.linspace(0, 2 * np.pi)
y = solver.geo.span_arr[0] * np.cos(theta) * solver.r_rotor
z = solver.geo.span_arr[0] * np.sin(theta) * solver.r_rotor

ax1.plot(y * 0, y, z, color='k')

# =============================================================================
# WAKE VIZ
# =============================================================================

ax1.plot_wireframe(rings[0, :nspan - 1, :ntheta + 1],
                   rings[1, :nspan - 1, :ntheta + 1],
                   rings[2, :nspan - 1, :ntheta + 1],
                   color='green', cstride=0)
ax1.plot_wireframe(rings[0, nspan - 1:2 * nspan - 2, :ntheta + 1],
                   rings[1, nspan - 1:2 * nspan - 2, :ntheta + 1],
                   rings[2, nspan - 1:2 * nspan - 2, :ntheta + 1],
                   color='blue', cstride=0)
ax1.plot_wireframe(rings[0, (2 * nspan - 2):, :ntheta + 1],
                   rings[1, (2 * nspan - 2):, :ntheta + 1],
                   rings[2, (2 * nspan - 2):, :ntheta + 1],
                   color='red', cstride=0)

ax1.view_init(0, 180)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
ax1.set_xlim(-5, 5)
ax1.set_ylim(-50, 50)
ax1.set_zlim(-50, 50)
ax1.axis("auto")

# =============================================================================
# BEM COMPARISON
# =============================================================================
plt.close('All')

[BEM_rR, BEM_alpha, BEM_phi, BEM_rho, BEM_Ax, BEM_Az, BEM_Gamma,
 BEM_CT, BEM_CP, BEM_a, BEM_aline, BEM_vinf, BEM_radius] = read_matlab_data()

# Radial distribution alpha and phi

plt.figure(figsize=(8, 6), dpi=150)
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_alpha, 0), BEM_rR.shape)[0, :] * 180 / np.pi, '--', color=c[2], label=r'$\alpha$ BEM')
plt.plot(data[2][:nspan - 1, 0], np.degrees(np.resize(data[6], data[2].shape)[:nspan - 1, 0]), '-', color=c[2], label=r'$\alpha$ LLM')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_phi, 0), BEM_rR.shape)[0, :] * 180 / np.pi, '--', color=c[3], label=r'$\phi$ BEM')
plt.plot(data[2][:nspan - 1, 0], np.degrees(np.resize(data[7], data[2].shape)[:nspan - 1, 0]), '-', color=c[3], label=r'$\phi$ LLM')
plt.xlabel('Radial location r/R (-)', fontsize=15)
plt.ylabel('Angle (deg)', fontsize=15)
plt.title(r'Radial distribution of $\alpha$ and $\phi$', fontsize=16)
plt.legend(fontsize=15)
plt.grid(True)
# plt.savefig('BEMcomp_figures/alpha_phi.eps', bbox_inches='tight', format='eps')
plt.savefig('BEMcomp_figures/alpha_phi.eps', format='eps')

# Radial distribution F_tan en F_ax

F_nondim_BEM = 0.5 * BEM_rho[0] * (BEM_vinf[0]**2) * BEM_radius[0]
F_nondim_LLM = 0.5 * (solver.geo.v_inf ** 2) * solver.geo.radius

plt.figure(figsize=(8, 6), dpi=150)
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Ax, 0), BEM_rR.shape)[0, :]/F_nondim_BEM, '--', color=c[2], label=r'$F_{ax}$ BEM')
plt.plot(data[2][:nspan - 1, 0], np.resize(data[3], data[2].shape)[:nspan - 1, 0]/F_nondim_LLM, '-', color=c[2], label=r'$F_{ax}$ LLM')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Az, 0), BEM_rR.shape)[0, :]/F_nondim_BEM, '--', color=c[3], label=r'$F_{tan}$ BEM')
plt.plot(data[2][:nspan - 1, 0], np.resize(data[4], data[2].shape)[:nspan - 1, 0]/F_nondim_LLM, '-', color=c[3], label=r'$F_{tan}$ LLM')
plt.xlabel('Radial location r/R (-)', fontsize=15)
plt.ylabel('Force $F$ (-)', fontsize=15)
plt.title(r'Radial distribution of $F_{ax}$ and $F_{az}$', fontsize=16)
plt.legend(fontsize=15)
plt.grid(True)
# plt.savefig('BEMcomp_figures/forces.eps', bbox_inches='tight', format='eps')
plt.savefig('BEMcomp_figures/forces.eps', format='eps')

# Radial distribution circulation

# made non-dimensional with (np.pi * Uinf**2) / (NBlades*Omega)
circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)

plt.figure(figsize=(8, 6), dpi=150)
plt.plot(BEM_rR[0, :], BEM_Gamma[0, :], '--', color=c[2], label=r'$\Gamma$ BEM')
plt.plot(data[2][:nspan - 1, 0], np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim, '-', color=c[2], label=r'$\Gamma$ LLM')
plt.xlabel('Radial location r/R (-)', fontsize=15)
plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
plt.legend(fontsize=15)
plt.grid(True)
# plt.savefig('BEMcomp_figures/circulation.eps', bbox_inches='tight', format='eps')
plt.savefig('BEMcomp_figures/circulation.eps', format='eps')

# Radial distribution of induction factors

# fig, ax = plt.subplots(1, 2, dpi=150)
# ax[0].plot(BEM_rR[0, :], BEM_a[0, :], label="BEM")
# ax[0].plot(data[2][:nspan - 1], data[0][:nspan - 1], label="Lifting Line")
# ax[1].plot(BEM_rR[0, :], BEM_aline[0, :], label="BEM")
# ax[1].plot(data[2][:nspan - 1], data[1][:nspan - 1], label="Lifting Line")
# ax[0].set_xlabel("r/R [-]")
# ax[0].set_ylabel("a [-]")
# ax[1].set_xlabel("r/R [-]")
# ax[1].set_ylabel("a' [-]")
# ax[0].grid()
# ax[1].grid()
# ax[0].legend()

plt.figure(figsize=(8, 6), dpi=150)
plt.plot(BEM_rR[0, :], BEM_a[0, :], '--', color=c[2], label="$a$ BEM")
plt.plot(data[2][:nspan - 1], data[0][:nspan - 1], '-', color=c[2], label="$a$ LLM")
plt.plot(BEM_rR[0, :], BEM_aline[0, :], '--', color=c[3], label="$a'$ BEM")
plt.plot(data[2][:nspan - 1], data[1][:nspan - 1], '-', color=c[3], label="$a'$ LLM")
plt.xlabel('Radial location r/R (-)', fontsize=15)
plt.ylabel('Induction factor a (-)', fontsize=15)
plt.title(r'Radial distribution of $a$ and $a^\prime$', fontsize=16)
plt.legend(fontsize=15)
plt.grid(True)
# plt.savefig('BEMcomp_figures/alpha_phi.eps', bbox_inches='tight', format='eps')
plt.savefig('BEMcomp_figures/inductionfactors.eps', format='eps')

# RADIAL DISTRIBUTION CT

print('CT Carlos:', np.sum(CT_LLM[:nspan-1]))
print('CT a:', CT_LLM2)

# plt.figure()
# plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_CT, 0), BEM_rR.shape)[0, :], '-r', label=r'$C_T$ BEM')
# plt.plot(data[2][:nspan-1, 0], CT_LLM[:nspan-1], '--r', label=r'$C_T$ LLM Carlos')
# plt.plot(data[2][:nspan-1, 0], CT_LLM2[:nspan-1], '--g', label=r'$C_T$ LLM 2')
# plt.xlabel('r/R (-)')
# plt.ylabel(r'$C_T$ (-)')
# plt.legend()
# plt.grid(True)

# Radial distribution CP

print('CP Carlos:', np.sum(CP_LLM[:nspan-1]))
print('CP a:', CP_LLM2)

# plt.figure()
# plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_CP, 0), BEM_rR.shape)[0, :], '-r', label=r'$C_P$ BEM')
# plt.plot(data[2][:nspan-2, 0], CP_LLM[:nspan-2], '--r', label=r'$C_P$ LLM Carlos')
# plt.plot(data[2][:nspan-1, 0], CP_LLM2[:nspan-1], '--g', label=r'$C_P$ LLM 2')
# plt.xlabel('r/R (-)')
# plt.ylabel('$C_P$ (-)')
# plt.legend()
# plt.grid(True)