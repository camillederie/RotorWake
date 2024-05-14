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
fig_circ, ax_circ = plt.subplots(1, 2, dpi=150)
fig_ind, ax_ind = plt.subplots(1, 2, dpi=150)
fig_aoa, ax_aoa = plt.subplots(1, 2, dpi=150)
fig_rotor, ax_rotor = plt.subplots(1, 1, dpi=150)
nspan = 20
ntheta = 200
nblades = 3
spacing = 'equal'
nrotor = 2

CP = []
CT = []

for color, phase_shift in zip(colors, phase_shifts):


    prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, phase_diff=phase_shift,
                             double_rotor=True,
                             n_theta=ntheta, spacing=spacing, a=0.33, xshift=0, yshift=100, zshift=0)
    blade = prop_geo.bladepanels
    rings = prop_geo.filaments
    solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6,
                               n_iter=100)
    data_double = solver.run_solver()

    # non-dimensionalise data

    # forces [3, 4]
    data_double[3] = data_double[3] / (0.5 * solver.u_inf ** 2 * solver.r_rotor)
    data_double[4] = data_double[4] / (0.5 * solver.u_inf ** 2 * solver.r_rotor)
    # circulation [5]

    omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius
    CPandCT = solver.CP_and_CT(np.resize(data_double[0], data_double[2].shape),
                                        np.resize(data_double[1], data_double[2].shape),
                                        data_double[2],
                                        np.resize(data_double[3], data_double[2].shape),
                                        np.resize(data_double[4], data_double[2].shape),
                                        solver.geo.v_inf, omega, solver.geo.radius, nblades)

    # =============================================================================
    # Double Rotor Plotting
    # =============================================================================

    def plotDoubleRotor():
        prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
                                 spacing='cosine', a=0, xshift=0, yshift=100, zshift=0, phase_diff=120,
                                 double_rotor=True)
        prop_geo.doubleRotor()
        rings = prop_geo.filaments

        fig = plt.figure(dpi=150)
        ax = Axes3D(fig)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        # ax.set_xlim3d(0, 50)
        # ax.set_ylim3d(-50, 200)
        # ax.set_zlim3d(-50, 50)

        c = ['green', 'blue', 'red', 'green', 'blue', 'red']
        for idx in range(nblades * nrotor):
            ax.plot_wireframe(rings[0, idx * (nspan - 1):(idx + 1) * (nspan - 1), :ntheta + 1],
                              rings[1, idx * (nspan - 1):(idx + 1) * (nspan - 1), :ntheta + 1],
                              rings[2, idx * (nspan - 1):(idx + 1) * (nspan - 1), :ntheta + 1],
                              color=c[idx], cstride=0)

    plotDoubleRotor()


    # =============================================================================
    # Rotor performance plots
    # =============================================================================
    def plottingFunction(solver, prop_geo, data):
        plt.plot(prop_geo.centerPoints, data[-1][:len(prop_geo.centerPoints)])
        plt.show()
        return


    # plottingFunction(solver,prop_geo,data)

    # =============================================================================
    # DOUBLE ROTOR RESULTS
    # =============================================================================
    plt.close('All')
    r_R = data_double[2][:nspan-1]

    [BEM_rR, BEM_alpha, BEM_phi, BEM_rho, BEM_Ax, BEM_Az, BEM_Gamma , BEM_CT, BEM_CP, BEM_a,
     BEM_aline, BEM_vinf, BEM_radius] = read_matlab_data()
    BEM_Ax = BEM_Ax / (0.5 * 10**2 * 50)
    BEM_Az = BEM_Az / (0.5 * 10**2 * 50)

    ax_aoa[0].plot(r_R, np.degrees(data_double[6][:nspan-1]), linestyle='-', c=color,
                   label=r"$\gamma$: {}$^\circ$".format(phase_shift))
    ax_aoa[1].plot(r_R, np.degrees(data_double[6][3*(nspan-1):4*(nspan-1)]), linestyle='--', c=color,
                                   label=r"$\gamma$: {}$^\circ$".format(phase_shift))

    # Radial distribution circulationd

    # made non-dimensional with (np.pi * Uinf**2) / (NBlades*Omega)
    circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)


    # ax_circ.plot(data_single[2][:nspan - 1, 0], np.resize(data_single[5], data_single[2].shape)[:nspan - 1,
    #                                         0] / circ_nondim)
    ax_circ[0].plot(r_R, data_double[5][:nspan - 1] / circ_nondim,
                    linestyle='-', c=color, label=r"$\gamma$: {}$^\circ$".format(phase_shift))
    ax_circ[1].plot(r_R, data_double[5][3*(nspan-1):4*(nspan-1)] / circ_nondim, linestyle='--', c=color)


    # RADIAL DISTRIBUTION INDUCTION FACTORS

    # INDUCTION FACTORS
    ax_ind[0].plot(r_R, data_double[0][:nspan - 1],
                   linestyle='-', c=color, label=r"$\gamma$: {}$^\circ$".format(phase_shift))
    ax_ind[1].plot(r_R, data_double[0][3*(nspan-1):4*(nspan-1)], linestyle='--', c=color)
    # ax_ind[0].plot(data_single[2][-(nspan - 1):], data_single[0][-(nspan - 1):], ':r')
    # ax_ind[1].plot(data_double[2][:nspan - 1], data_double[1][:nspan - 1], linestyle='-')
    # ax_ind[1].plot(data_double[2][-(nspan - 1):], data_double[1][-(nspan - 1):], linestyle='--')
    # ax_ind[1].plot(data_single[2][:nspan - 1], data_single[1][:nspan - 1], ':r', label="LLM - single rotor")

    # CP and CT
    CT.append(CPandCT[-1])
    CP.append(CPandCT[-2])

ax_circ[0].set_xlabel('r/R (-)', fontsize=14)
ax_circ[0].set_ylabel(r'$\Gamma$ (-)', fontsize=14)
ax_circ[0].legend(prop={"size": 10})
ax_circ[0].grid(True)
ax_circ[1].set_xlabel('r/R (-)', fontsize=14)
ax_circ[1].set_ylabel(r'$\Gamma$ (-)', fontsize=14)
ax_circ[1].grid(True)

ax_ind[0].set_xlabel("r/R (-)", fontsize=14)
ax_ind[0].set_ylabel("a (-)", fontsize=14)
ax_ind[1].set_xlabel("r/R (-)", fontsize=14)
ax_ind[1].set_ylabel("a' (-)", fontsize=14)
ax_ind[0].legend(prop={"size": 10})
ax_ind[0].grid()
ax_ind[1].grid()

ax_aoa[0].set_xlabel('r/R (-)', fontsize=14)
ax_aoa[0].set_ylabel(r'$\alpha$ ($^\circ$)', fontsize=14)
ax_aoa[0].legend(prop={"size": 10})
ax_aoa[0].grid(True)
ax_aoa[1].set_xlabel('r/R (-)', fontsize=14)
ax_aoa[1].set_ylabel(r'$\alpha$ ($^\circ$)', fontsize=14)
ax_aoa[1].grid(True)
ax_aoa.set_title("Blue")