"""
Sensitivity Plots
"""
from create_geometry import BladeGeometry
from create_geometry_varwake import BladeGeometry2
from lifting_line_solver import LiftingLineSolver
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from read_BEMdata_into_Python import read_matlab_data
import numpy as np

plot_Spacing = False
plot_DiscrAzimuthal = False
plot_ConvecSpeed = False
plot_ConvecSpeed2 = False
plot_WakeLength = False

[BEM_rR, BEM_alpha, BEM_phi, BEM_rho, BEM_Ax, BEM_Az, BEM_Gamma,
 BEM_CT, BEM_CP, BEM_a, BEM_aline, BEM_vinf, BEM_radius] = read_matlab_data()

c=["lawngreen", "deepskyblue", "orangered", "darkviolet"]
# =============================================================================
#  Spacing Discretization
# =============================================================================
if plot_Spacing:
    nspan = 20
    ntheta = 200
    nblades = 3
    spacingList = ['equal','cosine']
    c=["lawngreen", "deepskyblue", "orangered", "darkviolet"]
    nrotor = 2
    
    plt.figure(figsize=(8, 6), dpi=150)
    idx=0
    plt.plot(BEM_rR[0, :], BEM_Gamma[0, :], '--',c=c[idx], label=r'$\Gamma$ BEM equal spacing')
    for spacing in spacingList:
        idx+=1
        prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
            spacing=spacing, a=0, xshift=0, yshift=100, zshift=0, phase_diff=0, double_rotor=False)
        solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6,
                                   n_iter=200)
        data = solver.run_solver()
        omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius    
        circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)
        
        plt.plot(data[2][:nspan - 1, 0], np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim, label='Spacing: '+spacing,c=c[idx])
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    [CP_LLM, CT_LLM, CP_LLM2Spacing, CT_LLM2Spacing] = solver.CP_and_CT(np.resize(data[0], data[2].shape), np.resize(data[1], data[2].shape), data[2],
                                    np.resize(data[3], data[2].shape), np.resize(data[4], data[2].shape),
                                    solver.geo.v_inf, omega, solver.geo.radius, nblades)
    plt.savefig('sensitivity_spacing.pdf')
# =============================================================================
#  Wake Discretization
# =============================================================================
if plot_DiscrAzimuthal:
    nspan = 20
    nthetaList = [60,120,240]
    nblades = 3
    spacing = 'equal'

    nrotor = 2
    c=["lawngreen", "deepskyblue", "orangered", "darkviolet"]
    
    plt.figure(figsize=(8, 6), dpi=150)
    idx=0
    plt.plot(BEM_rR[0, :], BEM_Gamma[0, :], '--',c=c[idx], label=r'$\Gamma$ BEM')
    for ntheta in nthetaList:
        idx+=1
        prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
            spacing=spacing, a=0, xshift=0, yshift=100, zshift=0, phase_diff=0, double_rotor=False)
        solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6,
                                   n_iter=200)

        data = solver.run_solver()
        omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius
        circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)
        

        plt.plot(data[2][:nspan - 1, 0], np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim, label='$Segments/rot$= '+str(int(ntheta/15)),c=c[idx])
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    [CP_LLM, CT_LLM, CP_LLM2Azi, CT_LLM2Azi] = solver.CP_and_CT(np.resize(data[0], data[2].shape), np.resize(data[1], data[2].shape), data[2],
                                    np.resize(data[3], data[2].shape), np.resize(data[4], data[2].shape),
                                    solver.geo.v_inf, omega, solver.geo.radius, nblades)
    plt.savefig('sensitivity_azimuthal.pdf')
# =============================================================================
# ConvecSpeed
# =============================================================================
if plot_ConvecSpeed:

    nspan = 20
    ntheta = 200
    nblades = 3
    spacing = 'equal'

    tsr_list = [4, 6, 8]
    circ_list = []
    Fax_list = []
    Faz_list = []

    for tsr in tsr_list:

        prop_geo = BladeGeometry(radius=50.0, tsr=tsr, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
                                 spacing=spacing, a=0, xshift=0, yshift=100, zshift=0, phase_diff=0, double_rotor=False)
        # prop_geo.doubleRotor()
        solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6, n_iter=200)

        data = solver.run_solver()
        omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius

        [CP_LLM, CT_LLM, CP_LLM2, CT_LLM2] = solver.CP_and_CT(np.resize(data[0], data[2].shape),
                                                              np.resize(data[1], data[2].shape), data[2],
                                                              np.resize(data[3], data[2].shape),
                                                              np.resize(data[4], data[2].shape),
                                                              solver.geo.v_inf, omega, solver.geo.radius, nblades)
        print('-------------')
        print(CP_LLM2)
        print(CT_LLM2)

        circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)
        circ_list.append(np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim)

        F_nondim_LLM = 0.5 * (solver.geo.v_inf ** 2) * solver.geo.radius
        Fax_list.append(np.resize(data[3], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)
        Faz_list.append(np.resize(data[4], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], circ_list[0], '-', color=c[0], label='$\lambda$ = ' + str(tsr_list[0]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[1], '-', color=c[1], label='$\lambda$ = ' + str(tsr_list[1]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[2], '-', color=c[2], label='$\lambda$ = ' + str(tsr_list[2]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/convec_circ.eps', format='eps')

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], Fax_list[0], '-', color=c[0], label=r'$F_{ax}$ | $\lambda$ = ' + str(tsr_list[0]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[1], '-', color=c[1], label=r'$F_{ax}$ | $\lambda$ = ' + str(tsr_list[1]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[2], '-', color=c[2], label=r'$F_{ax}$ | $\lambda$ = ' + str(tsr_list[2]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[0], '--', color=c[0], label=r'$F_{az}$ | $\lambda$ = ' + str(tsr_list[0]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[1], '--', color=c[1], label=r'$F_{az}$ | $\lambda$ = ' + str(tsr_list[1]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[2], '--', color=c[2], label=r'$F_{az}$ | $\lambda$ = ' + str(tsr_list[2]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel('Force $F$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $F_{ax}$ and $F_{az}$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/convec_forces.eps', format='eps')

if plot_ConvecSpeed2:

    nspan = 20
    ntheta = 200
    nblades = 3
    spacing = 'equal'

    a_list = [0.2, 0.25, 0.3]
    circ_list = []
    Fax_list = []
    Faz_list = []

    for adf in a_list:

        prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
                                 spacing=spacing, a=adf, xshift=0, yshift=100, zshift=0, phase_diff=0, double_rotor=False)
        # prop_geo.doubleRotor()
        solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, afix=adf, weight=0.5, tol=1e-6, n_iter=200)

        data = solver.run_solver()
        omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius

        [CP_LLM, CT_LLM, CP_LLM2, CT_LLM2] = solver.CP_and_CT(np.resize(data[0], data[2].shape),
                                                              np.resize(data[1], data[2].shape), data[2],
                                                              np.resize(data[3], data[2].shape),
                                                              np.resize(data[4], data[2].shape),
                                                              solver.geo.v_inf, omega, solver.geo.radius, nblades)
        print('-------------')
        print(CP_LLM2)
        print(CT_LLM2)

        circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)
        circ_list.append(np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim)

        F_nondim_LLM = 0.5 * (solver.geo.v_inf ** 2) * solver.geo.radius
        Fax_list.append(np.resize(data[3], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)
        Faz_list.append(np.resize(data[4], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], circ_list[0], '-', color=c[0], label='$a$ = ' + str(a_list[0]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[1], '-', color=c[1], label='$a$ = ' + str(a_list[1]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[2], '-', color=c[2], label='$a$ = ' + str(a_list[2]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/convec_circ_a.eps', format='eps')

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], Fax_list[0], '-', color=c[0], label=r'$F_{ax}$ | $a$ = ' + str(a_list[0]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[1], '-', color=c[1], label=r'$F_{ax}$ | $a$ = ' + str(a_list[1]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[2], '-', color=c[2], label=r'$F_{ax}$ | $a$ = ' + str(a_list[2]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[0], '--', color=c[0], label=r'$F_{az}$ | $a$ = ' + str(a_list[0]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[1], '--', color=c[1], label=r'$F_{az}$ | $a$ = ' + str(a_list[1]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[2], '--', color=c[2], label=r'$F_{az}$ | $a$ = ' + str(a_list[2]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel('Force $F$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $F_{ax}$ and $F_{az}$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/convec_forces_a.eps', format='eps')

if plot_WakeLength:

    nspan = 20
    ntheta = 200
    nblades = 3
    spacing = 'equal'

    N_rotations = [2, 4, 8, 12]
    circ_list = []
    Fax_list = []
    Faz_list = []
    error_list = []

    fig = plt.figure(dpi=150)
    i = 0

    for rotation in N_rotations:

        prop_geo = BladeGeometry2(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
                                  spacing=spacing, a=0, n_rotations=rotation, xshift=0, yshift=100, zshift=0,
                                  phase_diff=0, double_rotor=False)

        blade = prop_geo.bladepanels
        rings = prop_geo.filaments

        solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6, n_iter=200)

        data = solver.run_solver()
        omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius

        [CP_LLM, CT_LLM, CP_LLM2, CT_LLM2] = solver.CP_and_CT(np.resize(data[0], data[2].shape),
                                                              np.resize(data[1], data[2].shape), data[2],
                                                              np.resize(data[3], data[2].shape),
                                                              np.resize(data[4], data[2].shape),
                                                              solver.geo.v_inf, omega, solver.geo.radius, nblades)
        print('-------------')
        print(CP_LLM2)
        print(CT_LLM2)

        circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)
        circ_list.append(np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim)

        F_nondim_LLM = 0.5 * (solver.geo.v_inf ** 2) * solver.geo.radius
        Fax_list.append(np.resize(data[3], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)
        Faz_list.append(np.resize(data[4], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)
        error_list.append(data[8])

        plot_wake = False

        if plot_wake:

            i += 1

            ax = fig.add_subplot(1, 3, i, projection="3d")

            # =============================================================================
            # WAKE VIZ
            # =============================================================================

            ax.plot_wireframe(rings[0, :nspan - 1, :ntheta + 1]/(40),
                               rings[1, :nspan - 1, :ntheta + 1],
                               rings[2, :nspan - 1, :ntheta + 1],
                               color='green', ccount=0, rcount=5)
            ax.plot_wireframe(rings[0, nspan - 1:2 * nspan - 2, :ntheta + 1]/(40),
                               rings[1, nspan - 1:2 * nspan - 2, :ntheta + 1],
                               rings[2, nspan - 1:2 * nspan - 2, :ntheta + 1],
                               color='blue', ccount=0, rcount=5)
            ax.plot_wireframe(rings[0, (2 * nspan - 2):, :ntheta + 1]/(40),
                               rings[1, (2 * nspan - 2):, :ntheta + 1],
                               rings[2, (2 * nspan - 2):, :ntheta + 1],
                               color='red', ccount=0, rcount=5)

            # split blade panels
            blades_split = np.split(blade, 3)

            for blade_i in (blades_split):
                x = blade_i[:, [0, 9]]
                y = blade_i[:, [1, 10]]
                z = blade_i[:, [2, 11]]

                ax.plot_surface(x, y, z, color='k')

            # =============================================================================
            # Draw rotor hub
            # =============================================================================
            theta = np.linspace(0, 2 * np.pi)
            y = solver.geo.span_arr[0] * np.cos(theta) * solver.r_rotor
            z = solver.geo.span_arr[0] * np.sin(theta) * solver.r_rotor

            ax.plot(y * 0, y, z, color='k')

            ax.view_init(0, 180)
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.set_xlim(0, N_rotations[-1])
            ax.set_xlabel(r"Wakelength $\frac{L}{2 \pi}$ (-)")
            ax.set_title('{} Rotations'.format(rotation))
            ax.view_init(elev=10, azim=-450)

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], circ_list[0], '-', color=c[0], label='# Rotations = ' + str(N_rotations[0]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[1], '-', color=c[1], label='# Rotations = ' + str(N_rotations[1]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[2], '-', color=c[2], label='# Rotations = ' + str(N_rotations[2]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[3], '-', color=c[3], label='# Rotations = ' + str(N_rotations[3]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/varwake_circ.eps', format='eps')

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], Fax_list[0], '-', color=c[0], label=r'$F_{ax}$ | # Rot = ' + str(N_rotations[0]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[1], '-', color=c[1], label=r'$F_{ax}$ | # Rot = ' + str(N_rotations[1]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[2], '-', color=c[2], label=r'$F_{ax}$ | # Rot = ' + str(N_rotations[2]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[3], '-', color=c[3], label=r'$F_{ax}$ | # Rot = ' + str(N_rotations[3]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[0], '--', color=c[0], label=r'$F_{az}$ | # Rot = ' + str(N_rotations[0]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[1], '--', color=c[1], label=r'$F_{az}$ | # Rot = ' + str(N_rotations[1]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[2], '--', color=c[2], label=r'$F_{az}$ | # Rot = ' + str(N_rotations[2]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[3], '--', color=c[3], label=r'$F_{az}$ | # Rot = ' + str(N_rotations[3]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel('Force $F$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $F_{ax}$ and $F_{az}$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/varwake_forces.eps', format='eps')

    plt.figure(figsize=(8, 6), dpi=150)
    plt.semilogy(np.arange(1, len(error_list[0]) + 1), error_list[0], '-', color=c[0], label='# Rotations = ' + str(N_rotations[0]))
    plt.semilogy(np.arange(1, len(error_list[1]) + 1), error_list[1], '-', color=c[1], label='# Rotations = ' + str(N_rotations[1]))
    plt.semilogy(np.arange(1, len(error_list[2]) + 1), error_list[2], '-', color=c[2], label='# Rotations = ' + str(N_rotations[2]))
    plt.semilogy(np.arange(1, len(error_list[3]) + 1), error_list[3], '-', color=c[3], label='# Rotations = ' + str(N_rotations[3]))
    # plt.plot(np.arange(1, len(error_list[0]) + 1), error_list[0], '-r', label='# Rotations = ' + str(N_rotations[0]))
    # plt.plot(np.arange(1, len(error_list[1]) + 1), error_list[1], '-b', label='# Rotations = ' + str(N_rotations[1]))
    # plt.plot(np.arange(1, len(error_list[2]) + 1), error_list[2], '-g', label='# Rotations = ' + str(N_rotations[2]))
    plt.xlabel('Iterations N', fontsize=15)
    plt.ylabel(r'Log error $\epsilon$ (-)', fontsize=15)
    plt.title('Error', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/varwake_convergence.eps', format='eps')

plt.show()