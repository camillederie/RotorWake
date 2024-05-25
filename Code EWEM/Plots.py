import matplotlib.pyplot as plt
#from Main import system_geom
from Variables import *
import numpy as np

def read_array_from_file(filename):
    #For every TSR, read the results from the file consisting of 3 collumns (one for TSR 6, 8 and 10) comma delimited, and return the results as a list of lists
    results = []
    with open(filename, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            results.append([float(value) for value in values])
    TSR_6 = [item[0] for item in results]
    TSR_8 = [item[1] for item in results]
    TSR_10 = [item[2] for item in results]
    return TSR_6, TSR_8, TSR_10

def plot_blade_geometry(system_geom):
    # Extract control points, vortex rings, and blade panels from rotor_wake_system
    controlpoints = system_geom['controlpoints']
    rings = system_geom['rings']
    bladepanels = system_geom['bladepanels']

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot control points
    for cp in controlpoints:
        x = cp['coordinates'][0]
        y = cp['coordinates'][1]
        z = cp['coordinates'][2]
        ax.scatter(x, y, z, color='red')

    # Plot vortex rings
    for ring in rings:
        filaments = ring['filaments']
        for filament in filaments:
            x1 = filament['x1']
            y1 = filament['y1']
            z1 = filament['z1']
            x2 = filament['x2']
            y2 = filament['y2']
            z2 = filament['z2']
            ax.plot([x1, x2], [y1, y2], [z1, z2], color='blue', linewidth=0.5)  # Make the lines thinner

    # Plot blade panels
    for panel in bladepanels:
        p1 = panel['p1']
        p2 = panel['p2']
        p3 = panel['p3']
        p4 = panel['p4']
        xs = [p1[0], p2[0], p3[0], p4[0], p1[0]]
        ys = [p1[1], p2[1], p3[1], p4[1], p1[1]]
        zs = [p1[2], p2[2], p3[2], p4[2], p1[2]]
        ax.plot(xs, ys, zs, color='green', linewidth=0.5)  # Make the lines thinner

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Blade Geometry')
    #plt.show()

#Run the plot functions

#plt.show()

def plot_results(results):
    r_BEM = np.linspace(0.2,1,401)
    # Extract the results from the dictionary
    results_TSR_6 = results[f"TSR_{TSR_list[0]}"][0]
    results_TSR_8 = results[f"TSR_{TSR_list[1]}"][0]
    results_TSR_10 = results[f"TSR_{TSR_list[2]}"][0]
    # print(results_TSR_6)
    indeces_b1 = results[f"TSR_{TSR_list[0]}"][7]
    # Plot 1: Gamma vs radial position
    plt.figure()
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[9][indeces_b1], label='TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[9][indeces_b1], label='TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[9][indeces_b1], label='TSR 10')
    # plt.title('Gamma vs radial position')
    plt.xlabel('Radial position [-]')
    plt.ylabel('Gamma [-]')
    plt.legend()
    plt.grid()
    # plt.axis('square')
    plt.savefig('Gamma_vs_radial_position.png')
    plt.close()

    # Plot 2: Alpha vs radial position
    alpha_TSR_6, alpha_TSR_8, alpha_TSR_10 = read_array_from_file('Arrays/Alpha_no_Prandtl.txt')
    plt.figure()
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[3][indeces_b1], color = 'tab:blue', label='LLM TSR 6')
    plt.plot(r_BEM, alpha_TSR_6, linestyle = 'dashed', color = 'tab:blue', label='BEM TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[3][indeces_b1], color = 'tab:orange', label='LLM TSR 8')
    plt.plot(r_BEM, alpha_TSR_8, linestyle = 'dashed', color = 'tab:orange', label='BEM TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[3][indeces_b1], color = 'tab:green', label='LLM TSR 10')
    plt.plot(r_BEM, alpha_TSR_10, color = 'tab:green', linestyle = 'dashed', label='BEM TSR 10')

    # ax1.set_title('Alpha vs radial position')
    plt.xlabel('Radial position [-]')
    plt.ylabel('Alpha [deg]')
    plt.legend()
    plt.grid()
    # ax1.axis('square')
    plt.savefig('Alpha_vs_radial_position.png')
    plt.close()

    # Plot 3: Phi vs radial position
    phi_TSR_6, phi_TSR_8, phi_TSR_10 = read_array_from_file('Arrays/Flow angle_no_Prandtl.txt')
    plt.figure()
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[4][indeces_b1], color = 'tab:blue', label='LLM TSR 6')
    plt.plot(r_BEM, phi_TSR_6, linestyle = 'dashed', color = 'tab:blue', label='BEM TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[4][indeces_b1], color = 'tab:orange', label='LLM TSR 8')
    plt.plot(r_BEM, phi_TSR_8, linestyle = 'dashed', color = 'tab:orange', label='BEM TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[4][indeces_b1], color = 'tab:green', label='LLM TSR 10')
    plt.plot(r_BEM, phi_TSR_10, color = 'tab:green', linestyle = 'dashed', label='BEM TSR 10')
    # ax.set_title('Phi vs radial position')
    plt.xlabel('Radial position [-]')
    plt.ylabel('Phi [deg]')
    plt.legend()
    plt.grid()
    plt.savefig('Phi_vs_radial_position.png')
    plt.close()

    # Plot 4: Induction factor vs radial position
    a_TSR_6, a_TSR_8, a_TSR_10 = read_array_from_file('Arrays/Axial_no_Prandtl.txt')
    plt.figure()
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[7][indeces_b1], color = 'tab:blue', label='LLM TSR 6')
    plt.plot(r_BEM, a_TSR_6, linestyle = 'dashed', color = 'tab:blue', label='BEM TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[7][indeces_b1], color = 'tab:orange', label='LLM TSR 8')
    plt.plot(r_BEM, a_TSR_8, linestyle = 'dashed', color = 'tab:orange', label='BEM TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[7][indeces_b1],color = 'tab:green', label='LLM TSR 10')
    plt.plot(r_BEM, a_TSR_10, linestyle = 'dashed', color = 'tab:green', label='BEM TSR 10')
    # plt.title('Induction factor vs radial position')
    plt.xlabel('Radial position [-]')
    plt.ylabel('Induction factor [-]')
    plt.legend()
    plt.grid()
    plt.savefig('Axial_Induction_factor_vs_radial_position.png')
    #plt.close()

    ap_TSR_6, ap_TSR_8, ap_TSR_10 = read_array_from_file('Arrays/Azimuthal_no_Prandtl.txt')
    plt.figure()
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[8][indeces_b1], color = 'tab:blue', label='LLM TSR 6')
    plt.plot(r_BEM, ap_TSR_6, linestyle = 'dashed', color = 'tab:blue', label='BEM TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[8][indeces_b1], color = 'tab:orange', label='LLM TSR 8')
    plt.plot(r_BEM, ap_TSR_8, linestyle = 'dashed', color = 'tab:orange', label='BEM TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[8][indeces_b1],color = 'tab:green', label='LLM TSR 10')
    plt.plot(r_BEM, ap_TSR_10, linestyle = 'dashed', color = 'tab:green', label='BEM TSR 10')
    # plt.title('Induction factor vs radial position')
    plt.xlabel('Radial position [-]')
    plt.ylabel('Azimuthal Induction factor [-]')
    plt.legend()
    plt.grid()
    plt.savefig('Azimuthal_Induction_factor_vs_radial_position.png')

    # Plot 5: Normal load coefficient vs radial position
    plt.figure()
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[11][indeces_b1], label='TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[11][indeces_b1], label='TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[11][indeces_b1], label='TSR 10')
    # plt.title('Normal load coefficient vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Normal load coefficient')
    plt.legend()
    plt.grid()
    plt.savefig('Normal_load_coefficient_vs_radial_position.png')

    # Plot 6: Tangential load coefficient vs radial position
    plt.figure()
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[10][indeces_b1], label='TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[10][indeces_b1], label='TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[10][indeces_b1], label='TSR 10')
    # plt.title('Tangential load coefficient vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Tangential load coefficient')
    plt.legend()
    plt.grid()
    plt.savefig('Tangential_load_coefficient_vs_radial_position.png')

    
