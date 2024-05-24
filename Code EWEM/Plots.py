import matplotlib.pyplot as plt
#from Main import system_geom
from Variables import *


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
            ax.plot([x1, x2], [y1, y2], [z1, z2], color='blue')

    # Plot blade panels
    for panel in bladepanels:
        p1 = panel['p1']
        p2 = panel['p2']
        p3 = panel['p3']
        p4 = panel['p4']
        xs = [p1[0], p2[0], p3[0], p4[0], p1[0]]
        ys = [p1[1], p2[1], p3[1], p4[1], p1[1]]
        zs = [p1[2], p2[2], p3[2], p4[2], p1[2]]
        ax.plot(xs, ys, zs, color='green')

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Blade Geometry')
    #plt.show()

#Run the plot functions

#plt.show()

def plot_results(results):
    
    # Extract the results from the dictionary
    results_TSR_6 = results[f"TSR_{TSR_list[0]}"][0]
    results_TSR_8 = results[f"TSR_{TSR_list[1]}"][0]
    results_TSR_10 = results[f"TSR_{TSR_list[2]}"][0]
    # print(results_TSR_6)
    indeces_b1 = results[f"TSR_{TSR_list[0]}"][7]
    # Plot 1: Gamma vs radial position
    plt.figure(figsize=(10, 10))
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[9][indeces_b1], label='TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[9][indeces_b1], label='TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[9][indeces_b1], label='TSR 10')
    # plt.title('Gamma vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Gamma')
    plt.legend()
    plt.grid()
    # plt.axis('square')
    plt.savefig('Gamma_vs_radial_position.png')
    plt.close()

    # Plot 2: Alpha vs radial position
    plt.figure(figsize=(10, 10))
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[3][indeces_b1], label='TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[3][indeces_b1], label='TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[3][indeces_b1], label='TSR 10')
    # ax1.set_title('Alpha vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Alpha')
    plt.legend()
    plt.grid()
    # ax1.axis('square')
    plt.savefig('Alpha_vs_radial_position.png')
    plt.close()

    # Plot 3: Phi vs radial position
    plt.figure(figsize=(10, 10))
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[4][indeces_b1], label='TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[4][indeces_b1], label='TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[4][indeces_b1], label='TSR 10')
    # ax.set_title('Phi vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Phi')
    plt.legend()
    plt.grid()
    plt.savefig('Phi_vs_radial_position.png')
    plt.close()

    # Plot 4: Induction factor vs radial position
    plt.figure(figsize=(10, 10))
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[7][indeces_b1], label='TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[7][indeces_b1], label='TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[7][indeces_b1], label='TSR 10')
    # plt.title('Induction factor vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Induction factor')
    plt.legend()
    plt.grid()
    plt.savefig('Induction_factor_vs_radial_position.png')
    #plt.close()

    # Plot 5: Normal load coefficient vs radial position
    plt.figure(figsize=(10, 10))
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
    plt.figure(figsize=(10, 10))
    plt.plot(results_TSR_6[6][indeces_b1], results_TSR_6[10][indeces_b1], label='TSR 6')
    plt.plot(results_TSR_8[6][indeces_b1], results_TSR_8[10][indeces_b1], label='TSR 8')
    plt.plot(results_TSR_10[6][indeces_b1], results_TSR_10[10][indeces_b1], label='TSR 10')
    # plt.title('Tangential load coefficient vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Tangential load coefficient')
    plt.legend()
    plt.grid()
    plt.savefig('Tangential_load_coefficient_vs_radial_position.png')

    
