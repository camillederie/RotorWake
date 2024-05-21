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
    plt.show()

#Run the plot functions

#plt.show()

def plot_results(results, indeces_b1, indeces_b2, indeces_b3):
    
    # plot gamma vs radial position
    plt.figure(figsize=(10, 5))
    plt.plot(results[6][indeces_b1], results[2][indeces_b1], label='Blade 1')
    plt.plot(results[6][indeces_b2], results[2][indeces_b2], label='Blade 2')
    plt.plot(results[6][indeces_b3], results[2][indeces_b3], label='Blade 3')
    plt.title('Gamma vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Gamma')
    plt.legend()
    plt.grid()
    plt.show()
    # plot alpha vs radial position and phi vs radial position

    fig, ax = plt.subplots(2, 1, figsize=(10, 10))
    ax[0].plot(results[6][indeces_b1], results[3][indeces_b1], label= 'Blade 1')
    ax[0].plot(results[6][indeces_b2], results[3][indeces_b2], label= 'Blade 2')
    ax[0].plot(results[6][indeces_b3], results[3][indeces_b3], label= 'Blade 3')

    ax[0].set_title('Alpha vs radial position') 
    ax[0].set_xlabel('Radial position')
    ax[0].set_ylabel('Alpha')
    ax[0].legend()
    ax[0].grid()

    ax[1].plot(results[6][indeces_b1], results[4][indeces_b1], label= 'Blade 1')
    ax[1].plot(results[6][indeces_b2], results[4][indeces_b2], label= 'Blade 2')
    ax[1].plot(results[6][indeces_b3], results[4][indeces_b3], label= 'Blade 3')
    ax[1].set_title('Phi vs radial position')
    ax[1].set_xlabel('Radial position')
    ax[1].set_ylabel('Phi')
    ax[1].legend()
    ax[1].grid()
    plt.show()

    plt.figure(figsize=(10, 5))
    #Plot induction factor vs radial position
    plt.plot(results[6][indeces_b1], results[7][indeces_b1], label='Blade 1')
    plt.plot(results[6][indeces_b2], results[7][indeces_b2], label='Blade 2')
    plt.plot(results[6][indeces_b3], results[7][indeces_b3], label='Blade 3')
    plt.title('Induction factor vs radial position')
    plt.xlabel('Radial position')
    plt.ylabel('Induction factor')
    plt.legend()
    plt.grid()
    plt.show()
    



