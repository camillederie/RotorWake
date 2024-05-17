import matplotlib.pyplot as plt
from Variables import *
from Main import rotor_wake_system

def plot_blade_geometry():
# Extract control points, vortex rings, and blade panels from rotor_wake_system
    controlpoints = rotor_wake_system['controlpoints']
    rings = rotor_wake_system['rings']
    bladepanels = rotor_wake_system['bladepanels']

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

#Run the plot functions
plot_blade_geometry()
plt.show()