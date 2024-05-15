#3D induced velocities, blade loads, run solver
import numpy as np


def LiftingLineSolver(system_geom, V_inf, Omega, R):
    # Inputs
    relax = 0.3
    n_iterations = 1200
    # system_geom: Contains the geometry of horseshoe vortex rings and control points at the blade

    control_points = system_geom['control_points']
    rings = system_geom['rings']
    
    matrix_u = []
    matrix_v = []
    matrix_w = []
    
    gamma_updated = np.zeros(len(control_points))
    gamma = np.zeros(len(control_points))

    for i in range(len(control_points)):
        matrix_u.append([])
        matrix_v.append([])
        matrix_w.append([])
        
        for j in range(len(rings)):
            for k in range(len(rings['filaments'])):
                rings[i]['filaments'][j]['gamma'] = 1 # Set ring strength to unity for calculating induced velocity at control point p
                v_ind = np.zeros(3)
                for r in rings: 
                    for f in r['filaments']:
                        gamma = f['gamma']
                        x1, y1, z1 = f['x1'],f['y1'],f['z1']	# Start point of filament
                        x2, y2, z2 = f['x2'],f['y2'],f['z2']
                        xP, yP, zP = control_points[i]['x'],control_points[i]['y'],control_points[i]['z']
                        r1 = np.sqrt((xP - x1)**2 + (yP - y1)**2 + (zP - z1)**2)
                        r2 = np.sqrt((xP - x2)**2 + (yP - y2)**2 + (zP - z2)**2)
                        r1xr2_x = (yP - y1) * (zP - z2) - (zP - z1) * (yP - y2)
                        r1xr2_y = -(xP - x1) * (zP - z2) + (zP - z1) * (xP - x2)
                        r1xr2_z = (xP - x1) * (yP - y2) - (yP - y1) * (xP - x2)
                        r1xr_sqr = r1xr2_x**2 + r1xr2_y**2 + r1xr2_z**2
                        if r1xr_sqr < 0.0001:
                            r1xr_sqr = 0.0001
                        if r1 < 0.0001:
                            r1 = 0.0001
                        if r2 < 0.0001:
                            r2 = 0.0001
                        r0r1 = (x2 - x1) * (xP - x1) + (y2 - y1) * (yP - y1) + (z2 - z1) * (zP - z1)
                        r0r2 = (x2 - x1) * (xP - x2) + (y2 - y1) * (yP - y2) + (z2 - z1) * (zP - z2)
                        k = gamma / (4 * np.pi * r1xr_sqr) * (r0r1 / r1 - r0r2 / r2)
                        v_ind[0] += k * r1xr2_x
                        v_ind[1] += k * r1xr2_y
                        v_ind[2] += k * r1xr2_z
                matrix_u[i].append(v_ind[0])
                matrix_v[i].append(v_ind[1])
                matrix_w[i].append(v_ind[2])
   
    

    
    for i in range(n_iterations):
        for j in range(len(control_points)):
            gamma_updated[j] = 0
            for k in range(len(rings)):
                gamma_updated[j] += matrix_u[j][k] + matrix_v[j][k] + matrix_w[j][k]
            gamma_updated[j] = gamma_updated[j] / 3
            gamma[j] = relax * gamma_updated[j] + (1 - relax) * gamma[j]    

                


