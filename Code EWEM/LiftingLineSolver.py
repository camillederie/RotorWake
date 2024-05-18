#3D induced velocities, blade loads, run solver
import numpy as np
from BEM import calculate_BEM


def LiftingLineSolver(system_geom, V_inf, Omega, R):
    # Inputs
    relax = 0.3
    n_iterations = 1 # 1200
    error_limit = 0.01
    
    # system_geom: Contains the geometry of horseshoe vortex rings and control points at the blade

    control_points = system_geom['controlpoints']
    rings = system_geom['rings']
    
    matrix_u = []
    matrix_v = []
    matrix_w = []
    
    gamma_updated = np.zeros(len(control_points))
    gamma = np.zeros(len(control_points))
    print('length rings =',len(rings))
    print('length filaments =',len(rings[0]['filaments']))
    print('length control points =',len(control_points))

    for i in range(len(control_points)):
        matrix_u.append([])
        matrix_v.append([])
        matrix_w.append([])
        #print('i =',i)
        # for j in range(len(rings)):
        #     for k in range(len(rings[0]['filaments'])):
        #         #rings[j]['filaments'][k]['Gamma'] = 1 # Set ring strength to unity for calculating induced velocity at control point p
        
        for r in rings: 
            v_ind = np.zeros(3)
            for f in r['filaments']:
                gamma = f['Gamma']
                x1, y1, z1 = f['x1'],f['y1'],f['z1']	# Start point of filament
                x2, y2, z2 = f['x2'],f['y2'],f['z2']
                xP, yP, zP = control_points[i]['coordinates'][0], control_points[i]['coordinates'][1], control_points[i]['coordinates'][2]
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
    print('Induced velocities calculated')

    
    F_norm_list = np.zeros(len(control_points))
    F_tan_list = np.zeros(len(control_points))

    for iter in range(n_iterations):
        gamma = gamma_updated.copy()

        for i in range(len(control_points)):
    
            pos_radial = np.linalg.norm(control_points[i]['coordinates'])
            
            u = v = w = 0
            for j in range(len(rings)):
                u += matrix_u[i][j] * gamma[j]
                v += matrix_v[i][j] * gamma[j]
                w += matrix_w[i][j] * gamma[j]

            # Calculate the velocity at the control point
           
            v_rotational = np.cross(np.array([-Omega, 0, 0]), np.array(control_points[i]['coordinates']))
            #print('v_rotational =',v_rotational)
            v_inflow = np.array([V_inf, 0, 0]) # Assume no  yaw in the inflow so only x - direcion is considered	
            v_total = v_inflow + v_rotational + np.array([u, v, w])
            #print('v_total =',v_total)
            # For BEM, the azimuthal and axial velocities are needed:
            azimuth = np.cross([-1/pos_radial, 0, 0], np.array(control_points[i]['coordinates']))
            v_azim = np.dot(azimuth, v_total)
            v_axial = np.dot([1, 0, 0], v_total)
            BEM = calculate_BEM(v_azim, v_axial, Omega, pos_radial/ R)
            F_norm_list[i] = BEM[0]
            F_tan_list[i] = BEM[1]
            gamma_updated[i] = BEM[2]


     

        error = max(abs(np.array(gamma_updated) - np.array(gamma)))
        if error < error_limit:
            break

        for i in range(len(control_points)):
            gamma_updated[i] = relax * gamma_updated[i] + (1 - relax) * gamma[i]

    return F_norm_list, F_tan_list, gamma_updated






           
    
                


