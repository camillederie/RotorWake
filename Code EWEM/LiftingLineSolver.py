#3D induced velocities, blade loads, run solver
import math

import matplotlib.pyplot as plt
import numpy as np
from BEM_for_LLM import calculate_BEM
from Geometry import geo_blade
from Variables import *


def LiftingLineSolver(system_geom, V_inf, Omega, R):
    # Inputs
    relax = 0.1
    n_iterations = 500 # 1200
    error_limit = 0.001
    
    # system_geom: Contains the geometry of horseshoe vortex rings and control points at the blade

    control_points = system_geom['controlpoints']
    rings = system_geom['rings']
    
    matrix_u = []
    matrix_v = []
    matrix_w = []
    
    gamma_updated = np.ones(len(control_points))
    gamma = np.ones(len(control_points))
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
                GAMMA = f['Gamma']
                X1, Y1, Z1 = f['x1'],f['y1'],f['z1']	# Start point of filament
                X2, Y2, Z2 = f['x2'],f['y2'],f['z2']
                XP, YP, ZP = control_points[i]['coordinates'][0], control_points[i]['coordinates'][1], control_points[i]['coordinates'][2]
                R1 = math.sqrt((XP - X1) ** 2 + (YP - Y1) ** 2 + (ZP - Z1) ** 2)
                R2 = math.sqrt((XP - X2) ** 2 + (YP - Y2) ** 2 + (ZP - Z2) ** 2)
                R1XR2_X = (YP - Y1) * (ZP - Z2) - (ZP - Z1) * (YP - Y2)
                R1XR2_Y = -(XP - X1) * (ZP - Z2) + (ZP - Z1) * (XP - X2)
                R1XR2_Z = (XP - X1) * (YP - Y2) - (YP - Y1) * (XP - X2)
                R1XR_SQR = R1XR2_X ** 2 + R1XR2_Y ** 2 + R1XR2_Z ** 2

                if R1XR_SQR < 0.0001:
                    R1XR_SQR = 0.0001
                if R1 < 0.0001:
                    R1 = 0.0001
                if R2 < 0.0001:
                    R2 = 0.0001

                R0R1 = (X2 - X1) * (XP - X1) + (Y2 - Y1) * (YP - Y1) + (Z2 - Z1) * (ZP - Z1)
                R0R2 = (X2 - X1) * (XP - X2) + (Y2 - Y1) * (YP - Y2) + (Z2 - Z1) * (ZP - Z2)

                K = GAMMA / (4 * math.pi * R1XR_SQR) * (R0R1 / R1 - R0R2 / R2)
                v_ind[0] += K * R1XR2_X
                v_ind[1] += K * R1XR2_Y
                v_ind[2] += K * R1XR2_Z
            matrix_u[i].append(v_ind[0])
            matrix_v[i].append(v_ind[1])
            matrix_w[i].append(v_ind[2])
    print('Induced velocities calculated')

    #print('matrix_u =',matrix_u)
    F_norm_list = np.zeros(len(control_points))
    F_tan_list = np.zeros(len(control_points))
    alpha_list = np.zeros(len(control_points))
    phi_list = np.zeros(len(control_points))
    r_R_list = np.zeros(len(control_points))
    pos_radial_list = np.zeros(len(control_points))
    a_list = np.zeros(len(control_points))
    a_line_list = np.zeros(len(control_points))
    gamma_nondim = np.ones(len(control_points))
    error_list = []
    iter_list = []

    for iter in range(n_iterations):
        gamma = gamma_updated.copy()    
        print('iter =',iter)
        for i in range(len(control_points)):
    
            pos_radial = np.linalg.norm(control_points[i]['coordinates'])
            
            u = v = w = 0
            for j in range(len(rings)):
                u += matrix_u[i][j] * gamma[j]
                v += matrix_v[i][j] * gamma[j]
                w += matrix_w[i][j] * gamma[j]
            # print('u =',u)
            # Calculate the velocity at the control point
           
            v_rotational = np.cross(np.array([-Omega, 0, 0]), np.array(control_points[i]['coordinates']))
            v_rotational_mag = np.linalg.norm(v_rotational)
            #print('v_rotational_mag =',v_rotational_mag, 'Omega =',Omega, 'pos_radial =',pos_radial, 'R =',R)
            #print('v_rotational =',v_rotational)
            v_inflow = np.array([V_inf, 0, 0]) # Assume no  yaw in the inflow so only x - direcion is considered	
            v_total = v_inflow + v_rotational + np.array([u, v, w])
        
            # For BEM, the azimuthal and axial velocities are needed:
            azimuth = np.cross([-1/pos_radial, 0, 0], np.array(control_points[i]['coordinates']))
            v_azim = np.dot(azimuth, v_total)
            v_azim_test = Omega * pos_radial + np.dot(v_inflow + np.array([u, v, w]), azimuth)
            # if i == 10 and iter == 40:
            #     print ('v_azim =',v_azim)
            #     print ('v_azim_test =',v_azim_test)
            v_axial = np.dot([1, 0, 0], v_total)

            # if i == 10:
            #     print ('v_azim =',v_azim)
            #     print ('v_azim_test =',v_azim_test)
            #     print ('v_axial =',v_axial)
            # print('v_azim =',v_azim)
            # print('v_axial =',v_axial)
            BEM = calculate_BEM(v_azim, v_axial, Omega, pos_radial/ R)
            F_norm_list[i] = BEM[0]
            F_tan_list[i] = BEM[1]
            gamma_updated[i] = BEM[2]
            alpha_list[i] = BEM[3]
            phi_list[i] = BEM[4]
            r_R_list[i] = pos_radial / R
            a_list[i] = -(u + v_rotational[0]) / V_inf
            a_line_list[i] = v_azim / (Omega * pos_radial) - 1
            pos_radial_list[i] = pos_radial

     

        error = max(abs(np.array(gamma_updated) - np.array(gamma)))
        error_list.append(error)
        iter_list.append(iter)
        if error < error_limit:
            break

        for i in range(len(control_points)):
            gamma_updated[i] = relax * gamma_updated[i] + (1 - relax) * gamma[i]
            gamma_nondim[i] = gamma_updated[i] * (n_blades * Omega) / (V_inf**2 * np.pi)


    #plot the error convergence
    # plt.figure()
    # plt.plot(iter_list, error_list)
    # plt.xlabel('Iteration')
    # plt.ylabel('Error')
    # plt.title('Convergence of the Lifting Line Solver')
    # plt.show()

    return [F_norm_list, F_tan_list, gamma_updated, alpha_list, phi_list, pos_radial_list, r_R_list, a_list, a_line_list, gamma_nondim]


def calculate_results(system_geom, V_inf, Omega, R):
    # Calculate the results on the blade elements
    # system_geom: Contains the geometry of horseshoe vortex rings and control points at the blade
    # V_inf: Freestream velocity
    # Omega: Rotational speed
    # R: Rotor radius

    # Calculate the results on the blade elements
    results = LiftingLineSolver(system_geom, V_inf, Omega, R)
    number_of_cp_per_blade = int(len(system_geom['controlpoints'])/ n_blades) # Number of control points per blade
    #print('number of cp per blade =',number_of_cp_per_blade)

    indeces_b1 = np.arange(0, number_of_cp_per_blade)
    indeces_b2 = np.arange(number_of_cp_per_blade, 2*number_of_cp_per_blade)
    indeces_b3 = np.arange(2*number_of_cp_per_blade, 3*number_of_cp_per_blade)
    
    ## TEST script ##
    T_B1_test = 0
    # dr = results[5][14]- results[5][13]
    for i in range(len(indeces_b1)-1):
        dr = results[5][i+1]- results[5][i]
        T_B1_test += results[0][i]*dr
    print('T_B1_test =',T_B1_test)
    # Calculate the total results on the blade
    T_B1 = np.trapz(results[0][indeces_b1], results[5][indeces_b1])
    print('T_B1 =',T_B1)
    T_B2 = np.trapz(results[0][indeces_b2], results[5][indeces_b2])
    T_B3 = np.trapz(results[0][indeces_b3], results[5][indeces_b3]) 
    T = T_B1 + T_B2 + T_B3
    P_B1 = np.trapz(results[1][indeces_b1] * results[5][indeces_b1], results[5][indeces_b1]) * Omega / (2*np.pi)
    P_B2 = np.trapz(results[1][indeces_b2] * results[5][indeces_b2], results[5][indeces_b2]) * Omega / (2*np.pi)
    P_B3 = np.trapz(results[1][indeces_b3] * results[5][indeces_b3], results[5][indeces_b3]) * Omega / (2*np.pi)
    P = P_B1 + P_B2 + P_B3

    # Calculate the power and thrust coefficients
    CP = P / (0.5 * rho * V_inf**3 * np.pi * R**2)
    CT = T / (0.5 * rho * V_inf**2 * np.pi * R**2)
    print('a = ', np.mean(results[7]))
    CP_a = 4 * np.mean(results[7]) * (1 - np.mean(results[7]))**2
    CT_a = 4 * np.mean(results[7]) * (1 - np.mean(results[7]))
    

    return [results, T, P, CP, CT, CP_a, CT_a, indeces_b1, indeces_b2, indeces_b3]



           
    
                


