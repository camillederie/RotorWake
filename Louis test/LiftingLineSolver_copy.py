#3D induced velocities, blade loads, run solver
import math

import matplotlib.pyplot as plt
import numpy as np
from BEM_for_LLM import calculate_BEM
from Geometry import geo_blade
from Variables import *


def LiftingLineSolver(system_geom, V_inf, Omega, R, TSR):
    # Inputs
    relax = 0.2
    n_iterations = 200 # 1200
    error_limit = 0.01
    
    # system_geom: Contains the geometry of horseshoe vortex rings and control points at the blade

    control_points = system_geom[0]
    rings = system_geom[3]
    
    matrix_u = []
    matrix_v = []
    matrix_w = []
    
    gamma_updated = np.ones(len(control_points))
    gamma = np.ones(len(control_points))
    print('length rings =',len(rings))
    print('length filaments =',len(rings[0]))
    print('length control points =',len(control_points))

    for i in range(len(control_points)):
        
        matrix_u.append([])
        matrix_v.append([])
        matrix_w.append([])
        #print('i =',i)
        # for j in range(len(rings)):
        #     for k in range(len(rings[0])):
        #         rings[j]['filaments'][k]['Gamma'] = 1 # Set ring strength to unity for calculating induced velocity at control point p
        for ring in rings:
            v_ind = np.zeros(3)
            for j, r in enumerate(rings): 

                # for f in r:
                    # GAMMA = f['Gamma']
                X1, Y1, Z1 = r[j][0], r[j][1], r[j][2]	# Start point of filament
                X2, Y2, Z2 = r[j+1][0], r[j+1][1], r[j+1][2]
                XP, YP, ZP = control_points[i][0], control_points[i][1], control_points[i][2]

                R1 = math.sqrt((XP - X1) ** 2 + (YP - Y1) ** 2 + (ZP - Z1) ** 2)
                R2 = math.sqrt((XP - X2) ** 2 + (YP - Y2) ** 2 + (ZP - Z2) ** 2)
                
                R1XR2_X = (YP - Y1) * (ZP - Z2) - (ZP - Z1) * (YP - Y2)
                R1XR2_Y = -(XP - X1) * (ZP - Z2) + (ZP - Z1) * (XP - X2)
                R1XR2_Z = (XP - X1) * (YP - Y2) - (YP - Y1) * (XP - X2)
                R1XR_SQR = R1XR2_X ** 2 + R1XR2_Y ** 2 + R1XR2_Z ** 2
                
                R0R1 = (X2 - X1) * (XP - X1) + (Y2 - Y1) * (YP - Y1) + (Z2 - Z1) * (ZP - Z1)
                R0R2 = (X2 - X1) * (XP - X2) + (Y2 - Y1) * (YP - Y2) + (Z2 - Z1) * (ZP - Z2)

                if R1XR_SQR < 0.0001 ** 2:
                    R1XR_SQR = 0.0001 ** 2
                if R1 < 0.0001:
                    R1 = 0.0001
                if R2 < 0.0001:
                    R2 = 0.0001

                K = (1 / (4 * math.pi * R1XR_SQR)) * ((R0R1 / R1) - (R0R2 / R2))
                
                v_ind[0] += K * R1XR2_X
                v_ind[1] += K * R1XR2_Y
                v_ind[2] += K * R1XR2_Z
            matrix_u[i].append(v_ind[0])
            matrix_v[i].append(v_ind[1])
            matrix_w[i].append(v_ind[2])

    print('Induced velocities calculated')
    # save U_matrix to txt file
    np.savetxt(f'U_matrix_{Omega}_us.txt', matrix_u, fmt='%1.4e')
    #matrix_u = np.loadtxt(f'U_matrix_{Omega}_pim.txt')
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
    ctan_list = np.zeros(len(control_points))
    cnormal_list = np.zeros(len(control_points))
    V_azim_list = np.zeros(len(control_points))
    u_list = np.zeros(len(control_points))
    v_total_mag_list = np.zeros(len(control_points))
    v_axial_list = np.zeros(len(control_points))
    v_azim_list = np.zeros(len(control_points))

    error_list = []
    iter_list = []

    for iter in range(n_iterations):
        gamma = gamma_updated.copy() 
        print('iter =',iter)
        for i in range(len(control_points)):
    
            pos_radial = np.linalg.norm(control_points[i])

            # u = v = w = 0
            # for j in range(len(rings)):
            #     u += matrix_u[i][j] * gamma[i]
            #     v += matrix_v[i][j] * gamma[i]
            #     w += matrix_w[i][j] * gamma[i]
            # matrix_u_flat = [item for sublist in matrix_u for item in sublist] #Flatten list
            # matrix_v_flat = [item for sublist in matrix_v for item in sublist]
            # matrix_w_flat = [item for sublist in matrix_w for item in sublist]

            u  = np.dot(matrix_u, gamma)
            v  = np.dot(matrix_v, gamma)
            w  = np.dot(matrix_w, gamma)
            # if i == 10 and iter == 40:
            #     print('u_test =',u_test)
            #     print('u =',u)
            
            # Calculate the velocity at the control point
           
            v_rotational = np.cross(np.array([-Omega, 0, 0]), np.array(control_points[i]))
            v_rotational_mag = np.linalg.norm(v_rotational)
            #print('v_rotational_mag =',v_rotational_mag, 'Omega =',Omega, 'pos_radial =',pos_radial, 'R =',R)
            #print('v_rotational =',v_rotational)
            v_inflow = np.array([V_inf, 0, 0]) # Assume no  yaw in the inflow so only x - direcion is considered	
            v_total = v_inflow + v_rotational + np.array([u, v, w])

            # For BEM, the azimuthal and axial velocities are needed:
            azimuth = np.cross([-1/pos_radial, 0, 0], np.array(control_points[i]))
            v_azim = np.dot(azimuth, v_total)
            v_azim_test = Omega * pos_radial + np.dot(v_inflow + np.array([u, v, w]), azimuth)
            # if i == 10 and iter == 40:
            #     print ('v_azim =',v_azim)
            #     print ('v_azim_test =',v_azim_test)
            v_axial = np.dot([1, 0, 0], v_total)
            if i == 10 and iter == 100:
                print ('v_tot 1 =',np.sqrt(v_azim**2 + v_axial**2))
                print ('v_tot 2',np.linalg.norm(v_total))
                
            
            # if (np.linalg.norm(v_total) - np.sqrt(v_azim**2 + v_axial**2)) > 0.01:
            #     print('Error in total velocity')

            BEM = calculate_BEM(v_azim, v_axial, Omega, pos_radial/ R)
            F_norm_list[i] = BEM[0]
            F_tan_list[i] = BEM[1]
            gamma_updated[i] = BEM[2]
            alpha_list[i] = BEM[3]
            phi_list[i] = BEM[4]
            ctan_list[i] = BEM[5]
            cnormal_list[i] = BEM[6]
            r_R_list[i] = pos_radial / R
            a_list[i] = -(u + v_rotational[0]) / V_inf
            a_line_list[i] = v_azim / (Omega * pos_radial) - 1
            pos_radial_list[i] = pos_radial
            V_azim_list[i] = v_azim
            u_list[i] = u
            v_total_mag_list[i] = np.linalg.norm(v_total)
            v_axial_list[i] = v_axial
            v_azim_list[i] = v_azim

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
    if TSR == 8:
        print(w)
    return [F_norm_list, F_tan_list, gamma_updated, alpha_list, phi_list, pos_radial_list, r_R_list, a_list, a_line_list, gamma_nondim, ctan_list, cnormal_list, V_azim_list, u_list, v_total_mag_list, v_axial_list, v_azim_list]

def calculate_results(system_geom, V_inf, Omega, R, TSR):
    # Calculate the results on the blade elements
    # system_geom: Contains the geometry of horseshoe vortex rings and control points at the blade
    # V_inf: Freestream velocity
    # Omega: Rotational speed
    # R: Rotor radius

    # Calculate the results on the blade elements
    results = LiftingLineSolver(system_geom, V_inf, Omega, R, TSR)
    number_of_cp_per_blade = int(len(system_geom[0])/ n_blades) # Number of control points per blade
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



           
    
                


