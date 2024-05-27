import math as m
import numpy as np
from Geometry import create_rotor_geometry, geo_blade, spanwise_discretisation
from LiftingLineSolver import *
from matplotlib import pyplot as plt
from Variables import *

results = {}
TSR = 8
Omega = TSR * v_inf / R  # Rotational speed in rad/s
r_BEM = np.linspace(0.2,1,401)

def convection_speed_analysis():
    print("Convection speed analysis")
    a_lst = [0.2, 0.25, 0.3]
    for a in a_lst:
        span_array, theta_array = spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations, n_theta)
        system_geom = create_rotor_geometry(span_array, R, TSR, v_inf, theta_array, n_blades, a)
        results[a] = calculate_results(system_geom, v_inf, Omega, R)

    # Check if all a values have been added to the results
    if all(a in results for a in a_lst):
        results_02 = results[a_lst[0]][0]
        results_025 = results[a_lst[1]][0]
        results_03 = results[a_lst[2]][0]
        indeces_b1 = results[a_lst[0]][7]

        r_BEM = np.linspace(0.2,1,401)
        plt.figure()
        plt.plot(results_02[6][indeces_b1], results_02[9][indeces_b1], label='a = 0.2')
        plt.plot(results_025[6][indeces_b1], results_025[9][indeces_b1], label='a = 0.25')
        plt.plot(results_03[6][indeces_b1], results_03[9][indeces_b1], label='a = 0.3')
        plt.xlabel('Radial position [-]')
        plt.ylabel('Gamma [-]')
        plt.legend()
        plt.grid()
        plt.savefig('Figures/Gamma_vs_radial_position_convection.png')
        plt.close()

        #Plot the error list in function of number of iterations (results[17] is the error list, results[18] is the iteration list) on a semilog y scale
        plt.figure()
        plt.semilogy(results_02[18], results_02[17], label='a = 0.2')
        plt.semilogy(results_025[18], results_025[17], label='a = 0.25')
        plt.semilogy(results_03[18], results_03[17], label='a = 0.3')
        plt.xlabel('Iteration')
        plt.ylabel('Error')
        plt.legend()
        plt.grid()
        plt.savefig('Figures/Error_vs_iteration_convection.png')
        plt.close()

        #Calculate Cp and Ct for the different a values and print them
        for a in a_lst:
            print(f"a = {a}")
            print(f"Cp: {results[a][5]}")
            print(f"Ct: {results[a][6]}")
            print()

def spacing_method_analysis():
    print("spacing_method_analysis")
    for Method in ['Equal', 'Cosine']:
        span_array, theta_array = spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations, n_theta)
        system_geom = create_rotor_geometry(span_array, R, TSR, v_inf, theta_array, n_blades, a)
        results[Method] = calculate_results(system_geom, v_inf, Omega, R)
    
    # Check if both methods have been added to the results
    if 'Equal' in results and 'Cosine' in results:
        results_equal = results['Equal'][0]
        results_cosine = results['Cosine'][0]
        indeces_b1_equal = results['Equal'][7]
        indeces_b1_cosine = results['Cosine'][7]

        plt.figure()
        plt.plot(results_equal[6][indeces_b1_equal], results_equal[9][indeces_b1_equal], label='Equal spacing')
        plt.plot(results_cosine[6][indeces_b1_cosine], results_cosine[9][indeces_b1_cosine], label='Cosine spacing')
        plt.xlabel('Radial position [-]')
        plt.ylabel('Gamma [-]')
        plt.legend()
        plt.grid()
        plt.savefig('Figures/Gamma_vs_radial_position_spacing.png')
        plt.close()

        #Plot the error list in function of number of iterations (results[17] is the error list, results[18] is the iteration list)  on a semilog y scale
        plt.figure()
        plt.semilogy(results_equal[18], results_equal[17], label='Equal spacing')
        plt.semilogy(results_cosine[18], results_cosine[17], label='Cosine spacing')
        plt.xlabel('Iteration')
        plt.ylabel('Error')
        plt.legend()
        plt.grid()
        plt.savefig('Figures/Error_vs_iteration_spacing.png')
        plt.close()

        #Calculate Cp and Ct for the different a values and print them
        for Method in ['Equal', 'Cosine']:
            print(f"Method = {Method}")
            print(f"Cp: {results[Method][5]}")
            print(f"Ct: {results[Method][6]}")
            print()
    else:
        print("Error: One or both methods did not produce results.")

def segments_rotation_analysis():
    print("Segments and rotation analysis")
    segments = [5, 10, 30]
    for n_theta in segments:
        span_array, theta_array = spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations, n_theta)
        system_geom = create_rotor_geometry(span_array, R, TSR, v_inf, theta_array, n_blades, a)
        results[n_theta] = calculate_results(system_geom, v_inf, Omega, R)

    # Check if all n_theta values have been added to the results
    if all(n_theta in results for n_theta in segments):
        results_5 = results[segments[0]][0]
        results_10 = results[segments[1]][0]
        results_30 = results[segments[2]][0]
        indeces_b1 = results[segments[0]][7]

        plt.figure()
        plt.plot(results_5[6][indeces_b1], results_5[9][indeces_b1], label='5 segments per rotation')
        plt.plot(results_10[6][indeces_b1], results_10[9][indeces_b1], label='10 segments per rotation')
        plt.plot(results_30[6][indeces_b1], results_30[9][indeces_b1], label='30 segments per rotation')
        plt.xlabel('Radial position [-]')
        plt.ylabel('Gamma [-]')
        plt.legend()
        plt.grid()
        plt.savefig('Figures/Gamma_vs_radial_position_segments.png')
        plt.close()

        #Plot the error list in function of number of iterations (results[17] is the error list, results[18] is the iteration list) on a semilog y scale
        plt.figure()
        plt.semilogy(results_5[18], results_5[17], label='5 segments per rotation')
        plt.semilogy(results_10[18], results_10[17], label='10 segments per rotation')
        plt.semilogy(results_30[18], results_30[17], label='30 segments per rotation')
        plt.xlabel('Iteration')
        plt.ylabel('Error')
        plt.legend()
        plt.grid()
        plt.savefig('Figures/Error_vs_iteration_segments.png')
        plt.close()

        #Calculate Cp and Ct for the different a values and print them
        for n_theta in segments:
            print(f"n_theta = {n_theta}")
            print(f"Cp: {results[n_theta][5]}")
            print(f"Ct: {results[n_theta][6]}")
            print()

    else:
        print("Error: One or more n_theta values did not produce results.")

def wake_length_analysis():
    print("Wake length analysis")
    rotations = [2, 5, 10, 20]
    for n_rotations in rotations:
        span_array, theta_array = spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations, n_theta)
        system_geom = create_rotor_geometry(span_array, R, TSR, v_inf, theta_array, n_blades, a)
        results[n_rotations] = calculate_results(system_geom, v_inf, Omega, R)

    # Check if all n_rotations values have been added to the results
    if all(n_rotations in results for n_rotations in rotations):
        results_2 = results[rotations[0]][0]
        results_5 = results[rotations[1]][0]
        results_10 = results[rotations[2]][0]
        results_20 = results[rotations[3]][0]
        indeces_b1 = results[rotations[0]][7]

        plt.figure()
        plt.plot(results_2[6][indeces_b1], results_2[9][indeces_b1], label='2 rotations')
        plt.plot(results_5[6][indeces_b1], results_5[9][indeces_b1], label='5 rotations')
        plt.plot(results_10[6][indeces_b1], results_10[9][indeces_b1], label='10 rotations')
        plt.plot(results_20[6][indeces_b1], results_20[9][indeces_b1], label='20 rotations')
        plt.xlabel('Radial position [-]')
        plt.ylabel('Gamma [-]')
        plt.legend()
        plt.grid()
        plt.savefig('Figures/Gamma_vs_radial_position_rotations.png')
        plt.close()

        #Plot the error list in function of number of iterations (results[17] is the error list, results[18] is the iteration list)  on a semilog y scale
        plt.figure()
        plt.semilogy(results_2[18], results_2[17], label='2 rotations')
        plt.semilogy(results_5[18], results_5[17], label='5 rotations')
        plt.semilogy(results_10[18], results_10[17], label='10 rotations')
        plt.semilogy(results_20[18], results_20[17], label='20 rotations')
        plt.xlabel('Iteration')
        plt.ylabel('Error')
        plt.legend()
        plt.grid()
        plt.savefig('Figures/Error_vs_iteration_rotations.png')
        plt.close()

        #Calculate Cp and Ct for the different a values and print them
        for n_rotations in rotations:
            print(f"n_rotations = {n_rotations}")
            print(f"Cp: {results[n_rotations][5]}")
            print(f"Ct: {results[n_rotations][6]}")
            print()
    else:
        print("Error: One or more n_rotations values did not produce results.")
        
convection_speed_analysis()
spacing_method_analysis()
segments_rotation_analysis()
wake_length_analysis()
