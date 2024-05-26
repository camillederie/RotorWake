import math as m

import numpy as np
from Geometry_copy import create_rotor_geometry, geo_blade, spanwise_discretisation
from LiftingLineSolver_copy import *
from matplotlib import pyplot as plt
from Plots_copy import plot_blade_geometry, plot_results
from Variables import *

#Main

#_____Single Rotor_____
#Inputs: Check variables file 
results = {}
for TSR in TSR_list: 
    Omega =  TSR * v_inf / R # Rotational speed in rad/s
    #Run the functions
    span_array, theta_array = spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations)
    system_geom = create_rotor_geometry(n_blades, TSR, n_rotations, span_array, R, v_inf, theta_array)

    #write code to save the system_geom dictionary to a .txt file
    # np.savetxt('system_geom.txt', np.array(list(system_geom.items())), fmt='%s')

    #Run the lifting line solver
    results[f"TSR_{TSR}"] = calculate_results(system_geom, v_inf, Omega, R, TSR)
    # print('CP, CT = ' ,results[3], results[4], results[5], results[6])

plot_results(results)
# plot_blade_geometry(system_geom)
# plt.show()
# Print all CP en CT values 

for TSR, result in results.items():
    print(f"TSR: {TSR}")
    print(f"CP: {result[3]}")
    print(f"CT: {result[4]}")
    print(f"CP_a: {result[5]}")
    print(f"CT_a: {result[6]}")
    print()

