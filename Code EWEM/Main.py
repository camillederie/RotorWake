import math as m

import numpy as np
from Geometry import create_rotor_geometry, geo_blade, spanwise_discretisation
from LiftingLineSolver import *
from matplotlib import pyplot as plt
from Plots import plot_blade_geometry, plot_results
from Variables import *

#Main

#_____Single Rotor_____
#Inputs: Check variables file 

#Run the functions
span_array, theta_array = spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations)
system_geom = create_rotor_geometry(span_array, R, TSR, v_inf, theta_array, n_blades)

#plot_blade_geometry(system_geom)
#write code to save the system_geom dictionary to a .txt file
# np.savetxt('system_geom.txt', np.array(list(system_geom.items())), fmt='%s')

#Run the lifting line solver
results = calculate_results(system_geom, v_inf, Omega, R)
print('CP, CT = ' ,results[3], results[4])
plot_results(results[0], results[-3],results[-2],results[-1])   


