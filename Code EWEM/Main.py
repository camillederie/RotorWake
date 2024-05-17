import numpy as np
import math as m
from matplotlib import pyplot as plt
from Variables import *
from Geometry import spanwise_discretisation, geo_blade, create_rotor_geometry

#Main

#_____Single Rotor_____

#Run the functions
span_array, theta_array = spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations)
rotor_wake_system = create_rotor_geometry(span_array, R, TSR, v_inf, theta_array, n_blades)

