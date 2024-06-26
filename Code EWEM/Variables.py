#List of all variables used in the code, to be imported in the main code
import math as m

R = 50.0
R_Root_Ratio = 0.2
TSR_list = [6,8,10]
v_inf = 10
wind = [10.0,0,0]
a = 0.2

n_span = 20
n_blades = 3
n_rotations = 10
n_theta = 20 # Per rotation!
#Omega =  TSR * v_inf / R # Rotational speed in rad/s
rho = 1.225
Method = "Cosine"