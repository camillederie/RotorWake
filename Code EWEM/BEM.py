import numpy as np
import pandas as pd
from Geometry import geo_blade

files = ['Code EWEM\polar_DU95W180.xlsx']
#'Code EWEM\polar_DU95W180.xlsx'
#Initializing tables    
cl_tab=np.zeros([61])
cd_tab=np.zeros([61])
cm_tab=np.zeros([61])
aoa_tab=np.zeros([61])
#Readin of tables. Only do this once at startup of simulation
aoa_tab[:],cl_tab[:],cd_tab[:],cm_tab[:] = pd.read_excel(files[0],header=None,skiprows=4).values.T

def force_coeffs(localalpha,aoa_tab,cl_tab,cd_tab,cm_tab):
    Cl=np.interp (localalpha,aoa_tab,cl_tab)
    Cd=np.interp (localalpha,aoa_tab,cd_tab)
    Cm=np.interp (localalpha,aoa_tab,cm_tab)
    return Cl, Cd, Cm 

def calculate_BEM(v_azim, v_axial,Omega, r_R):
    # Calculate the azimuthal and axial induction factors
    V_mag = np.sqrt(v_azim**2 + v_axial**2)
    phi = np.arctan(v_azim / v_axial)

    alpha = geo_blade(r_R)[1] + phi #  pitch + twist + inflow angle
    chord = geo_blade(r_R)[0]
    Cl, Cd, Cm = force_coeffs(alpha,aoa_tab,cl_tab,cd_tab,cm_tab)
    # Calculate the normal and tangential forces
    L = 0.5 *V_mag**2 * Cl * chord
    D = 0.5 *V_mag**2 * Cd * chord
    Fnorm = L * np.cos(phi) + D * np.sin(phi)
    Ftan = L * np.cos(phi) - D * np.sin(phi)
    Gamma = 0.5 * V_mag * chord * Cl
    return [Fnorm, Ftan, Gamma]