import math as m
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

files = ['polar DU95W180.xlsx']

#Constants______________
R = 50 #m
R_root = 0.2*R
r = np.arange(0.2*R, R+0.1, 0.1)
B = 3
rho = 1.225 #kg/m3
Vo = 10 #m/s
TSR = 8 #abc Camille Derie
pitch = -2
twist = 14*(1-r/R)
chord = 3*(1-r/R)+1

def velocity_3d_from_vortex_filament(gamma, xv1, xv2, xvp1, core):
    """
    Calculate the velocity induced by a straight 3D vortex filament with circulation gamma
    at a point xvp1. The geometry of the vortex filament is defined by its edges: the filament
    starts at xv1 and ends at xv2. The input 'core' defines a vortex core radius, inside which
    the velocity is defined as a solid body rotation. The function adapts an algorithm from:
    Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics. Vol. 13. Cambridge university press, 2001.
    """
    # Read coordinates that define the vortex filament
    x1, y1, z1 = xv1  # Start point of vortex filament
    x2, y2, z2 = xv2  # End point of vortex filament
    
    # Read coordinates of target point where the velocity is calculated
    xp, yp, zp = xvp1
    
    # Calculate geometric relations for integral of the velocity induced by filament
    r1 = np.sqrt((xp - x1)**2 + (yp - y1)**2 + (zp - z1)**2)
    r2 = np.sqrt((xp - x2)**2 + (yp - y2)**2 + (zp - z2)**2)
    r1xr2_x = (yp - y1) * (zp - z2) - (zp - z1) * (yp - y2)
    r1xr2_y = -(xp - x1) * (zp - z2) + (zp - z1) * (xp - x2)
    r1xr2_z = (xp - x1) * (yp - y2) - (yp - y1) * (xp - x2)
    r1xr_sqr = r1xr2_x**2 + r1xr2_y**2 + r1xr2_z**2
    
    # Check if target point is in the vortex filament core, and modify to solid body rotation
    if r1xr_sqr < core**2:
        r1xr_sqr = core**2
    if r1 < core:
        r1 = core
    if r2 < core:
        r2 = core
    
    # Determine scalar
    r0r1 = (x2 - x1) * (xp - x1) + (y2 - y1) * (yp - y1) + (z2 - z1) * (zp - z1)
    r0r2 = (x2 - x1) * (xp - x2) + (y2 - y1) * (yp - y2) + (z2 - z1) * (zp - z2)
    k = gamma / (4 * np.pi * r1xr_sqr) * (r0r1 / r1 - r0r2 / r2)
    
    # Determine the three velocity components
    u = k * r1xr2_x
    v = k * r1xr2_y
    w = k * r1xr2_z
    
    # Output results, vector with the three velocity components
    results = [u, v, w]
    return results


def solve_lifting_line_system_matrix_approach(rotor_wake_system, wind, omega, rotor_radius):
    """
    Solves a lifting line model of a horizontal axis rotor.
    Inputs:
    - rotor_wake_system: Contains the geometry of horseshoe vortex rings and control points at the blade
    - wind: Unperturbed wind velocity (U_infinity) in the form [u, v, w]
    - omega: Rotational velocity of the rotor 
    - rotor_radius: The radius of the rotor
    """
    control_points = rotor_wake_system['controlpoints']
    rings = rotor_wake_system['rings']
    
    matrix_u = []
    matrix_v = []
    matrix_w = []
    
    gamma_new = np.zeros(len(control_points))
    gamma = np.zeros(len(control_points))
    
    n_iterations = 1200
    error_limit = 0.01
    conv_weight = 0.3
    
    # Initialize matrices
    for icp in range(len(control_points)):
        matrix_u.append([])
        matrix_v.append([])
        matrix_w.append([])
        for jring in range(len(rings)):
            # Set ring strength to unity for calculating induced velocity at control point icp
            updated_ring = update_gamma_single_ring(rings[jring], 1, 1)
            velocity_induced = velocity_induced_single_ring(updated_ring, control_points[icp]['coordinates'])
            matrix_u[icp].append(velocity_induced[0])
            matrix_v[icp].append(velocity_induced[1])
            matrix_w[icp].append(velocity_induced[2])
    
    # Iterative solution process
    for kiter in range(n_iterations):
        gamma = gamma_new.copy()
        
        for icp in range(len(control_points)):
            u = v = w = 0
            radial_position = np.linalg.norm(control_points[icp]['coordinates'])
            
            for jring in range(len(rings)):
                u += matrix_u[icp][jring] * gamma[jring]
                v += matrix_v[icp][jring] * gamma[jring]
                w += matrix_w[icp][jring] * gamma[jring]
            
            v_rot = np.cross([-omega, 0, 0], control_points[icp]['coordinates'])
            vel1 = np.array([wind[0] + u + v_rot[0], wind[1] + v + v_rot[1], wind[2] + w + v_rot[2]])
            
            azim_dir = np.cross([-1 / radial_position, 0, 0], control_points[icp]['coordinates'])
            v_azim = np.dot(azim_dir, vel1)
            v_axial = vel1[0]
            
            temp_loads = load_blade_element(v_axial, v_azim, radial_position / rotor_radius)
            gamma_new[icp] = temp_loads[2]
        
        refer_error = max(abs(gamma_new))
        refer_error = max(refer_error, 0.001)
        error = max(abs(gamma_new - gamma))
        error = error / refer_error
        
        if error < error_limit:
            break
        
        gamma_new = (1 - conv_weight) * gamma + conv_weight * gamma_new
    
    a_temp = -(u + v_rot[0]) / wind[0]
    aline_temp = v_azim / (radial_position * omega) - 1
    r_r_temp = radial_position / rotor_radius
    fnorm_temp = temp_loads[0]
    ftan_temp = temp_loads[1]
    gamma_temp = temp_loads[2]
    
    return {
        'a': a_temp,
        'aline': aline_temp,
        'r_R': r_r_temp,
        'Fnorm': fnorm_temp,
        'Ftan': ftan_temp,
        'Gamma': gamma_temp
    }