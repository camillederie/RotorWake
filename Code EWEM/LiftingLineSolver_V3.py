import numpy as np
import matplotlib.pyplot as plt
import math
from Variables import *
from BEM_for_LLM import calculate_BEM


"""
This is a test and its not working well :(
"""

def solve_lifting_line_system_matrix_approach(rotor_wake_system, wind, Omega, rotorradius):
    controlpoints = rotor_wake_system['controlpoints']
    rings = rotor_wake_system['rings']
    velocity_induced = []
    up, vp, wp = [], [], []
    u, v, w = 0, 0, 0
    GammaNew = [0] * len(controlpoints)
    Gamma = [0] * len(controlpoints)
    GammaNondim = [0] * len(controlpoints)
    MatrixU, MatrixV, MatrixW = [], [], []


    Niterations = 500
    errorlimit = 0.01
    error = 1.0
    ConvWeight = 0.1

    for icp in range(len(controlpoints)):
        MatrixU.append([])
        MatrixV.append([])
        MatrixW.append([])
        for jring in range(len(rings)):
            rings[jring] = update_Gamma_single_ring(rings[jring], 1, 1)
            velocity_induced = velocity_induced_single_ring(rings[jring], controlpoints[icp]['coordinates'])
            MatrixU[icp].append(velocity_induced[0])
            MatrixV[icp].append(velocity_induced[1])
            MatrixW[icp].append(velocity_induced[2])

    for kiter in range(Niterations):
        a_temp, aline_temp, r_R_temp, Fnorm_temp, Ftan_temp, Gamma_temp, alpha_temp = [], [], [], [], [], [], []
        for ig in range(len(GammaNew)):
            Gamma[ig] = GammaNew[ig]
            GammaNondim[ig] = Gamma[ig] * (n_blades * Omega) / (v_inf**2 * np.pi)
        for icp in range(len(controlpoints)):
            radialposition = math.sqrt(np.dot(controlpoints[icp]['coordinates'], controlpoints[icp]['coordinates']))
            u, v, w = 0, 0, 0
            for jring in range(len(rings)):
                u += MatrixU[icp][jring] * Gamma[jring]
                v += MatrixV[icp][jring] * Gamma[jring]
                w += MatrixW[icp][jring] * Gamma[jring]

            vrot = np.cross([-Omega, 0, 0], controlpoints[icp]['coordinates'])
            vel1 = [wind[0] + u + vrot[0], wind[1] + v + vrot[1], wind[2] + w + vrot[2]]
            azimdir = np.cross([-1 / radialposition, 0, 0], controlpoints[icp]['coordinates'])
            vazim = np.dot(azimdir, vel1)
            vaxial = np.dot([1, 0, 0], vel1)
            temploads = calculate_BEM(vaxial, vazim, Omega, radialposition / rotorradius)
            GammaNew[icp] = temploads[2]

            a_temp.append(-(u + vrot[0]) / wind[0])
            aline_temp.append(vazim / (radialposition * Omega) - 1)
            r_R_temp.append(radialposition / rotorradius)
            Fnorm_temp.append(temploads[0])
            Ftan_temp.append(temploads[1])
            Gamma_temp.append(temploads[2])
            alpha_temp.append(temploads[3])

        refererror = max(abs(np.array(GammaNew)))
        refererror = max(refererror, 0.001)
        error = max(abs(np.array(GammaNew) - np.array(Gamma)))
        error /= refererror
        if error < errorlimit:
            break

        for ig in range(len(GammaNew)):
            GammaNew[ig] = (1 - ConvWeight) * Gamma[ig] + ConvWeight * GammaNew[ig] 
            
    return {
        'a': a_temp,
        'aline': aline_temp,
        'r_R': r_R_temp,
        'Fnorm': Fnorm_temp,
        'Ftan': Ftan_temp,
        'Gamma': GammaNondim,
        'alpha': alpha_temp
    }

def update_Gamma_single_ring(ring, GammaNew, WeightNew):
    for i in range(len(ring['filaments'])):
        ring['filaments'][i]['Gamma'] = ring['filaments'][i]['Gamma'] * (1 - WeightNew) + WeightNew * GammaNew
    return ring

def velocity_induced_rings(rings, controlpoints):
    velind = [0, 0, 0]
    for ring in rings:
        tempvel1 = velocity_induced_single_ring(ring, controlpoints)
        velind = np.add(velind, tempvel1)
    return velind

def velocity_induced_single_ring(ring, controlpoint):
    velind = [0, 0, 0]
    CORE = 0.00001
    for filament in ring['filaments']:
        GAMMA = filament['Gamma']
        XV1 = [filament['x1'], filament['y1'], filament['z1']]
        XV2 = [filament['x2'], filament['y2'], filament['z2']]
        tempvel1 = velocity_3D_from_vortex_filament(GAMMA, XV1, XV2, controlpoint, CORE)
        velind[0] += tempvel1[0]
        velind[1] += tempvel1[1]
        velind[2] += tempvel1[2]
    return velind

def velocity_3D_from_vortex_filament(GAMMA, XV1, XV2, XVP1, CORE):
    X1, Y1, Z1 = XV1
    X2, Y2, Z2 = XV2
    XP, YP, ZP = XVP1

    R1 = math.sqrt((XP - X1) ** 2 + (YP - Y1) ** 2 + (ZP - Z1) ** 2)
    R2 = math.sqrt((XP - X2) ** 2 + (YP - Y2) ** 2 + (ZP - Z2) ** 2)
    R1XR2_X = (YP - Y1) * (ZP - Z2) - (ZP - Z1) * (YP - Y2)
    R1XR2_Y = -(XP - X1) * (ZP - Z2) + (ZP - Z1) * (XP - X2)
    R1XR2_Z = (XP - X1) * (YP - Y2) - (YP - Y1) * (XP - X2)
    R1XR_SQR = R1XR2_X ** 2 + R1XR2_Y ** 2 + R1XR2_Z ** 2
    R0R1 = (X2 - X1) * (XP - X1) + (Y2 - Y1) * (YP - Y1) + (Z2 - Z1) * (ZP - Z1)
    R0R2 = (X2 - X1) * (XP - X2) + (Y2 - Y1) * (YP - Y2) + (Z2 - Z1) * (ZP - Z2)

    if R1XR_SQR < CORE ** 2:
        R1XR_SQR = CORE ** 2
    if R1 < CORE:
        R1 = CORE
    if R2 < CORE:
        R2 = CORE

    K = GAMMA / (4 * math.pi * R1XR_SQR) * (R0R1 / R1 - R0R2 / R2)
    U = K * R1XR2_X
    V = K * R1XR2_Y
    W = K * R1XR2_Z

    return [U, V, W]

def solve_wing_lifting_line_system_matrix_approach(rotor_wake_system, Alpha):
    controlpoints = rotor_wake_system['controlpoints']
    rings = rotor_wake_system['rings']
    velocity_induced = []
    up, vp, wp = [], [], []
    u, v, w = 0, 0, 0
    GammaNew = [0] * len(controlpoints)
    ClNew = [0] * len(controlpoints)
    rNew = [controlpoint['coordinates'][1] for controlpoint in controlpoints]
    Niterations = 20
    MatrixU, MatrixV, MatrixW = [], [], []
    errorlimit = 0.01
    error = 1.0
    ConvWeight = 0.1

    for icp in range(len(controlpoints)):
        for jring in range(len(rings)):
            rings[jring] = update_Gamma_single_ring(rings[jring], 1, 1)
            velocity_induced = velocity_induced_single_ring(rings[jring], controlpoints[icp]['coordinates'])
            up.append(velocity_induced[0])
            vp.append(velocity_induced[1])
            wp.append(velocity_induced[2])
            velocity_induced = []

        MatrixU.append(up)
        MatrixV.append(vp)
        MatrixW.append(wp)
        up, vp, wp = [], [], []

    for kiter in range(Niterations):
        Gamma = GammaNew.copy()
        for icp in range(len(controlpoints)):
            u, v, w = 0, 0, 0
            for jring in range(len(rings)):
                u += MatrixU[icp][jring] * Gamma[jring]
                v += MatrixV[icp][jring] * Gamma[jring]
                w += MatrixW[icp][jring] * Gamma[jring]

            vel1 = [1 + u, v, 1 * math.sin(Alpha * math.pi / 180) + w]
            angle1 = math.atan(vel1[2] / vel1[0])
            ClNew[icp] = 2 * math.pi * math.sin(angle1)
            vmag = math.sqrt(np.dot(vel1, vel1))
            GammaNew[icp] = 0.5 * 1 * vmag * ClNew[icp]

        refererror = max(abs(np.array(GammaNew)))
        refererror = max(refererror, 0.001)
        error = max(abs(np.array(GammaNew) - np.array(Gamma)))
        error /= refererror
        ConvWeight = max((1 - error) * 0.3, 0.1)

        if error < errorlimit:
            break

        for ig in range(len(GammaNew)):
            GammaNew[ig] = (1 - ConvWeight) * Gamma[ig] + ConvWeight * GammaNew[ig]

    return GammaNew

