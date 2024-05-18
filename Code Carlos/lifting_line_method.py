import math
import numpy as np

def create_array_sequence(start, step, end):
    return np.arange(start, end + step, step)

def solve_pBEM_lifting_line(TSR, NELEMENTS, Nrotations):
    s_Array = create_array_sequence(0., (math.pi) / NELEMENTS, math.pi).tolist()
    r_Array = [(-1 * (math.cos(s) - 1) / 2 * 0.8 + 0.2) for s in s_Array]

    for i in range(len(s_Array)):
        s_Array[i] = (-1 * (math.cos(s_Array[i]) - 1) / 2 * 0.8 + 0.2) * 50

    maxradius = max(s_Array)
    theta_Array = create_array_sequence(0., math.pi / 10, Nrotations * 2 * math.pi).tolist()
    rotor_wake_system = create_rotor_geometry(s_Array, maxradius, TSR / (1 - 0.2), 1, theta_Array, 3)
    
    resultsLL = solve_lifting_line_system_matrix_approach(rotor_wake_system, [1, 0, 0], TSR / 50, maxradius)
    
    results = solve_BEM_model(1, r_Array, TSR * 0.02, 50, 3)
    adim = math.pi / (3 * TSR / 50)

    results2 = calculate_CT_CProtor_CPflow(results['a'], results['aline'], results['Fnorm'], results['Ftan'], 1, r_Array, TSR * 0.02, 50, 3)
    CTCPliftline = calculate_CT_CProtor_CPflow(resultsLL['a'], resultsLL['aline'], resultsLL['Fnorm'], resultsLL['Ftan'], 1, r_Array, TSR * 0.02, 50, 3)

    return {
        'BEM': results,
        'BEMloads': results2,
        'LiftLine': resultsLL,
        'LiftLineloads': CTCPliftline,
        'rotor_wake_system': rotor_wake_system
    }

def solve_wing_lifting_line(Nelements, AspectRatio, Alpha):
    s_Array = create_array_sequence(0., (math.pi) / Nelements, math.pi).tolist()

    for i in range(len(s_Array)):
        s_Array[i] = (-1 * (math.cos(s_Array[i]) - 1) / 2) * AspectRatio

    WakeLength = 1000 * AspectRatio
    rotor_wake_system = create_straight_wing_geometry(s_Array, Alpha, WakeLength)

    results = solve_wing_lifting_line_system_matrix_approach(rotor_wake_system, Alpha)

    return {'Cl': results[0], 'span': results[1], 'rotor_wake_system': rotor_wake_system}

def create_straight_wing_geometry(span_array, Alpha, WakeLength):
    filaments = []
    ring = []
    controlpoints = []

    for i in range(len(span_array) - 1):
        controlpoints.append({
            'coordinates': [0, (span_array[i] + span_array[i + 1]) / 2, 0],
            'chord': 1,
            'normal': [0, 0, 1],
            'tangential': [1, 0, 0]
        })
        filaments.append({
            'x1': WakeLength,
            'y1': span_array[i],
            'z1': WakeLength * math.sin(Alpha * math.pi / 180),
            'x2': 1.25,
            'y2': span_array[i],
            'z2': 0,
            'Gamma': 0
        })
        filaments.append({
            'x1': 1.25,
            'y1': span_array[i],
            'z1': 0,
            'x2': 0,
            'y2': span_array[i],
            'z2': 0,
            'Gamma': 0
        })
        filaments.append({
            'x1': 0,
            'y1': span_array[i],
            'z1': 0,
            'x2': 0,
            'y2': span_array[i + 1],
            'z2': 0,
            'Gamma': 0
        })
        filaments.append({
            'x1': 0,
            'y1': span_array[i + 1],
            'z1': 0,
            'x2': 1.25,
            'y2': span_array[i + 1],
            'z2': 0,
            'Gamma': 0
        })
        filaments.append({
            'x1': 1.25,
            'y1': span_array[i + 1],
            'z1': 0,
            'x2': WakeLength,
            'y2': span_array[i + 1],
            'z2': WakeLength * math.sin(Alpha * math.pi / 180),
            'Gamma': 0
        })
        ring.append({'filaments': filaments})
        filaments = []

    return {'controlpoints': controlpoints, 'rings': ring}

def geo_blade(r_R):
    pitch = 2
    chord = 3 * (1 - r_R) + 1
    twist = -14 * (1 - r_R)
    return [chord, twist + pitch]

def create_rotor_geometry(span_array, radius, tipspeedratio, Uinf, theta_array, nblades):
    filaments = []
    ring = []
    controlpoints = []
    bladepanels = []

    for krot in range(nblades):
        angle_rotation = 2 * math.pi / nblades * krot
        cosrot = math.cos(angle_rotation)
        sinrot = math.sin(angle_rotation)

        for i in range(len(span_array) - 1):
            r = (span_array[i] + span_array[i + 1]) / 2
            geodef = geo_blade(r / radius)
            angle = geodef[1] * math.pi / 180
            
            temp1 = {
                'coordinates': [0, r, 0],
                'chord': geodef[0],
                'normal': [math.cos(angle), 0, -math.sin(angle)],
                'tangential': [-math.sin(angle), 0, -math.cos(angle)]
            }
            temp1['coordinates'] = [
                0,
                temp1['coordinates'][1] * cosrot - temp1['coordinates'][2] * sinrot,
                temp1['coordinates'][1] * sinrot + temp1['coordinates'][2] * cosrot
            ]
            temp1['normal'] = [
                temp1['normal'][0],
                temp1['normal'][1] * cosrot - temp1['normal'][2] * sinrot,
                temp1['normal'][1] * sinrot + temp1['normal'][2] * cosrot
            ]
            temp1['tangential'] = [
                temp1['tangential'][0],
                temp1['tangential'][1] * cosrot - temp1['tangential'][2] * sinrot,
                temp1['tangential'][1] * sinrot + temp1['tangential'][2] * cosrot
            ]

            controlpoints.append(temp1)
            temp1 = {
                'x1': 0,
                'y1': span_array[i],
                'z1': 0,
                'x2': 0,
                'y2': span_array[i + 1],
                'z2': 0,
                'Gamma': 0
            }
            filaments.append(temp1)
            geodef = geo_blade(span_array[i] / radius)
            angle = geodef[1] * math.pi / 180
            temp1 = {
                'x1': geodef[0] * math.sin(-angle),
                'y1': span_array[i],
                'z1': -geodef[0] * math.cos(angle),
                'x2': 0,
                'y2': span_array[i],
                'z2': 0,
                'Gamma': 0
            }
            filaments.append(temp1)

            for j in range(len(theta_array) - 1):
                xt = filaments[-1]['x1']
                yt = filaments[-1]['y1']
                zt = filaments[-1]['z1']
                dy = (math.cos(-theta_array[j + 1]) - math.cos(-theta_array[j])) * span_array[i]
                dz = (math.sin(-theta_array[j + 1]) - math.sin(-theta_array[j])) * span_array[i]
                dx = (theta_array[j + 1] - theta_array[j]) / tipspeedratio * radius

                temp1 = {
                    'x1': xt + dx,
                    'y1': yt + dy,
                    'z1': zt + dz,
                    'x2': xt,
                    'y2': yt,
                    'z2': zt,
                    'Gamma': 0
                }
                filaments.append(temp1)

            geodef = geo_blade(span_array[i + 1] / radius)
            angle = geodef[1] * math.pi / 180
            temp1 = {
                'x1': 0,
                'y1': span_array[i + 1],
                'z1': 0,
                'x2': geodef[0] * math.sin(-angle),
                'y2': span_array[i + 1],
                'z2': -geodef[0] * math.cos(angle),
                'Gamma': 0
            }
            filaments.append(temp1)

            for j in range(len(theta_array) - 1):
                xt = filaments[-1]['x2']
                yt = filaments[-1]['y2']
                zt = filaments[-1]['z2']
                dy = (math.cos(-theta_array[j + 1]) - math.cos(-theta_array[j])) * span_array[i + 1]
                dz = (math.sin(-theta_array[j + 1]) - math.sin(-theta_array[j])) * span_array[i + 1]
                dx = (theta_array[j + 1] - theta_array[j]) / tipspeedratio * radius

                temp1 = {
                    'x1': xt,
                    'y1': yt,
                    'z1': zt,
                    'x2': xt + dx,
                    'y2': yt + dy,
                    'z2': zt + dz,
                    'Gamma': 0
                }
                filaments.append(temp1)

            for ifil in range(len(filaments)):
                temp1 = filaments[ifil]
                temp2 = [
                    temp1['y1'] * cosrot - temp1['z1'] * sinrot,
                    temp1['y1'] * sinrot + temp1['z1'] * cosrot,
                    temp1['y2'] * cosrot - temp1['z2'] * sinrot,
                    temp1['y2'] * sinrot + temp1['z2'] * cosrot
                ]
                temp1['y1'] = temp2[0]
                temp1['z1'] = temp2[1]
                temp1['y2'] = temp2[2]
                temp1['z2'] = temp2[3]
                filaments[ifil] = temp1

            ring.append({'filaments': filaments})
            filaments = []

            geodef = geo_blade(span_array[i] / radius)
            angle = geodef[1] * math.pi / 180
            geodef2 = geo_blade(span_array[i + 1] / radius)
            angle2 = geodef2[1] * math.pi / 180

            temp1 = {
                'p1': [-0.25 * geodef[0] * math.sin(-angle), span_array[i], 0.25 * geodef[0] * math.cos(angle)],
                'p2': [-0.25 * geodef2[0] * math.sin(-angle2), span_array[i + 1], 0.25 * geodef2[0] * math.cos(angle2)],
                'p3': [0.75 * geodef2[0] * math.sin(-angle2), span_array[i + 1], -0.75 * geodef2[0] * math.cos(angle2)],
                'p4': [0.75 * geodef[0] * math.sin(-angle), span_array[i], -0.75 * geodef[0] * math.cos(angle)]
            }
            temp1['p1'] = [
                0,
                temp1['p1'][1] * cosrot - temp1['p1'][2] * sinrot,
                temp1['p1'][1] * sinrot + temp1['p1'][2] * cosrot
            ]
            temp1['p2'] = [
                0,
                temp1['p2'][1] * cosrot - temp1['p2'][2] * sinrot,
                temp1['p2'][1] * sinrot + temp1['p2'][2] * cosrot
            ]
            temp1['p3'] = [
                0,
                temp1['p3'][1] * cosrot - temp1['p3'][2] * sinrot,
                temp1['p3'][1] * sinrot + temp1['p3'][2] * cosrot
            ]
            temp1['p4'] = [
                0,
                temp1['p4'][1] * cosrot - temp1['p4'][2] * sinrot,
                temp1['p4'][1] * sinrot + temp1['p4'][2] * cosrot
            ]

            bladepanels.append(temp1)

    return {'controlpoints': controlpoints, 'rings': ring, 'bladepanels': bladepanels}

def solve_lifting_line_system_matrix_approach(rotor_wake_system, wind, Omega, rotorradius):
    controlpoints = rotor_wake_system['controlpoints']
    rings = rotor_wake_system['rings']
    velocity_induced = []
    up, vp, wp = [], [], []
    u, v, w = 0, 0, 0
    GammaNew = [0] * len(controlpoints)
    Gamma = [0] * len(controlpoints)
    MatrixU, MatrixV, MatrixW = [], [], []

    a_temp, aline_temp, r_R_temp, Fnorm_temp, Ftan_temp, Gamma_temp = [], [], [], [], [], []

    Niterations = 1200
    errorlimit = 0.01
    error = 1.0
    ConvWeight = 0.3

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
        for ig in range(len(GammaNew)):
            Gamma[ig] = GammaNew[ig]

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
            temploads = load_blade_element(vaxial, vazim, radialposition / rotorradius)
            GammaNew[icp] = temploads[2]

            a_temp.append(-(u + vrot[0]) / wind[0])
            aline_temp.append(vazim / (radialposition * Omega) - 1)
            r_R_temp.append(radialposition / rotorradius)
            Fnorm_temp.append(temploads[0])
            Ftan_temp.append(temploads[1])
            Gamma_temp.append(temploads[2])

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
        'Gamma': Gamma_temp
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
    Niterations = 340
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

    return [ClNew, rNew]
