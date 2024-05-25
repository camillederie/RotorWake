import numpy as np
import matplotlib.pyplot as plt
import math

# Placeholder for imported classes
from turbine2 import Turbine
from experimental_data import Data

# Case specifics
dr = 0.01
Tb = Turbine(dr)
U0 = 10  # [m/s]
rho = 1.225  # [kg/m^3] Using ISA

TSR = [6, 8, 10]
omega = U0 * TSR[1] / Tb.R

# Experimental data
C_L = lambda alpha: Data(Tb.type).lift(np.degrees(alpha))
C_D = lambda alpha: Data(Tb.type).drag(np.degrees(alpha))

# Utility functions
def lift(alpha, chord, velocity, rho):
    return C_L(alpha) * 0.5 * rho * chord * velocity ** 2

def drag(alpha, chord, velocity, rho):
    return C_D(alpha) * 0.5 * rho * chord * velocity ** 2

# Geometry and velocity calculations
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
            geodef = Tb.geometry(r / radius)
            angle = geodef[1] * math.pi / 180

            temp1 = {
                'coordinates': [0, r, 0],
                'chord': geodef[0],
                'normal': [math.cos(angle), 0, -math.sin(angle)],
                'tangential': [-math.sin(angle), 0, -math.cos(angle)]
            }

            temp1['coordinates'] = [0, temp1['coordinates'][1] * cosrot - temp1['coordinates'][2] * sinrot,
                                    temp1['coordinates'][1] * sinrot + temp1['coordinates'][2] * cosrot]
            temp1['normal'] = [temp1['normal'][0],
                               temp1['normal'][1] * cosrot - temp1['normal'][2] * sinrot,
                               temp1['normal'][1] * sinrot + temp1['normal'][2] * cosrot]
            temp1['tangential'] = [temp1['tangential'][0],
                                   temp1['tangential'][1] * cosrot - temp1['tangential'][2] * sinrot,
                                   temp1['tangential'][1] * sinrot + temp1['tangential'][2] * cosrot]

            controlpoints.append(temp1)

            temp1 = {'x1': 0, 'y1': span_array[i], 'z1': 0, 'x2': 0, 'y2': span_array[i + 1], 'z2': 0, 'Gamma': 0}
            filaments.append(temp1)

            geodef = Tb.geometry(span_array[i] / radius)
            angle = geodef[1] * math.pi / 180
            temp1 = {'x1': geodef[0] * math.sin(-angle), 'y1': span_array[i],
                     'z1': -geodef[0] * math.cos(angle), 'x2': 0, 'y2': span_array[i], 'z2': 0, 'Gamma': 0}
            filaments.append(temp1)

            for j in range(len(theta_array) - 1):
                xt = filaments[-1]['x1']
                yt = filaments[-1]['y1']
                zt = filaments[-1]['z1']
                dy = (math.cos(-theta_array[j + 1]) - math.cos(-theta_array[j])) * span_array[i]
                dz = (math.sin(-theta_array[j + 1]) - math.sin(-theta_array[j])) * span_array[i]
                dx = (theta_array[j + 1] - theta_array[j]) / tipspeedratio * radius

                temp1 = {'x1': xt + dx, 'y1': yt + dy, 'z1': zt + dz, 'x2': xt, 'y2': yt, 'z2': zt, 'Gamma': 0}
                filaments.append(temp1)

            for fil in filaments:
                rotated_y1 = fil['y1'] * cosrot - fil['z1'] * sinrot
                rotated_z1 = fil['y1'] * sinrot + fil['z1'] * cosrot
                rotated_y2 = fil['y2'] * cosrot - fil['z2'] * sinrot
                rotated_z2 = fil['y2'] * sinrot + fil['z2'] * cosrot
                fil.update({'y1': rotated_y1, 'z1': rotated_z1, 'y2': rotated_y2, 'z2': rotated_z2})

            ring.append({'filaments': filaments})
            filaments = []

            geodef = Tb.geometry(span_array[i] / radius)
            angle = geodef[1] * math.pi / 180
            geodef2 = Tb.geometry(span_array[i + 1] / radius)
            angle2 = geodef2[1] * math.pi / 180

            temp1 = {
                'p1': [-0.25 * geodef[0] * math.sin(-angle), span_array[i], 0.25 * geodef[0] * math.cos(angle)],
                'p2': [-0.25 * geodef2[0] * math.sin(-angle2), span_array[i + 1], 0.25 * geodef2[0] * math.cos(angle2)],
                'p3': [0.75 * geodef2[0] * math.sin(-angle2), span_array[i + 1], -0.75 * geodef2[0] * math.cos(angle2)],
                'p4': [0.75 * geodef[0] * math.sin(-angle), span_array[i], -0.75 * geodef[0] * math.cos(angle)]
            }
            temp1['p1'] = [0, temp1['p1'][1] * cosrot - temp1['p1'][2] * sinrot, temp1['p1'][1] * sinrot + temp1['p1'][2] * cosrot]
            temp1['p2'] = [0, temp1['p2'][1] * cosrot - temp1['p2'][2] * sinrot, temp1['p2'][1] * sinrot + temp1['p2'][2] * cosrot]
            temp1['p3'] = [0, temp1['p3'][1] * cosrot - temp1['p3'][2] * sinrot, temp1['p3'][1] * sinrot + temp1['p3'][2] * cosrot]
            temp1['p4'] = [0, temp1['p4'][1] * cosrot - temp1['p4'][2] * sinrot, temp1['p4'][1] * sinrot + temp1['p4'][2] * cosrot]

            bladepanels.append(temp1)

    return {'controlpoints': controlpoints, 'rings': ring, 'bladepanels': bladepanels}

def update_Gamma_single_ring(ring, GammaNew, WeightNew):
    for filament in ring['filaments']:
        filament['Gamma'] = filament['Gamma'] * (1 - WeightNew) + WeightNew * GammaNew
    return ring

def velocity_induced_rings(rings, controlpoint):
    vel_ind = np.array([0.0, 0.0, 0.0])
    for ring in rings:
        temp_vel1 = velocity_induced_single_ring(ring, controlpoint)
        vel_ind += temp_vel1
    return vel_ind

def velocity_induced_single_ring(ring, controlpoint):
    vel_ind = np.array([0.0, 0.0, 0.0])
    CORE = 0.00001
    for filament in ring['filaments']:
        GAMMA = filament['Gamma']
        XV1 = np.array([filament['x1'], filament['y1'], filament['z1']])
        XV2 = np.array([filament['x2'], filament['y2'], filament['z2']])
        XVP = np.array(controlpoint)
        temp_vel1 = velocity_3D_from_vortex_filament(GAMMA, XV1, XV2, XVP, CORE)
        vel_ind += temp_vel1
    return vel_ind

def velocity_3D_from_vortex_filament(GAMMA, XV1, XV2, XVP, CORE):
    r0 = XV2 - XV1
    rp = XVP - XV1
    rq = XVP - XV2

    r0_cross_rp = np.cross(r0, rp)
    norm_r0_cross_rp = np.linalg.norm(r0_cross_rp)

    if norm_r0_cross_rp < CORE:
        return np.array([0.0, 0.0, 0.0])

    induced_velocity = GAMMA / (4 * np.pi * norm_r0_cross_rp ** 2) * r0_cross_rp

    return induced_velocity

def loadBladeElement(Vnorm, Vtan, r_R):
    Vmag2 = Vnorm**2 + Vtan**2
    InflowAngle = math.atan2(Vnorm, Vtan)
    chord, twist = Tb.geometry(r_R)
    alpha = twist + InflowAngle * 180 / math.pi
    if alpha > np.radians(30):
        alpha = np.radians(30)
    if alpha < np.radians(-16):
        alpha = np.radians(-16)
    C_L = lambda alpha: Data(Tb.type).lift(alpha)
    C_D = lambda alpha: Data(Tb.type).drag(alpha)
    cl = C_L(alpha)
    cd = C_D(alpha)
    cd = 0 * cd

    Lift = 0.5 * Vmag2 * cl * chord
    Drag = 0.5 * Vmag2 * cd * chord

    Fnorm = Lift * math.cos(InflowAngle) + Drag * math.sin(InflowAngle)
    Ftan = Lift * math.sin(InflowAngle) - Drag * math.cos(InflowAngle)

    Gamma = 0.5 * math.sqrt(Vmag2) * cl * chord

    return [Fnorm, Ftan, Gamma]

def solve_lifting_line_system_matrix_approach(rotor_wake_system, wind, Omega, rotorradius):
    controlpoints = rotor_wake_system['controlpoints']
    rings = rotor_wake_system['rings']

    GammaNew = np.zeros(len(controlpoints))
    Gamma = np.zeros_like(GammaNew)
    MatrixU = []
    MatrixV = []
    MatrixW = []

    a_temp = []
    aline_temp = []
    r_R_temp = []
    Fnorm_temp = []
    Ftan_temp = []
    Gamma_temp = []

    Niterations = 1200
    errorlimit = 0.01
    ConvWeight = 0.3
    print(len(rings))
    for icp, cp in enumerate(controlpoints):
        MatrixU.append([])
        MatrixV.append([])
        MatrixW.append([])
        for jring, ring in enumerate(rings):
            updated_ring = update_Gamma_single_ring(ring, 1, 1)
            velocity_induced = velocity_induced_single_ring(updated_ring, cp['coordinates'])
            MatrixU[icp].append(velocity_induced[0])
            MatrixV[icp].append(velocity_induced[1])
            MatrixW[icp].append(velocity_induced[2])

    for kiter in range(Niterations):
        Gamma[:] = GammaNew[:]

        for icp, cp in enumerate(controlpoints):
            u = v = w = 0
            radialposition = np.linalg.norm(cp['coordinates'])
            vrot = np.cross([-Omega, 0, 0], cp['coordinates'])
            for jring, _ in enumerate(rings):
                u += MatrixU[icp][jring] * Gamma[jring]
                v += MatrixV[icp][jring] * Gamma[jring]
                w += MatrixW[icp][jring] * Gamma[jring]

            vel1 = np.array(wind) + np.array([u, v, w]) + vrot
            azimdir = np.cross([-1 / radialposition, 0, 0], cp['coordinates'])
            vazim = np.dot(azimdir, vel1)
            vaxial = np.dot([1, 0, 0], vel1)

            temploads = loadBladeElement(vaxial, vazim, radialposition / rotorradius)
            GammaNew[icp] = temploads[2]
            a_temp.append(-((u + vrot[0]) / wind[0]))
            aline_temp.append((vazim / (radialposition * Omega) - 1))
            r_R_temp.append(radialposition / rotorradius)
            Fnorm_temp.append(temploads[0])
            Ftan_temp.append(temploads[1])
            Gamma_temp.append(temploads[2])

        error = np.max(np.abs(GammaNew - Gamma) / np.maximum(np.abs(GammaNew), 0.001))

        if error < errorlimit:
            break

        GammaNew = (1 - ConvWeight) * Gamma + ConvWeight * GammaNew

    return {'a': a_temp, 'aline': aline_temp, 'r_R': r_R_temp, 'Fnorm': Fnorm_temp, 'Ftan': Ftan_temp, 'Gamma': Gamma_temp}

# Testing the functions
span_array = np.linspace(0.2 * Tb.R, Tb.R, 20)
radius = Tb.R
tipspeedratio = 8
Uinf = [10, 0, 0]
theta_array = np.arange(0, 2 * np.pi, 0.1)
nblades = Tb.n
Omega = 30

rotor_wake_system = create_rotor_geometry(span_array, radius, tipspeedratio, Uinf[0], theta_array, nblades)
solution = solve_lifting_line_system_matrix_approach(rotor_wake_system, Uinf, Omega, radius)

# Plotting the results
def plot_a_vs_radial_position(results):
    print(results['r_R'])
    plt.plot(results['r_R'][-19:-1], results['a'][-19:-1], label='Induction factor $a$')
    plt.xlabel('Radial Position (r/R)')
    plt.ylabel('Induction factor $a$')
    plt.title('Induction factor $a$ vs Radial Position')
    plt.grid(True)
    plt.legend()
    plt.show()

plot_a_vs_radial_position(solution)

if __name__ == "__main__":
    print("Lifting Line Method")

