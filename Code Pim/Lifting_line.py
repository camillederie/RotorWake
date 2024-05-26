import matplotlib.backends
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use('TkAgg')

from experimental_data import Data
# file import
from turbine2 import Turbine

""""
Experimental data
"""


def C_L(Tb, alpha):
    return Data(Tb.type).lift(np.degrees(alpha))


def C_D(Tb, alpha):
    return Data(Tb.type).drag(np.degrees(alpha))


def plot_cl_cd_vs_alpha(Tb, CL=C_L, CD=C_D):
    alpha_range_rad = np.linspace(np.radians(-16), np.radians(30), 100)

    # Calculate C_L and C_D for each alpha
    cl_values = [CL(Tb, alpha) for alpha in alpha_range_rad]
    cd_values = [CD(Tb, alpha) for alpha in alpha_range_rad]
    plt.plot(alpha_range_rad * 180 / np.pi, cl_values, label='$C_L$')
    plt.xlabel('Angle of Attack (degrees)')
    # Plot C_D vs. Alpha
    plt.plot(alpha_range_rad * 180 / np.pi, cd_values, label='$C_D$')
    plt.title(r'$C_L$ and $C_D$ vs. $\alpha$')
    plt.grid(True)
    plt.legend()
    # plt.show()


def lift(Tb, alpha, chord, velocity, rho):
    """Determine the sectional 2D lift force of the blade section"""
    return C_L(Tb, alpha) * 0.5  * chord * velocity ** 2


def drag(Tb, alpha, chord, velocity, rho):
    """Determine the sectional 2D lift force of the blade section"""
    return C_D(Tb, alpha) * 0.5  * chord * velocity ** 2

""""
Operational calculations
"""

def blade_geom(Tb, blade, theta, TSR):
    """"
    Tb: Turbine class
    blade: Blade number
    Theta: Previous angular positions of a single blade
    TSR: Tip Speed Ratio
    """
    # TSR/=0.8
    rotation_angle = blade / Tb.n * 2 * np.pi
    cosrot = np.cos(rotation_angle)
    sinrot = np.sin(rotation_angle)
    chord, beta = Tb.geometry(Tb.r_cp/Tb.R)

    cp = np.dstack((np.dstack((np.zeros((Tb.Nrad)), Tb.r_cp * cosrot)), Tb.r_cp * sinrot))[0]
    normal = np.dstack((np.dstack((np.cos(beta), np.sin(beta) * sinrot)), -np.sin(beta) * cosrot))[0]
    tangential = np.dstack((np.dstack((-np.sin(beta), np.cos(beta) * sinrot)), -np.cos(beta) * cosrot))[0]
    # Define the trailing locations
    bb, aa = Tb.r_ab[1:], Tb.r_ab[:-1]

    a = np.dstack((np.dstack((np.zeros((Tb.Nrad)), aa)), np.zeros((Tb.Nrad))))[0]
    b = np.dstack((np.dstack((np.zeros((Tb.Nrad)), bb)), np.zeros((Tb.Nrad))))[0]
    # a = np.dstack((np.dstack((np.zeros((Tb.Nrad)), aa * cosrot)), aa * sinrot))[0]
    # b = np.dstack((np.dstack((np.zeros((Tb.Nrad)), bb * cosrot)), bb * sinrot))[0]
    # Find trailing geometry parameters
    a_chord, a_beta = Tb.geometry(aa/Tb.R)
    b_chord, b_beta = Tb.geometry(bb/Tb.R)


    # create trailing filaments
    a_c = np.dstack((np.dstack((-a_chord * np.sin(a_beta), aa)), - a_chord * np.cos(a_beta)))[0]
    b_c = np.dstack((np.dstack((-b_chord * np.sin(b_beta), bb)), - b_chord * np.cos(b_beta)))[0]

    # rotate a_c and b_c
    # for point in range(Tb.Nrad):
    #     a_c[point][1] = a_c[point][1] * cosrot - a_c[point][2] * sinrot
    #     a_c[point][2] = a_c[point][2] * cosrot + a_c[point][1] * sinrot
    #     b_c[point][1] = b_c[point][1] * cosrot - b_c[point][2] * sinrot
    #     b_c[point][2] = b_c[point][2] * cosrot + b_c[point][1] * sinrot
    # create wake filaments
    filaments = np.array([[a_c[i], a[i], b[i], b_c[i]] for i in range(Tb.Nrad)])

    for i in range(len(theta) - 1):
        dx = (theta[i + 1] - theta[i]) / TSR * Tb.R
        dy = lambda radius: radius * (np.cos(-theta[i + 1]) - np.cos(-theta[i]))
        dz = lambda radius: radius * (np.sin(-theta[i + 1]) - np.sin(-theta[i]))

        a_last = np.array([filaments[j][0] for j in range(Tb.Nrad)])
        b_last = np.array([filaments[j][-1] for j in range(Tb.Nrad)])

        # update x, y and z values
        for k in range(Tb.Nrad):
            a_last[k][0] += dx
            a_last[k][1] += dy(aa[k])
            a_last[k][2] += dz(aa[k])
            b_last[k][0] += dx
            b_last[k][1] += dy(bb[k])
            b_last[k][2] += dz(bb[k])
        filaments = np.array([np.concatenate((a_last[i].reshape(1, 3), filaments[i], b_last[i].reshape(1, 3)))
                              for i in range(Tb.Nrad)])

    # rotate all points
    copy = filaments
    for r, ring in enumerate(filaments):
        for p, point in enumerate(ring):
            temp1 = [copy[r][p][0], copy[r][p][1], copy[r][p][2]]
            point[1] = temp1[1] * cosrot - temp1[2] * sinrot
            point[2] = temp1[2] * cosrot + temp1[1] * sinrot

    return cp, normal, tangential, filaments

def rotor_geom(Tb, TSR):
    # TSR/0.8
    cp = np.empty((0, 3))
    normal = np.empty((0, 3))
    tangential = np.empty((0, 3))
    filaments = np.empty((0, (len(Tb.wake_steps) + 1) * 2, 3))  # empty, horseshoe size, [x,y,z]

    for blade in range(Tb.n):
        blade_cp, blade_normal, blade_tangential, blade_filaments = blade_geom(Tb, blade, Tb.wake_steps, TSR)
        cp = np.vstack((cp, blade_cp))
        normal = np.vstack((normal, blade_normal))
        tangential = np.vstack((tangential, blade_tangential))
        filaments = np.vstack((filaments, blade_filaments))
    #store filaments to txt file
    if TSR == 8:
        np.savetxt("filaments_pim", cp.reshape((3,-1)), fmt="%s")
    return cp, normal, tangential, filaments

def velocity(U_inf, velocity, n_azim, radius, omg):
    """"
    Determine the velocity components:
        U_inf: freestream velocity
        velocity: Local wake velocity vector [u, v, w]
        n_azim: normal vector from the blade [x, y, z]
        radius: Rotor radius
        omega: rotational velocity
    """
    v_axial = U_inf * (1 - velocity[0])  # freestream + induced
    v_tan = omg * radius + np.cross((U_inf * np.array([1, 0, 0]) + velocity), n_azim)  # tangential velocity + induced
    v_p = np.sqrt(v_axial ** 2 + v_tan ** 2)

    return v_axial, v_tan, v_p


def angles(beta, U_inf, velocity, n_azim, radius, omega):
    """"
    Determine the associated angles:
        beta: geometric blade angle
        U_inf: freestream velocity
        velocity: Local wake velocity vector [u, v, w]
        n_azim: normal vector from the blade [x, y, z]
        radius: Rotor radius
    """
    v_axial, v_tan, v_p = velocity(U_inf, velocity, n_azim, radius, omega)
    phi = np.arctan2(v_axial, v_tan)  # incoming flow angle
    alpha = phi - beta  # Flow angle - Geometric angle
    # Ensure alpha does not exceed 30.0 or -16 degrees
    alpha[alpha > np.radians(30)] = np.radians(30)
    alpha[alpha < np.radians(-16)] = np.radians(-16)
    return phi, alpha


def forces(Tb, phi, alpha, chord, velocity, rho, CL=C_L):
    """"
    2D Azimuthal and axial forces. [N/m]
        phi: inflow angle angle
        alpha: angle of attack
        chord: local chord
        velocity: Local wake velocity vector [u, v, w]
        rho: density
        CL: Lift coefficient determined with alpha
    """
    angular = lift(Tb, alpha, chord, velocity, rho) * np.sin(phi) - \
              drag(Tb, alpha, chord, velocity, rho) * np.cos(phi)
    axial = lift(Tb, alpha, chord, velocity, rho) * np.cos(phi) + \
            drag(Tb, alpha, chord, velocity, rho) * np.sin(phi)
    gamma = CL(Tb, alpha) * velocity * chord / 2
    return axial, angular, gamma


def single_filament(a, b, p, core=0.0001):
    """"
    Calculate velocity at Blade segment
        gamma: vorticity
        a: start point of a filament [x, y, z]
        b: end point of a filament [x, y, z]
        p: control point
        core: filament radius
    """
    x1, y1, z1 = a  # Start vortex filament
    x2, y2, z2 = b  # Endpoint vortex filament
    xp, yp, zp = p  # location p
    # Cross product
    R1 = np.sqrt(np.power((xp - x1), 2) + np.power((yp - y1), 2) + np.power((zp - z1), 2))
    R2 = np.sqrt(np.power((xp - x2), 2) + np.power((yp - y2), 2) + np.power((zp - z2), 2))
    crossi = (yp - y1) * (zp - z2) - (zp - z1) * (yp - y2)
    crossj = (zp - z1) * (xp - x2) - (xp - x1) * (zp - z2)
    crossk = (xp - x1) * (yp - y2) - (yp - y1) * (xp - x2)
    product = np.power(crossi, 2) + np.power(crossj, 2) + np.power(crossk, 2)
    R01 = (x2 - x1) * (xp - x1) + (y2 - y1) * (yp - y1) + (z2 - z1) * (zp - z1)
    R02 = (x2 - x1) * (xp - x2) + (y2 - y1) * (yp - y2) + (z2 - z1) * (zp - z2)
    product = np.power(core, 2) if (product < np.power(core, 2)) else product
    R1 = core if (R1 < core) else R1
    R2 = core if (R2 < core) else R2
    K = 1 / (4 * np.pi * product) * (R01 / R1 - R02 / R2)
    u = K * crossi
    v = K * crossj
    w = K * crossk

    return u, v, w


def horseshoe(p, points):
    """"
    Calculate the total induced velocity of a horshoe ring
        p: control point
        points: start and end points of filaments
    """
    u, v, w = 0, 0, 0
    for i in range(len(points) - 1):
        a = points[i]
        b = points[i + 1]
        u_fil, v_fil, w_fil = single_filament(a, b, p)
        u += u_fil
        v += v_fil
        w += w_fil

    return u, v, w


def Lifting_Line_Method(Tb, TSR, Uinf, omega, rho, tol=1, iterations=1200):
    """
    Iterative process to calculate circulation and forces on the rotor blades.
    """
    cp, normal, tangential, filaments = rotor_geom(Tb, TSR)

    pos = 3
    for blade in np.arange(pos, 30, 10):
        fil = filaments[blade]
        x, y, z = [], [], []
        for coord in fil:
            x.append(coord[0])
            y.append(coord[1])
            z.append(coord[2])
        plt.plot(y, z, "b")

        # local blade horseshoe
        bound_loc = int(len(x)/2)
        plt.plot(y[bound_loc-2:bound_loc+2], z[bound_loc-2:bound_loc+2], "r")
        plt.plot(cp[blade][1], cp[blade][2], "k.")

    # plt.show()

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    # for ring in filaments:
    #     for point in range(1):
    #         x = [ring[point][0], ring[point+1][0]]
    #         y = [ring[point][1], ring[point+1][1]]
    #         z = [ring[point][2], ring[point+1][2]]
    #         # print(x, y, z)
    #         ax.plot(x, y, z)
    #
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    #
    # # plt.show()

    size = (Tb.n * Tb.Nrad, Tb.n * Tb.Nrad)
    U_matrix = np.zeros(size)
    V_matrix = np.zeros(size)
    W_matrix = np.zeros(size)
    for icp in range(Tb.n * Tb.Nrad):
        for jring in range(Tb.n * Tb.Nrad):
            u, v, w = horseshoe(p=cp[icp], points=filaments[jring])
            U_matrix[icp][jring] = u
            V_matrix[icp][jring] = v
            W_matrix[icp][jring] = w
    # save U_matrix to txt file
    np.savetxt(f'U_matrix_{omega}_pim.txt', U_matrix, fmt='%1.4e')
    error = 1.0
    ConvWeight = 0.2

    # initialize
    it = 0
    gamma_new = np.zeros(size[0])  # Setting the initial circulation to zero at each control point
    converging = True

    # Initialize a and aline
    a = np.zeros(size[0])
    a_prime = np.zeros(size[0])

    while converging is True:
        it += 1
        gamma = gamma_new  # Update the bound circulation with new estimate
        radial_position = np.tile(Tb.r_cp, Tb.n)
        u = np.dot(U_matrix, gamma)
        v = np.dot(V_matrix, gamma)
        w = np.dot(W_matrix, gamma)

        v_rot = np.cross(np.array([-omega, 0, 0]), cp)
        vel = np.column_stack((u, v, w))
        v_eff = vel + v_rot + np.array([Uinf, 0, 0])

        rad_coord = np.column_stack((-1 / radial_position, np.zeros(size[0]), np.zeros(size[0])))
        tangential_direction = np.array([np.cross(rad_coord[j], cp[j]) for j in range(size[0])])

        v_tan = np.array([np.dot(tangential_direction[j], v_eff[j]) for j in range(size[0])])
        v_norm = np.array([np.dot([1, 0, 0], v_eff[j]) for j in range(size[0])])

        inflow_angle = np.arctan(v_norm / v_tan)
        chord, twist = Tb.geometry(radial_position / Tb.R)
        AoA = twist + inflow_angle
        AoA[AoA < np.radians(-16)] = np.radians(-16)
        AoA[AoA > np.radians(30)] = np.radians(30)

        v_eff_new = np.sqrt(v_norm ** 2 + v_tan ** 2)
        F_norm, F_tan, gamma_check = forces(Tb, inflow_angle, AoA, chord, v_eff_new, rho)

        # Calculate a and aline
        for icp in range(size[0]):
            a[icp] = -((u[icp] + v_rot[icp][0]) / Uinf)
            a_prime[icp] = (v_tan[icp] / (radial_position[icp] * omega) - 1)

        if np.all(np.abs((gamma_check - gamma)) < tol):
            gamma_new = gamma_check
            print("converged after:", it, "iterations")
            converging = False
        else:
            gamma_new = (1 - ConvWeight) * gamma + ConvWeight * gamma_check

        # Visual checks
        if it >= iterations:
            raise RuntimeError(f"Too many iterations")
        if it % 50 == 0:
            print(it)

    return gamma_new, AoA, F_norm, F_tan, inflow_angle, a, a_prime, u , np.linalg.norm(v_eff, axis=1), v_norm, v_tan

def plots():
    T = Turbine(10, 10)
    TSR = [6, 8, 10]
    U_inf = 10
    rho = 1.225

    circulations = []
    angles_of_attack = []
    inflow_angles = []
    normal_forces = []
    tangential_forces = []
    a_list = []
    u_list = []
    V_tot_mag_list = []
    v_norm_list = []
    v_tan_list = [] 
    for tsr in TSR:
        omega = U_inf * tsr / T.R
        gamma, alpha, Normal_force, Tangential_force, phi, a, aprime, u, V_tot_mag, v_norm, v_tan = Lifting_Line_Method(Tb=T,
                                                                                           TSR=tsr,
                                                                                           Uinf=U_inf,
                                                                                           omega=omega,
                                                                                           rho=rho)


        # Normalize the results
        gamma_normalized = [g / (np.pi * U_inf ** 2 / (omega * T.n)) for g in gamma]
        alpha_degrees = [a * 180 / np.pi for a in alpha]
        phi_degrees = [p * 180 / np.pi for p in phi]
        Normal_force_normalized = [nf / (0.5 * rho * U_inf ** 2 * T.R) for nf in Normal_force]
        Tangential_force_normalized = [tf / (0.5 * rho * U_inf ** 2 * T.R) for tf in Tangential_force]

        circulations.append(gamma_normalized)
        angles_of_attack.append(alpha_degrees)
        inflow_angles.append(phi_degrees)
        normal_forces.append(Normal_force_normalized)
        tangential_forces.append(Tangential_force_normalized)
        a_list.append(a)
        u_list.append(u)
        V_tot_mag_list.append(V_tot_mag)
        v_norm_list.append(v_norm)
        v_tan_list.append(v_tan)
    r_cp_normalized = T.r_cp / T.R

    # Plot circulation over blade for all TSRs
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, circulations[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Circulation over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Circulation')
    # plt.show()

    # Plot angle of attack over blade for all TSRs
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, angles_of_attack[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Angle of Attack over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Angle of Attack (degrees)')
    # plt.show()

    # Plot inflow angle over blade for all TSRs
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, inflow_angles[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Inflow Angle over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Inflow Angle (degrees)')
    # plt.show()

    # Plot normal force over blade for all TSRs
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, normal_forces[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Normal Force over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Normal Force (normalized)')
    # plt.show()

    # Plot tangential force over blade for all TSRs
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, tangential_forces[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Tangential Force over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Tangential Force (normalized)')
    # plt.show()

    # plot a 
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, a_list[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Axial Induction over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Axial Induction Factor')
    # plt.show()

    # plot u
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, u_list[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Induction Velocity U over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Induction Velocity U [m/s]')
    # plt.show()
    
    # plot V_tot_mag
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, V_tot_mag_list[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Total Velocity over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Total Velocity [m/s]')
    # plt.show()

    # plot v_norm and v_tan
    for i, tsr in enumerate(TSR):
        plt.plot(r_cp_normalized, v_norm_list[i][:T.Nrad], label=f'TSR {tsr}')
        plt.plot(r_cp_normalized, v_tan_list[i][:T.Nrad], label=f'TSR {tsr}')
    plt.legend()
    plt.title("Normal and Tangential Velocity over blade")
    plt.grid()
    plt.xlabel('r/R')
    plt.ylabel('Velocity [m/s]')
    # # plt.show()
    

if __name__ == "__main__":
    print("Lifting Line Method")
    plots()