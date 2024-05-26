import math as m

import numpy as np
from Variables import *

#Discretisation of blade, vortex rings and control points

# X is normal
# Y along the span
# Z is tangential

def spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations):
    # Function to discretize the spanwise direction of the blade
    # Method: 'Equal' or 'Cosine' for equal or cosine spacing
    # R_Root_Ratio: ratio of root radius to total radius
    # n_span: number of spanwise sections
    # n_rotations: number of rotations for theta

    if Method == 'Equal':
        r = np.linspace(R_Root_Ratio * R, R, n_span)

    elif Method == 'Cosine':
        middle_point = (1 - R_Root_Ratio) * R / 2
        r = np.zeros(n_span)
        r[0] = R_Root_Ratio * R
        r[-1] = R

        angle_ratio = np.pi / (n_span - 1)
        angle = angle_ratio

        for i in range(1, n_span):
            r[i] = R_Root_Ratio * R + middle_point * (1 - np.cos(angle))
            angle += angle_ratio

    theta = np.arange(0., n_rotations * 2 * m.pi, m.pi / 10)
    return r, theta

def geo_blade(r_R):
    # Function to define the geometry of the blade
    # r_R: non-dimensional radial position

    pitch = 2
    chord = 3 * (1 - r_R) + 1
    twist = -14 * (1 - r_R)
    return [chord, twist + pitch]

def create_blade_geometry(span_array, R, TSR, Uinf, theta_array, n_blades, blade):
    """"
    n_blades: Turbine class
    blade: Blade number
    Theta: Previous angular positions of a single blade
    TSR: Tip Speed Ratio
    """
    # TSR/=0.8
    r_cp = np.zeros(len(span_array) - 1)
    beta = np.zeros(len(span_array) - 1)
    Nrad = len(span_array)-1
    # for krot in range(n_blades):
    angle_rotation = 2 * m.pi / n_blades * blade
    cosrot = m.cos(angle_rotation)
    sinrot = m.sin(angle_rotation)

    for i in range(len(span_array) - 1):
        r_cp[i] = (span_array[i] + span_array[i + 1]) / 2
        geodef = geo_blade(r_cp[i] / R)
        beta[i] = geodef[1] * m.pi / 180

    cp = np.dstack((np.dstack((np.zeros((Nrad)), r_cp * cosrot)), r_cp * sinrot))[0]
    normal = np.dstack((np.dstack((np.cos(beta), np.sin(beta) * sinrot)), -np.sin(beta) * cosrot))[0]
    tangential = np.dstack((np.dstack((-np.sin(beta), np.cos(beta) * sinrot)), -np.cos(beta) * cosrot))[0]
    # Define the trailing locations
    bb, aa = span_array[1:], span_array[:-1]

    a = np.dstack((np.dstack((np.zeros((Nrad)), aa)), np.zeros((Nrad))))[0]
    b = np.dstack((np.dstack((np.zeros((Nrad)), bb)), np.zeros((Nrad))))[0]
    # a = np.dstack((np.dstack((np.zeros((n_blades.Nrad)), aa * cosrot)), aa * sinrot))[0]
    # b = np.dstack((np.dstack((np.zeros((n_blades.Nrad)), bb * cosrot)), bb * sinrot))[0]
    # Find trailing geometry parameters
    a_chord, a_beta = geo_blade(aa/R)
    b_chord, b_beta = geo_blade(bb/R)
    a_beta, b_beta = np.radians(a_beta), np.radians(b_beta)


    # create trailing filaments
    a_c = np.dstack((np.dstack((-a_chord * np.sin(a_beta), aa)), - a_chord * np.cos(a_beta)))[0]
    b_c = np.dstack((np.dstack((-b_chord * np.sin(b_beta), bb)), - b_chord * np.cos(b_beta)))[0]

    # rotate a_c and b_c
    # for point in range(n_blades.Nrad):
    #     a_c[point][1] = a_c[point][1] * cosrot - a_c[point][2] * sinrot
    #     a_c[point][2] = a_c[point][2] * cosrot + a_c[point][1] * sinrot
    #     b_c[point][1] = b_c[point][1] * cosrot - b_c[point][2] * sinrot
    #     b_c[point][2] = b_c[point][2] * cosrot + b_c[point][1] * sinrot
    # create wake filaments
    filaments = np.array([[a_c[i], a[i], b[i], b_c[i]] for i in range(Nrad)])

    for i in range(len(theta_array) - 1):
        dx = (theta_array[i + 1] - theta_array[i]) / TSR * R
        dy = lambda radius: radius * (np.cos(-theta_array[i + 1]) - np.cos(-theta_array[i]))
        dz = lambda radius: radius * (np.sin(-theta_array[i + 1]) - np.sin(-theta_array[i]))

        a_last = np.array([filaments[j][0] for j in range(Nrad)])
        b_last = np.array([filaments[j][-1] for j in range(Nrad)])

        # update x, y and z values
        for k in range(Nrad):
            a_last[k][0] += dx
            a_last[k][1] += dy(aa[k])
            a_last[k][2] += dz(aa[k])
            b_last[k][0] += dx
            b_last[k][1] += dy(bb[k])
            b_last[k][2] += dz(bb[k])
        filaments = np.array([np.concatenate((a_last[i].reshape(1, 3), filaments[i], b_last[i].reshape(1, 3)))
                              for i in range(Nrad)])

    # rotate all points
    copy = filaments
    for r, ring in enumerate(filaments):
        for p, point in enumerate(ring):
            temp1 = [copy[r][p][0], copy[r][p][1], copy[r][p][2]]
            point[1] = temp1[1] * cosrot - temp1[2] * sinrot
            point[2] = temp1[2] * cosrot + temp1[1] * sinrot

    return cp, normal, tangential, filaments

def create_rotor_geometry(n_blades, TSR, n_rotations, span_array, R, Uinf, theta_array):
    # TSR/0.8
    wake_steps = np.arange(0, n_rotations*2*np.pi, np.pi/10)
    cp = np.empty((0, 3))
    normal = np.empty((0, 3))
    tangential = np.empty((0, 3))
    filaments = np.empty((0, (len(wake_steps) + 1) * 2, 3))  # empty, horseshoe size, [x,y,z]

    for blade in range(n_blades):
        blade_cp, blade_normal, blade_tangential, blade_filaments = create_blade_geometry(span_array, R, TSR, Uinf, theta_array, n_blades, blade)
        cp = np.vstack((cp, blade_cp))
        normal = np.vstack((normal, blade_normal))
        tangential = np.vstack((tangential, blade_tangential))
        filaments = np.vstack((filaments, blade_filaments))

    #store filaments to txt file
    if TSR == 8:
        np.savetxt("filaments_louis", filaments.reshape((3,-1)), fmt="%s")
    # print(cp)
    return cp, normal, tangential, filaments