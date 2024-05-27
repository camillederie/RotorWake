import math as m

import numpy as np
from Variables import *

#Discretisation of blade, vortex rings and control points

# X is normal
# Y along the span
# Z is tangential

def spanwise_discretisation(Method, R_Root_Ratio, n_span, n_rotations, n_theta):
    # Function to discretize the spanwise direction of the blade
    # Method: 'Equal' or 'Cosine' for equal or cosine spacing
    # R_Root_Ratio: ratio of root radius to total radius
    # n_span: number of spanwise sections
    # n_rotations: number of rotations for theta

    if Method == 'Equal':
        r = np.linspace(R_Root_Ratio * R, R, n_span + 1)

    elif Method == 'Cosine':
        middle_point = (1 - R_Root_Ratio) * R / 2
        r = np.zeros(n_span + 1)
        r[0] = R_Root_Ratio * R
        r[-1] = R

        angle_ratio = np.pi / (n_span)
        angle = angle_ratio

        for i in range(1, n_span + 1):
            r[i] = R_Root_Ratio * R + middle_point * (1 - np.cos(angle))
            angle += angle_ratio

    theta = np.linspace(0., n_rotations * 2 * m.pi, n_theta * n_rotations)
    return r, theta

def geo_blade(r_R):
    # Function to define the geometry of the blade
    # r_R: non-dimensional radial position

    pitch = 2
    chord = 3 * (1 - r_R) + 1
    twist = -14 * (1 - r_R)
    return [chord, twist + pitch]

def create_rotor_geometry(span_array, R, TSR, Uinf, theta_array, n_blades, a):
    # Function to create the rotor geometry
    # span_array: array of spanwise positions
    # R: rotor radius
    # TSR: tip-speed ratio
    # Uinf: freestream velocity
    # theta_array: array of azimuthal positions
    # n_blades: number of blades
    filaments = []
    ring = []
    controlpoints = []
    bladepanels = []
    TSR /= (1-a)
    for krot in range(n_blades):
        angle_rotation = 2 * m.pi / n_blades * krot
        cosrot = m.cos(angle_rotation)
        sinrot = m.sin(angle_rotation)

        for i in range(len(span_array)-1):
            r = (span_array[i] + span_array[i + 1]) / 2
            geodef = geo_blade(r / R)
            angle = geodef[1] * m.pi / 180

            # Create control points
            temp1 = {
                'coordinates': [0, r, 0],
                'chord': geodef[0],
                'normal': [m.cos(angle), 0, -m.sin(angle)],
                'tangential': [-m.sin(angle), 0, -m.cos(angle)]
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

            # Create vortex rings
            # Bound filaments
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
            # Trailing filament 1
            geodef = geo_blade(span_array[i] / R)
            angle = geodef[1] * m.pi / 180
            temp1 = {
                'x1': geodef[0] * m.sin(-angle),
                'y1': span_array[i],
                'z1': -geodef[0] * m.cos(angle),
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
                dy = (m.cos(-theta_array[j + 1]) - m.cos(-theta_array[j])) * span_array[i]
                dz = (m.sin(-theta_array[j + 1]) - m.sin(-theta_array[j])) * span_array[i]
                dx = (theta_array[j + 1] - theta_array[j]) / TSR * R

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
            #Trailing filament 2
            geodef = geo_blade(span_array[i + 1] / R)
            angle = geodef[1] * m.pi / 180
            temp1 = {
                'x1': 0,
                'y1': span_array[i + 1],
                'z1': 0,
                'x2': geodef[0] * m.sin(-angle),
                'y2': span_array[i + 1],
                'z2': -geodef[0] * m.cos(angle),
                'Gamma': 0
            }
            filaments.append(temp1)

            for j in range(len(theta_array) - 1):
                xt = filaments[-1]['x2']
                yt = filaments[-1]['y2']
                zt = filaments[-1]['z2']
                dy = (m.cos(-theta_array[j + 1]) - m.cos(-theta_array[j])) * span_array[i + 1]
                dz = (m.sin(-theta_array[j + 1]) - m.sin(-theta_array[j])) * span_array[i + 1]
                dx = (theta_array[j + 1] - theta_array[j]) / TSR * R

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

            # Create blade panels
            geodef = geo_blade(span_array[i] / R)
            angle = geodef[1] * m.pi / 180
            geodef2 = geo_blade(span_array[i + 1] / R)
            angle2 = geodef2[1] * m.pi / 180

            temp1 = {
                'p1': [-0.25 * geodef[0] * m.sin(-angle), span_array[i], 0.25 * geodef[0] * m.cos(angle)],
                'p2': [-0.25 * geodef2[0] * m.sin(-angle2), span_array[i + 1], 0.25 * geodef2[0] * m.cos(angle2)],
                'p3': [0.75 * geodef2[0] * m.sin(-angle2), span_array[i + 1], -0.75 * geodef2[0] * m.cos(angle2)],
                'p4': [0.75 * geodef[0] * m.sin(-angle), span_array[i], -0.75 * geodef[0] * m.cos(angle)]
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