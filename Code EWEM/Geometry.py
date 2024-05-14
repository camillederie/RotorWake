import numpy as np


#Discretisation of blade, vortex rings and control points

# X is normal
# Y along the span
# Z is tangential


def spanwise_discretisation(Method, R_Root_Ratio, span_points):

    if Method == 'Equal':
        r = np.linspace(R_Root_Ratio, 1, span_points)

    elif Method == 'Cosine':
       middle_point = (1 - R_Root_Ratio) / 2
       r = np.zeros(span_points)
       r[0] = R_Root_Ratio
       r[-1] = 1

       angle_ratio = np.pi / (span_points - 1)  # Calculate the ratio of current iteration
       angle = angle_ratio

       for i in range(1, span_points):
           r[i] = R_Root_Ratio + middle_point * (1 - np.cos(angle))
           angle += angle_ratio

    return r


def geo_blade(r_R):
    pitch = 2
    chord = 3 * (1 - r_R) + 1
    twist = -14 * (1 - r_R)
    return [chord, twist + pitch]


def blade_discretization(n_blades, spanpoints):

    panels = np.zeros((n_blades * (span_points - 1), 4 * 3))

    for blade in range(n_blades):
        geodef = geo_blade(span_array[i] / radius)
        angle = geodef[1] * math.pi / 180
        geodef2 = geo_blade(span_array[i + 1] / radius)
        angle2 = geodef2[1] * math.pi / 180


