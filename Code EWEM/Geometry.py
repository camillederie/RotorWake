import numpy as np


#Discretisation of blade, vortex rings and control points

# X is normal
# Y along the span
# Z is tangential


def spanwise_discretisation(Method, R_Root_Ratio, span_points):

    if Method == 'Equal':
        r_segment = np.linspace(R_Root_Ratio, 1, span_points)

    elif Method == 'Cosine':
       middle_point = (1 - R_Root_Ratio) / 2
       r_segment = np.zeros(span_points)
       r_segment[0] = R_Root_Ratio
       r_segment[-1] = 1

       angle_ratio = np.pi / (span_points - 1)  # Calculate the ratio of current iteration
       angle = angle_ratio

       for i in range(1, span_points):
           r_segment[i] = R_Root_Ratio + middle_point * (1 - np.cos(angle))
           angle += angle_ratio


    span_array = r_segment

    return span_array



#r = spanwise_discretisation(Method='Cosine', R_Root_Ratio=0.2, span_points=10)
#print(r)

def geo_blade(r_R):
    pitch = 2
    chord = 3 * (1 - r_R) + 1
    twist = -14 * (1 - r_R)
    return [chord, twist + pitch]


def blade_discretization(n_blades, spanpoints):

    panels = np.zeros((n_blades * (span_points - 1), 4 * 3))

    for k in range(n_blades):
        angle_rotation = 2 * np.pi / n_blades * k







        r = (span_array[i] + span_array[i + 1]) / 2
        geodef = geo_blade(r / radius)
        angle = geodef[1] * math.pi / 180


        geodef = geo_blade(span_array[i] / radius)
        angle = geodef[1] * math.pi / 180
        geodef2 = geo_blade(span_array[i + 1] / radius)
        angle2 = geodef2[1] * math.pi / 180

        # define the 4 corners
        p1 = [-0.25 * chord1 * np.sin(-angle1), self.span_arr[:-1], 0.25 * chord1 * np.cos(angle1)]
        p2 = [-0.25 * chord2 * np.sin(-angle2), self.span_arr[1:], 0.25 * chord2 * np.cos(angle2)]
        p3 = [0.75 * chord2 * np.sin(-angle2), self.span_arr[1:], -0.75 * chord2 * np.cos(angle2)]
        p4 = [0.75 * chord1 * np.sin(-angle1), self.span_arr[:-1], -0.75 * chord1 * np.cos(angle1)]




    return panels