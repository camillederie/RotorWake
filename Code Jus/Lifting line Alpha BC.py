import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rotor_wake_turbine_debug import main_bem
#from plots import plots
from scipy.optimize import fsolve, root

#set directory
#os.chdir('C:\Users\camil\OneDrive\Documenten\TU DELFT MSC\Rotor and Wake Aerodynamics\Assignment 2\Code Jus')


def equal(result_dict1, result_dict2):
    if result_dict1.all() == result_dict2.all():
        print('equal')
    else:
        print("not equal")

def get_induced_matrices(p1, p2, d,printen=False):
    X1 = np.zeros((len(d), len(p1))) + p1[:, 0].flatten()
    Y1 = np.zeros((len(d), len(p1))) + p1[:, 1].flatten()
    Z1 = np.zeros((len(d), len(p1))) + p1[:, 2].flatten()

    X2 = np.zeros((len(d), len(p2))) + p2[:, 0].flatten()
    Y2 = np.zeros((len(d), len(p2))) + p2[:, 1].flatten()
    Z2 = np.zeros((len(d), len(p2))) + p2[:, 2].flatten()

    XP = np.tile(d[:, 0].flatten(), (len(p2), 1)).T
    YP = np.tile(d[:, 1].flatten(), (len(p2), 1)).T
    ZP = np.tile(d[:, 2].flatten(), (len(p2), 1)).T
    R1 = np.sqrt((XP - X1) ** 2 + (YP - Y1) ** 2 + (ZP - Z1) ** 2)
    R2 = np.sqrt((XP - X2) ** 2 + (YP - Y2) ** 2 + (ZP - Z2) ** 2)
    R1XR2_X = (YP - Y1) * (ZP - Z2) - (ZP - Z1) * (YP - Y2)
    R1XR2_Y = -(XP - X1) * (ZP - Z2) + (ZP - Z1) * (XP - X2)
    R1XR2_Z = (XP - X1) * (YP - Y2) - (YP - Y1) * (XP - X2)
    R1XR_SQR = R1XR2_X ** 2 + R1XR2_Y ** 2 + R1XR2_Z ** 2
    R0R1 = (X2 - X1) * (XP - X1) + (Y2 - Y1) * (YP - Y1) + (Z2 - Z1) * (ZP - Z1)
    R0R2 = (X2 - X1) * (XP - X2) + (Y2 - Y1) * (YP - Y2) + (Z2 - Z1) * (ZP - Z2)
    R1XR_SQR = np.where(R1XR_SQR<10**-10,0,R1XR_SQR)
    K = 1 / (4 * np.pi * R1XR_SQR) * (R0R1 / R1 - R0R2 / R2)
    K = np.where(R1XR_SQR==0,0,K)
    K[np.isnan(K)] = 0
    if printen == True:
        #print(R1XR_SQR)
        #print(K)
        pass
    U = K * R1XR2_X
    V = K * R1XR2_Y
    W = K * R1XR2_Z
    return U, V, W

def get_bound_coordinates(r_R_list, R, location_mid, theta):
    location_end = np.array(
        [location_mid[0], location_mid[1] + np.cos(theta) * R, location_mid[2] + np.sin(theta) * R])

    delta_loc = location_end - location_mid
    p1_bound = np.zeros((len(r_R_list) - 1, 3))
    p2_bound = np.zeros((len(r_R_list) - 1, 3))

    p1_bound[:, 0] = ((r_R_list[:-1]).T) * delta_loc[0] + location_mid[0]
    p1_bound[:, 1] = ((r_R_list[:-1]).T) * delta_loc[1] + location_mid[1]
    p1_bound[:, 2] = ((r_R_list[:-1]).T) * delta_loc[2] + location_mid[2]

    p2_bound[:, 0] = ((r_R_list[1:]).T) * delta_loc[0] + location_mid[0]
    p2_bound[:, 1] = ((r_R_list[1:]).T) * delta_loc[1] + location_mid[1]
    p2_bound[:, 2] = ((r_R_list[1:]).T) * delta_loc[2] + location_mid[2]
    return p1_bound, p2_bound

def get_trailing_coordinates(r_R_list, R, chord_distribution, twist, location_mid, theta):

    p1_trailing = np.zeros((len(r_R_list), 3))
    p2_trailing = np.zeros((len(r_R_list), 3))

    # y coordinate
    p1_trailing[:, 1] = (r_R_list * R).T
    p2_trailing[:, 1] = (r_R_list * R).T

    # x coordinate
    p2_trailing[:, 0] = (chord_distribution * np.sin(twist * np.pi / 180)).T

    # z coordinate
    p2_trailing[:, 2] = (chord_distribution * - np.cos(twist * np.pi / 180)).T

    #### rotate and translate
    p1_trailing_moved = np.zeros((len(r_R_list), 3))
    p2_trailing_moved = np.zeros((len(r_R_list), 3))

    # x coordinate
    p1_trailing_moved[:, 0] = p1_trailing[:, 0] + location_mid[0]
    p2_trailing_moved[:, 0] = p2_trailing[:, 0] + location_mid[0]

    # y coordinate
    p1_trailing_moved[:, 1] = (p1_trailing[:, 1] * np.cos(theta)) + location_mid[1]
    p2_trailing_moved[:, 1] = (p2_trailing[:, 1] * np.cos(theta)) - (p2_trailing[:, 2] * np.
                                                                     sin(theta)) + location_mid[1]
    # z coordinate
    p1_trailing_moved[:, 2] = (p1_trailing[:, 1] * np.sin(theta)) + location_mid[2]
    p2_trailing_moved[:, 2] = (p2_trailing[:, 1] * np.sin(theta)) + (p2_trailing[:, 2] * np.cos(theta)) + location_mid[
        2]

    return p1_trailing_moved, p2_trailing_moved

def get_control_coordinates_permebility(r_R_list_control, R, chord_distribution_control,
    twist_control, location_mid, theta):

    p_control = np.zeros((len(r_R_list_control), 3))

    p_control[:, 1] = (r_R_list_control * R).T

    # x coordinate
    p_control[:, 0] = (chord_distribution_control * 0.5 * np.sin(twist_control * np.pi / 180)
                       ).T
    # z coordinate
    p_control[:, 2] = (chord_distribution_control * 0.5 * - np.cos(twist_control * np.pi / 180)).T
    #### rotate and translate
    p_control_moved = np.zeros((len(r_R_list_control), 3))

    # x coordinate
    p_control_moved[:, 0] = p_control[:, 0] + location_mid[0]
    # y coordinate
    p_control_moved[:, 1] = (p_control[:, 1] * np.cos(theta)) - (p_control[:, 2] * np.sin( theta)) + location_mid[1]
    # z coordinate
    p_control_moved[:, 2] = (p_control[:, 1] * np.sin(theta)) + (p_control[:, 2] * np.cos( theta)) + location_mid[2]
    return p_control_moved

def get_control_coordinates(r_R_list_control,R,location_mid,theta):
    location_end = np.array([location_mid[0],location_mid[1]+np.cos(theta) * R,location_mid [2] + np.sin(theta) * R])
    delta_loc = location_end-location_mid

    d = np.zeros((len(r_R_list_control),3))

    d[:, 0] = ((r_R_list_control).T) * delta_loc[0] + location_mid[0]
    d[:, 1] = ((r_R_list_control).T) * delta_loc[1] + location_mid[1]
    d[:, 2] = ((r_R_list_control).T) * delta_loc[2] + location_mid[2]
    return d
def get_control_orientation(d,p_control):
    diff = p_control - d
    diff_mag = np.reshape(np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2),(len(diff),1))
    orientation = diff / diff_mag
    return orientation

def bound_to_trailing_circulation(bound_circulation):
    N = len(bound_circulation)
    convergence_matrix = np.eye(N + 1, N, k=-1,dtype=int) - np.eye(N + 1, N, dtype=int)

    return convergence_matrix
def trailing_to_bound_circulation(trailing_circulation):
    N = len(trailing_circulation)
    convergence_matrix = -np.tri(N-1,N,dtype=int)
    return convergence_matrix

def get_initial_rotor_data(location,R,r_R_list,theta_base,n_blades,chord_distribution,twist, r_R_list_control,chord_distribution_control,twist_control,permebility=False):
    blade_pitch_angle = 2*np.pi/n_blades
    rotor_data = {"theta" : [],
                  "p1_bound": [],
                  "p2_bound": [],
                  "p1_trailing": [],
                  "p2_trailing": [], "d": [],
                  "p_control": [],
                  "orientation": []}
    for i in range(n_blades):
        theta = theta_base + blade_pitch_angle * i
        p1_bound, p2_bound = get_bound_coordinates(r_R_list,R,location,theta)
        p1_trailing, p2_trailing = get_trailing_coordinates(r_R_list, R, chord_distribution, twist, location, theta)
        d = get_control_coordinates(r_R_list_control,R,location,theta)
        if permebility == True:
            p_control = get_control_coordinates_permebility(r_R_list_control,R, chord_distribution_control,twist_control,location,theta)
            orientation = get_control_orientation(d,p_control)
            rotor_data["p_control"].append(p_control)
            rotor_data["orienttion"].append(orientation)
            rotor_data["d"].append(p_control)
        else:
            rotor_data["d"].append(d)
        rotor_data["theta"].append(theta)
        rotor_data["p1_bound"].append(p1_bound)
        rotor_data["p2_bound"].append(p2_bound)
        rotor_data["p1_trailing"].append(p1_trailing)
        rotor_data["p2_trailing"].append(p2_trailing)
    return rotor_data

def init_wake():
    wake_dict = {}
    wake_dict["position"] = []

    return wake_dict

def get_wake2(rotor_data,wake_dict,a_w,r_R_list,v_inf,dt,v_wake_factor):
    # Move wake
    avg_a_w = np.average(np.average(a_w ,weights = (0.5*(r_R_list[1:]-r_R_list[:-1])).flatten
    (),axis=2))
    v_wake = v_inf * (1 - avg_a_w) * v_wake_factor
    for i in range(len(wake_dict["position"])):
        for j in range(len(wake_dict["position"][0])):
            for k in range(len(wake_dict["position"][0][0])):
                wake_dict["position"][i][j][k][0] = rotor_data["p2_trailing"][j][k,0] + v_wake * dt * (len(wake_dict["position"])- i -1)
    return wake_dict

def append_wake(rotor_data,wake_dict):
    # append wake
    wake_dict["position"].append(rotor_data["p2_trailing"])
    return wake_dict

def get_induced_velocity_matrices(rotor_data,wake_data):
    #Each row contains data to get induced velocity at specific rotor
    data_tab = []
    for i in range(len(rotor_data)):
        sub_data_tab = []
        rotor_dict_d = rotor_data[i]
        for j in range(len(rotor_data)):
            rotor_dict_circulation = rotor_data[j]
            wake_dict = wake_data[j]
            sub_data_tab.append(get_induced_velocity_matrix(rotor_dict_d,wake_dict, rotor_dict_circulation))

        data_tab.append(sub_data_tab)
    return np.array(data_tab)

def get_induced_velocity_matrix(rotor_data_d,wake_dict,rotor_data_circulation):
    #print("test")
    data_tab = []
    convergence_matrix = trailing_to_bound_circulation(np.zeros(len(rotor_data_d["d"][0])+1))
    for i in range(len(rotor_data_d["d"])):
        sub_data_tab = []
        d = rotor_data_d["d"][i]
        #print("d",d)
        for i_blade in range(len(rotor_data_circulation["p1_trailing"])):
            u_wake, v_wake, w_wake = 0, 0, 0
            for k in range(len(wake_dict["position"]) - 1):
                p1 = wake_dict["position"][k + 1][i_blade]
                p2 = wake_dict["position"][k][i_blade]

                u_temp, v_temp, w_temp = get_induced_matrices(p1, p2, d)
                u_wake += u_temp
                v_wake += v_temp
                w_wake += w_temp

                u_trailing, v_trailing, w_trailing = get_induced_matrices(rotor_data_circulation[
                                                                              "p1_trailing"][i_blade],
                                                                          rotor_data_circulation["p2_trailing"][
                                                                              i_blade], d)
                # Bound vortices
                u_bound, v_bound, w_bound = get_induced_matrices(rotor_data_circulation["p1_bound"][i_blade],
                                                                 rotor_data_circulation["p2_bound"][i_blade], d,
                                                                 printen=False)
                u_bound = u_bound @ convergence_matrix
                v_bound = v_bound @ convergence_matrix
                w_bound = w_bound @ convergence_matrix

                sub_data_tab.append([u_wake + u_trailing + u_bound, v_wake + v_trailing + v_bound
                                        , w_wake + w_trailing + w_bound])
        data_tab.append(sub_data_tab)

    return np.array(data_tab)

def get_induced_velocity(UVW_matrix_data,trailing_circulation):
    n_rotors = len(UVW_matrix_data)
    n_blades = len(UVW_matrix_data[0,0])
    data_tab = []

    for i in range(n_rotors):
        sub_data_tab = []
        for j in range(n_blades):
            u,v,w = 0,0,0
            for k in range(n_rotors):
                for l in range(n_blades):
                    u += (UVW_matrix_data[i, k, j, l, 0] @ (trailing_circulation[k,l].T)).flatten()
                    v += (UVW_matrix_data[i, k, j, l, 1] @ (trailing_circulation[k,l].T)).flatten()
                    w += (UVW_matrix_data[i, k, j, l, 2] @ (trailing_circulation[k,l].T)).flatten()
            sub_data_tab.append([u,v,w])
        data_tab.append(sub_data_tab)
    return np.array(data_tab)

def get_cp(Ftan,r_R_list,R,v_inf,omega,r_R_list_control):
    Ftan = np.array(Ftan)
    n_rotors = len(Ftan)
    CP_list = np.zeros(n_rotors)
    dr = (r_R_list[1:] - r_R_list[:-1]) * R
    for i in range(n_rotors):
        for j in range(len(Ftan[0])):
            CP_list[i] += np.sum(dr * Ftan[i,j] * r_R_list_control * R * omega /(0.5 * v_inf **3 * np.pi * R**2))

    return CP_list

def get_ct(Faxial,r_R_list,R,v_inf):
    Faxial = np.array(Faxial)
    n_rotors = len(Faxial)
    CT_list = np.zeros(n_rotors)
    dr = (r_R_list[1:] - r_R_list[:-1]) * R
    for i in range(n_rotors):
        for j in range(len(Faxial[0])):
            CT_list[i] += np.sum(dr*Faxial[i,j]/(0.5 * v_inf**2 * np.pi * R**2))

    return CT_list

def plot_system(wake_data,rotor_data_tab,location_tab,R,one_blade=False):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    colors = ['b', 'orange', 'g', 'r', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    if len(wake_data)==1:
        label_list = [["Blade 1","Blade 2","Blade 3"]]
    else:
        label_list = [["Rotor R, Blade 1", "Rotor R, Blade 2", "Rotor R, Blade 3"],["Rotor L, Blade 1", "Rotor L, Blade 2", "Rotor L, Blade 3"]]

    for i in range(len(wake_data)):
        rotor_data = rotor_data_tab[i]
        wake_dict = wake_data[i]

    # plot wake
    wake_dict["position"] = np.array(wake_dict["position"])
    n_blades = len(wake_dict["position"][0])
    n_trailing = len(wake_dict["position"][0, 0])
    if one_blade == True:
        n_blades = 1
    for i_trailing in range(n_trailing):
        for i_blade in range(n_blades):
            x = wake_dict["position"][:, i_blade, i_trailing, 0]
            y = wake_dict["position"][:, i_blade, i_trailing, 1]
            z = wake_dict["position"][:, i_blade, i_trailing, 2]

            ax.plot(x, y, z, c='pink')

    # plot rotor
    rotor_data["p1_trailing"] = np.array(rotor_data["p1_trailing"])
    rotor_data["p2_trailing"] = np.array(rotor_data["p2_trailing"])

    for i_blade in range(n_blades):
        # quarterchord
        x = rotor_data["p1_trailing"][i_blade, :, 0]
        y = rotor_data["p1_trailing"][i_blade, :, 1]
        z = rotor_data["p1_trailing"][i_blade, :, 2]
        ax.plot(x, y, z, c=colors[i * n_blades + i_blade], label=label_list[i][i_blade])
        # end point
        x = rotor_data["p2_trailing"][i_blade, :, 0]
        y = rotor_data["p2_trailing"][i_blade, :, 1]
        z = rotor_data["p2_trailing"][i_blade, :, 2]
        ax.plot(x, y, z, c=colors[i * n_blades + i_blade])
    for i_blade in range(n_blades):
        # root chord
        x = [rotor_data["p1_trailing"][i_blade, 0, 0], rotor_data["p2_trailing"][i_blade,
        0, 0]]
        y = [rotor_data["p1_trailing"][i_blade, 0, 1], rotor_data["p2_trailing"][i_blade,
        0, 1]]
        z = [rotor_data["p1_trailing"][i_blade, 0, 2], rotor_data["p2_trailing"][i_blade,
        0, 2]]
        ax.plot(x, y, z, c=colors[i * n_blades + i_blade])
        # tip chord
        x = [rotor_data["p1_trailing"][i_blade, -1, 0], rotor_data["p2_trailing"][i_blade
        , -1, 0]]
        y = [rotor_data["p1_trailing"][i_blade, -1, 1], rotor_data["p2_trailing"][i_blade
        , -1, 1]]
        z = [rotor_data["p1_trailing"][i_blade, -1, 2], rotor_data["p2_trailing"][i_blade
        , -1, 2]]
        ax.plot(x, y, z, c=colors[i * n_blades + i_blade])

    # ax.set_xlim(-50, 50)
    y_min = np.min(location_tab[:, 1]) - R
    y_max = np.max(location_tab[:, 1]) + R
    d = abs(y_max - y_min)
    ax.set_ylim(y_min,
                y_max)  # Set the y-axis limits from 5 to 11 ax.set_zlim(-d/2, d/2) # Set the z-axis limits from 10 to 16
    ax.set_ylim(-75, 175)
    ax.set_zlim(-75, 75)
    plt.legend()
    plt.show()

def run_lifting_line_alpha(n_rotors,TSR,n_sections,time_steps,dt,location_tab = np.array([[0,0,0],[0,120,0]]),phase_tab = np.array([0,0]),linear_scale=False,v_wake_factor=1):
    result_dict = {}
    # Open data
    airfoil = 'Code Jus\polarDU95W180.csv'
    data1=pd.read_csv(airfoil, header=0, names = ["Alfa", "Cl", "Cd", "cm"])

    polar_alpha = np.array(data1['Alfa'][:])
    polar_cl = np.array(data1['Cl'][:])
    polar_cd = np.array(data1['Cd'][:])  # plt.plot(polar_alpha,polar_cl)
    # plt.plot(polar_alpha,2*np.pi*np.sin(np.array(polar_alpha * np.pi /180))) #plt.show()
    # Define parameters

    v_inf = 10
    r_R_start = 0.2
    n_blades = 3
    R = 50
    omega = TSR * v_inf / R

    if linear_scale == False:
        r_R_list = np.flip(r_R_start + (np.cos(np.linspace(0, np.pi, n_sections + 1)) * 0.5 + 0.5) * (1 - r_R_start))
        r_R_list_control = np.flip(r_R_start + (np.cos(0.5 * (np.linspace(0, np.pi, n_sections + 1)[1:] + np.linspace(0, np.pi, n_sections + 1)[:-1])) * 0.5 + 0.5) * (1 - r_R_start))
    else:
        r_R_list = np.linspace(r_R_start, 1, n_sections + 1)
        r_R_list_control = 0.5 * (r_R_list[1:] + r_R_list[:-1])

    result_dict["x_control"] = r_R_list_control * R
    result_dict["x"] = r_R_list * R

    chord_distribution = 3 * (1 - (r_R_list)) + 1
    pitch = -2  # degrees
    twist = 14 * (1 - r_R_list) + pitch  # degrees
    twist_mid = 14 * (1 - r_R_list_control) + pitch  # degrees7
    chord_distribution_mid = 3 * (1 - r_R_list_control) + 1
    a_w_tab = np.zeros((n_rotors, n_blades, n_sections))
    bound_to_trailing_convergence = bound_to_trailing_circulation(np.zeros(n_sections))

    wake_tab = []
    for i in range(n_rotors):
        wake_tab.append(init_wake())
    trailing_circulation = np.zeros((n_rotors, n_blades, n_sections + 1))

    rotor_data_pre = []
    for i in range(n_rotors):
        theta_global = (time_steps - 1) * dt * omega + phase_tab[i]
        rotor_data_pre.append(get_initial_rotor_data(location_tab[i], R, r_R_list, theta_global, n_blades, chord_distribution, twist, r_R_list_control, chord_distribution_mid, twist_mid))

    for time_step in range(time_steps):
        rotor_data = []
        wake_data_pre = []
        for i in range(n_rotors):
            theta_global = time_step * dt * omega + phase_tab[i]
            rotor_data.append(get_initial_rotor_data(location_tab[i], R, r_R_list,
                                                 theta_global, n_blades, chord_distribution, twist, r_R_list_control,
                                                 chord_distribution_mid, twist_mid))
            wake_data_pre.append(append_wake(rotor_data[-1], wake_tab[i]))

    converging = True
    while converging:
        wake_data = []
        for i in range(n_rotors):
            wake_data.append(get_wake2(rotor_data_pre[i], wake_data_pre[i], np.array(a_w_tab), r_R_list, v_inf, dt, v_wake_factor))

        UVW_matrix_data = get_induced_velocity_matrices(rotor_data_pre, wake_data)
        UVW_data = get_induced_velocity(UVW_matrix_data, trailing_circulation)
        last_trailing_circulation = trailing_circulation + 0
        trailing_circulation = []
        a_w_tab = []
        a_prime_tab = []
        alpha_tab = []
        inflow_angle_tab = []
        bound_circulation_tab = []
        F_azi_tab = []
        F_axial_tab = []
        CT_tab = []
        CP_tab = []
        C_azi_tab = []
        C_axial_tab = []
        cl_tab = []

        for i in range(n_rotors):
            theta_global = time_step * dt * omega + phase_tab[i]
            sub_trailing_circulation = []
            sub_a_w_tab = []
            sub_a_prime_tab = []
            sub_alpha_tab = []
            sub_inflow_angle_tab = []
            sub_bound_circulation_tab = []
            sub_F_azi_tab = []
            sub_F_axial_tab = []
            sub_CT_tab = []
            sub_CP_tab = []
            sub_C_azi_tab = []
            sub_C_axial_tab = []
            sub_cl_tab = []

            for j in range(n_blades):
                theta_local = theta_global + j * 2 * np.pi / n_blades
                u_induced = UVW_data[i, j, 0]
                v_induced = UVW_data[i, j, 1]
                w_induced = UVW_data[i, j, 2]
                w_rotated = - v_induced * np.sin(theta_local) + w_induced * np.cos(
                    theta_local)
                ############## plus/minus
                V_tan = r_R_list_control * R * omega - w_rotated
                # V_tan = 0.5 * (r) * R * omega - w_rotated
                V_axial = v_inf + u_induced
                # print(V_axial)
                # print(V_tan)
                V_p = np.sqrt(V_tan ** 2 + V_axial ** 2)
                alpha = np.arctan2(V_axial, V_tan) * 180 / np.pi - twist_mid
                phi_rad = np.arctan2(V_axial, V_tan)
                cl = np.interp(alpha, polar_alpha, polar_cl)
                cd = np.interp(alpha, polar_alpha, polar_cd)
                L = 0.5 * V_p ** 2 * cl * chord_distribution_mid
                D = 0.5 * V_p ** 2 * cd * chord_distribution_mid
                Fazim = L * np.sin(phi_rad) - D * np.cos(phi_rad)
                Faxial = L * np.cos(phi_rad) + D * np.sin(phi_rad)
                bound_circulation = 0.5 * chord_distribution_mid * V_p * cl
                ############## plus/minus
                trailing_circulation_new = 0.9 * last_trailing_circulation[
                    i, j] + 0.1 * bound_to_trailing_convergence @ bound_circulation.T
                sub_trailing_circulation.append(trailing_circulation_new)
                a_w = 1 - V_axial / v_inf
                a_prime = V_tan / (r_R_list_control * R * omega) - 1

                sub_a_prime_tab.append(a_prime)
                sub_a_w_tab.append(a_w)
                sub_alpha_tab.append(alpha)
                sub_inflow_angle_tab.append(phi_rad * 180 / np.pi)
                sub_bound_circulation_tab.append(bound_circulation)
                sub_F_azi_tab.append(Fazim)
                sub_F_axial_tab.append(Faxial)
                Cazim = cl * np.sin(phi_rad) - cd * np.cos(phi_rad)  # non dim
                Caxial = cl * np.cos(phi_rad) + cd * np.sin(phi_rad) # non dim
                sub_C_azi_tab.append(Cazim)
                sub_C_axial_tab.append(Caxial)
                sub_cl_tab.append(cl)

            trailing_circulation.append(sub_trailing_circulation)
            a_w_tab.append(sub_a_w_tab)
            alpha_tab.append(sub_alpha_tab)
            inflow_angle_tab.append(sub_inflow_angle_tab)
            bound_circulation_tab.append(sub_bound_circulation_tab)
            F_azi_tab.append(sub_F_azi_tab)
            F_axial_tab.append(sub_F_axial_tab)
            CT_tab.append(sub_CT_tab)
            CP_tab.append(sub_CP_tab)
            C_azi_tab.append(sub_C_azi_tab)
            C_axial_tab.append(sub_C_axial_tab)
            a_prime_tab.append(sub_a_prime_tab)
            cl_tab.append(sub_cl_tab)

        trailing_circulation = np.array(trailing_circulation)
        if np.max(abs(trailing_circulation - last_trailing_circulation)) < 0.0001:
            time_step += 1
            converging = False
            plot_system(wake_data, rotor_data_pre, location_tab, R, one_blade=False)
            plt.show()

            print(np.average(np.average(a_w_tab, weights=(0.5 * (r_R_list[1:] - r_R_list[:-1]))
                                        .flatten(), axis=2)))

        print(np.average(np.average(a_w_tab, weights=(0.5 * (r_R_list[1:] - r_R_list[:-1])).flatten(), axis=2)))

    result_dict["CT"] = get_ct(F_axial_tab, r_R_list, R, v_inf)
    result_dict["CP"] = get_cp(F_azi_tab, r_R_list, R, v_inf, omega, r_R_list_control)
    result_dict['r_R_list'] = result_dict['x_control'] / R
    result_dict["trailing_circulation"] = np.array(trailing_circulation) / (v_inf ** 2 * np.pi / (
                n_blades * omega))
    result_dict["a_w"] = a_w_tab
    result_dict["a_prime"] = a_prime_tab
    result_dict["alpha"] = alpha_tab
    result_dict["inflow_angle"] = inflow_angle_tab
    result_dict["bound_circulation"] = np.array(bound_circulation_tab) / (v_inf ** 2 * np.pi / (
                n_blades * omega))
    result_dict["F_azi"] = F_azi_tab
    result_dict["F_axial"] = F_axial_tab
    result_dict["location_tab"] = location_tab
    result_dict["R"] = R
    result_dict["wake_data"] = wake_data
    result_dict["rotor_data"] = rotor_data
    result_dict['nondimFaxial'] = C_axial_tab
    result_dict['nondimFazim'] = C_azi_tab
    result_dict["cl"] = cl_tab

    # plot_system(wake_data, rotor_data_pre, location_tab, R, one_blade=False)
    return result_dict

def plot_circulation(Circ_bem, Circ_LL,r_R_list_ll,r_R_list_bem, TSR, omega, v_inf):

    circulation_BEM = np.array(Circ_bem)
    circulation_LL = np.array(Circ_LL)
    plt.grid()
    plt.plot(r_R_list_bem, circulation_BEM / (np.pi * v_inf ** 2 / (n_blades * omega)), label='Circulationm BEM with TSR = ' + str(TSR), marker='x', linestyle='--')
    plt.plot(r_R_list_ll, circulation_LL / (np.pi * v_inf ** 2 / (n_blades * omega)), label='Circulation_LL with TSR = ' +str(TSR), marker='.')

    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Non dimensionalised circulation [-]')
    plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega * NBlades} $')
    #plt.savefig('Images / BEMcomparison / Circulation')
    return

def plot_aoa(alpha_BEM,alpha_LL,r_R_list_ll,r_R_list_bem, TSR):
    plt.grid()
    plt.plot(r_R_list_bem,alpha_BEM,label = 'Angle of Attack BEM with TSR = ' +str(TSR), marker='x', linestyle = '--')
    plt.plot(r_R_list_ll,alpha_LL,label = 'Angle of Attack LLT with TSR = ' +str(TSR),marker= '.')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel(r'$\alpha$ [deg]')
    plt.title("Angle of attack BEM VS LL")
    #plt.savefig('Images/BEMcomparison/aoa')
    return

def plot_phi(phi_BEM,phi_LL,r_R_list_ll,r_R_list_bem, TSR):
    plt.grid()
    plt.plot(r_R_list_bem,phi_BEM,label = 'Angle of Inflow BEM with TSR = ' +str(TSR),marker=
    'x', linestyle = '--')
    plt.plot(r_R_list_ll,phi_LL,label = 'Angle of Inlfow LLT with TSR = ' +str(TSR),marker='.')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel(r'$\phi$ [deg]')
    plt.title("Angle of Inflow BEM VS LL") #plt.savefig('Images/BEMcomparison/phi')
    return

def plot_loads(F_Bem,r_R_list_ll,r_R_list_bem,F_LLT, TSR):
    plt.grid()
    plt.plot(r_R_list_bem, np.array(F_Bem) / (0.5 * rho * v_inf ** 2 * R), label = 'Azimuthal load of BEM with TSR = ' +str(TSR),marker='x', linestyle = '--')
    plt.plot(r_R_list_ll, np.array(F_LLT) / (0.5 * rho * v_inf ** 2 * R), label = 'Azimuthal load of LL with TSR = ' +str(TSR),marker='.')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Normalized Load [-]')
    plt.title('Normalized Azimuthal Load BEM VS TSR') #plt.savefig('Images/BEMcomparison/Azimuthalload')

def plot_thrust(F_Bem,r_R_list_ll,r_R_list_bem,F_LLT, TSR):
    plt.grid()
    plt.plot(r_R_list_bem, np.array(F_Bem) / (0.5 * rho * v_inf ** 2 * R), label = 'Axial load of BEM with TSR = ' +str(TSR),marker='x', linestyle = '--')
    plt.plot(r_R_list_ll, np.array(F_LLT) / (0.5 * rho * v_inf ** 2 * R), label = 'Axial load of LL with TSR = ' +str(TSR),marker='.')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Normalized Load [-]')
    plt.title('Normalized Azxial Load BEM VS TSR')
    #plt.savefig('Images/BEMcomparison/Axialload')
    return

n_rotors = 1
n_blades = 3
n_sections = 15
time_steps = 10 * 15 + 1
v_inf = 10
r_R_start = 0.2
R = 50
blade_pitch = -2
rho = 1
location_tab = np.array([[0, 0, 0], [0, 100, 0]])
phase_tab = np.array([0,0])
TSR = 8
omega = TSR * v_inf / R
dt = 2 * np.pi / (omega) / 15
result_dict = run_lifting_line_alpha(n_rotors,TSR,n_sections,time_steps,dt,location_tab= location_tab,phase_tab=np.array([0,0]))