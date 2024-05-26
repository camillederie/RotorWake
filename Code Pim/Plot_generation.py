# Import necessary libraries for plotting
import BEM as bem
import Lifting_line as ll
import matplotlib.pyplot as plt
import numpy as np
from turbine import Turbine as Tb1
from turbine2 import Turbine as Tb2

# Initialize turbine instances
TB1 = Tb1(0.01)
TB2 = Tb2(10, 20)

# Function to calculate CT and CP
def calculate_CT_CProtor_CPflow(Fnorm, Ftan, Uinf, r_Rarray, Omega, Radius, NBlades):
    """
    Calculates the performance of the rotor, returning CT and CPs.

    Parameters:
    Fnorm (numpy array): Normal force.
    Ftan (numpy array): Tangential force.
    Uinf (float): Free stream velocity.
    r_Rarray (numpy array): Radial positions normalized by the radius.
    Omega (float): Rotational speed.
    Radius (float): Rotor radius.
    NBlades (int): Number of blades.

    Returns:
    dict: A dictionary with CTrotor, CProtor, and CPflow.
    """
    CTrotor = 0  # thrust coefficient
    CProtor = 0

    for i in range(len(r_Rarray) - 1):
        r_R_temp = (r_Rarray[i] + r_Rarray[i + 1]) / 2
        drtemp = (-r_Rarray[i] + r_Rarray[i + 1])

        CTrotor += (drtemp * Fnorm[i] * NBlades) / (0.5 * Uinf * Uinf * np.pi * Radius)
        CProtor += (drtemp * Ftan[i] * r_R_temp * Omega * NBlades) / (0.5 * Uinf * Uinf * Uinf * np.pi)

    results = {'CTrotor': CTrotor, 'CProtor': CProtor}
    return results

def generate_plots():
    # Call the function from BEM.py
    # Define TSR values
    TSR_values = [6, 8, 10]
    colors = ['b', 'g', 'r']  # Define colors for each TSR

    # Initialize lists to store data
    gamma_bem_list = []
    gamma_ll_list = []
    phi_bem_list = []
    phi_ll_list = []
    alpha_bem_list = []
    alpha_ll_list = []
    f_ax_bem_list = []
    f_ax_ll_list = []
    f_ang_bem_list = []
    f_ang_ll_list = []
    CT_CP_list = []  # List to store CT and CP results

    # Loop through TSR values to calculate data
    for TSR in TSR_values:
        bem.omega = bem.U0 * TSR / 50
        f_ax_bem, f_ang_bem, a_1_bem, a_p_1_bem, gamma_bem, it, phi_bem, alpha_bem = bem.blade_element(TB1, bem.U0,
                                                                                                       bem.a_i,
                                                                                                       bem.a_p_i,
                                                                                                       bem.tol,
                                                                                                       it=0)
        gamma_bem_list.append(gamma_bem)
        alpha_bem_list.append(alpha_bem)
        f_ax_bem_list.append(f_ax_bem)
        f_ang_bem_list.append(f_ang_bem)
        phi_bem_list.append(phi_bem)
        gamma_ll, alpha_ll, f_ax_ll, f_ang_ll, phi_ll, a_ll, a_prime_ll = ll.Lifting_Line_Method(TB2, TSR, bem.U0,
                                                                                                 bem.omega, 1.225,
                                                                                                 tol=0.5,
                                                                                                 iterations=12000)
        gamma_ll_list.append(gamma_ll)
        alpha_ll_list.append(alpha_ll)
        f_ax_ll_list.append(f_ax_ll)
        f_ang_ll_list.append(f_ang_ll)
        phi_ll_list.append(phi_ll)


        cpct = calculate_CT_CProtor_CPflow(f_ax_ll, f_ang_ll, bem.U0, TB2.r_cp/TB2.R, bem.omega, TB2.R, TB2.n)
        print('The results for tip speed ratio', TSR, cpct)

    # Plot for Circulation
    plt.figure(figsize=(8, 6))
    for i, TSR in enumerate(TSR_values):
        plt.plot(TB1.r_R, gamma_bem_list[i] / ((np.pi * bem.U0 ** 2) / (bem.omega * TB1.n)),
                 label=f'BEM: TSR={TSR}', color=colors[i])
        plt.plot(TB2.r_cp / TB2.R, gamma_ll_list[i][0:TB2.Nrad] / ((np.pi * bem.U0 ** 2) / (bem.omega * TB1.n)),
                 label=f'LL: TSR={TSR}', color=colors[i], linestyle='--')

    plt.title("Circulation")
    plt.ylabel('Circulation')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    # plt.show()

    # Plot for Angles
    plt.figure(figsize=(8, 6))
    for i, TSR in enumerate(TSR_values):
        plt.plot(TB1.r_R, np.degrees(phi_bem_list[i]), label=f'BEM: TSR={TSR}', color=colors[i])
        plt.plot(TB2.r_cp / TB2.R, np.degrees(phi_ll_list[i][0:TB2.Nrad]), label=f'LL: TSR={TSR}', color=colors[i],
                 linestyle='--')

    plt.title("Inflow Angle")
    plt.ylabel('Angle, degrees')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    # plt.show()

    # Plot for Angles of Attack
    plt.figure(figsize=(8, 6))
    for i, TSR in enumerate(TSR_values):
        plt.plot(TB1.r_R, np.degrees(alpha_bem_list[i]), label=f'BEM: TSR={TSR}', color=colors[i])
        plt.plot(TB2.r_cp / TB2.R, np.degrees(alpha_ll_list[i][0:TB2.Nrad]), label=f'LL: TSR={TSR}', color=colors[i],
                 linestyle='--')

    plt.title("Angle of Attack")
    plt.ylabel('Angle, degrees')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    # plt.show()

    # Plot for Forces
    plt.figure(figsize=(8, 6))
    for i, TSR in enumerate(TSR_values):
        plt.plot(TB1.r_R, f_ax_bem_list[i] / (0.5 * bem.U0 ** 2 * TB1.R), label=f'BEM: TSR={TSR}', color=colors[i])
        plt.plot(TB2.r_cp / TB2.R, f_ax_ll_list[i][0:TB2.Nrad] / (0.5 * bem.U0 ** 2 * TB1.R), label=f'LL: TSR={TSR}',
                 color=colors[i], linestyle='--')

    plt.title("Axial Forces")
    plt.ylabel('Force')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    # plt.show()

    # Plot for Tangential Forces
    plt.figure(figsize=(8, 6))
    for i, TSR in enumerate(TSR_values):
        plt.plot(TB1.r_R, f_ang_bem_list[i] / (0.5 * bem.U0 ** 2 * TB2.R), label=f'BEM: TSR={TSR}', color=colors[i])
        plt.plot(TB2.r_cp / TB2.R, f_ang_ll_list[i][0:TB2.Nrad] / (0.5 * bem.U0 ** 2 * TB2.R),
                 label=f'LL: TSR={TSR}', color=colors[i], linestyle='--')

    plt.title("Tangential Forces")
    plt.ylabel('Force')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    # plt.show()

    # Print CT and CP results
    print("CT and CP Results:")
    for i, TSR in enumerate(TSR_values):
        print(f"TSR {TSR}: CT = {CT_CP_list[i]['CTrotor']}, CP = {CT_CP_list[i]['CProtor']}")

if __name__ == "__main__":
    # Generate plots and calculate CT and CP
    generate_plots()