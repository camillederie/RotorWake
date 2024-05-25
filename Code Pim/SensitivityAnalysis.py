# Import necessary libraries for plotting
import matplotlib.pyplot as plt
import numpy as np
import Lifting_line as ll
import BEM as bem
from turbine2 import Turbine as Tb2
from turbine import Turbine as Tb1

# Initialize turbine instances
TB1 = Tb1(0.01)
TB2 = Tb2(10, 10)

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
    dict: A dictionary with CTrotor, CProtor.
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

# Function to generate different radial discretizations
def generate_radial_discretization(method, N):
    if method == 'constant':
        return np.linspace(0.2, 1, N)
    elif method == 'cosine':
        return (1 - np.cos(np.linspace(0, np.pi, N))) / 2 * (1 - 0.2) + 0.2
    else:
        raise ValueError("Unknown discretization method")

# Function to test different convection speeds
def test_convection_speeds(TB2, ll, speeds):
    results = []

    for speed in speeds:
        #DIT STAAT VOOR NU OP UINF MAAR MOET DENK EEN ANDERE INPUT ZIJN, DE CONVECTION SPEED, MAAR VOOR NU KLOPT HET
        U_inf = speed * 10

        # Call the function from Lifting_line.py
        gamma_ll, alpha_ll, f_ax_ll, f_ang_ll, phi_ll, a_ll, a_prime_ll = ll.Lifting_Line_Method(TB2, bem.TSR, U_inf, bem.omega, 1.225, tol=1, iterations=1200)

        # Normalize circulation and forces
        circulation_norm = (np.pi * U_inf**2) / (bem.omega * TB2.n)
        force_norm = 0.5 * bem.U0**2 * TB2.R

        # Store results
        results.append({
            'convection_speed': speed,
            'gamma_ll': gamma_ll / circulation_norm,
            'f_ax_ll': f_ax_ll / force_norm,
            'f_ang_ll': f_ang_ll / force_norm
        })

    return results

# Function to compare cosine and regular spacing
def compare_spacing(TB2, ll):
    results = []

    for spacing in ['constant', 'cosine']:
        TB2.r_cp = generate_radial_discretization(spacing, TB2.Nrad) * 50

        # Call the function from Lifting_line.py
        gamma_ll, alpha_ll, f_ax_ll, f_ang_ll, phi_ll, a_ll, a_prime_ll = ll.Lifting_Line_Method(TB2, bem.TSR, bem.U0, bem.omega, 1.225, tol=1, iterations=12000)

        # Normalize circulation and forces
        circulation_norm = (np.pi * bem.U0**2) / (bem.omega * TB2.n)
        force_norm = 0.5 * bem.U0**2 * TB2.R

        # Store results
        results.append({
            'spacing': spacing,
            'gamma_ll': gamma_ll / circulation_norm,
            'f_ax_ll': f_ax_ll / force_norm,
            'f_ang_ll': f_ang_ll / force_norm
        })

    return results

def test_wake_segments(TB2, ll, bem, segments_list):
    results = []
    # Calculate CP and CT for the case with the most segments
    TB2.wake_steps = np.arange(0, TB2.rotation * 2 * np.pi,  2 * np.pi / segments_list[-1])
    gamma_ll, alpha_ll, f_ax_ll, f_ang_ll, phi_ll, a_ll, a_prime_ll = ll.Lifting_Line_Method(TB2, bem.TSR, bem.U0, bem.omega,1.225, tol=1, iterations=1200)
    results_ll = calculate_CT_CProtor_CPflow(f_ax_ll, f_ang_ll, bem.U0, TB2.r_cp / TB2.R, bem.omega, TB2.R, TB2.n)
    CP_reference = results_ll['CProtor']
    CT_reference = results_ll['CTrotor']

    for segments in segments_list:
        TB2.wake_steps = np.arange(0, TB2.rotation * 2 * np.pi, 2* np.pi / segments)
        # Call the function from Lifting_line.py
        gamma_ll, alpha_ll, f_ax_ll, f_ang_ll, phi_ll, a_ll, a_prime_ll = ll.Lifting_Line_Method(TB2, bem.TSR, bem.U0, bem.omega,1.225, tol=1, iterations=1200)

        # Calculate CP and CT for Lifting Line
        results_ll = calculate_CT_CProtor_CPflow(f_ax_ll, f_ang_ll, bem.U0, TB2.r_cp / TB2.R, bem.omega, TB2.R, TB2.n)
        CP_ll = results_ll['CProtor']
        CT_ll = results_ll['CTrotor']
        delta_cp = abs(CP_reference - CP_ll)
        delta_ct = abs(CT_reference - CT_ll)

        # Store results
        results.append({
            'segments': segments,
            'delta_cp': delta_cp,
            'delta_ct': delta_ct
        })

    return results

def test_wake_length(TB2, ll, bem, lengths):
    results = []
    # Calculate CP and CT for the case with the longest wake length
    TB2.rotation = lengths[-1]
    gamma_ll, alpha_ll, f_ax_ll, f_ang_ll, phi_ll, a_ll, a_prime_ll = ll.Lifting_Line_Method(TB2, bem.TSR, bem.U0, bem.omega,1.225, tol=1, iterations=1200)
    results_ll = calculate_CT_CProtor_CPflow(f_ax_ll, f_ang_ll, bem.U0, TB2.r_cp / TB2.R, bem.omega, TB2.R, TB2.n)
    CP_reference = results_ll['CProtor']
    CT_reference = results_ll['CTrotor']

    for length in lengths:
        TB2.rotation = length
        # Call the function from Lifting_line.py
        gamma_ll, alpha_ll, f_ax_ll, f_ang_ll, phi_ll, a_ll, a_prime_ll = ll.Lifting_Line_Method(TB2, bem.TSR, bem.U0, bem.omega, 1.225, tol=1, iterations=1200)

        # Calculate CP and CT for Lifting Line
        results_ll = calculate_CT_CProtor_CPflow(f_ax_ll, f_ang_ll, bem.U0, TB2.r_cp / TB2.R, bem.omega, TB2.R, TB2.n)
        CP_ll = results_ll['CProtor']
        CT_ll = results_ll['CTrotor']
        delta_cp = abs(CP_reference - CP_ll)
        delta_ct = abs(CT_reference - CT_ll)

        # Store results
        results.append({
            'length': length,
            'delta_cp': delta_cp,
            'delta_ct': delta_ct
        })

    return results


# Plotting functions for each analysis
def plot_convection_speeds(results):
    plt.figure()
    for res in results:
        plt.plot(TB2.r_cp / TB2.R, res['gamma_ll'][0:TB2.Nrad], label=f"Speed: {res['convection_speed']}")
    plt.title("Circulation vs r/R for Different Convection Speeds")
    plt.xlabel('r/R')
    plt.ylabel('Normalized Circulation')
    plt.grid()
    plt.legend()
    plt.show()

    plt.figure()
    for res in results:
        plt.plot(TB2.r_cp / TB2.R, res['f_ax_ll'][0:TB2.Nrad], label=f"Speed: {res['convection_speed']} (Axial)")
        plt.plot(TB2.r_cp / TB2.R, res['f_ang_ll'][0:TB2.Nrad], '--', label=f"Speed: {res['convection_speed']} (Tangential)")
    plt.title("Forces vs r/R for Different Convection Speeds")
    plt.xlabel('r/R')
    plt.ylabel('Normalized Force')
    plt.grid()
    plt.legend()
    plt.show()

def plot_spacing(results):
    plt.figure()
    for res in results:
        plt.plot(TB2.r_cp / TB2.R, res['gamma_ll'][0:TB2.Nrad], label=f"Spacing: {res['spacing']}")
    plt.title("Circulation vs r/R for Different Blade Spacing Methods")
    plt.xlabel('r/R')
    plt.ylabel('Normalized Circulation')
    plt.grid()
    plt.legend()
    plt.show()

    plt.figure()
    for res in results:
        plt.plot(TB2.r_cp / TB2.R, res['f_ax_ll'][0:TB2.Nrad], label=f"Spacing: {res['spacing']} (Axial)")
        plt.plot(TB2.r_cp / TB2.R, res['f_ang_ll'][0:TB2.Nrad], '--', label=f"Spacing: {res['spacing']} (Tangential)")
    plt.title("Forces vs r/R for Different Blade Spacing Methods")
    plt.xlabel('r/R')
    plt.ylabel('Normalized Force')
    plt.grid()
    plt.legend()
    plt.show()


def plot_wake_segments(results):
    plt.figure()
    segments = [res['segments'] for res in results]
    delta_cp = [res['delta_cp'] for res in results]
    delta_ct = [res['delta_ct'] for res in results]

    plt.subplot(2, 1, 1)
    plt.plot(segments, delta_cp, 'o-')
    plt.title("$\Delta C_P$ and $\Delta C_T$ vs Wake Segments compared to converged value")
    plt.ylabel('$\Delta C_P$')
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(segments, delta_ct, 'o-')
    plt.xlabel('Wake Segments')
    plt.ylabel('$\Delta C_T$')
    plt.grid()

    plt.tight_layout()
    plt.show()


def plot_wake_length(results):
    plt.figure()
    lengths = [res['length'] for res in results]
    delta_cp = [res['delta_cp'] for res in results]
    delta_ct = [res['delta_ct'] for res in results]

    plt.subplot(2, 1, 1)
    plt.plot(lengths, delta_cp, 'o-')
    plt.title("$\Delta C_P$ and $\Delta C_T$ vs Wake Length compared to converged value")
    plt.ylabel('$\Delta C_P$')
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(lengths, delta_ct, 'o-')
    plt.xlabel('Wake Length (rotations)')
    plt.ylabel('$\Delta C_T$')
    plt.grid()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Define parameters for analysis
    convection_speeds = [0.8, 1.0, 1.2]
    wake_segments_list = np.arange(2, 21, 1)
    wake_lengths = np.arange(1, 20, 2)
    bem.TSR = 8
    # Test different convection speeds
    #convection_results = test_convection_speeds(TB2, ll, convection_speeds)
    #plot_convection_speeds(convection_results)

    # Compare cosine and regular spacing
    #spacing_results = compare_spacing(TB2, ll)
    #plot_spacing(spacing_results)

    # Test different wake segments per rotation
    #wake_segments_results = test_wake_segments(TB2, ll, bem, wake_segments_list)
    #plot_wake_segments(wake_segments_results)

    # Test different wake lengths (number of rotations)
    wake_length_results = test_wake_length(TB2, ll, bem, wake_lengths)
    plot_wake_length(wake_length_results)