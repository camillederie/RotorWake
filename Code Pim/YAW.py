import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends
matplotlib.use('TkAgg')

# file import
from corrections import prandtl, glauert2
from turbine import Turbine
from experimental_data import Data


""""
Case Specifics
"""
dr = 0.01
Tb = Turbine(dr)

a_i = 0.3      # initial induction ratio guess
a_p_i = 0.4    # initial prime induction ratio guess
dr = dr        # space step in the blade

U0 = 10        # [m/s]
rho = 1.225    # [kg/m^3] Using ISA

# (i, j) for selecting between ([6, 8, 10, 7] , [0, 15, 30])
TSR, yaw = Tb.operations(1, 0)
omega = U0 * TSR / Tb.R

""""
Experimental data
"""
C_L = lambda alpha: Data(Tb.type).lift(np.degrees(alpha))
C_D = lambda alpha: Data(Tb.type).drag(np.degrees(alpha))

azimutal_angle = np.linspace(0,2*np.pi,360)

def lift(alpha, chord, velocity):
    """Determine the sectional 2D lift force of the blade section"""
    return C_L(alpha) * 0.5 * chord * velocity**2

def drag(alpha, chord, velocity):
    """Determine the sectional 2D lift force of the blade section"""
    return C_D(alpha) * 0.5 * chord * velocity**2


""""
Operational calculations
"""
def velocity(U0, a, a_p, radius):
    """"Determine the velocity components:"""
    v_axial = U0 * (1 - a)
    v_tan = omega * radius * (1 + a_p)
    v_p = np.sqrt(v_axial**2 + v_tan**2)
    return v_axial, v_tan, v_p

def angles(twist, pitch, a, a_p, radius):
    """Determine the associated angles:"""
    v_axial, v_tan, v_p = velocity(U0, a, a_p, radius)
    phi = np.arctan2(v_axial, v_tan)    # incoming flow angle             # static orientation + dynamic orientation
    alpha = phi + twist +pitch        # Flow angle - geometric angle
    # Ensure alpha does not exceed 30.0 degrees in radians or -16
    #That would be out of the data range
    for i in range(len(alpha)):
        if alpha[i]>np.radians(30):
            alpha[i]=np.radians(30)
        if alpha[i]<np.radians(-16):
            alpha[i]=np.radians(-16)
    return phi, pitch, alpha

def forces(phi, alpha, chord, velocity):
    """Azimuthal and axial 2D force can be determined."""
    angular = lift(alpha, chord, velocity) * np.sin(phi) - \
              drag(alpha, chord, velocity) * np.cos(phi)
    axial = lift(alpha, chord, velocity) * np.cos(phi) + \
            drag(alpha, chord, velocity) * np.sin(phi)
    gamma = 0.5 * velocity * C_L(alpha) * chord
    return axial, angular, gamma


def yaw_modulation(Tb, a, U0, K, xsi, r_R, psi, theta):
    
    #u and f_axial are assumed to be perpendicular on the rotor plane.
    #This means the angle between U0 and u is equal to theta.
    
    u = glauert4(a, U0, K, xsi, r_R, psi)
    u_rel = np.sqrt((U0 * np.cos(theta) - u)**2 + (U0 * np.sin(theta))**2)
    f_axial = np.pi * Tb.R**2 * 2 * u * u_rel
    power = f_axial * (np.cos(theta) - a)
    return u, u_rel, f_axial

def vortex_cylinder(U0, a, theta):
    
    #This is the vortex cylinder model for yawed condition by Coleman.
    #An approximation for xsi is derived from the Biot-Savart law
    #This is done to require an averaged induced velocity.
    #Far in the wake u_ind needs to be multiplied with 2.
    
    xsi = (0.6 * a + 1) * theta  # approximation to the real relation given in the slides.
    u_ind = a * U0 / np.cos(xsi/2)  # sec = 1/cos
    c_t = 4 * a * (np.cos(theta) + np.sin(theta)*np.tan(xsi/2) -
                   a / np.cos(xsi/2)**2)
    c_p = c_t * (np.cos(theta) - a)
    return xsi, u_ind, c_t

def blade_element(Tb, U0, a_0, a_p_0, tol, it=0):
    """"Calculate the performance of a single angular disk."""
    phi, _, alpha = angles(Tb.twist, Tb.pitch, a_0, a_p_0, Tb.r)
    _, _, u = velocity(U0, a_0, a_p_0, Tb.r)
    f_ax, f_ang, gamma = forces(phi, alpha, Tb.chord, u)

    c_t = (f_ax * Tb.n * dr) / \
          (0.5 * U0**2 * Tb.angular_area(Tb.r))


    a_1 = glauert2(c_t, yaw)

    prandtls, _, _ = prandtl(Tb.r_R, Tb.n, TSR, Tb.blade_start, a_1)
    # prandtls[np.isnan(prandtls)] = 0.0001

    a_1 = a_1 / prandtls
    a_p_1 = (f_ang * Tb.n) / (2 * np.pi * U0 * (1 - a_1) * omega * 2 * Tb.r**2)
    a_p_1 = a_p_1 / prandtls

    if np.any(np.abs(a_1 - a_0) > tol):
        a_1 = 0.25 * a_1 + 0.75 * a_0
        a_p_1 = 0.25 * a_p_1 + 0.75 * a_p_0
        it += 1
        # print(a_1, a_p_1, gamma, it)
        return blade_element(Tb, U0, a_1, a_p_1, tol, it)
    else:
        return f_ax, f_ang, a_1, a_p_1, gamma, it, phi, alpha


""""
Plot data
"""
def plot_aap(Tb, tol):
    f_ax, f_ang, a_1, a_p_1, gamma, it, phi, alpha = blade_element(Tb, U0, a_i, a_p_i, tol)
    plt.plot(Tb.r_R, a, 'b-', label='a')
    plt.plot(Tb.r_R, a_p, 'r--', label='a_p')
    plt.title("Axial and tangential induction")
    plt.ylabel(r'Induction factor')
    plt.xlabel('r/R')
    plt.ylim(0, 0.6)
    plt.grid()
    plt.legend()
    plt.show()

def plot_f_ax(Tb, tol):
    f_ax, f_ang, a_1, a_p_1, gamma, it, phi, alpha = blade_element(Tb, U0, a_i, a_p_i, tol)
    plt.plot(Tb.r_R, f_ax/(0.5*U0**2*Tb.R), 'b-', label='$C_x$')
    plt.plot(Tb.r_R, f_ang/(0.5*U0**2*Tb.R), 'g-', label='$C_y$')
    plt.title("Normal and tangential force, normalized by $0.5*U0**2*Tb.R$")
    plt.ylabel(r'Force')
    plt.xlabel('r/R')
    plt.ylim(0, 1.35)
    plt.grid()
    plt.legend()
    plt.show()

def plot_circulation(Tb, tol, U0, omega):
    f_ax, f_ang, a_1, a_p_1, gamma, it, phi, alpha = blade_element(Tb, U0, a_i, a_p_i, tol)
    plt.plot(Tb.r_R, gamma/((np.pi*U0**2)/(omega*Tb.n)))
    plt.title("Circulation normalized by $(np.pi*U0**2)/(omega*Tb.n)$")
    plt.ylabel("$\Gamma$")
    plt.xlabel("r_R")
    #plt.ylim(0.4,0.72)
    plt.grid()
    plt.show()

def plot_angles(Tb, tol):
    f_ax, f_ang, a_1, a_p_1, gamma, it, phi, alpha = blade_element(Tb, U0, a_i, a_p_i, tol)
    plt.plot(Tb.r_R, phi*180/np.pi, label = 'inflow angle')
    plt.plot(Tb.r_R, alpha*180/np.pi, label = 'angle of attack')
    plt.ylabel('Angle (degrees)')
    plt.xlabel('r_R')

    plt.grid()
    plt.legend()
    plt.show()



tol = 0.0001
f_ax, f_ang, a, a_p, gamma, it, phi, alpha = blade_element(Tb, U0, a_i, a_p_i, tol)
plot_angles(Tb,tol)
plot_aap(Tb, tol)
plot_f_ax(Tb, tol)
plot_circulation(Tb, tol, U0, omega)
#print(f_ax, a, a_p, it)

