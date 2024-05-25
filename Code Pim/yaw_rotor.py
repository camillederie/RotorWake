import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends

matplotlib.use('TkAgg')

# file import
from corrections import prandtl, glauert, glauert2
from turbine import Turbine
from experimental_data import Data

cases = ['turbine', 'propeller']
dr = 0.01


# Tb = Turbine(dr)
# Pp = Propeller(dr)


class YAW(Turbine):
    """"
    This class...
    """

    def __init__(self, dr, case):
        super().__init__(dr)
        self.a_i = 0.3  # initial induction ratio guess
        self.a_p_i = 0.0  # initial prime induction ratio guess
        self.dr = dr  # space step in the blade

        self.U0 = 10  # [m/s]
        self.rho = 1.225  # [kg/m^3] Using ISA

        # (i, j) for selecting between ([6, 8, 10] , [0, 15, 30])
        self.TSR, self.theta = self.operations(1, 1)
        self.omega = lambda radius: (self.TSR * self.U0) / radius

        self.C_L = lambda alpha: Data(case).lift(alpha)
        self.C_D = lambda alpha: Data(case).drag(alpha)

    def velocity(self, a, a_p, radius):
        """"Determine the velocity components:"""
        v_axial = self.U0 * (1 - a)
        v_tan = self.omega(radius) * radius * (1 + a_p)
        v_p = np.sqrt(v_axial ** 2 + v_tan ** 2)
        return v_axial, v_tan, v_p

    def angles(self, twist, pitch, a, a_p, radius):
        """Determine the associated angles:"""
        v_axial, v_tan, v_p = self.velocity(a, a_p, radius)
        phi = np.arctan2(v_axial, v_tan)  # incomming flow angle
        beta = twist + pitch  # static orientation + dynamic orientation
        alpha = phi - beta  #
        return phi, beta, alpha

    def lift(self, alpha, chord, velocity):
        """Determine the sectional 2D lift force of the blade section"""
        return self.C_L(alpha) * 0.5 * chord * velocity ** 2

    def drag(self, alpha, chord, velocity):
        """Determine the sectional 2D lift force of the blade section"""
        return self.C_D(alpha) * 0.5 * chord * velocity ** 2

    def forces(self, phi, alpha, chord, velocity):
        """Azimuthal and axial 2D force can be determined."""
        angular = self.lift(alpha, chord, velocity) * np.sin(phi) - \
                  self.drag(alpha, chord, velocity) * np.cos(phi)
        axial = self.lift(alpha, chord, velocity) * np.cos(phi) + \
                self.drag(alpha, chord, velocity) * np.sin(phi)
        gamma = 0.5 * velocity * self.C_L(alpha) * chord
        return axial, angular, gamma

    def ct_cp(self, a, glauert_corr=True):
        """"Calculate thrust coefficient for given induction."""
        c_t = 4 * a * (np.cos(self.theta) - a)
        c_p = 4 * a * (np.cos(self.theta) - a)**2
        if glauert_corr:
            c_t = glauert(c_t, a)
        return c_t, c_p

    def u_glauert(self):
        """"
        uniform induced velocity(think about upwash in front of a wing)
        This is only true for 0 and 90 degrees, but we assume it to be true for all angles yaw.
        Note: at high angle yaw: F and u are aligned with e_z (rotor axis)
        """
        F =
        u_rel = np.sqrt((self.U0 * np.cos(self.theta) - u)**2 + (self.U0 * np.sin(self.theta))**2)
        u = 2 * F / (np.pi * (2*self.r)**2 * self.rho * u_rel) # upwash velocity of the rotor plane

        return u

    def blade_element(self, a_0, a_p_0, tol, it=0):
        """"Calculate the performance of a single angular disk."""
        phi, _, alpha = self.angles(self.twist, self.pitch, a_0, a_p_0, self.r)
        _, _, velocity = self.velocity(a_0, a_p_0, self.r)
        f_ax, f_ang, _ = self.forces(phi, alpha, self.chord, velocity)

        c_t = (f_ax * self.n * self.dr) / \
              (0.5 * self.U0 ** 2 * self.angular_area(self.r))

        a_1 = glauert2(c_t)
        prandtls, _, _ = prandtl(self.r_R, self.n, self.TSR, self.blade_start, a_1)
        prandtls[np.isnan(prandtls)] = 0.0001

        a_1 = a_1 / prandtls

        a_p_1 = (f_ang * self.n) / (2 * np.pi * self.U0 * (1 - a_1) *
                                    self.omega(self.r) * 2 * (self.r_R * self.r) ** 2)

        a_p_1 = a_p_1 / prandtls

        if np.any(np.abs(a_1 - a_0) > tol):
            a_1 = 0.25 * a_1 + 0.75 * a_0
            it += 1
            print(it)
            return self.angular_disk(a_1, a_p_1, tol, it)
        else:
            print(a_1, a_p_1, c_t, it)
            return a_1, a_p_1, c_t, it

    def plot_aap(self, tol):
        a, a_p, c_t, _ = self.angular_disk(self.a_i, self.a_p_i, tol)
        plt.plot(self.r_R, a, 'b-', label='a')
        plt.plot(self.r_R, a_p, 'r--', label='a_p')

        plt.ylabel(r'$C_T$ or $C_P$')
        plt.xlabel('a')
        plt.ylim(0, 1)
        plt.grid()
        plt.legend()
        plt.show()
        return 0

    def plot_ct(self, tol):
        a, a_p, c_t, _ = self.angular_disk(self.a_i, self.a_p_i, tol)
        plt.plot(self.r_R, c_t, 'b-', label='$C_T$')

        plt.ylabel(r'$C_T$ or $C_P$')
        plt.xlabel('a')
        plt.ylim(0, 2)
        plt.grid()
        plt.legend()
        plt.show()
        return 0


if __name__ == "__main__":
    dr = 0.1
    tol = 0.05
    solution = BEM(dr, 'turbine')
    solution.blade_element(0.3, 0.0, tol)
    # solution.plot_aap(tol)
    # solution.plot_ct(tol)

