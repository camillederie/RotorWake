"""
Runs frozen vortex wake solver // lifting line theory
"""

import numpy as np
from create_geometry import BladeGeometry
from scipy.interpolate import interp1d


class LiftingLineSolver:

    def __init__(self, geo, r_rotor, weight=0.3, tol=1e-6, n_iter=1000):
        """
        :param geo: BladeGeometry class
        :param r_rotor: rotor radius [m]
        :param weight: weighting factor for next step
        :param tol: stopping criteria
        :param n_iter: number of iterations
        """

        # rotor discretization
        self.geo = geo
        # rotor properties
        self.u_rot = np.array([self.geo.tsr / r_rotor * geo.v_inf, 0, 0])
        self.u_inf = geo.v_inf
        self.r_rotor = r_rotor
        self.double_rotor = self.geo.double_rotor

        # solver settings
        self.weight = weight
        self.tol = tol
        self.n_iter = n_iter

        # airfoil variables
        self._polarAirfoil()  # computes minimum/maximum alpha, polynomials for CL and CD alpha approximations

    def _compute_circ(self, gamma, weight):
        self.geo.filaments[-1] = self.geo.filaments[-1] * (1 - weight) + (weight * gamma)

    def _velocity_3D_from_vortex_filament(self, cp_i, core):
        """
        Computes the induced velocity at a control point due to ALL rings (vectorized part)
        :param cp_i: single control point
        :param core: minimum value for circulation for stability
        :return: matrix of size (3, (n_blades x (n_span-1))) containing the induced velocities on the control point
        in u, v, w.
        """

        r_gamma = self.geo.filaments[-1]  # vortex strength of each filament
        x1 = self.geo.filaments[0]
        y1 = self.geo.filaments[1]
        z1 = self.geo.filaments[2]
        x2 = self.geo.filaments[3]
        y2 = self.geo.filaments[4]
        z2 = self.geo.filaments[5]
        xc, yc, zc = cp_i

        R1 = np.sqrt((xc - x1) ** 2 + (yc - y1) ** 2 + (zc - z1) ** 2)
        R2 = np.sqrt((xc - x2) ** 2 + (yc - y2) ** 2 + (zc - z2) ** 2)

        R12_xx = ((yc - y1) * (zc - z2)) - ((zc - z1) * (yc - y2))
        R12_xy = -((xc - x1) * (zc - z2)) + ((zc - z1) * (xc - x2))
        R12_xz = ((xc - x1) * (yc - y2)) - ((yc - y1) * (xc - x2))

        R12_sq = (R12_xx ** 2) + (R12_xy ** 2) + (R12_xz ** 2)

        R01 = ((x2 - x1) * (xc - x1)) + ((y2 - y1) * (yc - y1)) + ((z2 - z1) * (zc - z1))
        R02 = ((x2 - x1) * (xc - x2)) + ((y2 - y1) * (yc - y2)) + ((z2 - z1) * (zc - z2))

        # check if target point is in the vortex filament core,
        # and modify to solid body rotation
        R12_sq = np.where(R12_sq < core ** 2, core ** 2, R12_sq)
        R1 = np.where(R1 < core, core, R1)
        R2 = np.where(R1 < core, core, R2)

        K = (r_gamma / (4 * np.pi * R12_sq)) * ((R01 / R1) - (R02 / R2))

        U = np.sum(K * R12_xx, axis=1)
        V = np.sum(K * R12_xy, axis=1)
        W = np.sum(K * R12_xz, axis=1)
        print(x1)
        return np.array([U, V, W])

    def _compute_induced_velocity(self):
        # create u, v, w matrix to store induced velocities. The matrix consists of elements in which the rows
        # represent the control points and columns the net effect of a single "ring". 
        uvw_mat = np.zeros((3, self.geo.cp.shape[0], self.geo.filaments.shape[1]))  # three square matrices
        core = 0.00001  # no idea what this is
        for i, cp_i in enumerate(self.geo.cp[:, :3]):  # only loop over the coordinates of the control points
            temp_v = self._velocity_3D_from_vortex_filament(cp_i, core)
            uvw_mat[:, i, :] = temp_v

        return uvw_mat

    @staticmethod
    def _geo_blade(r_R):
        pitch = 2  # in deg
        chord = 3 * (1 - r_R) + 1
        twist = -14 * (1 - r_R)  # in deg
        result = [chord.flatten(), np.radians(twist + pitch).flatten()]

        return result

    def _polarAirfoil(self):
        data = np.loadtxt("polar_DU95W180.csv", delimiter=';')
        data = data[3:]
        data = np.array(data, dtype=float)

        alphaRad = np.radians(data[:, 0])
        self.amax = max(alphaRad)
        self.amin = min(alphaRad)
        self.fcl = interp1d(alphaRad, data[:, 1], fill_value=(data[0, 1], data[-1, 1]), bounds_error=False,
                            kind='cubic')
        self.fcd = interp1d(alphaRad, data[:, 2], fill_value=(data[0, 2], data[-1, 2]), bounds_error=False,
                            kind='cubic')

    def _compute_loads_blade(self, v_norm, v_tan, r_R):
        V_mag2 = (v_norm ** 2 + v_tan ** 2)  # Velocity magnitude squared
        phi = np.arctan(v_norm / v_tan)  # Inflow angle

        # Get chord and twist
        [chord, twist] = self._geo_blade(r_R)

        alpha = twist + phi

        cl = self.fcl(alpha)
        cd = self.fcd(alpha)

        L = 0.5 * V_mag2 * cl * chord 
        D = 0.5 * V_mag2 * cd * chord

        F_norm = L * np.cos(phi) + D * np.sin(phi)
        F_tan = L * np.sin(phi) - D * np.cos(phi)

        Gamma = 0.5 * np.sqrt(V_mag2) * cl * chord

        return np.array([F_norm, F_tan, Gamma, alpha, phi])

    def _initialize_solver(self):

        # update Gamma given Gamma matrix, weight, and new Gamma
        self._compute_circ(gamma=1, weight=1)  # updates self.geo itself

        # compute [ui, vi, wi] based on vortex strength and distance
        # between control point and vortex
        v_induced = self._compute_induced_velocity()

        return v_induced

    def run_solver(self):
        # determine radial position of control point
        pos_radial = np.sqrt(np.sum(self.geo.cp[:(self.geo.n_blades * (self.geo.n_span - 1)), :3] ** 2,
                                    axis=1)).reshape(-1, 1)
        r_R = pos_radial / self.r_rotor

        if self.double_rotor:  # copy r/R for second rotor
            pos_radial = np.tile(pos_radial, 2).reshape((-1, 1), order='F')
            cp = self.geo.cp[:, :3].copy()
            cp[int(len(cp) / 2):, :3] = self.geo._compute_cp(self.geo.phase_diff)[:, :3]

        # initialize gamma vectors new and old
        gamma_new = np.ones((len(self.geo.cp), 1))
        # initialize output variables
        a = np.ones(len(self.geo.cp)) * 0.33
        aline = np.ones(len(self.geo.cp))
        # r_R = np.ones(len(self.geo.cp))
        f_norm = np.ones(len(self.geo.cp))
        f_tan = np.ones(len(self.geo.cp))
        gamma = np.ones(len(self.geo.cp))
        alpha = np.ones(len(self.geo.cp))
        phi = np.ones(len(self.geo.cp))
        # initial error
        err = 1.0
        error_log = []

        for i in range(self.n_iter):

            # re-discretise wake sheet based on new induction factor
            self.geo.a = np.mean(a[(self.geo.n_span - 1):])  # take only the values of the first rotor

            if self.double_rotor:  # shift filament coords of the 1st rotor
                self.geo.doubleRotorUpdate()
            else:
                self.geo.filaments = self.geo.compute_ring()

            # compute system of linear eqn. (influence of each filament)
            uvw_mat = self._initialize_solver()
            # update circulation
            gamma_curr = gamma_new

            # calculate velocity, circulation, control points
            # directly compute total velocity at each control point by mat. vec. product
            u = uvw_mat[0] @ gamma_curr
            v = uvw_mat[1] @ gamma_curr
            w = uvw_mat[2] @ gamma_curr

            if self.double_rotor:
                # compute perceived velocity by blade element
                vel_rot = np.cross(-self.u_rot, cp)
                vel_per = np.hstack([self.u_inf + u + vel_rot[:, 0].reshape(-1, 1),
                                     v + vel_rot[:, 1].reshape(-1, 1),
                                     w + vel_rot[:, 2].reshape(-1, 1)])

                # calculate azimuthal and axial velocity
                inv_pos_radial = np.hstack([-1 / pos_radial, np.zeros(pos_radial.shape), np.zeros(pos_radial.shape)])
                azim_dir = np.cross(inv_pos_radial, cp)
                u_azim = np.array([azim @ vel for azim, vel in zip(azim_dir, vel_per)])
                u_axial = vel_per @ np.array([1, 0, 0])  # should be the same as [1, 0, 0] @ vel_per (dot product)
            else:
                # compute perceived velocity by blade element
                vel_rot = np.cross(-self.u_rot, self.geo.cp[:, :3])
                vel_per = np.hstack([self.u_inf + u + vel_rot[:, 0].reshape(-1, 1),
                                     v + vel_rot[:, 1].reshape(-1, 1),
                                     w + vel_rot[:, 2].reshape(-1, 1)])

                # calculate azimuthal and axial velocity
                inv_pos_radial = np.hstack([-1 / pos_radial, np.zeros(pos_radial.shape), np.zeros(pos_radial.shape)])
                azim_dir = np.cross(inv_pos_radial, self.geo.cp[:, :3])
                u_azim = np.array([azim @ vel for azim, vel in zip(azim_dir, vel_per)])
                u_axial = vel_per @ np.array([1, 0, 0])  # should be the same as [1, 0, 0] @ vel_per (dot product)

            # calculate loads using BEM
            blade_loads = self._compute_loads_blade(u_axial, u_azim, pos_radial / self.r_rotor)

            # update loads and circulation
            gamma_new = blade_loads[2].reshape(-1, 1)
            a = -(u + vel_rot[:, 0].reshape(-1, 1)) / self.u_inf
            aline = u_azim / (pos_radial.flatten() * self.u_rot[0]) - 1
            f_norm = blade_loads[0]
            f_tan = blade_loads[1]
            gamma = blade_loads[2]
            alpha = blade_loads[3]
            phi = blade_loads[4]

            # check convergence
            err_ref = max(0.001, np.max(np.abs(gamma_new)))  # choose highest value for reference error
            err = np.max(np.abs(gamma_new - gamma_curr)) / err_ref
            error_log.append(err)
            print("iteration: {} | current error: {}".format(i, err))
            if err < self.tol:
                print("solution converged")
                break

            # set new estimate of bound circulation
            gamma_new = (1 - self.weight) * gamma_curr + self.weight * gamma_new
        print("solution unconverged error: {}".format(err)) if err > self.tol else None
        return [a, aline, r_R, f_norm, f_tan, gamma, alpha, phi, error_log]

    def CP_and_CT(self, a, aline, r_R, f_norm, f_tan, v_inf, omega, radius, nblades):

        CT_LLM = []
        CP_LLM = []
        CP_flow = []

        for i in range(len(r_R) - 1):
            r_R_temp = (r_R[i] + r_R[i + 1]) / 2
            drtemp = (-r_R[i] + r_R[i + 1])
            # Prof. Ferreira
            CT_LLM.append((drtemp * f_norm[i] * nblades) / (0.5 * (v_inf ** 2) * np.pi * radius))
            CP_LLM.append((drtemp * f_tan[i] * r_R_temp * omega * nblades) / (0.5 * (v_inf ** 3) * np.pi))
            # another method
            U_tan = r_R_temp * radius * omega * aline[i]
            U_norm = v_inf * (1 - a[i])
            CP_flow.append((drtemp * (f_norm[i] * U_norm - f_tan[i] * U_tan) * nblades) / (0.5 * (v_inf ** 3) *
                                                                                           np.pi * radius))

        # Based on induction factor
        CT_LLM2 = 4 * np.mean(a) * (1 - np.mean(a))
        CP_LLM2 = 4 * np.mean(a) * (1 - np.mean(a))**2

        return [CP_LLM, CT_LLM, CP_LLM2, CT_LLM2]
