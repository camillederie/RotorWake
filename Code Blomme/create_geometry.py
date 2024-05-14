"""
Class object which discretises the rotor blade(s) into bound and trailing vortex
filaments
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

class BladeGeometry:

    def __init__(self, radius, tsr, v_inf, n_blades, n_span, n_theta, spacing, phase_diff, a, double_rotor=False,
                 xshift=0, yshift=100, zshift=0):
        # todo: check whether non-dim of span_arr is needed
        self.radius = radius
        self.tsr = tsr
        self.v_inf = v_inf
        self.n_blades = n_blades
        self.n_span = n_span
        self.n_theta = n_theta
        self.a = a
        self.span_arr = np.linspace(0.2, 1.0, n_span)
        self.theta_arr = np.linspace(0, 30 * np.pi, n_theta)
        self.phase_diff = np.radians(phase_diff)

        self.xshift = xshift
        self.yshift = yshift
        self.zshift = zshift

        self.rRootRatio = 0.2
        self.rTipRatio = 1.0
        self.spacing = spacing
        self.double_rotor = double_rotor

        self.discretize_spanwise()
        # logic to have consistent filaments array size
        if double_rotor:
            # initialize arrays that are twice the normal size (twice the single rotor case)
            n_single_rotor = n_blades * (n_span - 1)
            self.filaments = np.zeros((7, 2 * n_single_rotor, 2 * n_theta + 1))  # rows are double
            self.cp = np.zeros((2 * n_single_rotor, 9))  # coord; normal; tangential
            self.bladepanels = np.zeros((2 * n_single_rotor, 4 * 3))  # empty dict to

            self.cp[:n_single_rotor, :] = self._compute_cp()
            self.filaments[:, :n_single_rotor] = self.compute_ring()  # second rotor rows are empty
            self.bladepanels[:n_single_rotor, :] = self.discretize_blade()

            self.doubleRotor()
            self.doubleRotorUpdate()

        else:
            self.cp = self._compute_cp()
            self.bladepanels = self.discretize_blade()
            self.filaments = self.compute_ring()



    def discretize_spanwise(self):
        if self.spacing == 'equal':
            rsegment = np.linspace(self.rRootRatio, 1, self.n_span)

        elif self.spacing == 'cosine':
            midpoint = (1 - self.rRootRatio) / 2
            angleInc = np.pi / (self.n_span - 1)
            curAngle = angleInc
            rsegment = np.zeros(self.n_span)
            rsegment[0] = self.rRootRatio
            rsegment[-1] = 1
            for idx in range(1, self.n_span):
                rsegment[idx] = self.rRootRatio + midpoint * (1 - np.cos(curAngle))
                curAngle += angleInc

        else:
            print("Choose equal or cosine spacing")
            sys.exit()
        self.span_arr = rsegment
        return

    def discretize_blade(self, phase_diff=0.0):

        blade_panels = np.zeros((self.n_blades * (self.n_span - 1), 4 * 3))
        for blade in range(self.n_blades):
            bladeRot = 2 * np.pi / self.n_blades * blade + phase_diff
            angle1 = np.deg2rad(- 14 * (1 - self.span_arr[:-1]) + 2)
            angle2 = np.deg2rad(- 14 * (1 - self.span_arr[1:]) + 2)
            chord1 = (3 * (1 - self.span_arr[:-1]) + 1) / self.radius
            chord2 = (3 * (1 - self.span_arr[1:]) + 1) / self.radius

            # define the 4 corners
            p1 = [-0.25 * chord1 * np.sin(-angle1), self.span_arr[:-1], 0.25 * chord1 * np.cos(angle1)]
            p2 = [-0.25 * chord2 * np.sin(-angle2), self.span_arr[1:], 0.25 * chord2 * np.cos(angle2)]
            p3 = [0.75 * chord2 * np.sin(-angle2), self.span_arr[1:], -0.75 * chord2 * np.cos(angle2)]
            p4 = [0.75 * chord1 * np.sin(-angle1), self.span_arr[:-1], -0.75 * chord1 * np.cos(angle1)]

            # rotate coordinates
            p1Rot = np.column_stack([np.zeros(self.n_span - 1), p1[1] * np.cos(bladeRot) - p1[2] * np.sin(bladeRot),
                                     p1[1] * np.sin(bladeRot) + p1[2] * np.cos(bladeRot)])
            p2Rot = np.column_stack([np.zeros(self.n_span - 1), p2[1] * np.cos(bladeRot) - p2[2] * np.sin(bladeRot),
                                     p2[1] * np.sin(bladeRot) + p2[2] * np.cos(bladeRot)])
            p3Rot = np.column_stack([np.zeros(self.n_span - 1), p3[1] * np.cos(bladeRot) - p3[2] * np.sin(bladeRot),
                                     p3[1] * np.sin(bladeRot) + p3[2] * np.cos(bladeRot)])
            p4Rot = np.column_stack([np.zeros(self.n_span - 1), p4[1] * np.cos(bladeRot) - p4[2] * np.sin(bladeRot),
                                     p4[1] * np.sin(bladeRot) + p4[2] * np.cos(bladeRot)])

            # write to bladepanels
            blade_panels[blade * (self.n_span - 1):blade * (self.n_span - 1) + (self.n_span - 1), :] = \
                np.column_stack((p1Rot, p2Rot, p3Rot, p4Rot)) * self.radius

        return blade_panels

    def compute_ring(self, phase_diff=0.0):
        # create container array for filaments
        filaments = np.zeros((7, self.n_blades * (self.n_span - 1), 2 * self.n_theta + 1))
        # loop over different blades
        for blade in range(self.n_blades):
            r = self.span_arr
            bladeRot = 2 * np.pi / self.n_blades * blade + phase_diff

            # loop             
            for idx, span in enumerate(self.span_arr[:-1]):
                data_arr = np.empty((0, 7))  # rows: thetas | columns
                chord1 = (3 * (1 - span) + 1) / self.radius
                angle1 = np.deg2rad(- 14 * (1 - span) + 2)

                # Define bound filament
                # [x1, y1, z1, x2 ,y2, z2, gamma]
                data_arr = np.vstack((data_arr, [0, r[idx], 0, 0, r[idx + 1], 0, 0]))

                # Define trailing filaments at x1
                data = [chord1 * np.sin(-angle1), r[idx], -chord1 * np.cos(angle1), 0, r[idx], 0, 0]
                data_arr = np.vstack((data_arr, data))

                # Loop over trailing filaments at x1
                xt = chord1 * np.sin(-angle1)
                yt = r[idx]
                zt = -chord1 * np.cos(angle1)
                dx = np.cumsum(np.ones(len(self.theta_arr[:-1])) * self.theta_arr[1] / (
                            self.tsr / (1 - self.a)))
                dy = np.cumsum(np.cos(-self.theta_arr[1:]) - np.cos(-self.theta_arr[:-1])) * r[idx]
                dz = np.cumsum(np.sin(-self.theta_arr[1:]) - np.sin(-self.theta_arr[:-1])) * r[idx]
                data = np.vstack([xt + dx, yt + dy, zt + dz, np.insert((xt + dx)[:-1], 0, xt),
                                  np.insert((yt + dy)[:-1], 0, yt), np.insert((zt + dz)[:-1], 0, zt),
                                  np.zeros(len(self.theta_arr[:-1]))]).T
                data_arr = np.vstack((data_arr, data))

                # trailing filaments at x2
                chord2 = (3 * (1 - r[idx + 1]) + 1) / self.radius
                angle2 = np.deg2rad(- 14 * (1 - r[idx + 1]) + 2)
                data = [0, r[idx + 1], 0, chord2 * np.sin(-angle2), r[idx + 1], -chord2 * np.cos(angle2), 0]
                data_arr = np.vstack((data_arr, data))

                xt = chord2 * np.sin(-angle2)
                yt = r[idx + 1]
                zt = -chord2 * np.cos(angle2)
                dx = np.cumsum(np.ones(len(self.theta_arr[:-1])) * self.theta_arr[1] / (
                            self.tsr / (1 - self.a)))
                dy = np.cumsum(np.cos(-self.theta_arr[1:]) - np.cos(-self.theta_arr[:-1])) * r[idx + 1]
                dz = np.cumsum(np.sin(-self.theta_arr[1:]) - np.sin(-self.theta_arr[:-1])) * r[idx + 1]
                data = np.vstack([np.insert((xt + dx)[:-1], 0, xt), np.insert((yt + dy)[:-1], 0, yt),
                                  np.insert((zt + dz)[:-1], 0, zt),
                                  xt + dx, yt + dy, zt + dz, np.zeros(len(self.theta_arr[:-1]))]).T
                data_arr = np.vstack((data_arr, data))

                # rotate the filaments to correspond with blade orientation
                for filament in data_arr:
                    temp1 = filament.copy()  # store non-rotated filament for subsequent rotation
                    filament[1] = temp1[1] * np.cos(bladeRot) - temp1[2] * np.sin(bladeRot)  # y1
                    filament[2] = temp1[1] * np.sin(bladeRot) + temp1[2] * np.cos(bladeRot)  # z1
                    filament[4] = temp1[4] * np.cos(bladeRot) - temp1[5] * np.sin(bladeRot)  # y2
                    filament[5] = temp1[4] * np.sin(bladeRot) + temp1[5] * np.cos(bladeRot)  # z2

                filaments[:, blade * (self.n_span - 1) + idx, :] = data_arr.T * self.radius

        return filaments

    def _compute_cp(self, phase_diff=0.0):

        cp = np.zeros((self.n_blades * (self.n_span - 1), 9))  # coord; normal; tangential
        segmentCenter = self.radius * (self.span_arr + (np.roll(self.span_arr, -1) - self.span_arr) / 2)[:-1]
        self.centerPoints = segmentCenter
        for blade in range(self.n_blades):
            bladeRot = 2 * np.pi / self.n_blades * blade + phase_diff
            angle = np.deg2rad(2 - 14 * (1 - segmentCenter))

            # bound edge
            boundEdge = np.column_stack(
                (np.zeros(self.n_span - 1), segmentCenter * np.cos(bladeRot), segmentCenter * np.sin(bladeRot)))
            tangVect = np.column_stack(
                (np.cos(angle), np.sin(angle) * np.sin(bladeRot), -np.sin(angle) * np.cos(bladeRot)))
            normVect = np.column_stack(
                (-np.sin(angle), np.cos(angle) * np.sin(bladeRot), -np.cos(angle) * np.cos(bladeRot)))
            # Assign to cp
            cp_iblade = np.column_stack((boundEdge, tangVect, normVect))
            # return [coord, norm, tang] x,y,z
            cp[blade * (self.n_span - 1):blade * (self.n_span - 1) + self.n_span - 1, :] = cp_iblade
        return cp

    def doubleRotor(self):
        # shift control points
        idxMid = int(np.shape(self.cp)[0] / 2)
        self.cp[idxMid:, :] = self._compute_cp(self.phase_diff)
        self.cp[idxMid:, 1] += self.yshift
        # shift blade
        self.bladepanels[idxMid:, :] = self.discretize_blade(self.phase_diff)
        self.bladepanels[idxMid:, (1, 4, 7, 10)] += self.yshift

    def doubleRotorUpdate(self):
        # shift filaments
        idxMid = int(np.shape(self.filaments)[1] / 2)
        # assume that size filaments [7, 2 * nspan * nblade,2 * (2 * theta)
        self.filaments[:, :idxMid] = self.compute_ring()
        self.filaments[:, idxMid:] =  self.compute_ring(self.phase_diff)
        self.filaments[(1, 4), idxMid:] += self.yshift
