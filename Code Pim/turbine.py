import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends
matplotlib.use('TkAgg')


class Turbine:

    def __init__(self, dr):
        """
        All the properties of the turbine are listed here. This can all be accessed by the following method:
        from turbine import Turb
        Turbine(dr).~variable~
        """
        self.type = 'turbine'       # case
        self.R = 50                 # [m]
        self.n = 3                  # count
        self.blade_start = 0.2      # r/R-ratio Before circular root section assumed without influence
        self.dr = dr                # annular thickness [m]
        self.r_R = np.arange(self.blade_start + self.dr / 2, 1, self.dr)  # center of blade cells
        self.r_edges = np.arange(self.blade_start, 1 + self.dr, self.dr)  # edges of blade cells
        self.r = self.r_R * self.R  # blade cell center coordinates
        self.twist, self.chord = self.geometry(self.r_R)
        self.pitch = np.radians(2)  # [rad]
        self.Nrad = (1-self.blade_start)/self.dr
    @staticmethod
    def operations(i, j):
        """"Operational conditions: tip speed ratio and yaw angle"""
        TSR = [6, 8, 10, 7]
        yaw = [0, 15, 30]
        return TSR[i], yaw[j]

    @staticmethod
    def geometry(r_R):
        """"Define the twist and chord length of every annular blade section."""
        twist = np.radians(-14 * (1 - r_R))
        chord = 3 * (1 - r_R) + 1
        return twist, chord

    def angular_area(self, radius, dpsi=1, yaw=False):
        """"Area of the blade segment"""
        #HIER MOET HET EEN HALVE ZIJN NU HEB JE 2X DE AREA
        A = np.pi * 2 *self.r *self.dr
        if yaw:
            A = A * dpsi / (2 * np.pi)
        return A


if __name__ == "__main__":
    print("Turbine")
    turbine = Turbine(0.01)

    # plt.plot(turbine.r_R, turbine.angular_area(turbine.r_R), label='twist')
    # plt.grid()
    # plt.ylabel('area')
    # plt.xlabel('r/R')
    # plt.legend()
    # plt.show()
