import matplotlib.backends
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use('TkAgg')


class Turbine:

    def __init__(self, Nrad=20, rotation=1/2):
        """
        All the properties of the turbine are listed here. This can all be accessed by the following method:
        from turbine import Turb
        Turbine(dr).~variable~
        """
        self.type = 'turbine'        # case
        self.R = 50                  # [m]
        self.n = 3                   # count
        self.blade_start = 0.2       # r/R-ratio Before circular root section assumed without influence
        self.pitch = np.radians(2)  # [rad]
        self.Nrad = Nrad             # Number of radial elements
        self.rotation = rotation             # Number of circular elements
        #self.r_ab = np.linspace(self.blade_start, 1, self.Nrad+1) * self.R  # linear distribution
        self.r_ab = ((-np.cos(np.linspace(0, np.pi, self.Nrad+1)) + 1)/2 * (1 - self.blade_start) + self.blade_start) * self.R  # cosine distribution
        self.r_cp = (self.r_ab[1:] + self.r_ab[:-1])/2
        self.c_cp = 0.25         # quarter chord point
        self.wake_steps = np.arange(0, self.rotation*2*np.pi, np.pi/10)  # steps in the wake for the filament distribution
        self.Nrot = len(self.wake_steps)
        self.beta, self.chord = self.geometry(self.r_cp)     # geometric angle, chord length in control points


    @staticmethod
    def geometry(r_R):
        """"Define the twist and chord length of every annular blade section."""
        beta = np.radians(14 * (1 - r_R) - 2)       # Geometric angle
        chord = 3 * (1 - r_R) + 1                   # local chord length
        return chord, -beta


if __name__ == "__main__":
    print("Turbine")
    T = Turbine(9, 9)
    print(T.r_ab)
