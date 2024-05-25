import matplotlib.backends
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate

matplotlib.use('TkAgg')


class Data:

    def __init__(self, case="turbine"):
        if case == 'turbine':
            self.airfoil = 'Code Pim\polarDU95W180.csv'
        if case == 'propellor':
            self.airfoil = 'no file'

        self.alpha, self.cl, self.cd = self.extract_data()

    def extract_data(self):
        # And reading the airfoil data from the excel file:
        if self.airfoil == 'no file':
            print("experimental data file not given...........")
        data1 = pd.read_csv(self.airfoil, delimiter=';', header=1, names=["Alfa", "Cl", "Cd", "Cm"])
        polar_alpha = np.array(data1['Alfa'][:])
        polar_cl = np.array(data1['Cl'][:])
        polar_cd = np.array(data1['Cd'][:])
        polar_cm = np.array(data1['Cm'][:])
        return polar_alpha, polar_cl, polar_cd

    def lift(self, alpha):
        "Interpolate Cl at given alpha"
        cl = interpolate.interp1d(x=self.alpha, y=self.cl)
        return cl(alpha)

    def drag(self, alpha):
        "Interpolate Cd at given alpha"
        cd = interpolate.interp1d(x=self.alpha, y=self.cd)
        return cd(alpha)


if __name__ == "__main__":
    d = Data("turbine")
    a, cl, cd = d.extract_data()
    plt.plot(a, cl)
    plt.plot(a, cd)
    #print(a, cl, cd)

