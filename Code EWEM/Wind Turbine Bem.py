#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import math as m
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

files = ['Code EWEM\polar_DU95W180.xlsx']

#Initializing tables    
cl_tab=np.zeros([61])
cd_tab=np.zeros([61])
cm_tab=np.zeros([61])
aoa_tab=np.zeros([61])
#Readin of tables. Only do this once at startup of simulation
aoa_tab[:],cl_tab[:],cd_tab[:],cm_tab[:] = pd.read_excel(files[0],header=None,skiprows=4).values.T


#Functions____________
def plots(i):
    plt.figure()
    plt.plot(r/R, flowAngle_lst[:,i],label='Inflow Angle')
    plt.plot(r/R, localalpha_lst[:,i],label='Local Angle of Attack')
    plt.xlabel('r/R')
    plt.ylabel('Angle [deg]')
    plt.legend()
    plt.title('Inflow and Local Angle of Attack TSR = '+str(TSR[i]))

    plt.figure()
    plt.plot(r/R, axialInduction_lst[:,i],label='Axial Induction')
    plt.plot(r/R, azimuthalInduction_lst[:,i],label='Azimuthal Induction')
    plt.xlabel('r/R')
    plt.ylabel('Induction Factor')
    plt.legend()
    plt.title('Axial and Azimuthal Induction TSR = '+str(TSR[i]))

def force_coeffs(localalpha,aoa_tab,cl_tab,cd_tab,cm_tab):
    Cl=np.interp (localalpha,aoa_tab,cl_tab)
    Cd=np.interp (localalpha,aoa_tab,cd_tab)
    Cm=np.interp (localalpha,aoa_tab,cm_tab)
    return Cl, Cd, Cm 

def BEM(TSR,pitch,r,c,twist,aoa_tab,cl_tab,cd_tab,cm_tab,Vo):
    omega = TSR*Vo/R
    a = 0
    aprime = 0
    convergenceFactor = 1e-10
    delta = 1
    deltaPrime = 1

    relax = 0.1
    solidity = (B*c)/(2*m.pi*r)
    count = 0

    while(delta > convergenceFactor and deltaPrime > convergenceFactor):
        count = count + 1
        if (count > 1e4):
            print("No convergence!")
            break

        flowAngle = m.atan(((1 - a) * R) / ((1 + aprime) * TSR * r))
        localalpha =  m.degrees(flowAngle) - (pitch + twist)

        Cl,Cd,Cm = force_coeffs(localalpha, aoa_tab, cl_tab, cd_tab, cm_tab)
        Ct = Cl*m.sin(flowAngle) - Cd*m.cos(flowAngle)
        Cn = Cl*m.cos(flowAngle) + Cd*m.sin(flowAngle)

        F = 2 / m.pi * m.acos(m.exp(-B * (R-r) / (2 * r *m.sin(abs(flowAngle)))))
        CT = ((1 - a)**2 * Cn * solidity) / m.sin(flowAngle)**2

        aold = a
        if(aold < 0.33):
            a = (solidity * Cn * (1 - aold)) / (4 * F * m.sin(flowAngle)**2)
        else:
            aStar = CT / (4 *F*(1 - 1/4 * (5 - 3 * aold) * aold))
            a = relax * aStar + (1 - relax) * aold

        aprimeOld  = aprime
        aprimeStar = (solidity * Ct * (1 + aprimeOld)) / (4 * F * m.sin(flowAngle) * m.cos(flowAngle))
        aprime = relax * aprimeStar + (1 - relax) * aprimeOld

        delta = abs(aprime - aprimeOld)
        deltaPrime = abs(aprime - aprimeOld)

    Vrel = m.sqrt(Vo**2+ (omega * r)**2)

    Pn = 0.5 * rho * Vrel**2 * c * Cn
    Pt = 0.5 * rho * Vrel**2 * c * Ct
    Gamma = 0.5 * Vrel * c * Cl
    Gamma_nondim = Gamma * (B* omega) / (Vo**2 * np.pi)
    print('V_rel = ', Vrel)
    if (m.isnan(Pt)|(m.isnan(Pn))):
        Pt, Pn = 0,0

    return Pn, Pt, flowAngle, localalpha, a, aprime, Vrel, Gamma_nondim

def single_BEM_loop(TSR,i):
    omega = TSR*Vo/R
    for k in range(len(r)):
        Pn, Pt, flowAngle, localalpha, a, aprime, Vrel, Gamma = BEM(TSR,pitch,r[k],chord[k],twist[k],aoa_tab,cl_tab,cd_tab,cm_tab,Vo)
        # print(r[k], Pn)
        Pn_lst[k] = Pn
        Pt_lst[k] = Pt
        flowAngle_lst[k,i] = m.degrees(flowAngle)
        localalpha_lst[k,i] = localalpha
        axialInduction_lst[k,i] = a
        azimuthalInduction_lst[k,i] = aprime
        Vrel_lst[k,i] = Vrel
        gamma_lst[k,i] = Gamma
    T = np.trapz(Pn_lst, r) * B
    P = np.trapz(Pt_lst * r, r) * omega * B

    Cp = P/(0.5*rho*Vo**3*m.pi*R**2)
    Ct = T/(0.5*rho*Vo**2*m.pi*R**2)
    

    return P, T, Cp, Ct, flowAngle_lst, localalpha_lst, a, aprime, Vrel_lst
#Constants______________

R = 50 #m
r = np.arange(0.2*R, R+0.1, 0.1)
B = 3
rho = 1.225 #kg/m3
Vo = 10 #m/s

TSR =[6,8,10]


pitch = -2
twist = 14*(1-r/R) 
chord = 3*(1-r/R)+1
yaw = 0

Pn_lst = np.zeros(len(r))
Pt_lst = np.zeros(len(r))

#Stored for plotting
flowAngle_lst = np.zeros([len(r),len(TSR)])
localalpha_lst = np.zeros([len(r),len(TSR)])
axialInduction_lst = np.zeros([len(r),len(TSR)])
azimuthalInduction_lst = np.zeros([len(r),len(TSR)])
Vrel_lst = np.zeros([len(r),len(TSR)])
gamma_lst = np.zeros([len(r),len(TSR)])

for i in range(len(TSR)):
    P, T, Cp, Ct, flowAngle_lst, localalpha_lst, a, aprime, Vrel_lst = single_BEM_loop(TSR[i],i)
    plots(i)
# save vrel_lst to a txt file
np.savetxt('vrel_lst.txt', Vrel_lst, delimiter=',')
# save gamma_lst to a txt file
np.savetxt('gamma_lst.txt', gamma_lst, delimiter=',')
plt.show()

# %%
