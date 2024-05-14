# ================== Bem Simulation==========================
# authors:Justing Brusche , Cristian Cortonavu
# -------- Pieter van de Ven
# Date: 27 March 2023
# Course: Rotor Wake
# =============================================================
# ------------------ Dependencies ----------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate

# -------------------- Constants ------------------------
rho = 1 # ISA density [kg/m^3]
v_inf = 10 # free stream velocity [m/s]
r_R_start = 0.2 # statring percentage of spa
n_sections = 15 # number of anuli
n_blades = 3 # number of blades
blade_pitch = -2 # pitch angle [deg]
R = 50 # radius
p_inf = 101325 # pressure at infinity [Pa]

# ------------- Airfoil data ---------------------------
# Open data
airfoil = 'Code Jus\polarDU95W180.csv'
data1=pd.read_csv(airfoil, header=0,
                    names = ["Alfa", "Cl", "Cd", "cm"])
polar_alpha = np.array(data1['Alfa'][:])
polar_cl = np.array(data1['Cl'][:])
polar_cd = np.array(data1['Cd'][:])


# --------------------- Script ------------------------
# Interpoolations for cl and cd
f_cl = interpolate.interp1d(polar_alpha, polar_cl)
f_cd = interpolate.interp1d(polar_alpha, polar_cd)

def find_a(CT):
    if CT>1:
        #a = 1-(CT*(1.816-4*(np.sqrt(1.816)-1)))
        a = (CT-1.816+4*(np.sqrt(1.816)-1))/(4*(np.sqrt(1.816)-1))
    else:
        a = (1 - np.sqrt(1 - CT)) / 2
        if 4*(np.sqrt(1.816)-1)-(4-8*a)>0:
            #a = 1 - (CT * (1.816 - 4 * (np.sqrt(1.816) - 1)))
            a = (CT - 1.816 + 4 * (np.sqrt(1.816) - 1)) / (4 * (np.sqrt(1.816) - 1))

    return a

CT_list = np.linspace(0,1.5,100)
yaw = 0
a = 0.3
yaw2 = 0.00001*np.pi/180
a_new_list1 = []
a_new_list2 = []
for i in range(len(CT_list)):
    a_new_list1.append(find_a(CT_list[i]))
    #a_new_list2.append(ainduction_Carlos(CT_list[i], yaw2, a))

r_R_list = np.linspace(r_R_start, 1, n_sections)

#print(r_R_list)

r_R_list_uneven_1 = np.geomspace(r_R_start, 1, n_sections)
r_R_list_uneven_1 = [round(num, 5) for num in r_R_list_uneven_1]
#print(type(r_R_list_uneven_1), type(r_R_list))
r_R_list_uneven_1 = np.array(r_R_list_uneven_1)
#print(type(r_R_list_uneven_1), type(r_R_list))
r_R_list_uneven_2 = [0.2]
b = []
r_R_list = np.flip(r_R_start+(np.cos(np.linspace(0,np.pi,n_sections+1))*0.5+0.5)*(1-r_R_start))


for i in range(len(r_R_list_uneven_1) - 1):
    dspace = r_R_list_uneven_1[i + 1] - r_R_list_uneven_1[i]
    b.append(dspace)

for i in range(len(b)):
    r_R_list_uneven_2.append(r_R_list_uneven_2[i] + b[-i - 1])

def loop1time(R, v_inf,blade_pitch,n_blades,rho,TSR, r_R_list, n_sections):
    min_error = 0.00001 # Convergence error
    # Setup functions
    omega = TSR * v_inf / R
    r_R_start = 0.2  # statring percentage of spa
    # r_R_list = np.linspace(r_R_start, 1, n_sections)
    # delta_R = ((1 - r_R_start) / (n_sections-1)) * R
    #delta_R = []
    #for i in range(len(r_R_list) - 1):
        #dr = (round(r_R_list[i+1], 8) - round(r_R_list[i], 8)) * R
        #delta_R.append(round(dr, 8))
    delta_R = R*(r_R_list[1:]-r_R_list[:-1])

    # Lists for plots
    alpha_list = []
    phi_list = []
    a_list = []
    Ct_list = []
    a_prime_list = []
    azimuthal_loads = []
    thrusts = []
    loads3DN = []
    loads3DT = []
    loads3DNiter = []
    F_root_vals = []
    F_tip_vals  = []
    prandtl_vals = []
    circulation = []
    R_mid_list = []
    iterlist = []
    #r_R_mid_list = np.flip(r_R_start + (np.cos(0.5 * (np.linspace(0, np.pi, n_sections + 1)[1:] + np.linspace(0, np.pi, n_sections + 1)[:-1])) * 0.5 + 0.5) * (1 - r_R_start))

    for i in range(len(r_R_list)-1):
        a = 0.3
        a_prime = 0

        R_inner = r_R_list[i]*R
        R_outer = r_R_list[i+1]*R
        #R_mid = r_R_mid_list[i] * R
        R_mid = 0.5*(R_inner+R_outer)
        print(R_inner,R_outer,R_mid,delta_R[i])
        A_section = (R_outer**2-R_inner**2)*np.pi
        twist = 14*(1-R_mid/R)
        chord = 3*(1-R_mid/R)+1
        iter = 0
        converging = True
        while converging:
            a_last = a+0
            v_local = v_inf*(1-a)
            v_r_local = omega*R_mid*(1+a_prime)

            W = np.sqrt(v_local**2+v_r_local**2)
            phi_rad = np.arctan2(v_local,v_r_local)
            phi_deg = phi_rad * 180/np.pi
            alpha = phi_deg - (twist + blade_pitch)
            print("alpha",alpha)

            L = 0.5 * rho * W**2 * f_cl(alpha) * chord
            D = 0.5 * rho * W**2 * f_cd(alpha) * chord


            N = L*np.cos(phi_rad) + D*np.sin(phi_rad)
            T = L*np.sin(phi_rad) - D*np.cos(phi_rad)
            Load3D_N = N * n_blades * delta_R[i] #Thrust 3D [N/m]
            Load3D_T = T * n_blades * delta_R[i] # azimuthal 3D load [N/m]

            # Load3D_N = Nmultiplied * n_blades
            # Load3D_T = Tmultiplied * n_blades


            Ct_N = Load3D_N / (0.5 * rho * A_section * v_inf ** 2)
            print("CTn",Ct_N)
            Ct_T = Load3D_T / (0.5 * rho * A_section * v_inf ** 2)
            a_new = find_a(Ct_N)
            print("a_new",a_new)
            # Correction factor
            mu = R_mid/R
            f_tip = (2/np.pi)*np.arccos(np.exp(-n_blades/2*((1-mu)/mu)*np.sqrt(1+((TSR*mu)**2)/((1-a_new)**2))))

            f_root= (2/np.pi)*np.arccos(np.exp(-n_blades/2*((mu-r_R_start)/mu)*np.sqrt(1+((TSR*mu)**2)/((1-a_new)**2))))
            a_new = a_new/f_root/f_tip
            print("acor",a_new)
            if R_outer/R<0.3 or R_outer/R>0.9:
                a = 0.99 * a + 0.01 * a_new
            else:
                a = 0.75 * a + 0.25 * a_new
            # calc gamma
            gamma = 0.5 * np.sqrt(v_local**2 + v_r_local**2) * f_cl(alpha) * chord



            a_prime_new = T * n_blades / (2 * np.pi * v_inf * (1 - a) * omega * 2 * (R_mid) ** 2)
            a_prime_new = a_prime_new / f_root / f_tip
            a_prime = 0.75*a_prime+0.25*a_prime_new
            loads3DNiter.append(Load3D_N)
            iter += 1
            if abs(a-a_last)<min_error:
                converging = False
                prandtl = f_tip * f_root
                prandtl_vals.append(prandtl)
                a_list.append(a)
                a_prime_list.append(a_prime)

                alpha_list.append(alpha)
                phi_list.append(phi_deg)
                azimuthal_loads.append(T)
                thrusts.append(N)
                #Ct_list.append(Ct_torque)
                loads3DT.append(Load3D_T)
                loads3DN.append(Load3D_N)
                F_tip_vals.append(f_tip)
                F_root_vals.append(f_root)
                circulation.append(gamma)
                R_mid_list.append(R_mid)
                iterlist.append(iter)


    return a_list, a_prime_list, alpha_list, \
           phi_list, azimuthal_loads, \
           thrusts, \
           loads3DN, loads3DT,\
           F_root_vals, F_tip_vals,prandtl_vals\
           ,circulation, omega,R_mid_list, loads3DNiter, iterlist





def plot_aoa(alpha,r_R_list,phi):

    plt.grid()
    plt.plot(r_R_list[:-1],phi, label = 'Inflow Angle')
    plt.plot(r_R_list[:-1],alpha,label = 'Angle of Attack')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Angle [deg]')
    plt.title('Change of Angle of Attack over span')
#     plt.savefig('AOA_change_rR')
    return plt.show()

def plot_inductions(a,r_R_list, a_prime):
    plt.grid()
    plt.plot(r_R_list[:-1], a, label='a')
    plt.plot(r_R_list[:-1], a_prime, label='a')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Induction factor [deg]')
    plt.title('axial and azimuthal inductions')
#     plt.savefig('inductions')
    return plt.show()

def plot_loads(thrusts,r_R_list,azimuthal_loads):
    plt.grid()
    plt.plot(r_R_list[:-1], np.array(thrusts) * n_blades, label='Thrust')
    plt.plot(r_R_list[:-1], np.array(azimuthal_loads)* n_blades, label='Azimuthal load')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Loads [N/m]')
    plt.title('Thrust and Azimuthal Load over the span')
#     plt.savefig('Thrust_load_span')
    return plt.show()

def plot_loadsvsTSR(R, v_inf, blade_pitch, n_blades, rho, r_R_list):
    TSR = np.arange(5,12,1)
    load_distributionN = []
    load_distributionT = []
    for i in range(len(TSR)):
        a_list, a_prime_list, alpha_list, \
        phi_list, azimuthal_loads, \
        thrusts, totalNload, totalTload, \
        F_root_vals, F_tip_vals, \
        prantl_vals, circulation, omega,R_mid_list, \
        loads3DTiter, iterlist= loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR[i], r_R_list, n_sections)
        load_distributionN.append(sum(totalNload))
        load_distributionT.append(sum(totalTload))
    plt.grid()
    plt.plot(TSR, load_distributionT, label='Total load in azimuthal direction', marker = '.')
    plt.plot(TSR,load_distributionN , label='Total load in thrust direction', marker = '.')
    plt.legend()
    plt.xlabel('tip-speed Ratio[-]')
    plt.ylabel('Loads [N]')
    plt.title('Total loads vs tip-speed ratio')
#     plt.savefig('LoadsvsTSR')
    return plt.show()

def plot_prandtl(F_tip_vals,F_root_vals,r_R_list,prandtl_vals):
    ax = plt.figure()
    ax.set_figwidth(10)
    ax.set_figheight(2.5)
    plt.grid()
    plt.plot(r_R_list[:-1], F_tip_vals, label='prandtl tip', marker ='.')
    plt.plot(r_R_list[:-1], F_root_vals, label='prandtl root', marker = '.')
    plt.plot(r_R_list[:-1], prandtl_vals, label = 'prandtl', marker = '.', markersize=2)
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Loads [N]')
    plt.title('Prandtl tip-root correction')
#     plt.savefig('Prandtl correection')
    return ax.show()

def plot_circulation(r_R_list, circulation,omega):
    ax = plt.figure()
    ax.set_figwidth(10)
    ax.set_figheight(2.5)
    circulation = np.array(circulation)
    plt.grid()
    plt.plot(r_R_list[:-1], circulation/(np.pi*v_inf**2/(n_blades*omega)), label='Circulation', marker ='.')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Non dimensionalised circulation [-]')
    plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega * NBlades } $')
#     plt.savefig('Circulation')
    return plt.show()

def plot_spacing(r_R_list, r_R_list_uneven_1, r_R_list_uneven_2, TSR):
    plt.grid()
    r_R_list_uneven_1 = np.array(r_R_list_uneven_1)
    r_R_list_uneven_2 = np.array(r_R_list_uneven_2)
    thrust1 = np.array(loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list, 81)[5])
    thrust2 = np.array(loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_uneven_1, 81)[5])
    thrust3 = np.array(loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_uneven_2, 81)[5])
    plt.plot((r_R_list[:-1]+r_R_list[1:])/2, thrust1*n_blades, label='Thrust evenly spaced')
    plt.plot((r_R_list_uneven_1[:-1]+r_R_list_uneven_1[1:])/2, thrust2*n_blades, label='Thrust concentrated begin')
    plt.plot((r_R_list_uneven_2[:-1]+r_R_list_uneven_2[1:])/2, thrust3*n_blades, label='Thrust concentrated end')

    #print('sums of thursses', sum(thrust1), sum(thrust2), sum(thrust3))
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('thrust [N/m]')
    plt.title('Thrust distribution for different spacing')
#     plt.savefig('Thrust_spacing')
    return plt.show()

def plot_annuli(TSR):
    plt.grid()
    r_R_list_7 = np.linspace(r_R_start, 1, 11)
    r_R_list_8 = np.linspace(r_R_start, 1, 21)
    r_R_list_9 = np.linspace(r_R_start, 1, 31)
    r_R_list_10 = np.linspace(r_R_start, 1, 41)
    r_R_list_1 = np.linspace(r_R_start, 1, 51)
    r_R_list_2 = np.linspace(r_R_start, 1, 61)
    r_R_list_3 = np.linspace(r_R_start, 1, 71)
    r_R_list_4 = np.linspace(r_R_start, 1, 81)
    r_R_list_5 = np.linspace(r_R_start, 1, 91)
    r_R_list_6 = np.linspace(r_R_start, 1, 101)
    thrust1 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_1, len(r_R_list_1))[5]
    thrust2 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_2, len(r_R_list_2))[5]
    thrust3 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_3, len(r_R_list_3))[5]
    thrust4 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_4, len(r_R_list_4))[5]
    thrust5 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_5, len(r_R_list_5))[5]
    thrust6 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_6, len(r_R_list_6))[5]
    thrust7 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_7, len(r_R_list_7))[5]
    thrust8 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_8, len(r_R_list_8))[5]
    thrust9 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_9, len(r_R_list_9))[5]
    thrust10 = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list_10, len(r_R_list_10))[5]

    #print('sum of thrusses', sum(thrust6), sum(thrust7), sum(thrust8), sum(thrust9), sum(thrust10))

    plt.plot((r_R_list_7[:-1]+r_R_list_7[1:])/2, thrust7, label='10 sections')
    plt.plot((r_R_list_8[:-1]+r_R_list_8[1:])/2, thrust8, label='20 sections')
    plt.plot((r_R_list_9[:-1]+r_R_list_9[1:])/2, thrust9, label='30 sections')
    plt.plot((r_R_list_10[:-1]+r_R_list_10[1:])/2, thrust10, label='40 sections')
    plt.plot((r_R_list_1[:-1]+r_R_list_1[1:])/2, thrust1, label='50 sections')
    plt.plot((r_R_list_2[:-1]+r_R_list_2[1:])/2, thrust2, label='60 sections')
    plt.plot((r_R_list_3[:-1]+r_R_list_3[1:])/2, thrust3, label='70 sections')
    plt.plot((r_R_list_4[:-1]+r_R_list_4[1:])/2, thrust4, label='80 sections')
    plt.plot((r_R_list_5[:-1]+r_R_list_5[1:])/2, thrust5, label='90 sections')
    plt.plot((r_R_list_6[:-1]+r_R_list_6[1:])/2, thrust6, label='100 sections')
    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('thrust [N/m]')
    plt.title('Thrust for different number of sections')
#     plt.savefig('Thrust_spacing')
    return plt.show()



def plot_Pressure(p_inf,r_R_list,a_list,R_mid,Load_3DN):
    p1 = np.zeros(len(r_R_list)-1)+p_inf+0.5*rho*v_inf**2
    p2 = p1+0
    A_list = ((((r_R_list*R)[1:])**2)-(((r_R_list*R)[:-1])**2))*np.pi
    p3 = p2-Load_3DN/A_list
    p4=p3+0
    r_mid_list = ((r_R_list*R)[1:]+(r_R_list*R)[:-1])/2
    plt.plot(r_mid_list, p1, label="p1")
    plt.plot(r_mid_list,p2,label="p2")
    plt.plot(r_mid_list,p3,label="p3")
    plt.plot(r_mid_list, p4, label="p4")

    plt.legend()
    plt.xlabel('r/R [-]')
    plt.ylabel('Stagnation pressure [Pa]')
    plt.title('Stagnation pressure distribution at the 4 stations')
#     plt.savefig('pressure at different annuli')
    return plt.show()


def plot_convergence(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list, n_sections):
    iterlist =     loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list, n_sections)[-1]
    loads3DTiter = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list, n_sections)[-2]
    a = iterlist[0]+iterlist[1]+iterlist[2]+iterlist[3]
    b = iterlist[0]+iterlist[1]+iterlist[2]+iterlist[3] + iterlist[4]
    thrust = loads3DTiter[a:b]
    y = np.linspace(1, len(thrust), len(thrust))
    plt.grid()
    plt.plot(y, thrust, label=' Thrust of one section')
    plt.legend()
    plt.xlabel('iteration [-]')
    plt.ylabel('Thrust [N]')
    plt.title('Convergence rate')
    return plt.show()

def main_bem(R, v_inf, blade_pitch, n_blades, rho, r_R_list, n_sections, TSR):
    # ----------------- loop ---------------------
    a_list, a_prime_list, alpha_list,\
    phi_list, azimuthal_loads,\
    thrusts, Load_3DN, Load_3DT, \
    F_root_vals, F_tip_vals, \
    prantl_vals, circulation,omega,R_mid_list, \
    loads3DTiter, iterlist = loop1time(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list, n_sections)
    result_bem_dict = {}
    result_bem_dict["circulation"] = circulation
    result_bem_dict["alpha"] = alpha_list
    result_bem_dict["phi"] = phi_list
    result_bem_dict["thrust"] = thrusts
    result_bem_dict["azimuthal_load"] =azimuthal_loads
    result_bem_dict["a_list"] = a_list
    result_bem_dict["a_prime_list"] = a_prime_list
    R_mid = 0.5*(r_R_list[1:]+r_R_list[:-1])
    # plt.plot(R_mid,np.array(thrusts)/(0.5 * v_inf**2 * 50),"-x")
    # plt.plot(R_mid, np.array(azimuthal_loads) / (0.5 * v_inf ** 2 * 50), "-x")
    #plt.show()


    #-------------- 1 TSR val plots ----------------
    # plot_aoa(alpha_list, r_R_list, phi_list)
    # plot_inductions(a_list, r_R_list, a_prime_list)
    # plot_loads(thrusts, r_R_list, azimuthal_loads)
    # plot_prandtl(F_tip_vals,F_root_vals,r_R_list,prantl_vals)
    # plot_circulation(r_R_list,circulation, omega)
    # plot_Pressure(p_inf,r_R_list,a_list,R_mid_list,Load_3DN)
    # plot_convergence(R, v_inf, blade_pitch, n_blades, rho, TSR, r_R_list, n_sections)
    # plot_annuli(TSR)
    # # plot_spacing(r_R_list, r_R_list_uneven_1, r_R_list_uneven_2, TSR)
    #
    #
    # # --------------- for different TSR vals -----------------
    # plot_loadsvsTSR(R, v_inf, blade_pitch, n_blades, rho, r_R_list)
    # plot_spacing(r_R_list, r_R_list_uneven_1, r_R_list_uneven_2, TSR)




    return result_bem_dict

if __name__ == '__main__':
    r_R_list = np.flip(r_R_start + (np.cos(np.linspace(0, np.pi, n_sections + 1)) * 0.5 + 0.5) * (1 - r_R_start))

    #main(R, 10, blade_pitch, n_blades, rho, r_R_list, n_sections,6)
    main_bem(50, 10, -2,3, 1, r_R_list, n_sections, 6)














