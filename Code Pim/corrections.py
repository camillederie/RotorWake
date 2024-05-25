import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends
matplotlib.use('TkAgg')


def prandtl(r_R, n, TSR, blade_start, a):
    """Prandtl's correction for finite number of blades:"""
    mu = r_R
    f_tip_mu = (2 / np.pi) * np.arccos(np.exp(-(n / 2) * ((1 - mu) / mu) *
                                              np.sqrt(1 + (TSR**2 * mu**2) /
                                                      (1 - a)**2)))
    f_root_mu = (2 / np.pi) * np.arccos(np.exp(-(n / 2) * ((mu - blade_start) / mu) *
                                               np.sqrt(1 + (TSR**2 * mu**2) /
                                                       (1 - a)**2)))
    f_total_mu = f_tip_mu * f_root_mu
    return f_total_mu, f_tip_mu, f_root_mu


def glauert(a, corr=False):
    """CT correction for heavily loaded rotors:"""
    c_t = 4 * a * (1 - a)
    if corr is True:
        c_t1 = 1.816
        a1 = (1 - np.sqrt(c_t1) / 2)
        c_t[a>a1] = c_t1 - 4 * (np.sqrt(c_t1) - 1) * (1 - a[a>a1])
    return c_t


def glauert2(c_t,yaw):
    #a corrected for heavily loaded rotors:
    a = np.zeros(np.shape(c_t))
    c_t1 = 1.816
    c_t2 = 2 * np.sqrt(c_t1) - c_t1
    if yaw == 0:
        # Check the C_T value
        a[c_t>=c_t2] = 1 + (c_t[c_t>=c_t2]-c_t1)/(4*(np.sqrt(c_t1)-1))
        a[c_t<c_t2] = 0.5-0.5*np.sqrt(1-c_t[c_t<c_t2])
    else:
        print("This is the yaw case")
        #here we need to code the correction for a for the yawed case

    return a



"""def glauert2(c_t):
    #a corrected for heavily loaded rotors:
    c_t1 = 1.816
    c_t2 = 2 * np.sqrt(c_t1) - c_t1
    for i in range(len(c_t)):
        if c_t[i]>c_t2:
            a = (c_t[i]-c_t1+4*(np.sqrt(c_t1)-1))/(4*(np.sqrt(c_t1)-1))
        else:
            a = (1 - np.sqrt(1 - c_t[i])) / 2
            if 4 * (np.sqrt(c_t1) - 1) - (4 - 8 * a) > 0:
                a = (c_t[i] - c_t1 + 4 * (np.sqrt(c_t1) - 1)) / (4 * (np.sqrt(c_t1) - 1))
    return a
"""
def glauert3(a, theta, glauert=False):
    """"
    Yaw-condition;
    correction for CT due to low angle of attack for the lifting line
    This means F and u_ind can be assumed to be normal to the rotor plane.
    """
    c_t = 4 * a * (np.cos(theta) - a)
    if glauert:
        c_t = 4 * a * np.sqrt(1 - a * (2 * np.cos(theta) - a))
    c_p = c_t * (np.cos(theta) - a)
    return c_t, c_p

def glauert4(a, u0, K, xsi, r_R, psi):
    """"
    Induced velocity on the rotor is dependent on the azimuthal position,
    whether there is a positive or negative addition from the wake by u1.
    Note: u0 is not the free-stream velocity, this is the average induced velocity by the wake.
    Normally this is set to (2*)a * U0, but the vortex cylinder model tells different.
    Add colemanse model to get K(xsi)
    """
    u1 = u0 * 2 * np.tan(xsi/2) * r_R * np.sin(psi)
    return u0 + u1

if __name__ == "__main__":
    print("corrections")
    # Prandtl
    r_R = np.arange(0.2, 1, 0.01)
    pr_tot, pr_tip, pr_root = prandtl(r_R,
                                      3,
                                      8,
                                      0.2,
                                      0.3)
    plt.plot(r_R, pr_tot, 'b-', label='Total Prandtl correction')
    plt.plot(r_R, pr_tip, 'r.',label='Tip Prandtl correction')
    plt.plot(r_R, pr_root, 'g.',label='Root Prandtl correction')
    plt.ylabel('f(r/R)')
    plt.xlabel('r/R')

    plt.grid()
    plt.legend()
    plt.show()

    # Glauert
    a = np.arange(-0.5, 1.0, 0.01)
    Ct = glauert(a, False)
    Ct_corr = glauert(a, True)

    plt.plot(a, Ct, 'b-', label='$C_T$')
    plt.plot(a, Ct_corr, 'r--', label='$C_T$ Glauert correction')
    plt.plot(a, Ct_corr*(1-a), 'g--', label='$C_P$ Glauert correction')

    plt.ylabel(r'$C_T$ or $C_P$')
    plt.xlabel('a')
    plt.grid()
    plt.legend()
    plt.show()
