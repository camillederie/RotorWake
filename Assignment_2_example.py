import numpy as np


def solve_lifting_line_system_matrix_approach(rotor_wake_system, wind, Omega, rotorradius):
    # get controlpoints data structure
    controlpoints = rotor_wake_system["controlpoints"]
    # get horseshoe vortex rings data structure
    rings = rotor_wake_system["rings"]

    # initialize variables that we will use during the calculation
    velocity_induced = np.zeros(3)  # velocity induced by a horse vortex ring at a control point
    u, v, w = 0, 0, 0  # total velocity induced at one control point
    alpha = 0  # angle of attack
    GammaNew = np.zeros(len(controlpoints))  # new estimate of bound circulation
    Gamma = np.zeros(len(controlpoints))     # current solution of bound circulation
    MatrixU = []  # matrix of induction, for velocity component in x-direction
    MatrixV = []  # matrix of induction, for velocity component in y-direction
    MatrixW = []  # matrix of induction, for velocity component in z-direction

    # output variables
    a_temp = []  # output vector for axial induction
    aline_temp = []  # output vector for azimuthal induction
    r_R_temp = []  # output vector for radial position
    Fnorm_temp = []  # output vector for axial force
    Ftan_temp = []  # output vector for tangential force
    Gamma_temp = []  # output vector for circulation

    # the variables below are to setup the maximum number of iterations and convergence criteria
    Niterations = 1200
    errorlimit = 0.01
    error = 1.0
    ConvWeight = 0.3

    # initialize and calculate matrices for velocity induced by horseshoe vortex rings
    # two nested loops, each iterating over controlpoints and horseshoe vortex rings
    for icp in range(len(controlpoints)):
        MatrixU.append([])
        MatrixV.append([])
        MatrixW.append([])
        for jring in range(len(rings)):
            # set ring strength to unity, to calculate velocity induced by horseshoe vortex ring "jring"
            # at controlpoint "icp"
            rings[jring] = update_Gamma_sinle_ring(rings[jring], 1, 1)
            velocity_induced = velocity_induced_single_ring(rings[jring], controlpoints[icp]["coordinates"])
            # add component of velocity per unit strength of circulation to induction matrix
            MatrixU[icp].append(velocity_induced[0])
            MatrixV[icp].append(velocity_induced[1])
            MatrixW[icp].append(velocity_induced[2])

    # calculate solution through an iterative process
    for kiter in range(Niterations):
        Gamma = np.copy(GammaNew)  # update current bound circulation with new estimate

        # calculate velocity, circulation, and loads at the control points
        for icp in range(len(controlpoints)):
            # determine radial position of the control point
            radialposition = np.sqrt(np.dot(controlpoints[icp]["coordinates"], controlpoints[icp]["coordinates"]))
            u, v, w = 0, 0, 0  # initialize velocity
            # multiply icp line of Matrix with vector of circulation Gamma to calculate velocity at controlpoint
            for jring in range(len(rings)):
                u += MatrixU[icp][jring] * Gamma[jring]  # axial component of velocity
                v += MatrixV[icp][jring] * Gamma[jring]  # y-component of velocity
                w += MatrixW[icp][jring] * Gamma[jring]  # z-component of velocity

            # calculate total perceived velocity
            vrot = np.cross([-Omega, 0, 0], controlpoints[icp]["coordinates"])  # rotational velocity
            vel1 = [wind[0] + u + vrot[0], wind[1] + v + vrot[1], wind[2] + w + vrot[2]]  # total perceived velocity at section
            # calculate azimuthal and axial velocity
            azimdir = np.cross([-1 / radialposition, 0, 0], controlpoints[icp]["coordinates"])  # rotational direction
            vazim = np.dot(azimdir, vel1)  # azimuthal direction
            vaxial = np.dot([1, 0, 0], vel1)  # axial velocity
            # calculate loads using blade element theory
            temploads = loadBladeElement(vaxial, vazim, radialposition / rotorradius)
            # new point of new estimate of circulation for the blade section
            GammaNew[icp] = temploads[2]
            # update output vector
            a_temp.append(-(u + vrot[0]) / wind[0])
            aline_temp.append(vazim / (radialposition * Omega) - 1)
            r_R_temp.append(radialposition / rotorradius)
            Fnorm_temp.append(temploads[0])
            Ftan_temp.append(temploads[1])
            Gamma_temp.append(temploads[2])

        # check convergence of solution
        refererror = max(abs(GammaNew))
        refererror = max(refererror, 0.001)  # define scale of bound circulation
        error = max(abs(np.subtract(GammaNew, Gamma)))  # difference between iterations
        error = error / refererror  # relative error
        if error < errorlimit:
            # if error smaller than limit, stop iteration cycle
            break

        # set new estimate of bound circulation
        GammaNew = (1 - ConvWeight) * Gamma + ConvWeight * GammaNew

    # output results of converged solution
    return {"a": a_temp, "aline": aline_temp, "r_R": r_R_temp, "Fnorm": Fnorm_temp, "Ftan": Ftan_temp, "Gamma": Gamma_temp}

