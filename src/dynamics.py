import numpy as np

def evaluate_LMT(a_shoulder, b_shoulder, c_shoulder, a_elbow, b_elbow, c_elbow, l_base, l_multiplier, theta_shoulder, theta_elbow):
    # Returns the NORMALIZED muscle fiber length!
    # Note that in this model the tendons have been ignored and that to compute
    # the fiber length we do not need LMT-lt_slack

    l_full = a_shoulder*theta_shoulder + b_shoulder*np.sin(c_shoulder*theta_shoulder)/c_shoulder + a_elbow*theta_elbow + b_elbow*np.sin(c_elbow*theta_elbow)/c_elbow
    LMT = l_full*l_multiplier + l_base
    return LMT

def evaluate_LMT_vector(dM_coefficients, LMT_coefficients, theta_shoulder, theta_elbow):
    LMT_vector = evaluate_LMT(dM_coefficients[:,0], dM_coefficients[:,1], dM_coefficients[:,2], dM_coefficients[:,3], dM_coefficients[:,4], dM_coefficients[:,5], LMT_coefficients[:,0], LMT_coefficients[:,1], theta_shoulder, theta_elbow)
    return LMT_vector

def evaluate_VMT(a_shoulder, b_shoulder, c_shoulder, a_elbow, b_elbow, c_elbow, lM_multiplier, theta_shoulder, theta_elbow, dtheta_shoulder, dtheta_elbow):
    # Returns the fiber velocity NORMALIZED by the optimal fiber length. The
    # units are # optimal fiber lengths per second.
    nCoeff = a_shoulder.shape[0]
    v_full = a_shoulder*dtheta_shoulder + b_shoulder*np.cos(c_shoulder*theta_shoulder)*np.tile(dtheta_shoulder, (nCoeff, 1)) + a_elbow*dtheta_elbow + b_elbow*np.cos(c_elbow*theta_elbow)*np.tile(dtheta_elbow, (nCoeff, 1))
    VMT = lM_multiplier*v_full
    return VMT

def evaluate_VMT_vector(dM_coefficients, LMT_coefficients, theta_shoulder, theta_elbow, dtheta_shoulder, dtheta_elbow):
    VMT_vector = evaluate_VMT(dM_coefficients[:,0], dM_coefficients[:,1], dM_coefficients[:,2], dM_coefficients[:,3], dM_coefficients[:,4], dM_coefficients[:,5], LMT_coefficients[:,1], theta_shoulder, theta_elbow, dtheta_shoulder, dtheta_elbow)
    return VMT_vector

def getMuscleForce(q, qdot, auxdata):
    # Compute muscle fiber length and muscle fiber velocity 
    lMtilde = evaluate_LMT_vector(auxdata['dM_coefficients'], auxdata['LMT_coefficients'], q[0,:], q[1,:]) 
    vMtilde = evaluate_VMT_vector(auxdata['dM_coefficients'], auxdata['LMT_coefficients'], q[0,:], q[1,:], qdot[0,:], qdot[1,:]) 
    vMtilde_normalizedToMaxVelocity = vMtilde / auxdata['vMtilde_max'] 

    FMo = auxdata['Fiso']

    # Active muscle force-length characteristic
    Faparam = auxdata['Faparam']
    b11, b21, b31, b41, b12, b22, b32, b42 = Faparam

    b13 = 0.1
    b23 = 1
    b33 = 0.5*np.sqrt(0.5)
    b43 = 0
    num3 = lMtilde-b23
    den3 = b33+b43*lMtilde
    FMtilde3 = b13*np.exp(-0.5*num3**2/den3**2)

    num1 = lMtilde-b21
    den1 = b31+b41*lMtilde
    FMtilde1 = b11*np.exp(-0.5*num1**2/den1**2)

    num2 = lMtilde-b22
    den2 = b32+b42*lMtilde
    FMtilde2 = b12*np.exp(-0.5*num2**2/den2**2)

    FMltilde = FMtilde1+FMtilde2+FMtilde3

    # Active muscle force-velocity characteristic
    Fvparam = auxdata['Fvparam']
    e1, e2, e3, e4 = Fvparam

    FMvtilde = e1*np.log((e2*vMtilde_normalizedToMaxVelocity+e3)+np.sqrt((e2*vMtilde_normalizedToMaxVelocity+e3)**2+1))+e4 

    # Active muscle force
    Fce = FMltilde*FMvtilde

    # Passive muscle force-length characteristic
    Fpparam = auxdata['Fpparam']

    e0 = 0.6
    kpe = 4
    t5 = np.exp(kpe * (lMtilde - 1) / e0)
    Fpe = ((t5 - 1) - Fpparam[0]) / Fpparam[1]

    # Muscle force + damping 
    Fpv = auxdata['muscleDampingCoefficient']*vMtilde_normalizedToMaxVelocity 
    Fa = FMo*Fce
    Fp = FMo*(Fpe+Fpv)

    return Fa, Fp, lMtilde, vMtilde, FMltilde, FMvtilde, Fce, Fpe, Fpv

def evaluate_dM(a, b, c, theta):
    # Returns the muscle moment arm
    dM = a + b * np.cos(c * theta)
    return dM

def evaluate_dM_matrix(dM_coefficients, theta_shoulder, theta_elbow):
    dM_matrix = np.array([evaluate_dM(dM_coefficients[:,0], dM_coefficients[:,1], dM_coefficients[:,2], theta_shoulder), 
                          evaluate_dM(dM_coefficients[:,3], dM_coefficients[:,4], dM_coefficients[:,5], theta_elbow)]).T
    return dM_matrix

def TorqueForceRelation(Fm, q, auxdata):
    theta_shoulder = q[0,:]
    theta_elbow = q[1,:]
    dM_matrix = evaluate_dM_matrix(auxdata['dM_coefficients'], theta_shoulder, theta_elbow)

    T = dM_matrix @ Fm
    T = np.array([np.diag(T[:q.shape[1],:]).T, np.diag(T[q.shape[1]:,:]).T])
    return T

def armForwardDynamics(T, theta_elbow, qdot, T_EXT, auxdata):
    ddtheta = T_EXT.copy()
    dtheta_shoulder = qdot[0,:]
    dtheta_elbow = qdot[1,:]
    a1 = auxdata['I1'] + auxdata['I2'] + auxdata['m2'] * auxdata['l1'] ** 2
    a2 = auxdata['m2'] * auxdata['l1'] * auxdata['lc2']
    a3 = auxdata['I2']

    for i in range(T.shape[1]):
        M = np.array([[a1 + 2 * a2 * np.cos(theta_elbow[i]), a3 + a2 * np.cos(theta_elbow[i])],
                      [a3 + a2 * np.cos(theta_elbow[i]), a3]])

        C = a2 * np.sin(theta_elbow[i]) * np.array([-dtheta_elbow[i] * (2 * dtheta_shoulder[i] + dtheta_elbow[i]),
                                                    dtheta_shoulder[i] ** 2])

        # Joint friction matrix
        B = np.array([[0.05, 0.025],
                      [0.025, 0.05]])

        ddtheta[:, i] = np.linalg.solve(M, (T[:, i] + T_EXT[:, i] - C - B @ qdot[:, i]))

    return ddtheta

def forwardMusculoskeletalDynamics_motorNoise(X, u, T_EXT, wM, auxdata):
    q = X[:2]
    qdot = X[2:4]

    Fa, Fp, _, _, _, _, _, _, _ = getMuscleForce(q, qdot, auxdata)

    # u is vector of muscle activations
    u_bar = u + u * wM
    Fm = u_bar * Fa + Fp

    # wM is motor signal noise
    T = TorqueForceRelation(Fm, q, auxdata)

    F_forceField = auxdata['forceField'] * (auxdata['l1'] * np.cos(q[0,:]) + auxdata['l2'] * np.cos(q[0,:] + q[1,:]))
    T_forceField = -F_forceField * np.array([auxdata['l2'] * np.sin(q[0,:] + q[1,:]) + auxdata['l1'] * np.sin(q[0,:]), 
                                             auxdata['l2'] * np.sin(q[0,:] + q[1,:])])

    ddtheta = armForwardDynamics(T, q[1], qdot, T_EXT + T_forceField, auxdata)

    dX = np.concatenate([qdot, ddtheta])
    return dX