import numpy as np
import scipy.io as sio

def initializeModelParameters():
    auxdata = {}
    auxdata['nDOF'] = 2
    auxdata['nMuscles'] = 6
    auxdata['nMotorNoises'] = 6  # Dimension of the MotorNoises matrix

    auxdata['forceField'] = 0

    # Parameters of skeletal arm model
    auxdata['m1'] = 1.4
    auxdata['m2'] = 1
    auxdata['l1'] = 0.3
    auxdata['l2'] = 0.33
    auxdata['lc1'] = 0.11
    auxdata['lc2'] = 0.16
    auxdata['I1'] = 0.025
    auxdata['I2'] = 0.045

    # Load and store coefficients for muscle tendon lengths and moment arms
    auxdata['dM_coefficients'] = sio.loadmat('Muscle_LMT_dM/dM_coefficients.mat')['dM_coefficients']  # First shoulder muscles, then elbow muscles
    auxdata['LMT_coefficients'] = sio.loadmat('Muscle_LMT_dM/LMT_coefficients.mat')['LMT_coefficients']  # First shoulder muscles, then elbow muscles

    # Muscle parameters (Hill-type)
    auxdata['Fiso'] = 31.8 * np.array([18, 14, 22, 12, 5, 10])
    auxdata['vMtilde_max'] = np.array([10, 10, 10, 10, 10, 10])
    auxdata['muscleDampingCoefficient'] = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])

    # Parameters of active muscle force-velocity characteristic
    ActiveFVParameters = sio.loadmat('MuscleModel/ActiveFVParameters.mat')['ActiveFVParameters']
    Fvparam = np.zeros(4)
    Fvparam[0] = 1.475 * ActiveFVParameters[0, 0]
    Fvparam[1] = 0.25 * ActiveFVParameters[0, 1]
    Fvparam[2] = ActiveFVParameters[0, 2] + 0.75
    Fvparam[3] = ActiveFVParameters[0, 3] - 0.027
    auxdata['Fvparam'] = Fvparam

    # Parameters of active muscle force-length characteristic
    auxdata['Faparam'] = sio.loadmat('MuscleModel/Faparam.mat')['Faparam']

    # Parameters of passive muscle force-length characteristic
    e0 = 0.6
    kpe = 4
    t50 = np.exp(kpe * (0.2 - 1) / e0)
    pp1 = (t50 - 1)
    t7 = np.exp(kpe)
    pp2 = (t7 - 1)
    auxdata['Fpparam'] = np.array([pp1, pp2])

    return auxdata