import numpy as np
from parameters import initializeModelParameters
from dynamics import forwardMusculoskeletalDynamics_motorNoise as dynamics_function

auxdata = initializeModelParameters()

f = lambda x, u, w: dynamics_function(x, u, 0, w, auxdata)

x_init = np.array([0, 0, 0, 0]).reshape(-1, 1)
u = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

print(f(x_init, u, 0))