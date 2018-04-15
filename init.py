import numpy as np

# define parameters
m1 = 1
m2 = 1
l1 = 1
l2 = 1
I1 = 1
I2 = 1

# unitless parameters
M0 = m1/m2 + 1
L0 = l1/l2
I01 = I1/(l1*l2*m2)
I02 = I2/(l1*l2*m2)

# initial condition
init_cond = np.array([0, 0, 0, 0])
