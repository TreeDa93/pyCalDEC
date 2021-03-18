import math as ma
import numpy as np




Hy = 10 * 1e-3
Hp = 10 * 1e-3
dz = 10 * 1e-3
d_se = 10 * 1e-3
Bz = 10 * 1e-3
Bp = 10 * 1e-3
Hz = 10 * 1e-3
Q = 4
tz = Bp + Bz
Bi = tz * Q + Bz
marg = 0.1

phi_a = 0
phi_z = - ma.pi / 6
phi_b = - ma.pi*2/3
phi_x = - ma.pi
phi_c = ma.pi * 2 / 3
phi_y = ma.pi / 6

amplitude_current = 3e7
current_a = amplitude_current * np.exp(1j*phi_a)
current_b = amplitude_current * np.exp(1j*phi_b)
current_x = amplitude_current * np.exp(1j*phi_x)
current_c = amplitude_current * np.exp(1j*phi_c)
current_y = amplitude_current * np.exp(1j*phi_y)
current_z = amplitude_current * np.exp(1j*phi_z)

