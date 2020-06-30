

import math as ma
import numpy as np
"""The file contain general data and information about considering model"""


"""Geometry parameters of the model"""

Hi = 0.13       # Height of tooth of inductor [m]
Bi = 1          # Width of inductor [m]
Bse = 1         # Width of secondary element [m]
Hp1 = 0.07      # Height of slots [m]
Bp1 = 0.06      # Width of slots [m]
Bz1 = 0.033     # Width of tooth [m]
deltaz = 0.01   # The gap thickness between inductor and secondary element [m]
deltaSE = 8 * 10 ** (-3)    # Thickness of conducting layer of secondary element [m]
deltaEkv = deltaz + deltaSE  # Equivalent non-magnetic gap  [m]

"""Discrete parameters of model"""

Q = 24      # The number of slots
Qkp = 24    # The number of marginal slots
Qz = 14     # The number of layers of model along thickness

"""Coil parameters"""

If = 900   # phase current
IfA = np.sqrt(2) * If * np.exp(1j * 0)  # phase current A
IfB = np.sqrt(2) * If * np.exp(1j * (-2 / 3) * np.pi)   # Phase current B
IfC = np.sqrt(2) * If * np.exp(1j * (2 / 3) * np.pi)    # Phase current C
F1 = 28     # Frequency
m = 3   # The number of phase
q = 1           # The number of slots per phase and pole
p = Q / (2 * q * m)     # The number of pole pair
Kz = 0.5                # Winding fill factor of slots
Up = 29                 # The number of turns
omega = 2 * ma.pi * F1  # Angular frequency угловая частота
Spr = Hp1 * Bp1 * Kz / Up  # Cross section of coil conductor of inductor
Up = 29     # The number of turns in coil (slot)
L1 = 0.395      # Length of endwinding (winding faces)
Kz = 1      # Wire fill factor of slots
betta = 1   # Winding pitch shortening factor
gamma1 = 0.46 * 10 ** 6     # Electrical conductivity of coils of inductor [S/m]
Kmu = 1     # Saturation factor
ky = ma.sin(ma.pi * betta / 2)      # Shorting factor
kp = (ma.sin(ma.pi / 2 * m)) / (q * ma.sin(ma.pi / 2 * m * q))      # Distributing factor
Kw = ky * kp        # Winding coefficient of inductor
tz = (Bp1 + Bz1)    # Slot-tooth division
tau = tz * m * q    # Pole pitch


"Some physical properties"
gammaSE = 33 * 10 ** 6     # Electrical conductivity of conducting part of secondary element
muS = 1000      # Relative permeability
mu0 = 4 * ma.pi * 10 ** (-7)    # Absolute permeability


kr = 0.4        # Relaxation coefficient for calculating iteration process with saturation

"""Descriptiom of motion of secondary element"""
sk = 1          # Slip of magnetic field
Vsinhr = 2 * tau * F1   # Sinhronus velocity
v = (1 - sk) * Vsinhr   # Operating velocity
Flag = 'sk'

