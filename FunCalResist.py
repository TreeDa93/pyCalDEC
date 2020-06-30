"""Here we describe general function to calculate resistance of
magnetic system"""

"""Calculating resistance of inductor system"""
import numpy as np
def CalculateR1(Bi, L1, Up, gamma1, Spr):
    """Function is intended to calculate value of resistance
    of inductor winding
    R1 is value of the resistance [Ohm]
    Bi is the Width of inductor [m]
    L1 is the Length of endwinding (winding faces)
    Up is the The number of turns in coil (slot)
    gamma1 is the Electrical conductivity of coils of inductor [S/m]
    Spr is the Cross section of coil conductor of inductor
    """
    R1 = (2 * Bi + 2 * L1) * Up / (gamma1 * Spr)  # Resistance of winding
    return R1

def CalculateRf1(Q, m, Bi, L1, Up, gamma1, Spr):
    """Function is intended to calculate value of equivalent resistance
       of a phase of inductor
       R1 is the value of winding resistance consisting in the phase [Ohm]
       Rf1 is the value of the resistance [Ohm]
       Q is the  number of slots
       m is the  number of phase
       """
    R1 = (2 * Bi + 2 * L1) * Up / (gamma1 * Spr)  # Resistance of winding
    Rf1 = R1 * Q / (m * 2)  # Equivalent resistance of phase
    Rf = Rf1 * np.eye(3)
    return Rf1, Rf


def CulculateLf1(Hp1, Bp1, q, Bi, L1, betta, tau, F1, Up, Q, m,
                       p, omega):
    """The function is intended to calculate values of reactance
    and inductance of a phase of inductor.
    The function outputs Xdelta1 and Lf1.
     Lambdap1 is specific factor calculating by Kopylov
     labmdal1 is specific factor calculating by Kopylov
     Xdelta1 is the value of the reactance [H]
     Lf1 is value of the inductance, F1 is the supply frequency [Ohm],
     Hp1 is the Height of slots, Bp1 is the width of slots,
     q is the number of slots per phase and pole
     L1 is the Length of endwinding (winding faces)
      betta is the Winding pitch shortening factor
       tau is the Pole pitch, p is the number of poles,
       Up is the The number of turns in coil (slot)
       m is the number of phase, Bi is the Width of inductor [m],
       """
    lambdap1 = 1 / 3 * Hp1 / Bp1
    lambdal1 = (0.34 * q / Bi) * (L1 - 0.64 * betta * tau)
    Xdelta1 = 15.8 * F1 / 100 * (Up * Q / (200 * m)) ** 2 * (Bi / (p * q)) \
              * (lambdap1 + lambdal1)
    Lf1 = Xdelta1 / omega
    Lf = Lf1 *  np.eye(3)
    return Xdelta1, Lf1, Lf

