#from Data import *


def Bolton(epselon0, sk, tau, Bse, Bi):
    import numpy as np
    PP = np.sqrt(1 + 1j * epselon0 * sk)
    QQ = PP * pi / tau
    PR = pi * (Bse - Bi) / (2 * tau)
    SS = 1 / (1 + PP * np.tanh(0.5 * QQ * Bse) * np.tanh(PR))
    TT = SS * np.tanh(0.5 * QQ * Bse) / (0.5 * QQ * Bse)
    u = np.real(TT)
    v = np.imag(TT)
    kq = (1 - u - sk * epselon0 * v) / (1 - sk * epselon0 * v + sk ** 2 * epselon0 ** 2 * u)
    return kq
#test = Bolton(epselon0, sk, tau, Bse, Bi)
#print(test)