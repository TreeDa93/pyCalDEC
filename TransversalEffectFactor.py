import numpy as np


class TransversalEffectFactor:
    """Это класс отвечающий за расчет коэффициентов учета поперечного краевого эффекта.
    Пока в тестовом режиме он называется Boltoncalculate. Он расчитывает только коэффициент Болтона
    """

    def __init__(self, epselon0 = 25, sk = 1, tau = 0.279, Bse = 0.7, Bi = 0.5):
        self.epselon0 = epselon0 #electromagnetic Q factor
        self.sk = sk # Slip
        self.tau = tau # # Pole pitch
        self.Bse = Bse # the width of secondary element
        self.Bi = Bi # the width of inductor

    def Bolton(self):
        """The function is intended for calculating values of Bollton's coefficient.
        H. Bolton, "Transverse edge effect in sheet-rotor induction motors",
        Proceedings of the Institution of Electrical Engineers, vol. 116. # 5, pp. 725, 1969
         Parameters:
         sk is the slip; tau is pole pitch; Bse is the width of secondary element;
         Bi is the width of inductor; epselon0 is the electromagnetic Q factor;
         """
        PP = np.sqrt(1 + 1j * self.epselon0 * self.sk)
        QQ = PP * np.pi / self.tau
        PR = np.pi * (self.Bse - self.Bi) / (2 * self.tau)
        SS = 1 / (1 + PP * np.tanh(0.5 * QQ * self.Bse) * np.tanh(PR))
        TT = SS * np.tanh(0.5 * QQ * self.Bse) / (0.5 * QQ * self.Bse)
        u = np.real(TT)
        v = np.imag(TT)
        kq = (1 - u - self.sk * self.epselon0 * v) / (1 - self.sk * self.epselon0 * \
                                           v + self.sk ** 2 * self.epselon0 ** 2 * u)
        return kq
