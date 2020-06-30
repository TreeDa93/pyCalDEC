
import math as ma

def CalculateEpselon0(Bp1, deltaEkv, tz, tau, omega, mu0, deltaSE,
                      gammaSE, Kmu):
    """Recalculating factor to evaluate Karter's factor"""
    gammaZ1 = ((Bp1 / (3 * deltaEkv)) ** 2) / (5 + (Bp1 / (3 * deltaEkv)))
    Kdelta = tz / (tz - gammaZ1 * deltaEkv)  #
    KdeltaR = (2 * tau) * (ma.sinh(ma.pi * deltaEkv / (2 * tau))) \
              / (ma.pi * deltaEkv)
    deltaekv = deltaEkv * Kdelta * Kmu  # Equivalent gap Эквивалентный зазор

    epselon0 = omega * mu0 * deltaSE * gammaSE * tau ** 2 \
               / (ma.pi ** 2 * deltaekv)

    return epselon0



