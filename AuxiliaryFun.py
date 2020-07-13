
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


"""Функции для расчета времени выполнения фрагмента кода"""

import timeit as tit

fun1 = """
from FindMatrixResistance import TestFun
TestFun(10000)
"""
tt1 = tit.timeit(fun1, number=100) / 100
fun2 = """
from FindMatrixResistance import TestFun2
TestFun2(10000)
"""
tt2 = tit.timeit(fun2, number=100) / 100
fun3 = """
from FindMatrixResistance import TestFun3
TestFun3(10000)
"""
tt3 = tit.timeit(fun3, number=100) / 100
fun4 = """
from FindMatrixResistance import TestFun4
TestFun4(10000)
"""
tt4 = tit.timeit(fun4, number=100) / 100
print(tt1, tt2, tt3, tt4)



