#from Data import *
import numpy as np
"""Kfk = np.zeros((3, Q))
#RM = np.array([[1, 0, 0, -1, 0, 0], [0, 0, 1, 0, 0, -1], [0, -1, 0, 0, 1, 0]])
#RMK1 = np.array([[0.5, 0, 0, -1, 0, 0], [0, 0, 0.5, 0, 0, -1], [0, -0.5, 0, 0, 1, 0]])
#RMK2 = np.array([[1, 0, 0, -0.5, 0, 0], [0, 0, 1, 0, 0, -0.5], [0, -1, 0, 0, 0.5, 0]])
RM = np.hstack((Kfk, RMK1, RM, RM, RMK2, Kfk))"""

""" Тестовы Class для создания обмоточной матрицы ***пока что в разработке"""
class CoilFunction:
    def __init__(self, q=1, m=3, Q=24, Qkp=24):
        self.q = q
        self.m = m
        self.Q = Q
        self.Qkp = Qkp
    def Coil(self, q, m, Q, Qkp):
        Kf = np.zeros((m, Q))
        Qall = Q + 2 * Qkp
        k = 0
        while k < m:
            for j in range(0, q):
                for i in range(0, Q, 2 * q * m):
                    if (i < Qkp):
                        Kf[k][i + j + k] = 0.5
                    elif (i > Qkp) and (i < Qall - Qkp):
                        Kf[k][i + j + k] = 1
                    else:
                        Kf[k][i + j + k] = 0.5
            k += 1
        return Kf

    def CoilT(self):
        Kfk = np.zeros((self.m, self.Q))
        RM = np.array([[1, 0, 0, -1, 0, 0], [0, 0, 1, 0, 0, -1], [0, -1, 0, 0, 1, 0]])
        RMK1 = np.array([[0.5, 0, 0, -1, 0, 0], [0, 0, 0.5, 0, 0, -1], [0, -0.5, 0, 0, 1, 0]])
        RMK2 = np.array([[1, 0, 0, -0.5, 0, 0], [0, 0, 1, 0, 0, -0.5], [0, -1, 0, 0, 0.5, 0]])
        RM = np.hstack((Kfk, RMK1, RM, RM, RMK2, Kfk))
        return RM

