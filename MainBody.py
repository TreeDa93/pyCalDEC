

from Data import *  # Import the general data of the model

# Import file with function to calculate resistance of inductor
from FunCalResist import *
# The file consist of auxiliary function
from AuxiliaryFun import CalculateEpselon0
# The file consist function to calculate Bolton coefficient
from TransversalEffectFactor import TransversalEffectFactor #


# file below to be intended to recalculae and set values of velocityy and slip
#from Velocity import VelocityData4

"Calculating of inductor resistance"


Rf1, Rf = CalculateRf1(Q, m, Bi, L1, Up, gamma1, Spr)
Xdelta, Lf1, Lf = CulculateLf1(Hp1, Bp1, q, Bi, L1, betta, tau, F1, Up, Q, m,
                       p, omega)

# Calculate electromagnetic Q factor
epselon0 = CalculateEpselon0(Bp1, deltaEkv, tz, tau, omega, mu0, deltaSE,
                      gammaSE, Kmu)



# The file consist function to calculate Bolton coefficient
from TransversalEffectFactor import TransversalEffectFactor
"Calculating Bollton's factor"

kq = TransversalEffectFactor(epselon0, sk, tau, Bse, Bi)
kq1 = TransversalEffectFactor()
kq = kq.Bolton()
kq1 = kq1.Bolton()

print('Это работает и равно = ' + str(kq))
print('Это работает и равно = ' + str(kq1))


"We set the distribution of coils in space"





"""The part of model is intended to descreate domains
h - column of height size
mut , mun is tangential and normal components of permeability of each layers
 gamma is el. conductivity for each layer """






"""

from FindMatrixResistance import TestDiaFun, TestDiaFun2

import numpy as np
test1 = TestFun2()
test2 = TestFun5()
x = np.diag(test2)


from scipy.sparse import dia_matrix
x = dia_matrix(([rs1,s], [0,1]), shape=(10,10)).toarray()

x = dia_matrix((rs,rs1, rs2, [0,1,-1], shape(10,10)).toarray()

x = dia_matrix(([rs,rs1,rs2, [0,1,-1], shape=(10,10)).toarray()

rs, rs2, rs3, x = TestDiaFun2()"
"""
