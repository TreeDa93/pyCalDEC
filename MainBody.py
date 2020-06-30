

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

print('Это работает и равно = '+ str(kq))
print('Это работает и равно = '+ str(kq1))


"We set the distribution of coils in space"
from CoilMatrix import *

Kf = CoilFunction(q, m, Q, Qkp).CoilT() # Winding matrix
print(Kf)



"""The part of model is intended to descreate domains
h - column of height size
mut , mun is tangential and normal components of permeability of each layers
 gamma is el. conductivity for each layer """

from Discretization import *


h, mut, mun, gamma, gammaotn = Discrete().discretka(gammaSE, kq)


"""In the part of the program we calculate matrices of resistance"""
from FindMatrixResistance import *

Rn0, Rt0, g0, GX0, gb = FindResist(gammaSE,mu0,deltaz,deltaSE,kq,h,mun,tz,mut,gammaotn)
Ra0, Rk0 = FindResist2(deltaz, muS, Hi, Hp1)
Rtbaz, gbas = FindRtabazAndGbas(tz,deltaz,Bi,mu0,gb)
Vse = FindVse(Qkp,Q,v,tz,omega)

Z = FindZ(Qz,Qkp,Q,Rn0,Rt0,GX0,omega,v,tz,Rk0,Ra0)

R, Z0, r0 = FindResistMassive(Qkp,Qz,Rt0,Q, Z)

a = FindaABS(Qz,R,Z0,r0,Z, Q, Qkp)

b = Findb(Qz,Z0,r0,Z,R, Q, Qkp)

F0s = FindF0s(IfA, IfB, IfC, Up, Rtbaz, Kf)

FXX = FindFXX(b,a,Qz,Rtbaz,F0s, Q, Qkp)

D = FindD(Qkp,Q)

BnX, BtX = FindBtBn(Qz,Bi,h,tz,D,FXX,Q, Qkp)


