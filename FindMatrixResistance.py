

"""Here we are creating array of resitance"""

import numpy as np

def FindResist(gammaSE,mu0,deltaz,deltaSE,kq,h,mun,tz,mut,gammaotn):
    """The function is intended to calcualte y-component (normal) Rn0, x-component (tangential) Rt0
    realative conductivity g0, absolute conductivity GX0, and gb is the basis f conductivity.

    The name of functions are

    Rn0 is the matrix of magnetic resistance of  normal component to inductor plane or y-component;
     Rt0 is the matrix of magnetic resistance of  tangential component to inductor plane or y-component;
     gb is the value of basis electrical resistance take in account transversal effect (electrical conductivity of
     conducting layer of secondary element)
     g0 is the column of relative electrical conductivity for each layer of the model;
      GX0 is the column of absolute electrical conductivity for each layer of the model

      The variables are
      gammaSE is the electrical conductivity of conducting part of secondary element
      mu0 is the absolute permeability
      deltaz is the gap thickness between inductor and secondary element [m]
      deltaSE is the Thickness of conducting layer of secondary element [m]
      kq is the coefficient taken in account transversal edge effect
      h is the column of sizes of layers height
      mun is the normal component of relative permeability for each layer
      tz is the slot-tooth division
      mut is the tangential component of relative permeability for each laye
      gammaotn is the relative electrical conductivity for each layer of the model to the base vale.
      The base value is the electrical conductivity value of conducting layer of secondary element """

    gb=gammaSE * 0.5* mu0 * deltaz * deltaSE * kq # базисная электропроводность на 0.5deltaSE
    Rn0 = deltaz*h / (mun*tz**2) # нормальное магнитное сопротивление
    Rt0 = deltaz / (h * mut)
    g0 = gammaotn * h / deltaSE
    GX0 = g0 * gb
    flag = 'No print'
    if flag == 'print':
        print(gb, Rn0, Rt0, g0, GX0)
    return Rn0, Rt0, g0, GX0, gb

def FindResist2(deltaz, muS, Hi, Hp1):
    """The function FindResis2 is intended to calculate magnetic resistance of
    inductor teeth Ra0, and Rk0 - ****

    Parameers:
    deltaz is the gap thickness between inductor and secondary element [m]
    muS is the relative permeability of magnetic circuit (inductor);
    Hi is the height of inductor; Hp1 is the height of the slots
    """
    Ra0 = deltaz / (muS * (Hi-Hp1))
    Rk0 = deltaz / Hi
    flag = 'No print'
    if flag == 'print':
        print(Ra0, Rk0)
    return Ra0, Rk0


def FindRtabazAndGbas(tz,deltaz,Bi,mu0,gb):
    """The function calculate basis magnetic resistance Rtbaz
    and  magnetic resistance for each layer of model"""
    Rtbaz = tz / (mu0 * deltaz * Bi)
    gbas = gb * Rtbaz
    flag = 'No print'
    if flag == 'print':
        print(Rtbaz, gbas)
    return Rtbaz, gbas


def FindVse(Qkp,Q,v,tz,omega):
    Vse = np.zeros((Q+2*Qkp, Q+2*Qkp), dtype=complex)
    for x in range(Q+2*Qkp):
        if x > 1:
            Vse[x, x - 1] = - v/(2 * tz)
        else:
            Vse[x, x + 1] = v / (2 * tz)
        Vse[x,x] = -1j * omega
    flag = 'No print'
    if flag == 'print':
        print(Vse)
    return Vse

def FindLc(Qkp, Q):
    lc = 0
    Lc = np.zeros((2*Qkp+Q))
    for x in range(2*Qkp+Q):
        Lc[x, x] = lc
    flag = 'No print'
    if flag == 'print':
        print(Lc)
    return Lc

def FindZ(Qz,Qkp,Q,Rn0,Rt0,GX0,omega,v,tz,Rk0,Ra0):
    Z = np.zeros((Q+2*Qkp, Q+2*Qkp, Qz), dtype=complex)
    for k in range(Qz):
        for i in range(Q+2*Qkp):
            if k == 0:
                if i < Qkp:
                    Z[i, i, k] =  Rn0[k] + Rk0 + Rt0[k] + 1j*omega * GX0[k]
                elif i > Qkp & i<= (Qkp+Q):
                    Z[i,i,k] = Rn0[k] + Rt0[k] + Ra0 + 1j * omega * GX0[k]
                elif i > (Qkp+Q):
                    Z[i,i,k] = Rn0[k] + Rk0 + Rt0[k] + 1j*omega * GX0[k]
                elif i>0:
                    Z[i,i-1,k] = -0.5 * Rn0[k]- GX0[k]*v/(2 * tz)
                elif i<(2*Qkp+Q):
                    Z[i,i+1] = -0.5*Rn0[k] +GX0[k] * v/(2*tz)
            else:
                if i > 0:
                    Z[i,i,k] = Rn0[k]+Rn0[k-1]+Rt0[k]+Rt0[k-1]+1j*(GX0[k]+GX0[k-1])*omega
                    Z[i,i-1,k] = -0.5*Rn0[k]-0.5*Rn0[k-1]-v/(2*tz)*(GX0[k]+GX0[k-1])
                elif i < 2*Qkp+Q:
                    Z[i,i+1,k] = -0.5*Rn0[k]-0.5*Rn0[k-1]+v/(2*tz)*(GX0[k]+GX0[k-1])
    flag = 'No print'
    if flag == 'print':
        print(Z)
    return Z

def FindResistMassive(Qkp,Qz,Rt0,Q, Z):
    R = np.zeros((2 * Qkp + Q, 2 * Qkp + Q, Qz), dtype=complex)
    Z0 = np.zeros((2 * Qkp + Q, 2 * Qkp + Q, Qz), dtype=complex)
    r0  = np.zeros((2 * Qkp + Q, 2 * Qkp + Q, Qz), dtype=complex)
    for k in range(Qz):
        for i in range(2*Qkp+Q):
            R[i,i,k] = Rt0[k]
            Z0[i,i,k] = Rt0[k] ** -1 * Z[i,i,k]
            r0[i,i,k] = Rt0[k] ** -1 * R[i,i,k-1]
    flag = 'No print'
    if flag == 'print':
        print(R, Z0, r0)
    return R, Z0, r0

def FindaABS(Qz,R,Z0,r0,Z, Q, Qkp):
    a = np.zeros((Q+2*Qkp, Q+2*Qkp, Qz+1), dtype=complex)
    for k in range(1, Qz+1):
        for i in range(Q+2*Qkp):
            if k == 1:
                a[i,i, k] = R[i,i,k-1] ** -1
            elif k == 2:
                a[:, :, k] = Z0[:,:,k-1]*a[:,:,k-1]
            elif k == 14:
                a[:, :, Qz] = Z[:, :, Qz-1] * a[:, :, Qz-1] - R[:, :, Qz - 2] * a[:, :, Qz - 2]
            else:
                a[:,:,k] = Z0[:,:,k-1]*a[:,:,k-1]-r0[:,:,k-1]*a[:,:,k-2]
    flag = 'No print'
    if flag == 'print':
        print(a)
    return a

def Findb(Qz,Z0,r0,Z,R, Q, Qkp):
    b = np.zeros((Q + 2 * Qkp, Q + 2 * Qkp, Qz + 1), dtype=complex)
    for k in range(1, Qz + 1):
        for i in range(Q + 2 * Qkp):
            if k == 1:
                b[:, :, k] = Z0[:, :, k-1]
            elif k == 2:
                b[:, :, k] = Z0[:, :, k - 1] * b[:, :, k - 1]-r0[:,:,k-1]
            elif k == 14:
                b[:, :, Qz] = Z[:, :, Qz - 1] * b[:, :, Qz - 1] - R[:, :, Qz - 2] * b[:, :, Qz - 2]
            else:
                b[:, :, k] = Z0[:, :, k - 1] * b[:, :, k - 1] - r0[:, :, k - 1] * b[:, :, k - 2]
    flag = 'No print'
    if flag == 'print':
        print(b)
    return b

def FindFXX(b,a,Qz,Rtbaz,F0s, Q, Qkp):
    FXX = np.zeros((Q+2*Qkp,Q+2*Qkp, Qz), dtype=complex)
    for k in range(Qz):
        for i in range(Q+2*Qkp):
            if k == 0:
                FXX [i,i,k] = -(b[i,i,k+1]*Rtbaz**-2) * a[i,i,k+1]*F0s[i]
            else:
                FXX[i,i,k]=a[i,i,k]*F0s[i] + b[i,i,k]*FXX[i,i,0]
    flag = 'No print'
    if flag == 'print':
        print(FXX)
    return FXX


def FindF0s(IfA, IfB, IfC, Up, Rtbaz, Kf):
    Ifp = np.array([IfA, IfB, IfC])
    Ifm = np.dot(Kf.T, Ifp)
    Fs = Up * Ifm
    F0s = Fs/Rtbaz
    flag = 'No print'
    if flag == 'print':
        print(F0s)
    return F0s

def FindD(Qkp,Q):
    D = np.zeros((Q+2*Qkp,Q+2*Qkp))
    for i in range(Q+2*Qkp):
        if i == 0:
            D[i, i + 1] = 1
        elif i == 2*Qkp+Q-1:
            D[i, i - 1] = -1
        else:
            D[i, i + 1] = 1
            D[i, i - 1] = -1
    flag = 'No print'
    if flag == 'print':
        print(D)
    return D


def FindBtBn(Qz,Bi,h,tz,D,FXX,Q, Qkp):
    BnX = np.zeros((2*Qkp+Q, 2*Qkp+Q, Qz), dtype=complex)
    BtX = np.zeros((2 * Qkp + Q, 2 * Qkp + Q, Qz), dtype=complex)
    for k in range(Qz):
        BnX[:,:,k] = 1/(2 * tz * Bi) * D[:,:] * FXX[:,:,k]
    for k in range(Qz-1):
        BtX[:,:, k]=(FXX[:,:, k+1] - FXX[:,:, k]) / (Bi * h[k])
    flag = 'No print'
    if flag == 'print':
        print(BnX, BtX)
    return BnX, BtX

"""def CalculateEM(omega, Up, Kf, b, Qz, Rtbaz, a):
    Ifpm = np.eye((3,3))
    EM = -1j * omega * Up * Kf.T * ((b([:,:,Qz+1)] / Rtbaz ** -2) * a[:,:,Qz+1] / Rtbaz)*Up*Kf*Ifpm"""









