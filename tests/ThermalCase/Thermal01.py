# Use ODEINT to solve the differential equations defined by the vector field
from scipy.integrate import odeint
import numpy as np
from scipy.integrate import odeint
T0 = 100

x = np.linspace(0, 0.001, 9)
T1 = np.array([T0]*9)
T1[-1]=1000
T = np.array([T0]*7)
# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 10.0
numpoints = 250
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

def expression(T, x, i):
    #
    f =  T[i-1]*(x[i+1]-x[i]) / (x[i+1]-x[i-1]) + T[i+1]*(x[i]-x[i-1]) / (x[i+1]-x[i-1])-T[i]
    return f


def temperatureField(T, t, x, expression, T1):
    """
        Defines the differential equations for the coupled spring-mass system.

        Arguments:
            w :  vector of the state variables:
                      w = [x1,y1,x2,y2]
            t :  time
            p :  vector of the parameters:
                      p = [m1,m2,k1,k2,L1,L2,b1,b2]
        """

    # Create f = (dT1dt, dT2dt, dTndt)
    #f = [expression(T1, x, i+1) for i in range(len(T))]
    f = [
        T1[0]*(x[2]-x[1]) / (x[2]-x[0]) + T[1]*(x[1]-x[0]) / (x[2]-x[0])-T[0],
        T[0]*(x[3]-x[2]) / (x[3]-x[1]) + T[2]*(x[2]-x[1]) / (x[3]-x[1])-T[1],
        T[1]*(x[4]-x[3]) / (x[4]-x[2]) + T[3]*(x[3]-x[2]) / (x[4]-x[2])-T[2],
        T[2]*(x[5]-x[4]) / (x[5]-x[3]) + T[4]*(x[4]-x[3]) / (x[5]-x[3])-T[3],
        T[3]*(x[6]-x[5]) / (x[6]-x[4]) + T[5]*(x[5]-x[4]) / (x[6]-x[4])-T[4],
        T[4]*(x[7]-x[6]) / (x[7]-x[5]) + T[6]*(x[6]-x[5]) / (x[7]-x[5])-T[5],
        T[5]*(x[8]-x[7]) / (x[8]-x[6]) + T1[8]*(x[7]-x[6]) / (x[8]-x[6])-T[6]
        ]

    return f





wsol = odeint(temperatureField, T, t, args=(x, expression, T1),
              atol=abserr, rtol=relerr)

with open('temperature.dat', 'w') as f:
    # Print & save the solution.
    splitter = '   '
    for t1, T in zip(t, wsol):
        f.write(f"{t1} {splitter} {T[0]} {splitter} {T[1]} {splitter} {T[2]} {splitter} {T[3]} {splitter} {T[4]} "
                f"{splitter} {T[5]} {splitter} {T[6]}\n")

from numpy import loadtxt
from pylab import figure, plot, xlabel, grid, legend, title, savefig
from matplotlib.font_manager import FontProperties


T = loadtxt('temperature.dat', unpack=True)

figure(1, figsize=(6, 4.5))

xlabel('t')
grid(True)
lw = 1

plot(T[0], T[2], 'b', linewidth=lw)
plot(T[0], T[7], 'r', linewidth=lw)




