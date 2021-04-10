import numpy as np

import numpy as np
from scipy.integrate import odeint
T0 = 20
x = np.linspace(0, 0.01, 9) # 0 , 0.001 ... 9 точек
Tbc = [20, 1000]
T = np.array([T0]*7)

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 50
numpoints = 2000
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]




def spatialDerivativeties(T, x, i):
    f =  T[i-1]*(x[i+1]-x[i]) / (x[i+1]-x[i-1]) + T[i+1]*(x[i]-x[i-1]) / (x[i+1]-x[i-1])-T[i]
    return f

def spatialDerivativetiesBC(T, x, i, bc=None, Tbc=T0):
    if bc == 'left':
        f = Tbc*(x[i+1]-x[i]) / (x[i+1] - x[i-1]) + T[i+1] * (x[i] - x[i-1]) / (
                x[i + 1] - x[i - 1]) - T[i]
    elif bc == 'right':
        f = T[i-1] * (x[i+1] - x[i]) / (x[i+1] - x[i-1]) + Tbc * (x[i]-x[i-1]) / (x[i + 1] - x[i-1])-T[i]
    return f

def fun(T, t, x, Tbc):
    temperature = np.array([])
    for i in range(len(T)):
        if i == 0:
            expression = spatialDerivativetiesBC(T, x, i, bc='left', Tbc=Tbc[0])
            temperature = np.append(temperature, expression)
        elif i == len(x[1:-1])-1:
            expression = spatialDerivativetiesBC(T, x, i, bc='right', Tbc=Tbc[1])
            temperature = np.append(temperature, expression)
        else:
            expression = spatialDerivativeties(T, x, i)
            temperature = np.append(temperature, expression)
    return temperature


wsol = odeint(fun, T, t, args=(x, Tbc), atol=abserr, rtol=relerr)




with open('temperature2.dat', 'w') as f:
    # Print & save the solution.
    splitter = '   '
    for t1, T in zip(t, wsol):
        f.write(f"{t1} {splitter} {T[0]} {splitter} {T[1]} {splitter} {T[2]} {splitter} {T[3]} {splitter} {T[4]} "
                f"{splitter} {T[5]} {splitter} {T[6]}\n")

from numpy import loadtxt
from pylab import figure, plot, xlabel, grid, legend, title, savefig
from matplotlib.font_manager import FontProperties


Tsol = loadtxt('temperature2.dat', unpack=True)

figure(1, figsize=(6, 4.5))

xlabel('t')
grid(True)
lw = 1

plot(Tsol[0], Tsol[1], 'black', linewidth=lw)
plot(Tsol[0], Tsol[2], 'b', linewidth=lw)
plot(Tsol[0], Tsol[3], 'g', linewidth=lw)
plot(Tsol[0], Tsol[4], 'grey', linewidth=lw)
plot(Tsol[0], Tsol[5], 'brown', linewidth=lw)
plot(Tsol[0], Tsol[6], 'yellow', linewidth=lw)
plot(Tsol[0], Tsol[7], 'r', linewidth=lw)

figure(2, figsize=(6, 4.5))
plot(Tsol[:,-1], 'black', linewidth=lw)
