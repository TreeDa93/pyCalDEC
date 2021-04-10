"""Temperature evolution in a rod, computed by a ForwardEuler method."""

from numpy import linspace, zeros, linspace
import time
from Modules.ode_systemFE import ode_FE

def rhs(u, t):
    N = len(u) - 1
    rhs = zeros(N+1)
    rhs[0] = dsdt(t)
    for i in range(1, N):
        rhs[i] = (beta/dx**2)*(u[i+1] - 2*u[i] + u[i-1]) + g(x[i], t) # betta * uxx +g(x,t)
    i = N
    rhs[i] = (beta/dx**2)*(2*u[i-1] + 2*dx*dudx(t) - 2*u[i]) + g(x[N], t)
    return rhs

def dudx(t):
    return 0

def s(t):
    return 323

def dsdt(t):
    return 0

def g(x, t):
    return 0


L = 0.001
beta = 1
N = 10
x = linspace(0, L, N+1)
dx = x[1] - x[0]
u = zeros(N+1)

U_0 = zeros(N+1)
U_0[0] = s(0)
U_0[1:] = 283
dt = dx**2/(2*beta)
print('stability limit:', dt)
#dt = 0.00034375

t0 = time.clock()


u, t = ode_FE(rhs, U_0, dt, T=250*dt)
t1 = time.clock()
print('CPU time: %.1fs' % (t1 - t0))

from numpy import loadtxt
from pylab import figure, plot, xlabel, grid, legend, title, savefig
from matplotlib.font_manager import FontProperties

figure(1, figsize=(6, 4.5))

xlabel('t')
grid(True)
lw = 1

plot(u, 'black', linewidth=lw)

figure(2, figsize=(6, 4.5))
plot(u[-1], 'red', linewidth=lw)
plot(u[0], 'black', linewidth=lw)