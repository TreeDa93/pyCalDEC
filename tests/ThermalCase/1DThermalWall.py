from Modules.Geometry import Geometry
from scipy.integrate import odeint
h1=1
w1=0.001

geometry1 = Geometry(label='ThermalWall') # creat main node of geometry
geometry1.rect(h1, w1, intial=[0, 0], label='wall')

from Modules.DiscretizationNew import *

mesh1 = Mesh(label='Mesh1')
# разбили неактивные слои x and y
mesh1.discretRegion(startPoint=0, length=w1, axis='x', numberElem=11, lengthDisctr=None,label='wallX')
mesh1.discretRegion(startPoint=0, length=h1, axis='y', numberElem=2, lengthDisctr=None,label='wallX')
mesh1.builtOneDgrid(axis='x')
mesh1.builtOneDgrid(axis='y')
mesh1.createMesh()
mesh1.create1DMesh()



for i in mesh1.mesh1D[1:-1]:
    i.addTypeCell(cellType='General')
    i.addTypeApproximation()
    i.setBCtype(bc=None)

mesh1.mesh1D[0].addTypeApproximation(addTypeApproximation='right')
mesh1.mesh1D[-1].addTypeApproximation(addTypeApproximation='left')
mesh1.mesh1D[0].setBCtype(bc='bc1')
mesh1.mesh1D[-1].setBCtype(bc='bc2')
for i in mesh1.mesh1D[:]:
    i.setTemperatureInitial(20)




from Modules.PhysicHT import HeatTransfer

ht1 = HeatTransfer(mesh1, label='ht1')
T = ht1.initTemperature(0)
fun = ht1.temperatureFieldPDE
fun2 = ht1.temperatureFieldPDE2
fun3 = ht1.temperatureFieldPDE3
fun4 = ht1.temperatureFieldPDE4
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 100*5e-09
numpoints = 2500
t = ht1.timeRange(stoptime, numpoints)
wsol = odeint(fun, T, t, atol=abserr, rtol=relerr)
wsol2 = odeint(fun4, T, t, atol=abserr, rtol=relerr)

from numpy import loadtxt
from pylab import figure, plot, xlabel, grid, legend, title, savefig
from matplotlib.font_manager import FontProperties

figure(1, figsize=(6, 4.5))
xlabel('t')
grid(True)
lw = 1
plot(wsol, linewidth=lw)
plot(wsol2, linewidth=3)
figure(2, figsize=(6, 4.5))
plot(wsol[0], 'black', linewidth=lw)
plot(wsol[-1], 'red', linewidth=lw)
plot(wsol2[0], 'black', linewidth=lw)
plot(wsol2[-1], 'red', linewidth=lw)