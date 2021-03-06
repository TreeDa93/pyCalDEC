from data import *

# create geometry elements
from Modules.Geometry import Geometry

geometry1 = Geometry(label='LIM') # creat main node of geometry

geometry1.inductor(Hy, Hp, Bp, Bz, 6, intialCoord=(marg, marg),labelSetings=None, label='inductor1')
geometry1.coils(Hp, Bp, intial=(marg+Bz, Hy+marg), step=(0.1, 0), amount=1, label='coil1')
geometry1.coils(Hp, Bp, intial=(marg+Bz+tz, Hy+marg), step=[Bz, 0], amount=1, label='coil2')
geometry1.coils(Hp, Bp, intial=(marg+Bz+2*tz, Hy+marg), step=[Bz, 0], amount=1, label='coil3')
geometry1.coils(Hp, Bp, intial=(marg+Bz+3*tz, Hy+marg), step=[Bz, 0], amount=1, label='coil4')
geometry1.coils(Hp, Bp, intial=(marg+Bz+4*tz, Hy+marg), step=[Bz, 0], amount=1, label='coil5')
geometry1.coils(Hp, Bp, intial=(marg+Bz+5*tz, Hy+marg), step=[Bz, 0], amount=1, label='coil6')
geometry1.rect(height=d_se, width=0.1, intial=(marg, marg+Hy+Hp+dz), label='SE1')
geometry1.rect(height=d_se, width=0.1, intial=(marg, marg+Hy+Hp+dz+d_se), label='SE2')

#  create mesh

from Modules.DiscretizationNew import *

mesh1 = Mesh(label='MeshTest')
# разбили неактивные слои x and y
mesh1.discretRegion(startPoint=0, length=marg, axis='x', numberElem=4, lengthDisctr=None,label='margXLeft')
mesh1.discretRegion(startPoint=0, length=marg, axis='y', numberElem=4, lengthDisctr=None,label='margdownY')

"""Словарь с настройками
label - название узла сетки
length - длина региона
type - как делить область, numberElem - по количеству элементов, lengthDisctr - по длине элемента
pDiscret - значение параметра дискритизации если по количеству элементов делим, то это значение количества элементов
если делим по длине элемента то значение элемента
"""

body1 = {'label' : 'Tooth',
        'length' : Bz,
        'type' : 'numberElem',
        'pDiscret' : 4}

body2 = {'label' : 'Slot',
        'length' : Bp,
        'type' : 'numberElem',
        'pDiscret' : 4}

mesh1.discretMultiRegion(Listlengthes=[body1, body2],startPoint=mesh1.lastPoint(axis='x'),axis='x')
# copy nodes of Tooth and Slot
mesh1.copyMeshNodes(listNodesMesh=['Tooth', 'Slot'], CopyStep=tz, numberCopyTimes=5, axis='x')
# создадим последний зубец
mesh1.discretRegion(startPoint=mesh1.lastPoint(axis='x'), length=Bz, axis='x', numberElem=4,
                    lengthDisctr=None,label='LastTooth')
mesh1.discretRegion(startPoint=mesh1.lastPoint(axis='x'), length=marg, axis='x', numberElem=4,
                    lengthDisctr=None,label='margXRight')

mesh1.discretRegion(startPoint=mesh1.lastPoint(axis='y'), length=Hy, axis='y', numberElem=3,
                    lengthDisctr=None,label='yokeY')

mesh1.discretRegion(startPoint=mesh1.lastPoint(axis='y'), length=Hp, axis='y', numberElem=3,
                    lengthDisctr=None,label='slotY')

mesh1.discretRegion(startPoint=mesh1.lastPoint(axis='y'), length=dz, axis='y', numberElem=3,
                    lengthDisctr=None,label='gap')

mesh1.discretRegion(startPoint=mesh1.lastPoint(axis='y'), length=d_se, axis='y', numberElem=4,
                    lengthDisctr=None,label='se1')
mesh1.discretRegion(startPoint=mesh1.lastPoint(axis='y'), length=d_se, axis='y', numberElem=4,
                    lengthDisctr=None,label='se2')
mesh1.discretRegion(startPoint=mesh1.lastPoint(axis='y'), length=marg, axis='y', numberElem=4,
                    lengthDisctr=None,label='margYup')



mesh1.builtOneDgrid(axis='x')
mesh1.builtOneDgrid(axis='y')
mesh1.createMesh()

Bodies().defineBodies(geometry1.gNodes, mesh1)

#Определяем свойства материалов

from Modules.Materials import Material
#Создадим три основных шаблона
matAir = Material(conductivity=0, permeability=1, labelMat='air')
matCooper = Material(conductivity=5e7, permeability=1, labelMat='cooper')
matIron = Material(conductivity=1e4, permeability=1000, labelMat='iron')

#Определим свойства материала для каждого боди
Material().setMat(mesh1, material=matIron, body='inductor1')
Material().setMat(mesh1, material=matCooper, body='coil1')
Material().setMat(mesh1, material=matCooper, body='coil2')
Material().setMat(mesh1, material=matCooper, body='coil3')
Material().setMat(mesh1, material=matCooper, body='coil4')
Material().setMat(mesh1, material=matCooper, body='coil5')
Material().setMat(mesh1, material=matCooper, body='coil6')

Material().setMat(mesh1, material=matIron, body='SE1')
Material().setMat(mesh1, material=matIron, body='SE2')

# для области где тела не определены задаем свойства воздуха

Material().setMat(mesh1, material=matAir, body=None)

# Определим параметры магнитной системы


from Modules.PhysicMF4Cells import MagneticField

mf2 = MagneticField(mesh1, omega=314, label='mf2')
mf2.definiceCurrent(current=current_a, body='coil1')
mf2.definiceCurrent(current=current_z, body='coil2')
mf2.definiceCurrent(current=current_b, body='coil3')
mf2.definiceCurrent(current=current_x, body='coil4')
mf2.definiceCurrent(current=current_c, body='coil5')
mf2.definiceCurrent(current=current_y, body='coil6')

mf2.set_matrix_complex()

import Modules.Solvers as sol

solution = sol.solve_it_ling(mf2.weightMatrix, mf2.rightMatrix)

import Modules.PostProcessing4Cells as pp

postProc = pp.PostProcessing(solution, mf2)

postProc.dataCounturMatrix()
postProc.dataMagFluxMatrix()
from matplotlib.pyplot import  pcolormesh, pcolor, streamplot
pcolor(postProc.data['dataCounturMatrix'].real)
pcolor(postProc.data['magneticFluxDensity']['abs'].real)

omegat = ma.pi/2#ma.pi/4
multiplyer = np.exp(omegat*1j)
#z = postProc.data['magneticFlux']['abs']
#zfinal = (z * multiplyer).real
#pcolor(x, y, zfinal, shading='flat')


x=postProc.data['magneticFluxDensity']['axisX']
y=postProc.data['magneticFluxDensity']['axisY']
u=postProc.data['magneticFluxDensity']['x'].real
v=postProc.data['magneticFluxDensity']['y'].real

X, Y = np.meshgrid(x, y)
X = X.reshape(1 ,748)
Y = Y.reshape(1, 748)
U = u.reshape(1, 748)
V = v.reshape(1, 748)
from scipy import interpolate
#streamplot(X[0], Y[0], U[0], V[0])
f = interpolate.interp2d(X, Y, U)
xNew = np.linspace(x.min(), x.max(), 2000)
yNew = np.linspace(y.min(), y.max(), 2000)
zNew = f(xNew, yNew)
pcolor(xNew, yNew, zNew, shading='flat')




