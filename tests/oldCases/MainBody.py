from Preprocessing import *
from Modules.PhysicsMF import MagneticField
from Modules.PostProcessing import PostProcessing
import Modules.Solvers as sol

import time



def main():
    pass


# create properties of materials

from Modules.Materials import Material
matAir = Material(conductivity=0, permeability=1, labelMat='air')
matCooper = Material(conductivity=5e7, permeability=1, labelMat='cooper')
matIron = Material(conductivity=1e4, permeability=1000, labelMat='iron')

# create geometry elements
from Modules.Geometry import Geometry
from Modules.DiscretizationNew import *

a = Geometry(label='test') # creat main node of geometry
#a.rect(height=Hy, width=Bz, intialX=0, intialY=0, label='rect1')
a.coils(Hp, Bp, start=Coord(marg+Bz, Hy+marg), step=[Bz, 0], amount=1, label='coil1')
a.coils(Hp, Bp, start=Coord(marg+Bz+tz, Hy+marg), step=[Bz, 0], amount=1, label='coil2')
a.coils(Hp, Bp, start=Coord(marg+Bz+2*tz, Hy+marg), step=[Bz, 0], amount=1, label='coil3')
a.coils(Hp, Bp, start=Coord(marg+Bz+3*tz, Hy+marg), step=[Bz, 0], amount=1, label='coil4')
a.inductor(Hy, Bi, Hp, Bz, tz, 6, Coord(marg, marg), label='inductor1')
a.rect(height=d_se, width=0.1, intialX=marg, intialY=marg+Hy+Hp+dz, label='rect1')
#a.copy(a, gNodesLabel='coil1', step=[3*tz,0], number=1, label='coil4')


#  create mesh

k = 2
start = time.time()
mesh1 = Mesh(label='MeshTest')

mesh1.discretX((marg, k), ((Bz, k), (Bp, k)), True, 4)
mesh1.discretY(((marg, k), (Hy, k), (Hp, k), (dz, k), (d_se, 2), (marg, k)))
mesh1.createMesh()
# Определим какому телу (body) пренадлежит каждый элемент сетки
Bodies().defineBodies(a.gNodes, mesh1)
stop = time.time()
print(stop-start)


Material().setMat(mesh1, material=matIron, body='inductor1')
Material().setMat(mesh1, material=matCooper, body='coil1')
Material().setMat(mesh1, material=matCooper, body='coil2')
Material().setMat(mesh1, material=matCooper, body='coil3')
Material().setMat(mesh1, material=matCooper, body='coil4')
Material().setMat(mesh1, material=matIron, body='rect1')
Material().setMat(mesh1, material=matAir, body=None)

mf = MagneticField(mesh1, omega=314, label='mf1')
mf.definiceCurrent(current=current_a, body='coil1')
mf.definiceCurrent(current=current_z, body='coil2')
mf.definiceCurrent(current=current_b, body='coil3')
mf.definiceCurrent(current=current_x, body='coil4')

r, mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down = mf.set_matrix_complex_2()
mmf = mf.mmf(mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down)

import Modules.Solvers as sol

solution = sol.solve_it_ling(r, mmf)


from Modules.PostProcessing import PostProcessing
pp_class = PostProcessing(solution, mf)

reshape_sol = pp_class.reshape_data(solution)
magnetic_flux_x = pp_class.calculate_magnetic_flux_x(reshape_sol)
magnetic_flux_y = pp_class.calculate_magnetic_flux_y(reshape_sol)
reshape_mmf = pp_class.reshape_data(mmf)

pp_class.create_pcolor(magnetic_flux_y)
pp_class.create_pcolor(magnetic_flux_x)
pp_class.create_pcolor(reshape_mmf)


if __name__ == '__main__':
    main()


