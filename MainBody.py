

#from Data import *  # Import the general data of the model

from Discretization import *

Hy = 0.005
Hp = 6 * 0.001
dz = 1 * 0.001
d_se = 2 * 0.001
Bi = 17 * 0.001
Bz = 2 * 0.001
Bp = 3 * 0.001
Q = 3
marg = 10 * 0.001
tz = Bp + Bz

air = Material(0, 0)
copper = Material(100, 1)
steel = Material(0, 500)

Coord.scale_factor = 1

inductor = Body(steel).inductor(Hy, Bi, Hp, Bz, tz, Q, Coord(marg, 0.0))
coil = Body(copper).coils(Hp, Bp, Coord(marg+Bz, Hy), Coord(Bz, 0), Q, "phaseA")

secondary = Body(copper).rect(d_se, Bi, Coord(marg, Hy+Hp+dz))

axis_x = discret_X((marg, 1), ((Bz, 1), (Bp, 1)), True, Q)
axis_y = discret_Y(((Hy, 1), (Hp, 1), (dz, 1), (d_se, 1)))

data = body_grid(axis_x, axis_y, Body.bodies)


from Physics import MagneticField
a = MagneticField(data)
r, mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down = a.set_matrix()
mmf = a.mmf(mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down)

import Solvers as sl


solution = sl.solve_it(r, mmf)
bc_fluxes_up = a.set_bc(type_bc='up')
bc_fluxes_down = a.set_bc(type_bc='down')
bc_fluxes_right = a.set_bc(type_bc='right')
bc_fluxes_left = a.set_bc(type_bc='left')
test2 = a.set_matrix_complex_velocity()
# #test = a.set_matrix_complex(bc_fluxes_up=bc_fluxes_up, bc_fluxes_down=bc_fluxes_down, bc_fluxes_right=bc_fluxes_right,
#                             bc_fluxes_left=bc_fluxes_left)



