

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

air = Material(0, 1)
copper = Material(5e7, 1)
steel = Material(0, 1000)

Coord.scale_factor = 1

inductor = Body(steel).inductor(Hy, Bi, Hp, Bz, tz, Q, Coord(marg, 0.0))
coil = Body(copper).coils(Hp, Bp, Coord(marg+Bz, Hy), Coord(Bz, 0), Q, "phaseA", current=1e9)

secondary = Body(copper).rect(d_se, Bi, Coord(marg, Hy+Hp+dz))

axis_x = discret_X((marg, 10), ((Bz, 2), (Bp, 1)), True, Q)
axis_y = discret_Y(((Hy, 3), (Hp, 3), (dz, 3), (d_se, 2)))

data = body_grid(axis_x, axis_y, Body.bodies)


from Physics import MagneticField
a = MagneticField(data, omega=10)
r, mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down = a.set_matrix_complex()
mmf = a.mmf(mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down)



import Solvers as sol


solution = sol.solve_it(r, mmf)

from PostProcessing import PostProcessing

pp_class = PostProcessing(solution, a)
reshape_sol = pp_class.reshape_data(solution)
reshape_mmf = pp_class.reshape_data(mmf)

magnetic_flux_x = pp_class.calculate_magnetic_flux_x(reshape_sol)
magnetic_flux_y = pp_class.calculate_magnetic_flux_y(reshape_sol)

pp_class.create_simple_plot(magnetic_flux_x, layer=7)
pp_class.create_simple_plot(magnetic_flux_y, layer=7)
