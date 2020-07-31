

#from Data import *  # Import the general data of the model

from Discretization import *
import math as ma
import numpy as np
Hy = 10 * 1e-3
Hp = 10 * 1e-3
dz = 10 * 1e-3
d_se = 10 * 1e-3
Bz = 10 * 1e-3
Bp = 10 * 1e-3
Q = 3
tz = Bp + Bz
Bi = tz * Q + Bz
marg = 10

air = Material(0, 1)
# copper = Material(4.528e7, 1)
steel = Material(0, 1000)

phi_a = 0
phi_z = - ma.pi / 6
phi_b = - ma.pi*2/3
phi_x = - ma.pi
phi_c = ma.pi * 2 / 3
phi_y = ma.pi / 6

amplitude_current = 3e7
current_a = amplitude_current * np.exp(1j*phi_a)
current_b = amplitude_current * np.exp(1j*phi_b)
current_x = amplitude_current * np.exp(1j*phi_x)
current_c = amplitude_current * np.exp(1j*phi_c)
current_y = amplitude_current * np.exp(1j*phi_y)
current_z = amplitude_current * np.exp(1j*phi_z)


inductor = Body(steel).inductor(Hy, Bi, Hp, Bz, tz, 3, Coord(marg, marg))
coil_1 = Body(air).coils(Hp, Bp, start=Coord(marg+Bz, Hy+marg), amount=1, index="phaseA", current=current_a)
coil_2 = Body(air).coils(Hp, Bp, start=Coord(marg+Bz+tz, Hy+marg), amount=1, index="phaseZ", current=current_z)
coil_3 = Body(air).coils(Hp, Bp, start=Coord(marg+Bz+2*tz, Hy+marg), amount=1, index="phaseB", current=current_b)
# coil_4 = Body(copper).coils(Hp, Bp, start=Coord(marg+Bz+3*tz, Hy), amount=1, index="phaseX", current=current_x)
# coil_5 = Body(copper).coils(Hp, Bp, start=Coord(marg+Bz+4*tz, Hy), amount=1, index="phaseC", current=current_c)
# coil_6 = Body(copper).coils(Hp, Bp, start=Coord(marg+Bz+5*tz, Hy), amount=1, index="phaseY", current=current_y)
# secondary = Body(air).rect(50, 0.12, Coord(0, 0.07))

back_iron = Body(steel).rect(0.05, Bi, Coord(marg, marg+Hy+Hp+dz))

k = 8
axis_x = discret_X((marg, k), ((Bz, k), (Bp, k)), True, Q)
axis_y = discret_Y(((marg, k), (Hy, k), (Hp, k), (dz, k), (0.05, k)))

data = body_grid(axis_x, axis_y, Body.bodies)



# Hp = 20 * 0.001
# Bp = 20 * 0.001
# marg = 20 * 0.001
# copper = Material(0, 1)
# coil = Body(copper).coils(Hp, Bp, Coord(marg, marg), index="phaseA", current=3e7)
# back_iron = Body(steel).rect(0.05, Bp+2*marg, Coord(0, 2*marg+Hp))
# axis_x = discret_X(domains=((marg, 4), (Bp, 4), (marg, 4)))
# axis_y = discret_Y(((marg, 4), (Hp, 4), (marg, 4), (0.05, 4)))
#
# data = body_grid(axis_x, axis_y, Body.bodies)
#
#
# phi_a = 0
# phi_z = - ma.pi / 6
# phi_b = - ma.pi*2/3
# phi_x = - ma.pi
# phi_c = ma.pi * 2 / 3
# phi_y = ma.pi / 6

# amplitude_current = 3e7
# current_a = amplitude_current * np.exp(1j*phi_a)
# current_b = amplitude_current * np.exp(1j*phi_b)
# current_x = amplitude_current * np.exp(1j*phi_x)
# current_c = amplitude_current * np.exp(1j*phi_c)
# current_y = amplitude_current * np.exp(1j*phi_y)
# current_z = amplitude_current * np.exp(1j*phi_z)
#
#
# inductor = Body(air).inductor(0.05, 0.12, 0.02, 0.05, 0.07, 1, Coord(0.0, 0.0))
# coil_1 = Body(air).coils(0.02, 0.02, start=Coord(0.05, 0.05), amount=1, index="phaseA", current=current_a)
# # coil_2 = Body(copper).coils(Hp, Bp, start=Coord(marg+Bz+tz, Hy), amount=1, index="phaseZ", current=current_z)
# # coil_3 = Body(copper).coils(Hp, Bp, start=Coord(marg+Bz+2*tz, Hy), amount=1, index="phaseB", current=current_b)
# # coil_4 = Body(copper).coils(Hp, Bp, start=Coord(marg+Bz+3*tz, Hy), amount=1, index="phaseX", current=current_x)
# # coil_5 = Body(copper).coils(Hp, Bp, start=Coord(marg+Bz+4*tz, Hy), amount=1, index="phaseC", current=current_c)
# # coil_6 = Body(copper).coils(Hp, Bp, start=Coord(marg+Bz+5*tz, Hy), amount=1, index="phaseY", current=current_y)
# secondary = Body(air).rect(50, 0.12, Coord(0, 0.07))
# #back_iron = Body(steel).rect(0.05, Bi, Coord(marg, Hy+Hp+dz+d_se))
#
# axis_x = discret_X((0.05, 1), ((Bz, 1), (Bp, 1)), True, Q)
# axis_y = discret_Y(((0.05, 1), (Hp, 1), (Hy, 1)))
#
# data = body_grid(axis_x, axis_y, Body.bodies)




from Physics import MagneticField
a = MagneticField(data, omega=314)

r, mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down = a.set_matrix_complex_2()
mmf = a.mmf(mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down)
#r = a.test_fun()



import Solvers as sol


solution = sol.solve_it_ling(r, mmf)

solution_2 = sol.solve_it(r, mmf)




from PostProcessing import PostProcessing

pp_class = PostProcessing(solution, a)
pp_class_2 = PostProcessing(solution_2, a)
reshape_sol = pp_class.reshape_data(solution)
reshape_sol_2 = pp_class.reshape_data(solution_2)
reshape_mmf = pp_class.reshape_data(mmf)

magnetic_flux_x = pp_class.calculate_magnetic_flux_x(reshape_sol)
magnetic_flux_y = pp_class.calculate_magnetic_flux_y(reshape_sol)
magnetic_flux_x_2 = pp_class.calculate_magnetic_flux_x(reshape_sol_2)
magnetic_flux_y_2 = pp_class.calculate_magnetic_flux_y(reshape_sol_2)

pp_class.simple_plot(magnetic_flux_y, layer=1)

pp_class.create_pcolor(reshape_mmf)
pp_class.create_pcolor(reshape_sol_2)
pp_class.create_pcolor(magnetic_flux_y)
pp_class.create_pcolor(magnetic_flux_x)

# pp_class.create_simple_plot(magnetic_flux_x, layer=1)
# pp_class.create_simple_plot(magnetic_flux_y, layer=1)
# pp_class.create_simple_plot(magnetic_flux_x_2, layer=1)
# pp_class.create_simple_plot(magnetic_flux_y_2, layer=1)

# pp_class_2.create_2d_plot(magnetic_flux_y)
# pp_class_2.create_2d_plot(reshape_mmf)

# tang_up_dia = r.diagonal(offset=a.size_j)
# tang_down_dia = r.diagonal(offset=-a.size_j)
# major_dia = r.diagonal()
# up_dia = r.diagonal(offset=1)
# down_dia = r.diagonal(offset=-1)

