"""The part of the program is to set up generals features of the model"""
air = Material(0, 1)
# copper = Material(4.528e7, 1)
steel = Material(0, 1000)

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

from PhysicsMF import MagneticField
a = MagneticField(data, omega=314)

r, mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down = a.set_matrix_complex_2()
mmf = a.mmf(mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down)
#r = a.test_fun()

"""The part of the program is to run simulation"""

import Solvers as sol

solution = sol.solve_it_ling(r, mmf)

solution_2 = sol.solve_it(r, mmf)


"""The part of the program is to present the results"""

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


