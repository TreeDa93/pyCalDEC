

#from Data import *  # Import the general data of the model

from Discretization import *

air = Material(0, 1)
copper = Material(100, 1)
axis_x = discret_x(((10, 5), (20, 10), (10, 5)))
axis_y, mat_y = discret_y(((10, 1, air), (15, 2, copper), (30, 4, air)))
data = grid(axis_x, axis_y, mat_y, 10)

from Physics import Solver1


a = Solver1(data)

rmn_left, rmn_right, rmt_up, rmt_down, rmn, rmt = a.calculate_resistance()


import Solvers as sl

z = sl.create_dia_matrix(rmn, 0)

f = sl.create_f(rmn)

solution = sl.solve_it(z, f)




