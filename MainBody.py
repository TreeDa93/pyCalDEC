

#from Data import *  # Import the general data of the model

from Discretization import *

air = Material(0, 1)
copper = Material(100, 1)
axis_x = discret_x(((10, 200), (20, 200), (10, 10)))
axis_y, mat_y = discret_y(((10, 10, air), (15, 10, copper), (30, 20, air)))
data = grid(axis_x, axis_y, mat_y, 10)

from Physics import MagneticField
a = MagneticField(data)
r, mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down = a.set_matrix()
mmf = a.mmf_test()

import Solvers as sl


solution = sl.solve_it(r, mmf)




