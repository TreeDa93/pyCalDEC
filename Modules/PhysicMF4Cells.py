import numpy as np


class MagneticField:
    """"This class is intended for creating matrix of linear equations based on
    physical equation of magnetic filed and circuit theory
    Main parameters of the class is the data of models consisting to "a" massiv
    a is the input massive of data
    size_j is the size or number of layers along y component minus one. Other word is the number of resistances along
    y component
    size_i is the number of elements at each layer along x component minus one. Other word is the number of resistances
    along x component.
    L is width of model directing to plane of the model
    size is the number of unknown of linear equation system.

    """

    def __init__(self, mesh, omega=0, L=0.5, label='mf'):
        self.mesh = mesh
        self.cells = mesh.mesh  # Массив с описанием данных
        self.L = L
        self.omega = omega
        self.size_i = mesh.sizeX - 1
        self.size_j = mesh.sizeY - 1
        self.size = (mesh.sizeX - 1)*(mesh.sizeY - 1)
        self.i = mesh.sizeX
        self.j = mesh.sizeY
        self.size_cell = self.i * self.j


    def definiceCurrent(self, current=1000, body='coil'):
        for x in range(self.mesh.sizeX):
            for y in range(self.mesh.sizeY):
                if self.cells[x][y].body == body:
                    self.cells[x][y].defineCurrent(current)

    def set_matrix_complex(self, bc_fluxes_up=None, bc_fluxes_down=None,
                                 bc_fluxes_right=None, bc_fluxes_left=None):
        """The function is intended to build the stiffness matrix of linear equations system.
        Matrix filling occurs from left to right and to the top. Namely, the first equation is being built for
        the first cell with indexes (0, 0). The function write the self resistance at general diagonal (0, 0) velue of
        right normal mutual resistance at (0 , 1) because the resistance is mutual with second loop flux. and R_up is
        written to (0, j_size) because the resistance is mutual with loop upper flux, that we can obtain index of
        the flux we should add number of cells at the layer. Next equation well be written for sencod loop flux with
        index (0, 1) and so on.

                0       1       2        3       j
         0   |R_self  R-right    -       -       R_up   |
         1   |R_l     R_self   R-right                  |
         j   |-         R_l     R_self   R-right        |
         3   |-                  R_l     R_self         |
             |                                          |
         j   |R_d                                 R_self|

        r is the variables contained values of stiffness matrix elements. It has size the number of cells minus one
        at layer and the number of layers minus one.
        counter is special advanced variable to understand current cell number.
        """
        r = np.zeros((self.size, self.size), dtype=complex)
        mmf_bc_up = np.zeros(self.size, dtype=complex)
        mmf_bc_down = np.zeros(self.size, dtype=complex)
        mmf_bc_right = np.zeros(self.size, dtype=complex)
        mmf_bc_left = np.zeros(self.size, dtype=complex)
        if bc_fluxes_up is None:
            bc_fluxes_up = np.zeros(self.size_i, dtype=complex)
        if bc_fluxes_down is None:
            bc_fluxes_down = np.zeros(self.size_i, dtype=complex)
        if bc_fluxes_right is None:
            bc_fluxes_right = np.zeros(self.size_j, dtype=complex)
        if bc_fluxes_left is None:
            bc_fluxes_left = np.zeros(self.size_j, dtype=complex)
        counter = 0

        for j in range(self.size_j):
            for i in range(self.size_i):
                r[counter, counter] = self.formula_self_resistance(i, j) + self.formula_inductance_term(i, j)
                counter += 1

        counter = 0
        counter_bc = 0
        for j in range(self.size_j):
            for i in range(self.size_i):
                if counter >= (self.size - self.size_i):
                    mmf_bc_up[counter] = -self.formula_resistance_rmt_up(i, j) * bc_fluxes_up[counter_bc]
                    counter += 1
                    counter_bc += 1
                else:
                    r[counter, counter + self.size_i] = -self.formula_resistance_rmt_up(i, j)
                    counter += 1

        counter = 0
        counter_bc = 0
        for j in range(self.size_j):
            for i in range(self.size_i):
                if counter < self.size_i:
                    mmf_bc_down[counter] = -self.formula_resistance_rmt_down(i, j) * bc_fluxes_down[counter_bc]
                    counter += 1
                    counter_bc += 1
                else:
                    r[counter, counter-self.size_i] = -self.formula_resistance_rmt_down(i, j)
                    counter += 1

        counter = 0
        counter_bc = 0
        stopper = 0
        for j in range(self.size_j):
            for i in range(self.size_i):
                if counter == stopper:
                    mmf_bc_left[counter] = -self.formula_resistance_rmn_left(i, j) * bc_fluxes_left[counter_bc]
                    counter += 1
                    counter_bc += 1
                    stopper += self.size_i
                else:
                    r[counter, counter-1] = -self.formula_resistance_rmn_left(i, j)
                    counter += 1

        counter = 0
        counter_bc = 0
        stopper = self.size_i - 1
        for j in range(self.size_j):
            for i in range(self.size_i):
                if counter == stopper:
                    mmf_bc_right[counter] = -self.formula_resistance_rmn_right(i, j) * bc_fluxes_right[counter_bc]
                    counter += 1
                    counter_bc += 1
                    stopper += self.size_i
                else:
                    r[counter, counter+1] = -self.formula_resistance_rmn_right(i, j)
                    counter += 1
            mmf = self.mmf(mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down)
            self.weightMatrix = r
            self.rightMatrix = mmf

    def formula_resistance_rmn_left(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to lefter
        R = h / (2 * w * L * mu) + h / (2 * w * L * mu)
        a[i][j].center.y is value y coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].width() is the width of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                """
        rmn_left = self.cells[i][j].height() / (2 * self.cells[i][j].width() * self.L * self.cells[i][j].mu()) + \
                   self.cells[i][j + 1].height() / (2 * self.cells[i][j + 1].width() * self.L * self.cells[i][j + 1].mu())

        return rmn_left

    def formula_resistance_rmn_right(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to righter

                a[i][j].center.y is the value y coordianat of the center of  the i-th cell at the j-th layer
                a[i][j].width() is the width of the i-th cell at the j-th layer
                a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell
                located at the j-th layer
                        """

        rmn_right = self.cells[i + 1][j].height() / (2 * self.cells[i + 1][j].width() * self.L * self.cells[i + 1][j].mu()) + \
                    self.cells[i + 1][j + 1].height() / (2 * self.cells[i + 1][j + 1].width() * self.L * self.cells[i + 1][j + 1].mu())

        return rmn_right

    def formula_resistance_rmt_up(self, i, j):
        """ The function describes formula of  mutual tangent magnetic resistance of circuit loop with the loop located
        to upper
        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                        """
        rmt_up = self.cells[i][j + 1].width() / (2 * self.cells[i][j + 1].height() * self.L * self.cells[i][j + 1].mu()) + \
                 self.cells[i + 1][j + 1].width() / (2 * self.cells[i + 1][j + 1].height() * self.L * self.cells[i + 1][j + 1].mu())
        return rmt_up


    def formula_resistance_rmt_down(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to downer

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                               """
        rmt_down = self.cells[i][j].width() / (2 * self.cells[i][j].height() * self.L * self.cells[i][j].mu()) + \
                   self.cells[i + 1][j].width() / (2 * self.cells[i + 1][j].height() * self.L * self.cells[i + 1][j].mu())
        return rmt_down

    def formula_self_resistance(self, i, j):
        r_self = (self.formula_resistance_rmn_left(i, j) + self.formula_resistance_rmn_right(i, j) +
                  self.formula_resistance_rmt_up(i, j) + self.formula_resistance_rmt_down(i, j))
        return r_self

    def formula_mmf_coil(self, i, j):
        mmf1 = 0.25 * self.cells[i][j].current * self.cells[i][j].calculateSquare()
        mmf2 = 0.25 * self.cells[i+1][j].current * self.cells[i+1][j].calculateSquare()
        mmf3 = 0.25 * self.cells[i][j+1].current * self.cells[i][j+1].calculateSquare()
        mmf4 = 0.25 * self.cells[i+1][j+1].current * self.cells[i+1][j+1].calculateSquare()
        mmf = mmf1 + mmf2 + mmf3 + mmf4
        return mmf

    def formula_inductance_term(self, i, j):
        square = 0.25 * (self.cells[i][j].calculateSquare() + self.cells[i+1][j].calculateSquare() +
                         self.cells[i][j+1].calculateSquare() + self.cells[i+1][j+1].calculateSquare())
        sigma = 0.25 * (self.cells[i][j].sigma() + self.cells[i+1][j].sigma() + self.cells[i][j+1].sigma() +
                        self.cells[i+1][j+1].sigma())

        ind_term = 1j * self.omega * square * sigma / self.L
        return ind_term

  

    def mmf(self, mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down):
        """The function is intended to build the vector-column of right free terms of linear equations system """
        counter = 0
        mmf = np.zeros(self.size, dtype=complex)
        for j in range(self.size_j):
            for i in range(self.size_i):
                mmf[counter] = self.formula_mmf_coil(i, j) + mmf_bc_right[counter] + mmf_bc_left[counter] + \
                               mmf_bc_up[counter] + mmf_bc_down[counter]
                counter += 1
        return mmf

    def set_bc(self, bc_fluxes=None, number_element=None, values=0, type_bc='right'):
        if bc_fluxes is None:
            if type_bc == 'right' or 'left':
                bc_fluxes = np.zeros(self.size_j)
                if number_element is None:
                    number_element = np.arange(self.size_j)
                    bc_fluxes[number_element] = values
                else:
                    bc_fluxes[number_element] = values
            else:
                bc_fluxes = np.zeros(self.size_i)
                if number_element is None:
                    number_element = np.arange(self.size_i)
                    bc_fluxes[number_element] = values
                else:
                    bc_fluxes[number_element] = values
        else:
            if type_bc == 'right' or 'left':
                if number_element is None:
                    number_element = np.arange(self.size_j)
                    bc_fluxes[number_element] = values
                else:
                    bc_fluxes[number_element] = values
            else:
                if number_element is None:
                    number_element = np.arange(self.size_i)
                    bc_fluxes[number_element] = values
                else:
                    bc_fluxes[number_element] = values
        return bc_fluxes

    

