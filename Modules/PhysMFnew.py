import numpy as np

class MagneticField2:
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
        self.a = mesh.mesh  # Массив с описанием данных
        self.L = L
        self.size_i = mesh.sizeX - 1
        self.size_j = mesh.sizeY - 1
        self.size = self.size_i * self.size_j
        self.omega = omega
        self.i = mesh.sizeX
        self.j = mesh.sizeY
        self.size_cell = self.i * self.j

    def definiceCurrent(self, current=1000, body='coil'):
        for x in range(self.i):
            for y in range(self.j):
                if self.a[x][y].body == body:
                    self.a[x][y].defineCurrent(current)



    def formula_resistance_rmn_left(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to lefter
        R = h / (2 * w * L * mu) + h / (2 * w * L * mu)
        a[i][j].center.y is value y coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].width() is the width of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                """
        rmn_left = self.a[i][j].height() / (2 * self.a[i][j].width() * self.L * self.a[i][j].mu()) + \
                   self.a[i][j + 1].height() / (2 * self.a[i][j + 1].width() * self.L * self.a[i][j + 1].mu())

        return rmn_left

    def formula_resistance_rmn_right(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to righter

                a[i][j].center.y is the value y coordianat of the center of  the i-th cell at the j-th layer
                a[i][j].width() is the width of the i-th cell at the j-th layer
                a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell
                located at the j-th layer
                        """

        rmn_right = self.a[i + 1][j].height() / (2 * self.a[i + 1][j].width() * self.L * self.a[i + 1][j].mu()) + \
                    self.a[i + 1][j + 1].height() / (2 * self.a[i + 1][j + 1].width() * self.L * self.a[i + 1][j + 1].mu())

        return rmn_right

    def formula_resistance_rmt_up(self, i, j):
        """ The function describes formula of  mutual tangent magnetic resistance of circuit loop with the loop located
        to upper
        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                        """
        rmt_up = self.a[i][j + 1].width() / (2 * self.a[i][j + 1].height() * self.L * self.a[i][j + 1].mu()) + \
                 self.a[i + 1][j + 1].width() / (2 * self.a[i + 1][j + 1].height() * self.L * self.a[i + 1][j + 1].mu())
        return rmt_up


    def formula_resistance_rmt_down(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to downer

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                               """
        rmt_down = self.a[i][j].width() / (2 * self.a[i][j].height() * self.L * self.a[i][j].mu()) + \
                   self.a[i + 1][j].width() / (2 * self.a[i + 1][j].height() * self.L * self.a[i + 1][j].mu())
        return rmt_down


    def formula_resistance_rmt(self, i, j):
        """ The function describes formula of down tangent self magnetic resistance of circuit loop

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                                       """
        rmt = self.formula_resistance_rmt_up(i, j) + self.formula_resistance_rmt_down(i, j)
        return rmt

    def formula_resistance_rmn(self, i, j):
        """ The function describes formula of down normal self magnetic resistance of circuit loop

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                                               """
        rmn = self.formula_resistance_rmn_left(i, j) + self.formula_resistance_rmn_right(i, j)
        return rmn

    def formula_self_resistance(self, i, j):
        r_self = (self.formula_resistance_rmn_left(i, j) + self.formula_resistance_rmn_right(i, j) +
                  self.formula_resistance_rmt_up(i, j) + self.formula_resistance_rmt_down(i, j))
        return r_self

    def formula_mmf_coil(self, i, j):

        mmf_coil = (self.a[i][j].current * 1 / 4 * self.a[i][j].width() * self.a[i][j].height()
                    + self.a[i + 1][j + 1].current * 1 / 4 * self.a[i+1][j+1].width() * self.a[i+1][j+1].height()
                    + self.a[i][j + 1].current * 1 / 4 * self.a[i][j+1].width() * self.a[i][j+1].height()
                    + self.a[i+1][j].current * 1 / 4 * self.a[i+1][j].width() * self.a[i+1][j].height())

        # mmf_coil = self.a[i][j].width() * self.a[i][j].height() / 4 * \
        #            (self.a[i][j].current + self.a[i + 1][j + 1].current
        #             + self.a[i][j + 1].current + self.a[i+1][j].current)
        return mmf_coil

    def formula_inductance_term(self, i, j):
        square = (self.a[i + 1][j].center.x - self.a[i][j].center.x) * \
                 (self.a[i][j + 1].center.y - self.a[i][j].center.y)
        ind_term = 1j * square / 4 * (self.a[i][j].sigma() + self.a[i + 1][j].sigma() + self.a[i][j + 1].sigma() +
                                      self.a[i + 1][j + 1].sigma())
        return ind_term


    def formula_resistance_rmn_left_velocity(self, i, j, velocity_x=0):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to lefter

        a[i][j].center.y is value y coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].width() is the width of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                """
        delta_y = (self.a[i][j + 1].center.y - self.a[i][j].center.y) / 2
        rmn_left = delta_y / (2 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu() + 1 / self.a[i][j + 1].mu()) \
                   - (velocity_x * delta_y) / (4 * self.L) * (self.a[i][j].sigma() + self.a[i][j + 1].sigma())

        return rmn_left

    def formula_resistance_rmn_right_velocity(self, i, j, velocity_x=0):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to righter

                a[i][j].center.y is the value y coordianat of the center of  the i-th cell at the j-th layer
                a[i][j].width() is the width of the i-th cell at the j-th layer
                a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell
                located at the j-th layer
                        """
        delta_y = (self.a[i + 1][j + 1].center.y - self.a[i + 1][j].center.y) / 2
        rmn_right = delta_y / (2 * self.a[i + 1][j].width() * self.L) * (1 / self.a[i + 1][j].mu() +
                                                                         1 / self.a[i + 1][j + 1].mu()) + \
                    (velocity_x * delta_y) / (4 * self.L) * (self.a[i + 1][j].sigma() + self.a[i + 1][j + 1].sigma())
        return rmn_right

    def formula_resistance_rmt_up_velocity(self, i, j, velocity_y=0):
        """ The function describes formula of  mutual tangent magnetic resistance of circuit loop with the loop located
        to upper

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                        """
        delta_x = (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) / 2
        rmt_up = delta_x / (2 * self.a[i][j + 1].height() * self.L) * \
                 (1 / self.a[i][j + 1].mu() + 1 / self.a[i + 1][j + 1].mu()) - \
                 (velocity_y * delta_x) / (4 * self.L) * (self.a[i][j + 1].sigma() + self.a[i + 1][j + 1].sigma())
        return rmt_up

    def formula_resistance_rmt_down_velocity(self, i, j, velocity_y=0):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to downer

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                               """
        delta_x = (self.a[i + 1][j].center.x - self.a[i][j].center.x) / 2
        rmt_down = delta_x / (2 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu() +
                                                                     1 / self.a[i + 1][j].mu()) + \
                   (velocity_y * delta_x) / (4 * self.L) * (self.a[i][j].sigma() + self.a[i + 1][j].sigma())
        return rmt_down


    def initMatrixResistence(self):
        self.matrixMagResistence = np.zeros((self.size, self.size), dtype=complex)


    def calculateSelfResistence(self):
        counter = 0
        for j in range(self.size_j):
            for i in range(self.size_i):
                self.matrixMagResistence[counter, counter] = self.formula_self_resistance(i, j) \
                                                             + self.formula_inductance_term(i, j)
                counter += 1
        return self.matrixMagResistence



    def runner(self):
        for indexX, valueX in enumerate(self.a):
            for indexY, value, in enumerate(valueX):
                return value, indexX, indexY

    def defineMatrix(self):

        #cell, i, j = self.runner()
        for i, valueX in enumerate(self.a):
            for j, cell, in enumerate(valueX):
                self.matrixMagResistence[cell.index, cell.index] = self.formula_self_resistance(i, j) \
                                                                     + self.formula_inductance_term(i, j)
        return self.matrixMagResistence
