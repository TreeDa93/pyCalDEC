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


    def __init__(self, data):
        self.a = data # Массив с описанием данных
        self.L = 1
        self.size_j = len(data[0])-1
        self.size_i = len(data)-1
        self.size = self.size_i*self.size_j

    def formula_resistance_rmn_left(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to lefter

        a[i][j].center.y is value y coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].width() is the width of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                """
        rmn_left = (self.a[i][j + 1].center.y - self.a[i][j].center.y) / \
                            (4 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu() + 1 /
                                                                   self.a[i][j + 1].mu())
        return rmn_left

    def formula_resistance_rmn_right(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to righter

                a[i][j].center.y is the value y coordianat of the center of  the i-th cell at the j-th layer
                a[i][j].width() is the width of the i-th cell at the j-th layer
                a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell
                located at the j-th layer
                        """
        rmn_right = (self.a[i + 1][j + 1].center.y - self.a[i + 1][j].center.y) / \
                             (4 * self.a[i + 1][j].width() * self.L) * (1 / self.a[i + 1][j].mu() + 1
                                                                        / self.a[i + 1][j + 1].mu())
        return rmn_right

    def formula_resistance_rmt_up(self, i, j):
        """ The function describes formula of  mutual tangent magnetic resistance of circuit loop with the loop located
        to upper

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                        """
        rmt_up = (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) / \
                          (4 * self.a[i][j + 1].height() * self.L) * \
                          (1 / self.a[i][j + 1].mu() + 1 / self.a[i][j].mu())
        return rmt_up

    def formula_resistance_rmt_down(self, i, j):
        """ The function describes formula of  mutual normal magnetic resistance of circuit loop with the loop located
        to downer

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                               """
        rmt_down = (self.a[i + 1][j].center.x - self.a[i][j].center.x) / \
                            (4 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu() + 1
                                                                    / self.a[i][j + 1].mu())
        return rmt_down

    def formula_resistance_rmt(self, i, j):
        """ The function describes formula of down tangent self magnetic resistance of circuit loop

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                                       """
        rmt = ((self.a[i + 1][j].center.x - self.a[i][j].center.x) /
                        (4 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu()
                                                                + 1 / self.a[i][j + 1].mu())
                        + (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) /
                        (4 * self.a[i][j + 1].height() * self.L) * (1 / self.a[i][j + 1].mu()
                                                                    + 1 / self.a[i][j].mu()))
        return rmt

    def formula_resistance_rmn(self, i, j):
        """ The function describes formula of down normal self magnetic resistance of circuit loop

        a[i][j].center.x is the value x coordianat of the center of  the i-th cell at the j-th layer
        a[i][j].height() is the height of the i-th cell at the j-th layer
        a[i][j].mu() is the method calling values of magnetic permeability at the i-th cell located at the j-th layer
                                               """
        rmn = ((self.a[i][j + 1].center.y - self.a[i][j].center.y) /
                        (4 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu()
                                                               + 1 / self.a[i][j + 1].mu())
                        + (self.a[i + 1][j + 1].center.y - self.a[i + 1][j].center.y) /
                        (4 * self.a[i + 1][j].width() * self.L) * (1 / self.a[i + 1][j].mu() + 1 /
                                                                   self.a[i + 1][j + 1].mu()))
        return rmn

    def set_matrix(self):
        """The function is intended to build the stiffness matrix of linear equations system.
        Matrix filling occurs from left to right and to the top. Namely, the first equation is being built for
        the first cell with indexes (0, 0). The function write the self resistance at general diagonal (0, 0) velue of
        right normal mutual resistance at (0 , 1) because the resistance is mutual with second loop flux. and R_up is
        written to (0, j_size) because the resistance is mutual with loop upper flux, that we can obtain index of
        the flux we should add number of cells at the layer. Next equation well be written for sencod loop flux with
        index (0, 1) and so on.

                0       1       2        3       j      |
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
        r = np.zeros((self.size, self.size))
        i_1 = np.arange(self.size)
        j_1 = np.arange(self.size)
        mmf_bc_up = np.zeros(self.size)
        mmf_bc_down = np.zeros(self.size)
        mmf_bc_right = np.zeros(self.size)
        mmf_bc_left = np.zeros(self.size)
        bc_fluxes_up = np.zeros(self.size_i)
        bc_fluxes_down = np.zeros(self.size_i)
        bc_fluxes_right = np.zeros(self.size_j)
        bc_fluxes_left = np.zeros(self.size_j)
        counter = 0

        for j in range(self.size_j):
            for i in range(self.size_i):
                r[i_1[counter], j_1[counter]] = 2 * (self.formula_resistance_rmt(i, j)
                                                     + self.formula_resistance_rmn(i, j))
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
                    r[i_1[counter], j_1[counter] + self.size_j] = -self.formula_resistance_rmt_up(i, j)
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
                    r[i_1[counter], j_1[counter] - self.size_j] = \
                        -self.formula_resistance_rmt_down(i, j)
                    counter += 1

        counter = 0
        counter_bc = 0
        stopper = self.size_i - 1
        for j in range(self.size_j):
            for i in range(self.size_i):
                if counter == stopper:
                    mmf_bc_left[counter] = -self.formula_resistance_rmn_left(i, j) * bc_fluxes_left[counter_bc]
                    counter += 1
                    counter_bc += 1
                    stopper += self.size_i
                else:
                    r[i_1[counter] + 1, j_1[counter]] = -self.formula_resistance_rmn_left(i, j)
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
                    r[i_1[counter], j_1[counter] + 1] = -self.formula_resistance_rmn_right(i, j)
                    counter += 1
        return r, mmf_bc_right, mmf_bc_left, mmf_bc_up, mmf_bc_down

    def mmf(self):
        """The function is intended to build the vector-column of right free terms of linear equations system """
        counter = 0
        mmf = np.zeros((self.size))
        for j in range(self.size_j):
            for i in range(self.size_i):
                if self.a[i][j].type == 'coil':
                    mmf[counter] = 100000
                    counter += 1
                else:
                    counter += 1

    def mmf_test(self):

        counter = 0
        mmf = np.zeros(self.size)
        for j in range(self.size_j):
            for i in range(self.size_i):
                if j == 3:
                    mmf[counter] = 100000
                    counter += 1
                else:
                    counter += 1
        return mmf

