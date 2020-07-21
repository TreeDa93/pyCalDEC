import numpy as np












import numpy as np
data = np.array([[3, 8, 3], [6, 2, 3], [11, 2, 3]])


###################################

import numpy as np
data = np.array([[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]])
i = data[0,:]
j = data[1,:]

def test_fun(i, j):
    y = i ** 2 + j
    return y


class Test:


    def __init__(self, data):
        self.a = data

    def test_fun(self):
        i = self.a[0, :]
        j = self.a[1, :]
        y1 = i ** 2 + j
        y2 = 2 * i ** 2 + j
        return y1, y2

    from scipy.sparse import dia_matrix
    def create_dia_matrix(self,y,k):
        size = len(y[0])
        x = dia_matrix((y, [0, 1]), shape=(size, size)).toarray()
        return x


#####################################

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
        i_1 = np.arange(self.size)
        j_1 = np.arange(self.size)
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
                r[i_1[counter], j_1[counter]] = 2 * (self.formula_resistance_rmt(i, j)
                                                     + self.formula_resistance_rmn(i, j)) \
                                                + self.formula_inductance_term(i, j)
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




