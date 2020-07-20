

import numpy as np
class MagneticField:
    """Это класс отвечающий за расчет коэффициентов учета поперечного краевого эффекта.
    Пока в тестовом режиме он называется Boltoncalculate. Он расчитывает только коэффициент Болтона
    """


    def __init__(self, data):
        self.a = data # Массив с описанием данных
        self.L = 1
        self.size_j = len(data[0])-1
        self.size_i = len(data)-1
        self.size = self.size_i*self.size_j

    def formula_resistance_rmn_left(self, i, j):
        rmn_left = (self.a[i][j + 1].center.y - self.a[i][j].center.y) / \
                            (4 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu() + 1 /
                                                                   self.a[i][j + 1].mu())
        return rmn_left

    def formula_resistance_rmn_right(self, i, j):
        rmn_right = (self.a[i + 1][j + 1].center.y - self.a[i + 1][j].center.y) / \
                             (4 * self.a[i + 1][j].width() * self.L) * (1 / self.a[i + 1][j].mu() + 1
                                                                        / self.a[i + 1][j + 1].mu())
        return rmn_right

    def formula_resistance_rmt_up(self, i, j):
        rmt_up = (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) / \
                          (4 * self.a[i][j + 1].height() * self.L) * \
                          (1 / self.a[i][j + 1].mu() + 1 / self.a[i][j].mu())
        return rmt_up

    def formula_resistance_rmt_down(self, i, j):
        rmt_down = (self.a[i + 1][j].center.x - self.a[i][j].center.x) / \
                            (4 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu() + 1
                                                                    / self.a[i][j + 1].mu())
        return rmt_down

    def formula_resistance_rmt(self, i, j):
        rmt = ((self.a[i + 1][j].center.x - self.a[i][j].center.x) /
                        (4 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu()
                                                                + 1 / self.a[i][j + 1].mu())
                        + (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) /
                        (4 * self.a[i][j + 1].height() * self.L) * (1 / self.a[i][j + 1].mu()
                                                                    + 1 / self.a[i][j].mu()))
        return rmt

    def formula_resistance_rmn(self, i, j):
        rmn = ((self.a[i][j + 1].center.y - self.a[i][j].center.y) /
                        (4 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu()
                                                               + 1 / self.a[i][j + 1].mu())
                        + (self.a[i + 1][j + 1].center.y - self.a[i + 1][j].center.y) /
                        (4 * self.a[i + 1][j].width() * self.L) * (1 / self.a[i + 1][j].mu() + 1 /
                                                                   self.a[i + 1][j + 1].mu()))
        return rmn

    def calculate_resistance(self):
        """Рассчитываем взамное магнитное сопротивление вевей
        """
        rmn_left = np.zeros(self.size)
        rmn_right = np.zeros(self.size)
        rmt_up = np.zeros(self.size)
        rmt_down = np.zeros(self.size)
        rmn = np.zeros(self.size)
        rmt = np.zeros(self.size)
        counter = 0
        for i in range(self.size_i):
            for j in range(self.size_j):
                rmn_left[counter] = (self.a[i][j+1].center.y - self.a[i][j].center.y) / \
                                    (4 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu() + 1 /
                                                                           self.a[i][j+1].mu())

                rmn_right[counter] = (self.a[i+1][j+1].center.y - self.a[i+1][j].center.y) / \
                                     (4 * self.a[i+1][j].width() * self.L) * (1 / self.a[i+1][j].mu() + 1
                                                                              / self.a[i+1][j+1].mu())

                rmt_down[counter] = (self.a[i+1][j].center.x - self.a[i][j].center.x) / \
                                    (4 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu() + 1
                                                                            / self.a[i][j+1].mu())

                rmt_up[counter] = (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) / \
                                  (4 * self.a[i][j + 1].height() * self.L) * \
                                  (1 / self.a[i][j + 1].mu() + 1 / self.a[i][j].mu())

                rmn[counter] = ((self.a[i][j + 1].center.y - self.a[i][j].center.y) /
                                (4 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu()
                                                                       + 1 / self.a[i][j+1].mu())
                                + (self.a[i + 1][j + 1].center.y - self.a[i + 1][j].center.y) /
                                (4 * self.a[i + 1][j].width() * self.L) * (1 / self.a[i + 1][j].mu() + 1 /
                                                                           self.a[i+1][j+ 1].mu()))

                rmt[counter] = ((self.a[i + 1][j].center.x - self.a[i][j].center.x) /
                                (4 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu()
                                + 1 / self.a[i][j + 1].mu())
                                + (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) /
                                (4 * self.a[i][j + 1].height() * self.L) * (1 / self.a[i][j + 1].mu()
                                                                            + 1 / self.a[i][j].mu()))
                counter += 1
        return rmn_left, rmn_right, rmt_up, rmt_down, rmn, rmt

    def bd_flux(self):
        flux_down = np.zeros(self.size_i)
        flux_up = np.zeros(self.size_i)
        flux_left = np.zeros(self.size_j)
        flux_right = np.zeros(self.size_j)
        return flux_down, flux_up, flux_left, flux_right

    def create_mmf(self):
        if a.type == 'coil'
            if a.name_body == 'phaseA'
                J
        mmf = np.zeros((self.size))
        return mmf

    def set_matrix(self, mmf, flux_bc_right, flux_bc_left, flux_bc_up, flux_bc_down):
        right_term = np.zeros((self.size)) # изменить размер
        r_self = np.array([])
        index_r_self_i = np.array([])
        index_r_self_j = np.array([])
        rmn_left = np.array([])
        index_rmn_left_i = np.array([])
        index_rmn_left_j = np.array([])
        rmn_right = np.array([])
        index_rmn_right_i = np.array([])
        index_rmn_right_j = np.array([])
        rmt_up = np.array([])
        index_rmt_up_i = np.array([])
        index_rmt_up_j = np.array([])
        rmt_down = np.array([])
        index_rmt_down_i = np.array([])
        index_rmt_down_j = np.array([])
        counter = 0
        for j in range(self.size_j):
            for i in range(self.size_i):
                if i == 0 and j == 0:
                    right_term[counter] = mmf[counter] + self.formula_resistance_rmt_down(i,j) * flux_bc_down[i]\
                                          + self.formula_resistance_rmn_left(i,j) * flux_bc_left[j]

                    rmt_up = np.append(rmt_up, self.formula_resistance_rmt_up(i,j))
                    index_rmt_up_i = np.append(index_rmt_up_i, counter)
                    index_rmt_up_j = np.append(index_rmt_up_j, self.size_i+counter)

                    rmn_right = np.append(rmn_right, self.formula_resistance_rmn_right(i,j))
                    index_rmn_right_i = np.append(index_rmn_right_i, counter)
                    index_rmn_right_j = np.append(index_rmn_right_j, counter+1)

                    r_self = np.append(r_self, 2 * (self.formula_resistance_rmt(i, j))
                                       + self.formula_resistance_rmn(i, j))
                    index_r_self_i = np.append(index_r_self_i, counter)
                    index_r_self_j = np.append(index_r_self_j, counter)

                elif i == self.size_i and j == self.size_j:
                    right_term[counter] = mmf[counter] + self.formula_resistance_rmt_up(i,j) * flux_bc_up[i]\
                                          + self.formula_resistance_rmn_right(i,j) * flux_bc_right[j]

                    rmt_down = np.append(rmt_down, self.formula_resistance_rmt_down(i,j))
                    index_rmt_down_i = np.append(index_rmt_down_i, counter)
                    index_rmt_down_j = np.append(index_rmt_down_j, self.size-self.size_i)

                    rmn_left = np.append(rmn_left, self.formula_resistance_rmn_left(i,j))
                    index_rmn_left_i = np.append(index_rmn_left_i, counter)
                    index_rmn_left_j = np.append(index_rmn_left_j, self.size-1)

                    r_self = np.append(r_self, 2 * (self.formula_resistance_rmt(i, j))
                                       + self.formula_resistance_rmn(i, j))
                    index_r_self_i = np.append(index_r_self_i, counter)
                    index_r_self_j = np.append(index_r_self_j, counter)
                elif i == 0:
                    right_term[counter] = mmf[counter] + self.formula_resistance_rmt_down(i, j) * flux_bc_down[i]

                    rmt_up = np.append(rmt_up, self.formula_resistance_rmt_up(i, j))
                    index_rmt_up_i = np.append(index_rmt_up_i, counter)
                    index_rmt_up_j = np.append(index_rmt_up_j, self.size_i + counter)

                    rmn_right = np.append(rmn_right, self.formula_resistance_rmn_right(i, j))
                    index_rmn_right_i = np.append(index_rmn_right_i, counter)
                    index_rmn_right_j = np.append(index_rmn_right_j, counter+1)

                    rmn_left = np.append(rmn_left, self.formula_resistance_rmn_left(i, j))
                    index_rmn_left_i = np.append(index_rmn_left_i, counter)
                    index_rmn_left_j = np.append(index_rmn_left_j, counter-1)

                    r_self = np.append(r_self, 2 * (self.formula_resistance_rmt(i, j))
                                       + self.formula_resistance_rmn(i, j))
                    index_r_self_i = np.append(index_r_self_i, counter)
                    index_r_self_j = np.append(index_r_self_j, counter)
                elif j == 0:
                    right_term[counter] = mmf[counter] + self.formula_resistance_rmn_left(i, j) * flux_bc_left[j]

                    rmt_up = np.append(rmt_up, self.formula_resistance_rmt_up(i, j))
                    index_rmt_up_i = np.append(index_rmt_up_i, counter)
                    index_rmt_up_j = np.append(index_rmt_up_j, self.size_i + counter)

                    rmn_right = np.append(rmn_right, self.formula_resistance_rmn_right(i, j))
                    index_rmn_right_i = np.append(index_rmn_right_i, counter)
                    index_rmn_right_j = np.append(index_rmn_right_j, counter+1)

                    rmt_down = np.append(rmt_down, self.formula_resistance_rmt_down(i, j))
                    index_rmt_down_i = np.append(index_rmt_down_i, counter)
                    index_rmt_down_j = np.append(index_rmt_down_j, counter-self.size_i)

                    r_self = np.append(r_self, 2 * (self.formula_resistance_rmt(i, j))
                                       + self.formula_resistance_rmn(i, j))
                    index_r_self_i = np.append(index_r_self_i, counter)
                    index_r_self_j = np.append(index_r_self_j, counter)
                elif j == self.size_j: # возможно это не так!
                    right_term[counter] = mmf[counter] + self.formula_resistance_rmn_right(i, j) * flux_bc_right[j]

                    rmt_up = np.append(rmt_up, self.formula_resistance_rmt_up(i, j))
                    index_rmt_up_i = np.append(index_rmt_up_i, counter)
                    index_rmt_up_j = np.append(index_rmt_up_j, self.size_i + counter)

                    rmn_left = np.append(rmn_left, self.formula_resistance_rmn_left(i, j))
                    index_rmn_left_i = np.append(index_rmt_down_i, counter)
                    index_rmn_left_j = np.append(index_rmt_down_j, counter-1)

                    rmt_down = np.append(rmt_down, self.formula_resistance_rmt_down(i, j))
                    index_rmt_down_i = np.append(index_rmt_down_i, counter)
                    index_rmt_down_j = np.append(index_rmt_down_j, counter - self.size_i)

                    r_self = np.append(r_self, 2 * (self.formula_resistance_rmt(i, j))
                                       + self.formula_resistance_rmn(i, j))
                    index_r_self_i = np.append(index_r_self_i, counter)
                    index_r_self_j = np.append(index_r_self_j, counter)
                elif i == self.size_i:
                    right_term[counter] = mmf[counter] + self.formula_resistance_rmt_up(i, j) * flux_bc_up[i]

                    rmn_left = np.append(rmn_left, self.formula_resistance_rmn_left(i, j))
                    index_rmn_left_i = np.append(index_rmn_left_i, counter)
                    index_rmn_left_j = np.append(index_rmn_left_j, counter-1)

                    rmn_right = np.append(rmn_right, self.formula_resistance_rmn_right(i, j))
                    index_rmn_right_i = np.append(index_rmn_right_i, counter)
                    index_rmn_right_j = np.append(index_rmn_right_j, counter+1)

                    rmt_down = np.append(rmt_down, self.formula_resistance_rmt_down(i, j))
                    index_rmt_down_i = np.append(index_rmt_down_i, counter)
                    index_rmt_down_j = np.append(index_rmt_down_j, counter-self.size_i)

                    r_self = np.append(r_self, 2 * (self.formula_resistance_rmt(i, j))
                                       + self.formula_resistance_rmn(i, j))
                    index_r_self_i = np.append(index_r_self_i, counter)
                    index_r_self_j = np.append(index_r_self_j, counter)
                else:
                    right_term[counter] = mmf[counter]

                    rmt_up = np.append(rmt_up, self.formula_resistance_rmt_up(i, j))
                    index_rmt_up_i = np.append(index_rmt_up_i, counter)
                    index_rmt_up_j = np.append(index_rmt_up_j, self.size_i + counter)

                    rmn_left = np.append(rmn_left, self.formula_resistance_rmn_left(i, j))
                    index_rmn_left_i = np.append(index_rmn_left_i, counter)
                    index_rmn_left_j = np.append(index_rmn_left_j, counter-1)

                    rmn_right = np.append(rmn_right, self.formula_resistance_rmn_right(i, j))
                    index_rmn_right_i = np.append(index_rmn_right_i, counter)
                    index_rmn_right_j = np.append(index_rmn_right_j, counter+1)

                    rmt_down = np.append(rmt_down, self.formula_resistance_rmt_down(i, j))
                    index_rmt_down_i = np.append(index_rmt_down_i, counter)
                    index_rmt_down_j = np.append(index_rmt_down_j, counter)

                    r_self = np.append(r_self, 2 * (self.formula_resistance_rmt(i, j))
                                       + self.formula_resistance_rmn(i, j))
                    index_r_self_i = np.append(index_r_self_i, counter)
                    index_r_self_j = np.append(index_r_self_j, counter)
                counter += 1
        values_matrix = (right_term, rmt_up, rmn_left, rmn_right, rmt_down, r_self)
        indexes = (index_rmt_up_i, index_rmt_up_j, index_rmn_left_i, index_rmn_left_j,
                   index_rmn_right_i, index_rmn_right_j, index_rmt_down_i, index_rmt_down_j,
                   index_r_self_i, index_r_self_j)
        return values_matrix, indexes






