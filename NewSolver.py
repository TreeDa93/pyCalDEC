import numpy as np
from numpy.random import rand
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import spsolve

class Solver1:
    """Это класс отвечающий за расчет коэффициентов учета поперечного краевого эффекта.
    Пока в тестовом режиме он называется Boltoncalculate. Он расчитывает только коэффициент Болтона
    """


    def __init__(self, data):
        self.a = data # Массив с описанием данных
        self.L = 1

    def calculate_resistance(self):
        """Рассчитываем взамное магнитное сопротивление вевей
        """
        i = self.a[0, :]
        j = self.a[1, :]
        rmn_left = {(self.a[i,j+1].centerY - self.a[i,j].centerY) /
                    (4 * self.a[i,j].width * self.L) * (1 / self.a[i, j].mat.mu + 1 / self.a[i, j + 1].mat.mu)
                    }
        rmn_right = {(self.a[i+1,j+1].centerY - self.a[i+1,j].centerY) /
                     (4 * self.a[i+1,j].width * self.L) * (1 / self.a[i + 1, j].mat.mu + 1
                                                           / self.a[i + 1, j + 1].mat.mu)
                     }
        self.a = a
        rmt_down = {(self.a[i+1,j].centerX - self.a[i, j].centerX) /
                    (4 * self.a[i, j].hieght * self.L) * (1 / self.a[i, j].mat.mu + 1 / self.a[i, j + 1].mat.mu)
                    }
        rmt_up = {(self.a[i + 1, j + 1].centerX - self.a[i, j + 1].centerX) /
                  (4 * self.a[i, j + 1].hieght * self.L) * (1 / self.a[i, j + 1].mat.mu + 1 / self.a[i, j].mat.mu)
                  }
        rmn = {(self.a[i, j + 1].centerY - self.a[i, j].centerY) /
               (4 * self.a[i, j].width * self.L) * (1 / self.a[i, j].mat.mu + 1 / self.a[i, j + 1].mat.mu) +
               (self.a[i + 1, j + 1].centerY - self.a[i + 1, j].centerY) /
               (4 * self.a[i + 1, j].width * self.L) * (1 / self.a[i + 1, j].mat.mu + 1 / self.a[i + 1, j + 1].mat.mu)
               }
        rmt = {(self.a[i + 1, j].centerX - self.a[i, j].centerX) /
               (4 * self.a[i, j].hieght * self.L) * (1 / self.a[i, j].mat.mu + 1 / self.a[i, j + 1].mat.mu) +
               (self.a[i + 1, j + 1].centerX - self.a[i, j + 1].centerX) /
               (4 * self.a[i, j + 1].hieght * self.L) * (1 / self.a[i, j + 1].mat.mu + 1 / self.a[i, j].mat.mu)
               }
        return rmn_left, rmn_right, rmt_up, rmt_down, rmn, rmt


    def rmn_left(a, i, j, L=1):
        """Рассчитываем взамное магнитное сопротивление вевей
        """
        rmn_left = {(a[i,j+1].centerY - a[i,j].centerY) /
                    (4 * a[i,j].width * L) * (1/a[i,j].mat.mu + 1/a[i,j+1].mat.mu)
                    }
        return rmn_left

    def Rmn_right(a, i, j, L=1):
        rmn_right = {(a[i+1,j+1].centerY - a[i+1,j].centerY) /
                    (4 * a[i+1,j].width * L) * (1/a[i+1,j].mat.mu + 1/a[i+1,j+1].mat.mu)
                     }
        return rmn_right

    def Rmt_down(a, i, j, L=1):
        rmt_down = {(a[i + 1, j].centerX - a[i, j].centerX) /
                    (4 * a[i, j].hieght * L) * (1 / a[i, j].mat.mu + 1 / a[i, j + 1].mat.mu)
                    }
        return rmt_down

    def Rmt_up(a, i, j, L=1):
        rmt_up = {(a[i + 1, j + 1].centerX - a[i, j + 1].centerX) /
                  (4 * a[i, j + 1].hieght * L) * (1 / a[i, j + 1].mat.mu + 1 / a[i, j].mat.mu)
                  }
        return rmt_up
    def Rmn(a, i, j, L=1):
        rmn = {(a[i, j + 1].centerY - a[i, j].centerY) /
               (4 * a[i, j].width * L) * (1 / a[i, j].mat.mu + 1 / a[i, j + 1].mat.mu) +
               (a[i + 1, j + 1].centerY - a[i + 1, j].centerY) /
               (4 * a[i + 1, j].width * L) * (1 / a[i + 1, j].mat.mu + 1 / a[i + 1, j + 1].mat.mu)
               }
        return rmn
    def Rmt(a, i, j, L=1):
        rmt = {(a[i + 1, j].centerX - a[i, j].centerX) /
               (4 * a[i, j].hieght * L) * (1 / a[i, j].mat.mu + 1 / a[i, j + 1].mat.mu) +
               (a[i + 1, j + 1].centerX - a[i, j + 1].centerX) /
               (4 * a[i, j + 1].hieght * L) * (1 / a[i, j + 1].mat.mu + 1 / a[i, j].mat.mu)
               }
        return rmt

    def test_fun(self):
        i = self.a[0, :]
        j = self.a[1, :]
        y1 = i ** 2 + j
        y2 = 2 * i ** 2 + j
        return y1, y2


    def create_dia_matrix(self, y, k):
        size = len(y[0])
        x = dia_matrix((y, k), shape=(size, size))
        #x = dia_matrix((y, k), shape=(size, size)).toarray()
        return x

    def create_f(self, y):
        size = len(y[0])
        f = rand(size)
        return f


    def solve_it(self, z, f):
        z = z.tocsr()
        result = spsolve(z, f)
        return result


"""
Instruction to run the class and to obtain sparse matrix
we should create set of data
import numpy as np 
data = np.array([[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]])

импортируем класс Solver1
from NewSolver import Solver1 

a = Solver1(data)
y = a.test_fun()
k = [0, 1]
z = a.create_dia_matrix(y,k)
f = a.create_f(y)
solution = a.solve_it(z, f)

"""






