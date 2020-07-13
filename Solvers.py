

import numpy as np
from scipy.sparse import dia_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from numpy.random import rand


def create_dia_matrix(y, k):
        size = len(y)
        y = np.asarray(y)
        x = dia_matrix((y, k), shape=(size, size))
        # x = dia_matrix((y, k), shape=(size, size)).toarray()
        return x


def create_csr_marix(data):
        row = np.arange(len(data))
        x = csr_matrix(data, (row, row))
        return x


def create_f(y):
        size = len(y)
        f = rand(size)
        return f


def solve_it(z, f):
        z = z.tocsr()
        result = spsolve(z, f)
        return result
