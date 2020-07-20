

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


def create_csr_matrix(data, row, column):
        x = csr_matrix((data, (row, column)))
        return x

def solve_it(r, mmf):
        r = csr_matrix(r)
        result = spsolve(r, mmf)
        return result
