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






