import numpy as np

h = np.array([])

matN = 5
layerN1 = 5

hcell = 2
for i in range(matN):
    for N in range(layerN1):
        h = 1


grid = []
#Create a 8x8 grid
for row in range(8):
    grid.append([])
    for col in range(8):
        grid[row].append("0")


"""
SOLID принципы
1. принцип единственной ответсвенности класса
KISS =  keep it simple stupid

иногда не надо классы, лучше простота.
"""

class Ball:
    def __init__(self):
        R = randint(20, 50)
        x =

import numpy as np

a = np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])

b = np.array([[[1, 2, 3, 4],
              [1, 2, 3, 4],
              [1, 2, 3, 4],
              [1, 2, 3, 4]],

             [[1, 2, 3, 4],
              [1, 2, 3, 4],
              [1, 2, 3, 4],
              [1, 2, 3, 4]]])

"двухмерный массив"
for j in range(5):
    d2 = []
    for i in range(5):
        d2.append(0)
    d1.append(d2)

K = []
K1 = np.empty()
K1 = np.zeros((3,3), dtype = np.complex)
K3 = np.eye(5)
K3 = np.eye(5, dtype = np.complex)
K4 = np.empty(5 , dtype = np.complex)

K5 = np.zeros((10,10))
for i in range(0, 10):
    for j in range(0,10):
        K5[i][j] = 1

h1 = np.array([1, 2 ,3, 5, 6])
R1 = np.array([], dtype= complex)
R2 = np.array([])
R3 = np.zeros((2,2), dtype=complex)
R3 = np.zeros((10,1))
for i in range(len(R3)):
    R3[i] = 1+R3[i-1]*10

R4 = np.zeros((10,10))
for i in range(len(R4)):
    for j in range(len(R4)):
        if i == j:
            R4[i][j] = 1
"двухмерный массив"
R5 = np.zeros((10,10))
for i in range(len(R5)):
    for j in range(len(R5)):
            R5[j][j] = 1

"расширение массива"
ret = np.array([])
for i in range(100000):
    tmp = gi(i)
    ret = np.append(ret, np.zeros(len(tmp)))
    ret = np.append(ret, np.ones(fixed_length))
"трехмерный массив"
R5 = np.zeros((10,10,10))
for i in range(len(R5)):
    for j in range(len(R5)):
            R5[i][j][j] = 1

"трехмерный массив"
R5 = np.zeros((10,10,10))
for i in range(len(R5)):
    for j in range(len(R5)):
            R5[i,j,j] = 1


T1 = np.array([2,2])
T2 = np.array((2,2))
#############
A = np.array( [[1,1,0], [0,1,1], [1,2,3]])
y = np.array([1,4,6])
A_inv = np.linalg.inv(A)
x = np.matmul(A_inv, y)

print(x)

##############
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt

np.random.seed(0)

S0 = 100.0
N_days = 30
r =0.01
sigma = 0.30
r = r/ 250.
dt = 1.0
sigma = sigma/np.sqrt(252.0)

epsilon = np.random.normal( size = (N_days) )
Lambda = r* dt * sigma * np.sqrt(dt) * epsilon

ones = -np.ones( (N_days+1))
ones[0] = 1

d = [Lambda + 1, ones]
K = [-1, 0]

M = scipy.sparse.diags(d, K, format = 'csc')
p = np.zeros( (N_days + 1, 1) )
p[0] = S0

s = scipy.sparse.linalg.spsolve(M, p)

plt.spy(M)

b = np.array([[2, 3], [5, 7]])
np.append(b,[[11, 13]], axis=0)
np.append(b, [[11], [13]], axis = 1)
np.append(b,[[11, 13]])

T1 = np.array([[1, 2, 3], [3, 4, 5], [6, 7, 8]])
T2 = np.array([[1, 2, 3], [3, 4, 5], [6, 7, 8]])
T3 = np.append(T1,T2, axis = 0)

