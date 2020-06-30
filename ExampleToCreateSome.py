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

