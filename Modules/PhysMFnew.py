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