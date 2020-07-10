import numpy as np
import math


class Coord:
    """
    Class that contains x and y coordinate
    """
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"X={self.x};Y={self.y}"

    def __str__(self):
        return f"X={self.x};Y={self.y}"


class Material:
    """
    Class that contains electrical conductivity and magnetic permeability
    """
    def __init__(self, conductivity=0, permeability=0):
        self.sigma = conductivity
        self.mur = permeability

    def __repr__(self):
        return f"Sigma={self.sigma},mur={self.mur}"

    def __str__(self):
        return f"Sigma={self.sigma},mur={self.mur}"


class Cell:
    """
    Class describes the properties of cell in mesh.
    Atributes:
        mat - Material
        center - Cell center
        type - Cell type
        init - Cell initial coordinates
    Methods:
        height - Calculates the cell height
        width - Calculates the cell width
        mu - Calculates the cell absolute magnetic permeability
        sigma - Return the cell electrical conductivity
    """
    def __init__(self, material=Material(0, 0), center=Coord(0, 0), initial=Coord(0, 0), tp=0):
        self.mat = material
        self.center = center
        self.type = tp
        self.init = initial

    def height(self):
        height = (self.center.y - self.init.y)*2
        return height

    def width(self):
        width = (self.center.x - self.init.x) * 2
        return width

    def mu(self):
        mu = self.mat.mur * math.pi * 4 * 10 ** -7
        return mu

    def sigma(self):
        sigma = self.mat.sigma
        return sigma

    def __repr__(self):
        return f"Material({self.mat}) Center({self.center}) Initial=({self.init})"

    def __str__(self):
        return f"Material({self.mat}) Center({self.center}) Initial=({self.init})"


def discret_x(domains=((10, 1), (15, 2), (30, 4))):
    """
    The function for discretization of x coordinate.
    :param domains - Tuple with patern ((the first domain width, number of divisions), (the second domain width,
     number of divisions),..., (the Nth domain width, number of divisions))
    :return divided axis
    """
    # Function for discretization
    axis_init = list(np.zeros(len(domains)))
    counter = 0
    start = 0
    for i in domains:
        axis_init[counter] = list(np.linspace(start, i[0]+start, i[1], endpoint=False))
        start += i[0]
        counter += 1
    axis_init[-1].append(start)
    axis = []
    for i in axis_init:
        for j in i:
            axis.append(j)
    #axis = np.asarray(axis)
    return axis


def discret_y(domains=((10, 1, "Air"), (15, 2, "Copper"), (30, 4, "Air"))):
    """
    The function for discretization of y coordinate.
    :param domains - Tuple with patern ((the first domain width, number of divisions, the first domain material),
     (the second domain width, number of divisions, the second domain material),..., (the Nth domain width,
      number of divisions, the Nth domain material))
    :return divided axis and material for every point
    """
    axis_init = [[0 for j in range(len(domains))] for i in range(2)]
    counter = 0
    start = 0
    for i in domains:
        axis_init[0][counter] = list(np.linspace(start, i[0]+start, i[1], endpoint=False))
        axis_init[1][counter] = [i[2] for j in range(i[1])]
        start += i[0]
        counter += 1
    axis_init[0][-1].append(start)
    axis_init[1][-1].append(domains[-1][-1])
    element_sum = 7
    axis = []
    for i in axis_init[0]:
        for j in i:
            axis.append(j)
    mat = []
    for i in axis_init[1]:
        for j in i:
            mat.append(j)
    return axis, mat


def center(axis):
    """
    Function for calculation of cell centers
    :param axis - divided axis
    :return: cell centers
    """
    center = [0 for j in range(len(axis)-1)]
    for i in range(1, len(axis)):
        center[i-1] = (axis[i] - axis[i-1]) / 2 + axis[i-1]
    return center


def grid(axis_x, axis_y, material_y, marginal_width):
    """
    function for forming a grid consisting of classes Cell.
    :param axis_x: divided x axis
    :param axis_y: divided y axis
    :param material_y: materials for every point in y axis
    :param marginal_width: marginal width of computational domain
    :return: grid - list of classes Cell
    """
    mesh = [[0 for j in range(len(axis_y)-1)] for i in range(len(axis_x)-1)]
    center_x = center(axis_x)
    center_y = center(axis_y)
    for i in range(len(axis_x)-1):
        for j in range(len(axis_y)-1):
            if axis_x[i] >= marginal_width and axis_x[-1]-marginal_width:
                mesh[i][j] = Cell(material_y[j], Coord(center_x[i], center_y[j]), Coord(axis_x[i], axis_y[j]))
            else:
                mesh[i][j] = Cell(Material(0, 1), Coord(center_x[i], center_y[j]), Coord(axis_x[i], axis_y[j]))
    return mesh


# class Discrete:
#     """The class is intended to create general elements for descreatization of model.
#     It’s crutches so that the model is considered. A detailed study of this class is required."""
#
#
#     def __init__(self, matN =1 , layerN=0, sizeLayer=0):
#         self.matN = matN
#         self.layer = layerN # лист где по очереди каждый элемент количество слоев
#         self.sizeLayer = sizeLayer # пока что константа размер слоя
#
#     def discretka(self, gammaSE, kq):
#         """The function is intended to create clasic winding matrix"""
#         h = 10 ** -3 * np.array([5, 5, 4, 4, 2, 2, 2, 2, 2, 3, 8, 9, 5, 5], dtype=float)
#         mut = np.array([1, 1, 1, 1, 20, 25, 30, 40, 60, 100, 100, 100, 1, 1], dtype=float)
#         mun = mut
#         gammaotn = np.array([0, 0, 1, 1, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125,
#                              0.0125, 0.0125, 0.0125, 0, 0], dtype=float)
#         gamma = gammaotn * gammaSE * kq
#         return h, mut, mun, gamma, gammaotn
#
#
#
# #mut = mut.shape(1, len(mut))



