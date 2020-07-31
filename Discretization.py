import numpy as np
import math


class Coord:
    """
    Class that contains x and y coordinate
    x - x coordinate
    y - y coordinate
    scale_factor - multiplicator of coordinates
    """
    scale_factor = 1.0
    def __init__(self, x=0.0, y=0.0):
        self.x = x * self.scale_factor
        self.y = y * self.scale_factor

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
        vel - value of Cell velocity
        index - Cell text index
    Methods:
        height - Calculates the cell height
        width - Calculates the cell width
        mu - Calculates the cell absolute magnetic permeability
        sigma - Return the cell electrical conductivity
    """
    scale_factor = 1

    def __init__(self, material=Material(0, 0), center=Coord(0, 0), initial=Coord(0, 0), tp="", velocity=0.0, index="",
                 current=0.0):
        self.mat = material
        self.center = center
        self.type = tp
        self.init = initial
        self.vel = velocity
        self.index = index
        self.current = current

    def height(self):
        height = (self.center.y - self.init.y) * 2 * self.scale_factor
        return height

    def width(self):
        width = (self.center.x - self.init.x) * 2 * self.scale_factor
        return width

    def mu(self):
        mu = self.mat.mur * math.pi * 4 * 10 ** -7
        return mu

    def sigma(self):
        sigma = self.mat.sigma
        return sigma

    def init_x(self):
        return self.init.x * self.scale_factor

    def init_y(self):
        return self.init.y * self.scale_factor

    def center_x(self):
        return self.center.x * self.scale_factor

    def center_y(self):
        return self.center.y * self.scale_factor

    def __repr__(self):
        return f"Material({self.mat}) Center({self.center}) Initial=({self.init} Type=({self.type}) " \
               f"Velocity=({self.vel})"

    def __str__(self):
        return f"Material({self.mat}) Center({self.center}) Initial=({self.init} Type=({self.type}) " \
               f"Velocity=({self.vel})"


class Body:
    """
    Class describes the properties of cell in mesh.
    Atributes:
        mat - Material
        vel - value of Body velocity
        type - Body type
        distr - array of initial and end coordinates of the Body parts
        bodies - array of created Bodies
        index - text Body index
    Methods:
        coils - creates coils
        inductor - creates inductor
        rect - creates rectangular
    """
    bodies = []

    def __init__(self, material=Material(0, 0), velocity=0.0):
        self.mat = material
        self.vel = velocity
        self.type = ""
        self.distr = 0
        self.bodies.append(self)
        self.index = ""
        self.current = 0.0

    def __repr__(self):
        return f"Material({self.mat}) Type={self.type} Velocity={self.vel}"

    def __str__(self):
        return f"Material({self.mat}) Type={self.type} Velocity={self.vel}"

    def coils(self, height, width, start=Coord(0, 0), step=Coord(0, 0), amount=1, index="", current=0.0):
        self.type = "coil"
        initial = [start, Coord(start.x+width, start.y+height)]
        if amount == 1:
            self.distr = [initial]
        else:
            self.distr = [0 for j in range(amount)]
            self.distr[0] = initial
            for i in range(1, amount):
                self.distr[i] = [Coord(self.distr[i-1][1].x+step.x, self.distr[i-1][0].y+step.y),
                                 Coord(self.distr[i-1][1].x+step.x+width*bool(step.x),
                                       self.distr[i-1][1].y+step.y+height*bool(step.y))]
        self.index = index
        self.current = current
        return self

    def inductor(self, yoke_height, yoke_width, tooth_height, tooth_width, slot_pitch, slots_number, start=Coord(0, 0),
                 orient=0):
        self.type = "core"
        yoke = [start, Coord(start.x+yoke_width, start.y+yoke_height)]
        initial_tooth = [Coord(start.x, start.y+yoke_height),
                         Coord(start.x+tooth_width, start.y+yoke_height+tooth_height)]
        amount = slots_number+2
        self.distr = [0 for j in range(amount)]
        self.distr[0] = yoke
        self.distr[1] = initial_tooth
        for i in range(2, amount):
            self.distr[i] = [Coord(self.distr[i-1][0].x+slot_pitch, self.distr[i-1][0].y),
                             Coord(self.distr[i-1][0].x+slot_pitch+tooth_width, self.distr[i-1][1].y)]
        return self

    def rect(self, height, width, start=Coord(0, 0)):
        self.type = "other"
        self.distr = [[start, Coord(start.x+width, start.y+height)]]
        return self


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


def discret_X(marginal=(0, 0), domains=((15, 2), (30, 4)), periodic=False, period_amount=0):
    """
    The function for discretization of x coordinate.
    :param marginal - Tuple with patern (marginal width, number of divisions)
    :param domains - Tuple with patern ((the first domain width, number of divisions), (the second domain width,
    :param number of divisions),..., (the Nth domain width, number of divisions))
    :param periodic - key word of domains period
    :param period_amount - number of repetitions
    :return divided axis
    """
    if marginal[0]:
        marginal_size = 2
    else:
        marginal_size = 0

    if periodic:
        domain_size = len(domains) * period_amount + 1
    else:
        domain_size = len(domains)

    axis_init = list(np.zeros(domain_size + marginal_size))

    if marginal[0] and periodic:
        counter = 0
        start = 0
        total_domains = [marginal]
        for i in range(period_amount):
            for k in domains:
                total_domains.append(k)
        total_domains.append(domains[0])
        total_domains.append(marginal)
        for i in total_domains:
            axis_init[counter] = list(np.linspace(start, i[0] + start, i[1], endpoint=False))
            start += i[0]
            counter += 1
        axis_init[-1].append(start)

    if not(marginal[0]) and periodic:
        counter = 0
        start = 0
        total_domains = []
        for i in range(period_amount):
            for k in domains:
                total_domains.append(k)
        total_domains.append(domains[0])
        for i in total_domains:
            axis_init[counter] = list(np.linspace(start, i[0] + start, i[1], endpoint=False))
            start += i[0]
            counter += 1
        axis_init[-1].append(start)

    if marginal[0] and not periodic:
        start = 0
        total_domains = [marginal, domains, marginal]
        counter = 0
        for i in total_domains:
            axis_init[counter] = list(np.linspace(start, i[0] + start, i[1], endpoint=False))
            start += i[0]
            counter += 1
        axis_init[-1].append(start)

    if not marginal[0] and not periodic:
        start = 0
        counter = 0
        for i in domains:
            axis_init[counter] = list(np.linspace(start, i[0] + start, i[1], endpoint=False))
            start += i[0]
            counter += 1
        axis_init[-1].append(start)

    axis = []
    for i in axis_init:
        for j in i:
            axis.append(j)
    return axis


def discret_Y(domains=((10, 1), (15, 2), (30, 4))):
    """
    The function for discretization of y coordinate.
    :param domains - Tuple with patern ((the first domain width, number of divisions),
     (the second domain width, number of divisions),..., (the Nth domain width,
      number of divisions))
    :return divided axis and material for every point
    """
    axis_size = int(np.asarray(domains).sum(0)[1]) + 1
    axis_init = [0 for j in range(len(domains))]
    counter = 0
    start = 0
    for i in domains:
        axis_init[counter] = list(np.linspace(start, i[0]+start, i[1], endpoint=False))
        start += i[0]
        counter += 1
    axis_init[-1].append(start)
    axis = [0 for j in range(axis_size)]
    counter = 0
    for i in axis_init:
        for j in i:
            axis[counter] = j
            counter += 1
    return axis


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


def body_grid(axis_x, axis_y, bodies=()):
    """
    function for forming a grid consisting of classes Cell.
    :param axis_x: divided x axis
    :param axis_y: divided y axis
    :param bodies: array of created bodies
    :return: grid - list of classes Cell
    """
    size_x = len(axis_x)
    size_y = len(axis_y)
    mesh = [[0 for j in range(size_y-1)] for i in range(size_x-1)]
    center_x = center(axis_x)
    center_y = center(axis_y)
    for i in range(size_x-1):
        for j in range(size_y-1):
            mesh[i][j] = Cell(Material(0, 1), Coord(center_x[i], center_y[j]), Coord(axis_x[i], axis_y[j]),
                              "other", velocity=0.0, current=0.0)
            for k in bodies:
                size_distr = len(k.distr)
                for counter in range(size_distr):
                    if ((center_x[i] >= k.distr[counter][0].x) and (center_y[j] >= k.distr[counter][0].y) and
                        (center_x[i] < k.distr[counter][1].x) and (center_y[j] < k.distr[counter][1].y)
                            and (k.type == "coil")):
                        mesh[i][j] = Cell(k.mat, Coord(center_x[i], center_y[j]),
                                          Coord(axis_x[i], axis_y[j]), k.type, k.vel, k.index, current=k.current)
                    elif ((center_x[i] >= k.distr[counter][0].x) and (center_y[j] >= k.distr[counter][0].y) and
                        (center_x[i] < k.distr[counter][1].x) and (center_y[j] < k.distr[counter][1].y)):
                        mesh[i][j] = Cell(k.mat, Coord(center_x[i], center_y[j]),
                                          Coord(axis_x[i], axis_y[j]), k.type, k.vel, current=k.current)
    return mesh


# # Example of code
# # """
# Hy = 5
# Hp = 6
# dz = 1
# d_se = 2
# Bi = 17
# Bz = 2
# Bp = 3
# Q = 3
# marg = 10
# tz = Bp + Bz
#
# air = Material(0, 1)
# copper = Material(100, 1)
# steel = Material(0, 500)
#
# Coord.scale_factor = 1
#
# inductor = Body(steel).inductor(Hy, Bi, Hp, Bz, tz, Q, Coord(marg, 0.0))
# coil = Body(copper).coils(Hp, Bp, Coord(marg+Bz, Hy), Coord(Bz, 0), Q, "phaseA", current=100)
# secondary = Body(copper).rect(d_se, Bi, Coord(marg, Hy+Hp+dz))
#
# axis_x = discret_X((marg, 1), ((Bz, 1), (Bp, 1)), True, Q)
# axis_y = discret_Y(((Hy, 1), (Hp, 1), (dz, 1), (d_se, 1)))
#
# Cell.scale_factor = 0.001
#
# data = body_grid(axis_x, axis_y, Body.bodies)
# """


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



