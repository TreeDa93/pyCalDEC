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

    def rewrite(self, xNew, yNew):
        """Метод позволяет переписать значение координат"""
        self.x = xNew * self.scale_factor
        self.y = yNew * self.scale_factor
        return self

class Bodies:

    def __init__(self):
        self.geometries = None


    def logicalExpr(self, centerX, centerY, elem):
        """Логическое выражение для определения какому телу пренадлежит клетка"""
        logicalExpr = ((centerX >= elem[0].x) and
                      (centerY >= elem[0].y) and
                      (centerX < elem[1].x) and
                      (centerY < elem[1].y))
        return logicalExpr

    def center(self, axis):
        """
        Function for calculation of cell centers
        :param axis - divided axis
        :return: cell centers
        """
        center = [0 for j in range(len(axis) - 1)]
        for i in range(1, len(axis)):
            center[i - 1] = (axis[i] - axis[i - 1]) / 2 + axis[i - 1]
        return center


    def defineBodies(self, gNodes, mesh):
        """Метод задает имена боди в клетках сетки на основе элементов геометрии и их лейблов"""
        listCenterX = self.center(mesh.axisX)
        listCenterY = self.center(mesh.axisY)
        for indexX, centerX in enumerate(listCenterX):
            for indexY, centerY in enumerate(listCenterY):
                for nodeLabel in gNodes:
                    for elem in gNodes[nodeLabel]:
                        if self.logicalExpr(centerX, centerY, elem):
                            mesh.mesh[indexX][indexY].defineBody(body=nodeLabel)
        return mesh


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


    def __init__(self,body=None, material=None, center=None, initial=None, tp=None, velocity=None, index=None,
                 current=0.0):
        self.mat = material  # словарь с набором физических свойств клетки
        self.center = center  # координаты центра клетки (x, y)
        self.type = tp  # тип клетки  string
        self.init = initial  # координаты нижнего левого угла клетки (x, y)
        self.vel = velocity  # значение скорости в клетке
        self.index = index  # значение индекса клетки
        self.current = current  # значение тока в клетке
        self.body = body   # имя тела в этой  клетке string

    def height(self):
        """Метод вычисляет высоту клетки"""
        height = (self.center.y - self.init.y) * 2
        return height

    def width(self):
        """Метод вычисляет ширину клетки"""
        width = (self.center.x - self.init.x) * 2
        return width

    def mu(self):
        """Метод вычисляет абсолютное значение магнитной проницаемости"""
        mu = self.mat.mur * math.pi * 4 * 10 ** -7
        return mu

    def sigma(self):
        """Метод выдает значение электропроводности клетки"""
        sigma = self.mat.sigma
        return sigma

    def init_x(self):
        """Вызывает координаты x нижнего левого угла клеток  list[][]"""
        return self.init.x

    def init_y(self):
        """Вызывает координаты y нижнего левого угла клеток  list[][]"""
        return self.init.y

    def center_x(self):
        """Вызывает координаты x центра клеток  list[][]"""
        return self.center.x

    def center_y(self):
        """Вызывает координаты y центра клеток  list[][]"""
        return self.center.y

    def defineMat(self, material=''):
        """Вызывает значения метериала по ключу"""
        self.mat = material

    def defineBody(self, body=''):
        """Метод задает имя тела в этой клетке"""
        self.body = body

    def defineCenter(self, center=''):
        """Метод добавляет значение центра клетки"""
        self.center = center

    def defineType(self, type=''):
        "Метод позволяет определить тип клетки"
        self.type = type

    def defineVelocity(self, velocity=''):
        "Метод задает значение скорости"
        self.velocity = velocity

    def defineIndex(self, index=''):
        """Метод задает индекс клетки"""
        self.index = index

    def defineCurrent(self, current):
        """Метод задает значение тока в клетке"""
        self.current = current

    def calculateSquare(self):
        """метод вычисляет площадь клетки"""
        square = self.height()*self.width()
        return square

    def __repr__(self):
        return f"Body({self.body}) Material({self.mat}) Center({self.center}) Initial=({self.init} Type=({self.type}) " \
               f"Velocity=({self.vel})"

    def __str__(self):
        return f"Body({self.body}) Material({self.mat}) Center({self.center}) Initial=({self.init} Type=({self.type}) " \
               f"Velocity=({self.vel})"




class Mesh():


    def __init__(self, label='Mesh1', axisX={}, axisY={}, regionAxisX={}, regionAxisY={}):
        self.label = label  # element label
        self.nodes = {} # dictionary of mesh nodes
        self.axisX = axisX  # global discrete X axis
        self.axisY = axisY  # global discrete Y axis
        self.regionAxisX = regionAxisX  # X grid for a number of regions
        self.regionAxisY = regionAxisY  # Y grid for a number of regions


    def discretRegion(self, startPoint=0, length=0, axis='x', numberElem=None, lengthDisctr=None, label=''):
        """Функция создает узел сетки (словарь) в котором хранится одномерный массив с координатами сетки
        startPoint - начальная точка дискретизации
        length - длина дисритизируемой области
        axis - забает по какой оси проводить дискретизацию
        numberElem - описывает на сколько элементов делить length (если это значение не None то будет
        функция будет дискритизировать область по этмоу параметру)
        lengthDisctr - определяет длину одного элемента дискритизации длины (работает если numberElem = None)
        """
        if numberElem != None:
            grid = [i for i in np.linspace(startPoint, startPoint + length, numberElem, endpoint=True)]
        elif lengthDisctr != None:
            grid = [i for i in np.arange(startPoint, startPoint + length+lengthDisctr/100, lengthDisctr)]

        if axis == 'x':
            self.regionAxisX[label] = np.array(grid)
        elif axis =='y':
            self.regionAxisY[label] = np.array(grid)

    def lastPoint(self, axis='x'):
        """Считает крайнюю точку на заданной оси (x or y) axis"""
        if axis == 'x':
            reqAxis = self.regionAxisX
        elif axis == 'y':
            reqAxis = self.regionAxisY

        lastkey = list(reqAxis.keys())[-1]
        point = reqAxis[lastkey][-1]
        return point

    def lastPointBody(self, axis='x', labelBody=''):
        """Считает крайнюю точку на заданой оси в заданном элементе по labelBody"""
        if axis == 'x':
            reqAxis = self.regionAxisX
        elif axis == 'y':
            reqAxis = self.regionAxisY

        point = reqAxis[labelBody][-1]
        return point

    def center(self, axis):
        """
        Function for calculation of cell centers
        :param axis - divided axis
        :return: cell centers
        """
        center = [0 for j in range(len(axis) - 1)]
        for i in range(1, len(axis)):
            center[i - 1] = (axis[i] - axis[i - 1]) / 2 + axis[i - 1]
        return center

    def discretMultiRegion(self, Listlengthes=[], startPoint=0, axis='x'):
        """
        Функция предназначена для дискритизации заданного множества тел по выбранной оси
        Набор тел задается словарем
        Listlengthes={'название тела 1': {словарь с параметрами дискретизации 1}
                       'название тела 2':{словарь с параметрами дискретизации 2}
                       'название тела N':{словарь с параметрами дискретизации N}}

        словарь с параметрами дискретизации = {'label' : 'название узла сетки'
                                                'length' : значение длины,
                                                'type' : 'numberElem' or 'lengthDisctr'
                                                'pDiscret' : значение numberElem или lengthDisctr}
        startPoint - начало дискретизации
        axis - ооьс на который производится дискритизация
        """
        for index, body in enumerate(Listlengthes):

            if body['type'] == 'numberElem':
                grid = [i for i in np.linspace(startPoint, startPoint + body['length'],
                                               body['pDiscret'], endpoint=True)]
            elif body['type'] == 'lengthDisctr':
                grid = [i for i in np.arange(startPoint, startPoint + body['length']+
                                             body['pDiscret']/100, body['pDiscret'])]

            if axis == 'x':
                self.regionAxisX[body['label']] = np.array(grid)
            elif axis == 'y':
                self.regionAxisY[body['label']] = np.array(grid)

            startPoint += (index+1) * body['length']

    def copyMeshNodes(self, listNodesMesh=[''], CopyStep=1, numberCopyTimes=1, axis='x', prefixNewLabel='copy'):
        """
        Функция копирует существующие узлы сетки и сохраняет предыдущие
        listNodesMesh - лист с labels существующими узлами сетки
        CoptStep - значение расстояния на которое будут пересены копируемые узлы сетки
        numberCopyTimes - сколько раз выполнять эту операцию
        axis - ось вдоль которой происходит копирование
        prefixNewLabel - дополнительное имя, которое будет добавлено к скопированному узлу
        """
        if axis=='x':
            reqAxis = self.regionAxisX
        elif axis=='y':
            reqAxis = self.regionAxisY
        for indexCopy in range(numberCopyTimes):
            for label in listNodesMesh:
                newLabel = str(f'{label}{prefixNewLabel}{indexCopy}')
                reqAxis[newLabel] = reqAxis[label] + (indexCopy+1) * CopyStep

    def builtOneDgrid(self, axis='x'):
        """
        Построить глабальный узел сетки вдоль выбранной оси из локальных узло ( для регионов)

        """
        if axis=='x':
            reqAxis = self.regionAxisX

        if axis=='y':
            reqAxis = self.regionAxisY

        newReqAxis = np.array([])

        for values in reqAxis:
            newReqAxis = np.append(newReqAxis, reqAxis[values][:-1])


        if axis == 'x':
            self.axisX = np.append(newReqAxis, self.lastPoint(axis='x'))
        if axis=='y':
            self.axisY = np.append(newReqAxis, self.lastPoint(axis='y'))

        return newReqAxis

    def createMesh(self):
        """Создает сетку и определяет в ней центры элементов и начальные точки, а также индекс"""
        listCenterX = self.center(self.axisX)
        listCenterY = self.center(self.axisY)
        mesh = [[0 for j in range(len(listCenterY))] for i in range(len(listCenterX))]
        i = 0 # index of the element
        for indexX, centerX in enumerate(listCenterX):
            for indexY, centerY in enumerate(listCenterY):
                mesh[indexX][indexY] = Cell(center=Coord(centerX, centerY),
                    initial=Coord(self.axisX[indexX], self.axisY[indexY]),
                    index=i)
                i+=1
        self.mesh = mesh
        self.sizeX = len(listCenterX)
        self.sizeY = len(listCenterY)
