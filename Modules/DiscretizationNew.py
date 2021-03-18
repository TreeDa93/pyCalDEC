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
        self.x = xNew * self.scale_factor
        self.y = yNew * self.scale_factor
        return self

class Bodies:

    def __init__(self):
        self.geometries = None


    def logicalExpr(self, centerX, centerY, elem):

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
        self.mat = material
        self.center = center
        self.type = tp
        self.init = initial
        self.vel = velocity
        self.index = index
        self.current = current
        self.body = body

    def height(self):
        height = (self.center.y - self.init.y) * 2
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

    def init_x(self):
        return self.init.x

    def init_y(self):
        return self.init.y

    def center_x(self):
        return self.center.x

    def center_y(self):
        return self.center.y

    def defineMat(self, material=''):
        self.mat = material

    def defineBody(self, body=''):
        self.body = body

    def defineCenter(self, center=''):
        self.center = center

    def defineType(self, type=''):
        self.type = type

    def defineVelocity(self, velocity=''):
        self.velocity = velocity

    def defineIndex(self, index=''):
        self.index = index

    def defineCurrent(self, current):
        self.current = current

    def calculateSquare(self):
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
        self.label = label
        self.nodes = {}
        self.axisX = axisX
        self.axisY = axisY
        self.regionAxisX = regionAxisX
        self.regionAxisY = regionAxisY


    def discretRegion(self, startPoint=0, length=0, axis='x', numberElem=None, lengthDisctr=None, label=''):
        if numberElem != None:
            grid = [i for i in np.linspace(startPoint, startPoint + length, numberElem, endpoint=True)]
        elif lengthDisctr != None:
            grid = [i for i in np.arange(startPoint, startPoint + length+lengthDisctr/100, lengthDisctr)]

        if axis == 'x':
            self.regionAxisX[label] = np.array(grid)
        elif axis =='y':
            self.regionAxisY[label] = np.array(grid)

    def lastPoint(self, axis='x'):
        """Считает крайнюю точку во всей сетке по заданной оси (x or y) axis"""
        if axis == 'x':
            reqAxis = self.regionAxisX
        elif axis == 'y':
            reqAxis = self.regionAxisY

        lastkey = list(reqAxis.keys())[-1]
        point = reqAxis[lastkey][-1]
        return point

    def lastPointBody(self, axis='x', labelBody=''):
        """Считает крайнюю точку в заданном элементе по labelBody"""
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
        """Listlengthes={'body1':{dictionaries with parameters}
                        'body2':{dictionaries with parameters}}
        dictionaries with paramiters = {'label' : 'name of label'
                                        'length' : L,
                                        'type' : 'numberElem' or 'lengthDisctr'
                                        'pDiscret' : paramValue}
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

        if axis=='x':
            reqAxis = self.regionAxisX
        elif axis=='y':
            reqAxis = self.regionAxisY
        for indexCopy in range(numberCopyTimes):
            for label in listNodesMesh:
                newLabel = str(f'{label}{prefixNewLabel}{indexCopy}')
                reqAxis[newLabel] = reqAxis[label] + (indexCopy+1) * CopyStep

    def builtOneDgrid(self, axis='x'):
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
        """Создает сетку и определяет в ней центы элементов и начальные точки, а также индекс"""
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