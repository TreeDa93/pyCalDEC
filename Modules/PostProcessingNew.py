import numpy as np
import matplotlib
matplotlib.use('TkAgg')   # 'TkAgg' - to run in Tkinter GUI; GTK3Agg
import pandas as pd
from scipy.interpolate import interp2d, NearestNDInterpolator

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.tri as tri
import matplotlib.cm as cm
class PostProcessing:

    def __init__(self, solution, physic):
        self.solution = solution
        self.physic = physic
        self.data = {}
        self.variables = {}

    def dataCounturMatrix(self, label='dataCounturMatrix'):
        """ФУнкция прнвращает стобец решений в матрицу согласно координатам, x - i; y - j"""
        dataCounturMatrix = np.zeros((self.physic.i, self.physic.j), dtype=complex)
        for i, valueX in enumerate(self.physic.cells):
            for j, cell in enumerate(valueX):
                dataCounturMatrix[i][j] = self.solution[cell.index]

        self.data[label] = {'value': dataCounturMatrix,
                             'centerX': self.physic.mesh.centerX,
                             'centerY': self.physic.mesh.centerY}

        return self.data[label]

    def dataCounturColumn(self, label='dataCounterColumn'):
        """Создаем набор данных в виде стобца"""
        centerX = np.array([]) #np.zeros(self.mf.size_cell)
        centerY = np.array([]) #np.zeros(self.mf.size_cell)
        indexX  = np.array([]) #np.zeros(self.mf.size_cell)
        indexY = np.array([]) #np.zeros(self.mf.size_cell)
        index = np.array([]) #np.zeros(self.mf.size_cell)
        value = np.array([])

        for i, valueX in enumerate(self.physic.cells):
            for j, cell in enumerate(valueX):
                centerX = np.append(centerX, cell.center.x)
                centerY = np.append(centerY, cell.center.y)
                indexX = np.append(indexX, i)
                indexY = np.append(indexY, j)
                index = np.append(index, cell.index)
                value = np.append(value, self.solution[cell.index])

        self.data[label] = {'centerX': centerX,
                    'centerY': centerY,
                    'indexX': indexX,
                    'indexY':indexY,
                    'index': index,
                    'value': value,
                    }
        return self.data[label]

    def reshapeToColumn(self, axisMesh=('axisX','axisY'), dataMatrix=None, label='rehapeedDataToColumn'):
        """Функция переводит данные из матрицы в столбец"""
        if axisMesh == ('axisX','axisY'):
            axisX = self.physic.mesh.axisX
            axisY = self.physic.mesh.axisY
        else:
            axisX = axisMesh[0]
            axisY = axisMesh[1]

        self.data[label] ={'X':axisX,
                     'Y':axisY,
                     'value': dataMatrix.reshape((1, dataMatrix.size))}
        return self.data[label]


    def createGrid(self, label='grid'):
        """Функция создает сетку по центрам клеток относительно координат и индектосов"""

        gridXY = np.meshgrid(self.physic.mesh.centerX, self.physic.mesh.centerY)
        gridIndex = np.meshgrid(range(self.physic.mesh.sizeX), range(self.physic.mesh.sizeY))
        self.data[label] = {'gridXY': gridXY,
                            'gridIndex': gridIndex}
        return self.data[label]

    def calculateMfFluxMatrix(self, data='', label='MfFluxMatrix'):
        """Считаем магнитный поток"""
        mag_flux_x = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        mag_flux_y = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        mag_flux_abs = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        for i, valueX in enumerate(self.physic.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.physic.j - 1) or i == 0 or i == (self.physic.i - 1):  # первый слой по x (x - fixed)
                    mag_flux_x[i, j] = 0
                    mag_flux_y[i, j] = 0
                    mag_flux_abs[i, j] = 0
                else:
                    mag_flux_x[i, j] = data[i, j] - data[i, j + 1]
                    mag_flux_y[i, j] = data[i + 1, j] - data[i, j]
                    mag_flux_abs[i, j] = (mag_flux_x[i, j] ** 2 + mag_flux_y[i, j] ** 2) ** (1 / 2)

        self.data[label] = {'X': mag_flux_x,
                         'Y': mag_flux_y,
                         'absolute': mag_flux_abs,
                         'axisX': self.physic.mesh.axisX,
                         'axisY': self.physic.mesh.axisY}
        return  self.data[label]

    def calculateFluxDensity(self, data='', label='FluxDensity'):
        """Считаем магнитную индукцию data - matrix of MF counter"""
        mag_fluxDensity_x = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        mag_fluxDensity_y = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        mag_fluxDensity_abs = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        for i, valueX in enumerate(self.physic.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.physic.j - 1) or i == 0 or i == (self.physic.i - 1):
                    mag_fluxDensity_x[i, j] = 0
                    mag_fluxDensity_y[i, j] = 0
                    mag_fluxDensity_abs[i, j] = 0
                else:
                    mag_fluxDensity_x[i, j] = (data[i, j] - data[i, j + 1]) / \
                                              (self.physic.L * cell.height())
                    mag_fluxDensity_y[i, j] = data[i + 1, j] - data[i, j] / \
                                              (self.physic.L * cell.width())
                    mag_fluxDensity_abs[i, j] = (mag_fluxDensity_x[i, j] ** 2 + mag_fluxDensity_y[i, j] ** 2) ** (1 / 2)

        self.data[label] = {'X': mag_fluxDensity_x,
                            'Y': mag_fluxDensity_y,
                            'absolute': mag_fluxDensity_abs,
                            'axisX': self.physic.mesh.axisX,
                            'axisY': self.physic.mesh.axisY}

    def calculateCurrentDensity(self, data='', label='CurrentDensity'):
        """Вычислить плотность тока, data - MFFluxDensityMatrix"""
        currentDiensty_z = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        for i, valueX in enumerate(self.physic.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.physic.j - 1) or i == 0 or i == (self.physic.i - 1):  # первый слой по x (x - fixed)
                    currentDiensty_z[i, j] = 0
                else:
                    dBydx = (data['Y'][i + 1, j] - data['Y'][i, j]) / (cell.width())
                    dBxdy = (data['X'][i, j + 1] - data['Y'][i, j]) / (cell.height())

                    currentDiensty_z[i, j] = dBydx / cell.mu() - dBxdy / cell.mu()

        self.data[label] = {'Z' : currentDiensty_z,
                            'axisX': self.physic.mesh.axisX,
                            'axisY': self.physic.mesh.axisY}

    def calculateLorentzForce(self, dataCurrent='', dataFluxDensity='', label='MatrixLorentzForce'):

        avgLorentzForce_x = -0.5 * (dataCurrent['Z'] * dataFluxDensity['Y'].conjugate()).real
        avgLorentzForce_y = 0.5 * (dataCurrent['Z'] * dataFluxDensity['X'].conjugate()).real

        LorentzForce_x = dataCurrent['Z'] * dataFluxDensity['Y']
        LorentzForce_y = dataCurrent['Z'] * dataFluxDensity['X']

        self.data[label] = {'avgX': avgLorentzForce_x,
                             'avgY': avgLorentzForce_y,
                             'X': LorentzForce_x,
                             'Y': LorentzForce_y,
                            'axisX': self.physic.mesh.axisX,
                            'axisY': self.physic.mesh.axisY
                            }

    def calculateJouleLosses(self, dataCurrentDensity='', label='MatrixJoule'):
        JouleLosses = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        for i, valueX in enumerate(self.physic.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.physic.j - 1) or i == 0 or i == (self.physic.i - 1):
                    JouleLosses[i, j] = 0
                else:
                    small = 1e-6
                    JouleLosses[i, j] = abs(dataCurrentDensity['Z'][i, j]) ** 2 / (cell.sigma() + small)

        self.data[label] = {'value': JouleLosses,
                            'axisX': self.physic.mesh.axisX,
                            'axisY': self.physic.mesh.axisY
                            }

    def interpolateCounturMatrix(self, kind='linear'):
        """"‘linear’, ‘cubic’, ‘quintic’"""
        m = len(self.centerX)
        s = (m - (2 * m) ** 0.5, m + (2 * m) ** 0.5)
        self.funCountur_real1 = interp2d(self.centerX, self.centerX, self.dataCounturMatrix.real, kind=kind)
        self.funCountur_imag1 = interp2d(self.centerX, self.centerX, self.dataCounturMatrix.imag, kind=kind)
    ###################################################################
    ###############################################################3

    def plot2DIndex(self, data, intertype=None, figsize=(9, 6)):
        """
        Type interpolation
        [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
           'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
           'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']"""


        fig, ax= plt.subplots(figsize=figsize)
        ax.imshow(data.real, interpolation=intertype)

    def plotCounterSurface(self,y, x, z, levels=20):

        fig, ax = plt.subplots()
        ax.contour(y, x, z.real, levels=levels, linewidths=0.5, colors='k')
        cntr1 = ax.contourf(y, x, z.real, levels=levels, cmap="RdBu_r")
        fig.colorbar(cntr1, ax=ax)
        ax.set_title('')

    def plot2DCoord(self, coordList=['x', 'y'], valueList=[], numberPointsXY=['xNumber', 'yNumber'],
                       iterpType='linear'):

        f = interp2d(coordList[0], coordList[1], valueList.real, kind=iterpType)
        xNew = np.linspace(min(coordList[0]), max(coordList[0]), numberPointsXY[0])
        yNew = np.linspace(min(coordList[0]), max(coordList[0]), numberPointsXY[1])
        Z = f(xNew, yNew)
        fig = plt.imshow(Z, extent=[min(coordList[0]),max(coordList[0]),min(coordList[1]),max(coordList[1])],
           origin="lower")
        plt.scatter(coordList[0], coordList[1], 400, facecolors='none')





###################################################################
###############################################################3

class Variabls:

    def __init__(self):

        self.listVariabls = {'fluxMatrix': self.calculateMfFluxMatrix,
                              'densitFluxMatrix' : self.calculateFluxDensity,
                                'currentDensityMatrix' : self.calculateCurrentDensity,
                              'LorentzForce': self.calculateLorentzForce,
                             'JouleHeating': self.JouleLosses}

    def derivativeX(self, i, j, data=''):
        value = data[i]

        return

