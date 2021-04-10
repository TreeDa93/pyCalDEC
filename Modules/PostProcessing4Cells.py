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
        self.data[label] = self.solution.reshape(self.physic.size_j, self.physic.size_i)
        self.i = self.physic.size_i
        self.j = self.physic.size_j
        return self.data[label]

    def dataMagFluxMatrix(self, label='magneticFlux'):
        mag_flux_x = np.zeros((self.j, self.i), dtype=complex)
        mag_flux_y = np.zeros((self.j, self.i), dtype=complex)
        absFlux = np.zeros((self.j, self.i), dtype=complex)
        for j in range(self.j-1):
            for i in range(self.i):
                mag_flux_x[j, i] = self.data['dataCounturMatrix'][j + 1, i] - self.data['dataCounturMatrix'][j, i]
        for j in range(self.j):
            for i in range(self.i - 1):
                mag_flux_y[j, i] = self.data['dataCounturMatrix'][j, i + 1] - self.data['dataCounturMatrix'][j, i]

        absFlux = (mag_flux_x ** 2 + mag_flux_y ** 2) ** 0.5
        self.data[label] = {'x':mag_flux_x,
                                     'y':mag_flux_y,
                                     'abs':absFlux,
                            'axisX': self.physic.mesh.axisX[1:-1],
                            'axisY': self.physic.mesh.axisY[1:-1]}

        return self.data[label]

    def dataMagFluxMatrix(self, label='magneticFluxDensity'):
        mag_flux_x = np.zeros((self.j, self.i), dtype=complex)
        mag_flux_y = np.zeros((self.j, self.i), dtype=complex)
        absFlux = np.zeros((self.j, self.i), dtype=complex)
        for j in range(self.j-1):
            for i in range(self.i):
                mag_flux_x[j, i] = self.data['dataCounturMatrix'][j + 1, i] - self.data['dataCounturMatrix'][j, i]
        for j in range(self.j):
            for i in range(self.i - 1):
                mag_flux_y[j, i] = self.data['dataCounturMatrix'][j, i + 1] - self.data['dataCounturMatrix'][j, i]

        absFlux = (mag_flux_x ** 2 + mag_flux_y ** 2) ** 0.5
        self.data[label] = {'x':mag_flux_x,
                                     'y':mag_flux_y,
                                     'abs':absFlux,
                            'axisX': self.physic.mesh.axisX[1:-1],
                            'axisY': self.physic.mesh.axisY[1:-1]}

        return self.data[label]

    def calculateCurrentDensity(self, label='CurrentDensity'):
        """Вычислить плотность тока, data - MFFluxDensityMatrix"""
        currentDiensty_z = np.zeros((self.physic.i + 1, self.physic.j + 1), dtype=complex)
        phiX = self.data['magneticFlux']['x']
        phiY = self.data['magneticFlux']['y']
        for i, valueX in enumerate(self.physic.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.physic.j - 1) or i == 0 or i == (self.physic.i - 1):  # первый слой по x (x - fixed)
                    currentDiensty_z[i, j] = 0
                else:
                    dBydx = (phiY[i + 1, j] - phiY[i, j]) / (cell.width())
                    dBxdy = (phiX[i, j + 1] - phiX[i, j]) / (cell.height())

                    currentDiensty_z[i, j] = dBydx / cell.mu() - dBxdy / cell.mu()
