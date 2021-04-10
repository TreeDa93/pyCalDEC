import numpy as np



# Use ODEINT to solve the differential equations defined by the vector field
import numpy as np
from scipy.integrate import odeint


class HeatTransfer():

    def __init__(self, mesh, label='ht1'):
        self.mesh = mesh
        self.cells = mesh.mesh1D  # Массив с описанием данных
        self.i = mesh.sizeX
        self.j = mesh.sizeY
        self.size_cell = self.i * self.j


    def initTemperature(self, t):
        self.T0 = np.array([])
        for cell in self.cells:
            if cell.bcType == 'bc1':
                self.T0 = np.append(self.T0, self.T(t))
            else:
                self.T0 = np.append(self.T0, 283)
        return self.T0


    def temperatureFieldPDE(self, T, t):
        N = len(T) - 1
        Temperature = np.zeros([len(self.T0)])
        Temperature[0] = self.dTdt(t)
        for i in range(1, N):
            Temperature[i] = (T[i+1]- 2*T[i]+T[i-1])/(0.0001**2)

        Temperature[N] = (2*T[N-1]+ 0.0002*self.dTdx(t) - 2*T[N])/(0.0001**2)

        return Temperature

    def temperatureFieldPDE4(self, T, t):
        Temperature = np.zeros([len(self.T0)])
        Temperature[0] = self.dTdt(t)
        for cell in self.cells:
            print(cell.index)
            if cell.index == 0:
                pass
            elif cell.index != self.cells[-1].index:
                deltaX = self.cells[cell.index + 1].init.x - self.cells[cell.index - 1].init.x
                print(deltaX)
                Temperature[cell.index] = (T[cell.index + 1] - 2 * T[cell.index] + T[cell.index - 1]) / (deltaX**2)

            Temperature[self.cells[-1].index] = (2 * T[cell.index - 1] + 2*cell.width() * self.dTdx(t) - 2 * T[cell.index]) / (
                        cell.width() ** 2)

        return Temperature

    def temperatureFieldPDE2(self, T, t):
        Temperature = np.zeros([len(self.T0)])
        for cell in self.cells:
            if cell.typeApproximation == 'central':
                # central approximation
                if cell.bcType == 'bc1':
                    Temperature[cell.index] = self.dTdt(t)
                if cell.bcType == 'bc2':
                    Temperature[cell.index] = 2 * (self.dTdx(t))
                else:
                    deltaX = self.cells[cell.index + 1].init.x - self.cells[cell.index - 1].init.x
                    Temperature[cell.index] = (T[cell.index + 1] - 2 * T[cell.index] + T[cell.index - 1]) / cell.width()**2
            if cell.typeApproximation == 'left':
                # left approximation
                if cell.bcType == 'bc1':
                    Temperature[cell.index] = self.dTdt(t)
                if cell.bcType == 'bc2':
                    Temperature[cell.index] = 2 * (self.dTdx(t))
                else:
                    deltaX = self.cells[cell.index].init.x - self.cells[cell.index - 1].init.x
                    Temperature[cell.index] = 2 * (T[cell.index-1] - T[cell.index]) / cell.width()**2
            if cell.typeApproximation == 'right':
                # right approximation
                if cell.bcType == 'bc1':
                    Temperature[cell.index] = self.dTdt(t)
                if cell.bcType == 'bc2':
                    Temperature[cell.index] = 2 * (self.dTdx(t))
                else:
                    deltaX = self.cells[cell.index + 1].init.x - self.cells[cell.index].init.x
                    Temperature[cell.index] = 2 * (T[cell.index + 1] - T[cell.index]) / cell.width()**2
        return Temperature

    def temperatureFieldPDE3(self, T, t):
        Temperature = np.zeros([len(self.T0)])
        for cell in self.cells:
            if cell.bcType == 'bc1':
                Temperature[cell.index] = self.dTdt(t)
            if cell.bcType == 'bc2':
                Temperature[cell.index] = (2*T[cell.index-1]+ 0.0002*self.dTdx(t) - 2*T[cell.index])/(0.001**2) #2 * (self.dTdx(t))
            else:
                Temperature[cell.index] = (T[cell.index+1]- 2*T[cell.index]+T[cell.index-1])/(0.001**2)
        return Temperature


    def secondDerivative(self,fun, t, cell):
        """bc1 and bc2"""
        alpha = 0
        betta = 1
        if cell.typeApproximation == 'general':
            # central approximation
            if cell.bcType == 'bc1':
                value = self.dTdt(t)
            if cell.bcType == 'bc2':
                value = 2*(self.dTdx(t))
            else:
                deltaX = self.cells[cell.index + 1].init.x - self.cells[cell.index - 1].init.x
                value = (fun[cell.intex + 1] - 2 * fun[cell.index] + fun[cell.index - 1]) / deltaX

        if cell.typeApproximation == 'left':
            # left approximation
            if cell.bcType == 'bc1':
                value = self.dTdt(t)
            if cell.bcType == 'bc2':
                value = 2*(self.dTdx(t))
            else:
                deltaX = self.cells[cell.index].init.x - self.cells[cell.index - 1].init.x
                value = 2*(fun[cell.index-1] - fun[cell.index]) / deltaX

        if cell.typeApproximation == 'right':
            # right approximation
            if cell.bcType == 'bc1':
                value = self.dTdt(t)
            if cell.bcType == 'bc2':
                value = 2*(self.dTdx(t))
            else:
                deltaX = self.cells[cell.index+1].init.x - self.cells[cell.index].init.x
                value = 2 * (fun[cell.index + 1] - fun[cell.index]) / deltaX

        return value







    def derivative(self, fun, t, cell, bc=None):
        """bc1 and bc2"""

        if cell.typeApproximation == 'general':
            # central approximation
            deltaX = self.cells[cell.index+1].init.x-self.cells[cell.index-1].init.x
            derivative = (fun[cell.index+1]-fun[cell.index-1])/(deltaX)
        if cell.typeApproximation == 'left':
            #left approximation
            deltaX = self.cells[cell.index].init.x - self.cells[cell.index-1].init.x
            derivative = (fun[cell.index] - fun[cell.index-1]) / (deltaX)
        if cell.typeApproximation == 'right':
            # right approximation
            deltaX = self.cells[cell.index+1].init.x - self.cells[cell.index].init.x
            derivative = (fun[cell.index+1] - fun[cell.index]) / (deltaX)
        return derivative

    def dTdx(self, t):
        return 0        #dT/dx

    def T(self, t):
        return 373   #s(t) = T(0, t)

    def dTdt(self, t):
        return 0        #dsdt = dT/dt

    def g(self, x, t):
        return 0        #soure (x, t

    def timeRange(self, stoptime, numpoints):
        self.t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
        return self.t
