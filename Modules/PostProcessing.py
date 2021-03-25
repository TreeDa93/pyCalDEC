import numpy as np
import matplotlib
matplotlib.use('TkAgg')   # 'TkAgg' - to run in Tkinter GUI; GTK3Agg


import matplotlib.pyplot as plt
import matplotlib
import matplotlib.tri as tri
import matplotlib.cm as cm
class PostProcessing:

    def __init__(self, solution, a):
        self.solution = solution
        self.size = a.size
        self.size_i = a.size_i
        self.size_j = a.size_j
        self.i = a.size_i+1
        self.j = a.size_j+1
        self.data = a.cells
        self.data_2 = a


    def reshape_data(self, data):
        reshape_sol = data.reshape(self.size_j, self.size_i)
        self.reshape_sol = reshape_sol
        return reshape_sol



    def create_plot(self, reshape_sol):
        x = np.arange(len(reshape_sol))
        y = np.arange(len(reshape_sol[0]))
        z = reshape_sol.real
        fig, ax = plt.subplots()
        ax.pcolormesh(x, y, z)

    def create_simple_plot(self, data, layer=0):
        x = np.arange(len(data[0]))
        fig, ax = plt.subplots()
        ax.plot(x, data[layer].real)


    def create_2d_plot(self, reshape_sol):
        x = np.arange(len(reshape_sol))
        y = np.arange(len(reshape_sol[0]))
        Z = reshape_sol.real
        fig, ax = plt.subplots()
        im = ax.imshow(Z, interpolation='bilinear', cmap=cm.RdYlGn,
                       origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
                       vmax=abs(Z).max(), vmin=abs(Z).min())


    def calculate_magnetic_flux_x(self, reshape_sol):
        mag_flux_x = np.zeros((self.j, self.i), dtype=complex)
        for j in range(self.size_j-1):
            for i in range(self.size_i):
                mag_flux_x[j, i] = (reshape_sol[j+1, i] - reshape_sol[j, i]) / (self.data[i][j].height() * self.data_2.L)
        return mag_flux_x



    def calculate_magnetic_flux_y(self, reshape_sol):
        mag_flux_y = np.zeros((self.j, self.i), dtype=complex)
        for j in range(self.size_j):
            for i in range(self.size_i-1):
                mag_flux_y[j, i] = reshape_sol[j, i+1] - reshape_sol[j, i] / (self.data[i][j].width() * self.data_2.L)
        return mag_flux_y


    def simple_plot(self, reshape_sol, layer=0):
        fig = plt.figure()
        plt.title('layer=' + str(layer))
        plt.xlabel('Номер сектора')
        plt.ylabel("Значение")
        plt.grid()
        plt.plot(reshape_sol[layer].real, "r--")

    def simple_subplots(self, reshape_sol, layer=0):
        fig = plt.figure()
        plt.figure(figsize=(9, 9))
        plt.subplot(2, 1, 1)

        plt.title('layer=' + str(layer))
        plt.xlabel('Номер сектора')
        plt.ylabel("Значение")
        plt.grid()
        plt.plot(reshape_sol[layer].real, "r--")

    def create_pcolor(self, reshape_sol):
        fig = plt.figure()
        plt.ylim()
        plt.pcolor(reshape_sol.real)
        plt.colorbar()


class Variables:

    def __init__(self):
        self.dictVariables = {}

    def calculateMagFlufx(self):
        """Считаем магнитный поток"""
        mag_flux_x = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        mag_flux_y = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        mag_flux_abs = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.mf.j - 1) or i == 0 or i == (self.mf.i - 1):  # первый слой по x (x - fixed)
                    mag_flux_x[i, j] = 0
                    mag_flux_y[i, j] = 0
                    mag_flux_abs[i, j] = 0
                else:
                    mag_flux_x[i, j] = self.dataCounturMatrix[i, j] - self.dataCounturMatrix[i, j + 1]
                    mag_flux_y[i, j] = self.dataCounturMatrix[i + 1, j] - self.dataCounturMatrix[i, j]
                    mag_flux_abs[i, j] = (mag_flux_x[i, j] ** 2 + mag_flux_y[i, j] ** 2) ** (1 / 2)

        self.mag_flux = {'X': mag_flux_x,
                         'Y': mag_flux_y,
                         'absolute': mag_flux_abs}

    def calculateFluxDensity(self):
        """Считаем магнитную индукцию"""
        mag_fluxDensity_x = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        mag_fluxDensity_y = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        mag_fluxDensity_abs = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.mf.j - 1) or i == 0 or i == (self.mf.i - 1):  # первый слой по x (x - fixed)
                    mag_fluxDensity_x[i, j] = 0
                    mag_fluxDensity_y[i, j] = 0
                    mag_fluxDensity_abs[i, j] = 0
                else:
                    mag_fluxDensity_x[i, j] = (self.dataCounturMatrix[i, j] - self.dataCounturMatrix[i, j + 1]) / \
                                              (self.mf.L * cell.height())
                    mag_fluxDensity_y[i, j] = self.dataCounturMatrix[i + 1, j] - self.dataCounturMatrix[i, j] / \
                                              (self.mf.L * cell.width())
                    mag_fluxDensity_abs[i, j] = (mag_fluxDensity_x[i, j] ** 2 + mag_fluxDensity_y[i, j] ** 2) ** (1 / 2)

        self.fluxDensity = {'x': mag_fluxDensity_x,
                            'y': mag_fluxDensity_y,
                            'abs:': mag_fluxDensity_abs}

    def calculateCurrentDensity(self):
        currentDiensty_z = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.mf.j - 1) or i == 0 or i == (self.mf.i - 1):  # первый слой по x (x - fixed)
                    currentDiensty_z[i, j] = 0
                else:
                    dBydx = (self.fluxDensity['y'][i + 1, j] - self.fluxDensity['y'][i, j]) / (cell.width())
                    dBxdy = (self.fluxDensity['x'][i, j + 1] - self.fluxDensity['x'][i, j]) / (cell.height())

                    currentDiensty_z[i, j] = dBydx / cell.mu() - dBxdy / cell.mu()

        self.currentDensity = currentDiensty_z

    def calculateLorentzForce(self):

        avgLorentzForce_x = -0.5 * (self.currentDensity * self.fluxDensity['y'].conjugate()).real
        avgLorentzForce_y = 0.5 * (self.currentDensity * self.fluxDensity['x'].conjugate()).real

        LorentzForce_x = self.currentDensity * self.fluxDensity['y']
        LorentzForce_y = self.currentDensity * self.fluxDensity['x']

        self.lorentzForce = {'avgX': avgLorentzForce_x,
                             'avgY': avgLorentzForce_y,
                             'x': LorentzForce_x,
                             'y': LorentzForce_y}

    def calculateJouleLosses(self):
        JouleLosses = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.mf.j - 1) or i == 0 or i == (self.mf.i - 1):  # первый слой по x (x - fixed)
                    JouleLosses[i, j] = 0
                else:
                    small = 1e-6
                    JouleLosses[i, j] = abs(self.currentDensity[i, j]) ** 2 / (cell.sigma() + small)
        self.JouleLosses = JouleLosses



    def calculate_magnetic_flux_x(self, reshape_sol):
        mag_flux_y = np.zeros((self.i, self.i), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell in enumerate(valueX):
                mag_flux_y[i, j] = self.dataMatrix[i + 1, j] - self.dataMatrix[i, j]
        return mag_flux_y


    def calculateMagFlufx(self):
        """Считаем магнитный поток"""
        mag_flux_x = np.zeros((self.mf.i +  1, self.mf.j + 1), dtype=complex)
        mag_flux_y = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        mag_flux_abs = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.mf.j - 1) or i == 0 or i == (self.mf.i - 1):  # первый слой по x (x - fixed)
                    mag_flux_x[i, j] = 0
                    mag_flux_y[i, j] = 0
                    mag_flux_abs[i, j] = 0
                else:
                    mag_flux_x[i, j] = self.dataCounturMatrix[i, j] - self.dataCounturMatrix[i, j + 1]
                    mag_flux_y[i, j] = self.dataCounturMatrix[i + 1, j] - self.dataCounturMatrix[i, j]
                    mag_flux_abs[i, j] = (mag_flux_x[i, j] ** 2 + mag_flux_y[i, j] ** 2) ** (1 / 2)

        self.mag_flux = {'magFluxX': mag_flux_x,
                         'magFluxY': mag_flux_y,
                         'magFluxABS': mag_flux_abs}

    def calculateFluxDensity(self):
        """Считаем магнитную индукцию"""
        mag_fluxDensity_x = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        mag_fluxDensity_y = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        mag_fluxDensity_abs = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.mf.j - 1) or i == 0 or i == (self.mf.i - 1):  # первый слой по x (x - fixed)
                    mag_fluxDensity_x[i, j] = 0
                    mag_fluxDensity_y[i, j] = 0
                    mag_fluxDensity_abs[i, j] = 0
                else:
                    mag_fluxDensity_x[i, j] = (self.dataCounturMatrix[i, j] - self.dataCounturMatrix[i, j + 1])/\
                                              (self.mf.L * cell.height())
                    mag_fluxDensity_y[i, j] = self.dataCounturMatrix[i + 1, j] - self.dataCounturMatrix[i, j]/\
                                              (self.mf.L * cell.width())
                    mag_fluxDensity_abs[i, j] = (mag_fluxDensity_x[i, j] ** 2 + mag_fluxDensity_y[i, j] ** 2) ** (1 / 2)

        self.fluxDensity = {'x':mag_fluxDensity_x,
                            'y':mag_fluxDensity_y,
                            'abs:':mag_fluxDensity_abs}


    def calculateCurrentDensity(self):
        currentDiensty_z = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.mf.j - 1) or i == 0 or i == (self.mf.i - 1):  # первый слой по x (x - fixed)
                    currentDiensty_z[i, j] = 0
                else:
                    dBydx = (self.fluxDensity['y'][i+1, j] - self.fluxDensity['y'][i, j])/(cell.width())
                    dBxdy = (self.fluxDensity['x'][i,j+1] - self.fluxDensity['x'][i, j])/(cell.height())

                    currentDiensty_z[i, j] = dBydx/cell.mu() - dBxdy/cell.mu()

        self.currentDensity = currentDiensty_z


    def calculateLorentzForce(self):

        avgLorentzForce_x = -0.5 * (self.currentDensity * self.fluxDensity['y'].conjugate()).real
        avgLorentzForce_y = 0.5 * (self.currentDensity * self.fluxDensity['x'].conjugate()).real

        LorentzForce_x = self.currentDensity * self.fluxDensity['y']
        LorentzForce_y = self.currentDensity * self.fluxDensity['x']

        self.lorentzForce = {'avgX':avgLorentzForce_x,
                             'avgY': avgLorentzForce_y,
                             'x':LorentzForce_x,
                             'y':LorentzForce_y}

    def calculateJouleLosses(self):
        JouleLosses = np.zeros((self.mf.i + 1, self.mf.j + 1), dtype=complex)
        for i, valueX in enumerate(self.mf.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0 or j == (self.mf.j - 1) or i == 0 or i == (self.mf.i - 1):  # первый слой по x (x - fixed)
                    JouleLosses[i, j] = 0
                else:
                    small = 1e-6
                    JouleLosses[i,j] = abs(self.currentDensity[i,j]) ** 2 / (cell.sigma() + small)
        self.JouleLosses = JouleLosses



