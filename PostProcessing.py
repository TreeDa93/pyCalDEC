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
        self.data = a.a
        self.data_2 = a


    def reshape_data(self, data):
        reshape_sol = data.reshape(self.size_j, self.size_i)
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







