import numpy as np
import matplotlib.pyplot as plt
import matplotlib
class PostProcessing:

    def __init__(self, solution, a):
        self.solution = solution
        self.size = a.size
        self.size_i = a.size_i
        self.size_j = a.size_j
        self.i = a.size_i+1
        self.j = a.size_j+1

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


    def calculate_magnetic_flux_x(self, reshape_sol):
        mag_flux_x = np.zeros((self.j, self.i), dtype=complex)
        for j in range(self.size_j):
            for i in range(self.size_i-1):
                mag_flux_x[j, i] = reshape_sol[j, i+1] - reshape_sol[j, i]
        return mag_flux_x



    def calculate_magnetic_flux_y(self, reshape_sol):
        mag_flux_y = np.zeros((self.j, self.i), dtype=complex)
        for j in range(self.size_j-1):
            for i in range(self.size_i):
                mag_flux_y[j, i] = reshape_sol[j+1, i] - reshape_sol[j, i]
        return mag_flux_y



