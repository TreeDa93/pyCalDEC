import numpy as np

class MagneticField2:
    """"This class is intended for creating matrix of linear equations based on
    physical equation of magnetic filed and circuit theory
    Main parameters of the class is the data of models consisting to "a" massiv
    a is the input massive of data
    size_j is the size or number of layers along y component minus one. Other word is the number of resistances along
    y component
    size_i is the number of elements at each layer along x component minus one. Other word is the number of resistances
    along x component.
    L is width of model directing to plane of the model
    size is the number of unknown of linear equation system.

    """

    def __init__(self, mesh, omega=0, L=0.5, label='mf'):
        self.mesh = mesh
        self.cells = mesh.mesh  # Массив с описанием данных
        self.paramertrs = ParametrsMF(omega, L)
        self.L = L
        self.omega = omega
        self.i = mesh.sizeX
        self.j = mesh.sizeY
        self.size_cell = self.i * self.j

    def definiceCurrent(self, current=1000, body='coil'):
        for x in range(self.mesh.sizeX):
            for y in range(self.mesh.sizeY):
                if self.cells[x][y].body == body:
                    self.cells[x][y].defineCurrent(current)

    def defineMatrixis(self):
        """Метод, который инициаоизирует матрицы жесткости и свободных членов и расчитывает их"""
        self.initMatrixResistenceCell()
        self.selfResistenceMatrixCell()
        self.mutualResistenceMatrixCell()

    def selfResistenceMatrixCell(self):
        """Метод расчитывает собственные сопротивления системы и определяет МДС по заданным токам"""
        for i, valueX in enumerate(self.cells):
            for j, cell, in enumerate(valueX):
                # расчитываем главную диагональ
                self.weightMatrix[cell.index, cell.index] = self.formulaIndTerCell(i, j) +\
                                                                  2 * self.formulaResistenceXCell(i, j) + \
                                                                  2 * self.formulaResistenceYCell(i, j)
                # расчитываем вектор свободных членов на основе МДС
                self.rightMatrix[cell.index] += self.formula_mmf_coil_cell(i,j)

        return self.weightMatrix

    def mutualResistenceMatrixCell(self):
        """Метод расчитывает взаимные сопротивления системы и определяет ГУ"""
        xnum = self.mesh.sizeX -1
        ynum = self.mesh.sizeY -1
        for i, valueX in enumerate(self.cells):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0:  # первый слой по x (x - fixed)
                    self.rightMatrix[cell.index] += 0
                else:
                    self.weightMatrix[cell.index, cell.index - 1] = -self.formulaResistenceXCell(i, j-1)

                # up resistence######
                if j == ynum:
                    self.rightMatrix[cell.index] += 0
                else:
                    self.weightMatrix[cell.index, cell.index + 1] = -self.formulaResistenceXCell(i, j+1)

                # left resistence ####
                if i == 0:  # первый слой по x (x - fixed)
                    self.rightMatrix[cell.index] += 0
                else:
                    self.weightMatrix[cell.index, cell.index-ynum] = -self.formulaResistenceYCell(i-1, j)

                # right resistence ####
                if i == xnum:  # последний слой по x (x - fixed)
                    self.rightMatrix[cell.index] += 0
                else:
                    self.weightMatrix[cell.index, cell.index+ynum] = -self.formulaResistenceYCell(i+1, j)



    #############################################################################################################
      ############################  !!!    Формулы для расчета сопротивлений   !!!     ####################
    #############################################################################################################



    def formulaResistenceYCell(self, i, j):
        r = self.cells[i][j].height() / (self.cells[i][j].mu() * self.cells[i][j].width() / 2 * self.L)
        return r

    def formulaResistenceXCell(self, i, j):
        r = self.cells[i][j].width() / (self.cells[i][j].mu() * self.cells[i][j].height() / 2 * self.L)
        return r

    def formulaIndTerCell(self, i, j):
        rm = -1j * self.omega * self.cells[i][j].calculateSquare() * self.cells[i][j].sigma() / self.L
        return rm

    def formula_mmf_coil_cell(self, i, j):
        mmf = self.cells[i][j].calculateSquare() * self.cells[i][j].current
        return mmf

    def initMatrixResistenceCell(self):
        self.weightMatrix = np.zeros((self.size_cell, self.size_cell), dtype=complex)
        self.rightMatrix = np.zeros(self.size_cell, dtype=complex)



class ParametrsMF:

    def __init__(self, L, omega):
        self.L = L
        self.omega = omega




