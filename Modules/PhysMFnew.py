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
        self.a = mesh.mesh  # Массив с описанием данных
        self.L = L
        self.size_i = mesh.sizeX - 1
        self.size_j = mesh.sizeY - 1
        self.size = self.size_i * self.size_j
        self.omega = omega
        self.i = mesh.sizeX
        self.j = mesh.sizeY
        self.size_cell = self.i * self.j

    def definiceCurrent(self, current=1000, body='coil'):
        for x in range(self.i):
            for y in range(self.j):
                if self.a[x][y].body == body:
                    self.a[x][y].defineCurrent(current)

    def defineMatrixis(self):
        """Метод, который инициаоизирует матрицы жесткости и свободных членов и расчитывает их"""
        self.initMatrixResistenceCell()
        self.selfResistenceMatrixCell()
        self.mutualResistenceMatrixCell()

    def selfResistenceMatrixCell(self):
        """Метод расчитывает собственные сопротивления системы и определяет МДС по заданным токам"""
        for i, valueX in enumerate(self.a):
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
        xnum = self.i -1
        ynum = self.j -1
        for i, valueX in enumerate(self.a):
            for j, cell, in enumerate(valueX):
                # down resistence####
                if j == 0:  # первый слой по x (x - fixed)
                    self.rightMatrix[cell.index] += 0
                else:
                    self.weightMatrix[cell.index, cell.index - 1] = self.formulaResistenceXCell(i, j)

                # up resistence######
                if j == ynum:
                    self.rightMatrix[cell.index] += 0
                else:
                    self.weightMatrix[cell.index, cell.index + 1] = self.formulaResistenceXCell(i, j)

                # left resistence ####
                if i == 0:  # первый слой по x (x - fixed)
                    self.rightMatrix[cell.index] += 0
                else:
                    self.weightMatrix[cell.index, cell.index-ynum] = self.formulaResistenceYCell(i, j)

                # right resistence ####
                if i == xnum:  # последний слой по x (x - fixed)
                    self.rightMatrix[cell.index] += 0
                else:
                    self.weightMatrix[cell.index, cell.index+ynum] = self.formulaResistenceYCell(i, j)



    #############################################################################################################
      ############################  !!!    Формулы для расчета сопротивлений   !!!     ####################
    #############################################################################################################



    def formulaResistenceYCell(self, i, j):
        r = self.a[i][j].height() / (self.a[i][j].mu() * self.a[i][j].width() / 2 * self.L)
        return r

    def formulaResistenceXCell(self, i, j):
        r = self.a[i][j].width() / (self.a[i][j].mu() * self.a[i][j].height() / 2 * self.L)
        return r

    def formulaIndTerCell(self, i, j):
        rm = -1j * self.omega * self.a[i][j].calculateSquare() * self.a[i][j].sigma() / self.L
        return rm

    def formula_mmf_coil_cell(self, i, j):
        mmf = self.a[i][j].calculateSquare() * self.a[i][j].current
        return mmf

    def initMatrixResistenceCell(self):
        self.weightMatrix = np.zeros((self.size_cell, self.size_cell), dtype=complex)
        self.rightMatrix = np.zeros(self.size_cell, dtype=complex)






