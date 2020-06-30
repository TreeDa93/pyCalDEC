import numpy as np


class Discrete:
    """The class is intended to create general elements for descreatization of model.
    It’s crutches so that the model is considered. A detailed study of this class is required."""


    def __init__(self, matN =1 , layerN=0, sizeLayer=0):
        self.matN = matN
        self.layer = layerN # лист где по очереди каждый элемент количество слоев
        self.sizeLayer = sizeLayer # пока что константа размер слоя

    def discretka(self, gammaSE, kq):
        """The function is intended to create clasic winding matrix"""
        h = 10 ** -3 * np.array([5, 5, 4, 4, 2, 2, 2, 2, 2, 3, 8, 9, 5, 5], dtype=float)
        mut = np.array([1, 1, 1, 1, 20, 25, 30, 40, 60, 100, 100, 100, 1, 1], dtype=float)
        mun = mut
        gammaotn = np.array([0, 0, 1, 1, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125,
                             0.0125, 0.0125, 0.0125, 0, 0], dtype=float)
        gamma = gammaotn * gammaSE * kq
        return h, mut, mun, gamma, gammaotn



#mut = mut.shape(1, len(mut))



