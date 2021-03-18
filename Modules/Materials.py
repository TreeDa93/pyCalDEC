
from Modules.DiscretizationNew import Cell
class Material:
    """
    Class that contains electrical conductivity and magnetic permeability
    """
    def __init__(self, conductivity=0, permeability=0, labelMat='mat1'):
        self.sigma = conductivity
        self.mur = permeability
        self.label = labelMat

    def __repr__(self):
        return f"Sigma={self.sigma},mur={self.mur}"

    def __str__(self):
        return f"Sigma={self.sigma},mur={self.mur}"

    def setMat(self, mesh, material='', body=''):
        for x in range(len(mesh.mesh)):
            for y in range(len(mesh.mesh[0])):
                if mesh.mesh[x][y].body == body:
                    mesh.mesh[x][y].defineMat(material=material)
        return mesh