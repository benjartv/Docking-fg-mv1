__author__ = 'Felipe'

# --------------------------- Gene Object ---------------------------
class Gene:

    def __init__(self):
        self.id = 0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.sph_theta = 0.0
        self.sph_phi = 0.0
        self.theta = 0.0
        self.score = 0.0
        self.probScore = 0.0
        self.rotateAtoms = []