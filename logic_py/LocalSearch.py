__author__ = 'Felipe_Gonzalez_Foncea'

from utils_py import Methods
import rosetta
import random
import math
import copy


# -------------------------- Local Search ---------------------------
class LocalSeach(object):

    # ------------------------------ Init -------------------------------
    def __init__(self, temp, tempMin, tempAlpha, innerLoop,  scoreFx, protein, spaceLen, spaceCenter, poseLigand):
        self.__temp = temp
        self.__tempMin = tempMin
        self.__tempAlpha = tempAlpha
        self.__innerIteration = innerLoop
        self.__scoreFx = scoreFx
        self.__protein = protein
        self.__spaceLen = spaceLen
        self.__spaceCenter = spaceCenter
        self.__ligand = None
        self.__pose_ligand = poseLigand


    # -------------------------- Init Process ---------------------------
    def initLocalProcess(self, selectedPop, ligand):
        cont = 0
        self.__ligand = ligand
        T = self.__temp
        oldSol = copy.deepcopy(selectedPop)
        oldSol.score = self.getScoring(oldSol)
        while( T > self.__tempMin):
            j = 1
            while(j <= self.__innerIteration):
                newSol = self.getNeighbor(oldSol)
                if(self.acceptance(oldSol, newSol, T) > random.random()):
                    oldSol = copy.deepcopy(newSol)
                j += 1
            T *= self.__tempAlpha
            cont= cont +1
        selectedPop = copy.deepcopy(oldSol)
        return selectedPop


    # --------------------------- Acceptance ----------------------------
    def acceptance(self, oldSol, newSol, temp):
        scoreOld = oldSol.score
        scoreNew = newSol.score
        if(scoreNew <= scoreOld):
            result = 1
        else:
            delthaE = abs(scoreOld - scoreNew)
            k = (1 / (1 + (delthaE / temp)))
            result = math.exp( - delthaE / (temp * k) )
        return result


    # ----------------------- Gets New Neighbor -------------------------
    def getNeighbor(self, oldSol):
        solution = copy.deepcopy(oldSol)
        solution = self.mutation(solution)
        solution.score = self.getScoring(solution)
        return solution


    #----------------------------- Mutation -----------------------------
    def mutation(self, selectedPop):
        rotation_bonds = len(selectedPop.rotateAtoms)
        if(rotation_bonds > 0):
            mutPos = random.randint(1, (6 + rotation_bonds))
        else:
             mutPos = random.randint(1, 6)
        if mutPos == 1:
            selectedPop.x = random.uniform(-self.__spaceLen, self.__spaceLen)
        elif mutPos == 2:
            selectedPop.y = random.uniform(-self.__spaceLen, self.__spaceLen)
        elif mutPos == 3:
            selectedPop.z = random.uniform(-self.__spaceLen, self.__spaceLen)
        elif mutPos == 4:
            selectedPop.sph_theta = math.pi * random.randint(0, 200) / 100.0
        elif mutPos == 5:
            selectedPop.sph_phi = math.pi * random.randint(0, 100) / 100.0
        elif mutPos == 6:
            selectedPop.theta = math.pi * random.randint(0, 200) / 100.0
        elif mutPos > 6:
            selectedPop.rotateAtoms[mutPos - 7] = math.pi * random.randint(0, 200) / 100.0

        return selectedPop


    # -------------------------- Gets Scoring ---------------------------
    def getScoring(self, solution):
        auxLigand = copy.deepcopy(self.__ligand)
        auxLigand.traslateToPoint([(self.__spaceCenter[0] + solution.x),
                                   (self.__spaceCenter[1] + solution.y),
                                   (self.__spaceCenter[2] + solution.z)])
        auxLigand.rotateByVector(Methods.spherePoint(1, solution.sph_theta, solution.sph_phi), solution.theta)
        for r in range(len(self.__ligand.rotationsPoints)):
            auxLigand.rotateSegment(r, solution.rotateAtoms[r])

        'Score'
        auxPose = rosetta.core.pose.Pose()
        auxPose.assign(self.__pose_ligand)
        for atom in xrange(auxPose.residue(auxPose.total_residue()).natoms()):
            atomV = rosetta.numeric.xyzVector_Real()
            atomV.x = round(auxLigand.x[atom],3)
            atomV.y = round(auxLigand.y[atom],3)
            atomV.z = round(auxLigand.z[atom],3)
            auxPose.residue(auxPose.total_residue()).set_xyz(auxLigand.atom[atom], atomV)
        result = self.__scoreFx.generateScoringByPose(auxPose, auxLigand, False)

        return result
