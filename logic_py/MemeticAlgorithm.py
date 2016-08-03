__author__ = 'Felipe'

from objects_py.Gene import *
from objects_py.TimeManager import *
from logic_py import LocalSearch, ROOTPATH
from utils_py import Methods, IO
import rosetta
import random
import math
import copy
import os


# ------------------------ Memetic Algorithm ------------------------
class MemeticAlgorithm(object):

    # ------------------------------ Init -------------------------------
    def __init__(self, params, scoreFx, protein, ligand):
        self.__searchSpaceSize = params.searchSpaceSize
        self.__spaceCenter = params.searchCenterPoint
        self.__newSpaceCenter = params.newCenterPoint
        self.__rotateAtoms = params.rotatesAtoms
        self.__fileDir = params.folderName
        self.__moleculeName = params.moleculeName
        self.__ScoreFx = scoreFx
        self.__Protein = protein
        self.__Ligand = ligand
        self.__LocalSearch = None
        self.__isLocalSearch = params.isLocalSearch
        self.__generations = params.generationNumber
        self.__pocketSize = params.populationSize
        self.__loop = 0
        self.__log = ""
        self.__dataLog = ""
        self.__bestLoop = 0
        self.__loopReset = params.loopReset
        self.__lastReset = 0
        self.__resetCount = 0
        self.__treeNodes = params.populationCast[0]
        self.__treeLevels = params.populationCast[1]
        self.__treeSize = self.calculatesTree()
        self.__testPdb = params.testPdb
        self.__tmpDir = os.path.join(ROOTPATH,'temp')
        self.__treePopulation = [None] * (self.__treeSize * self.__pocketSize)
        self.__treeCurrent = [None] * self.__treeSize
        self.__timeStop = params.timeStop

        'Timer start'
        self.__timer = TimeManager()

        'Add Rotate Atoms'
        self.__Ligand.addRotationsAtoms(self.__rotateAtoms)

        'Pose Rosetta'
        self.__pose_ligand = self.__ScoreFx.getPoseFromPdbStr(self.__Protein, self.__Ligand)

        'Init Local Search'
        if(self.__isLocalSearch):
            self.__LocalSearch = LocalSearch.LocalSeach(params.initTemperature,
                                                        params.minTemperature,
                                                        params.alphaTemperature,
                                                        params.innerGenerations,
                                                        copy.copy(self.__ScoreFx),
                                                        self.__Protein,
                                                        self.__searchSpaceSize,
                                                        self.__newSpaceCenter,
                                                        self.__pose_ligand)

    #--------------------------- Init Process ---------------------------
    def initProcess(self):
        self.__timer.start()
        self.initPopulation()

        while(self.validateEvolution()):
            self.generation()

            if(self.__ScoreFx.findBest()):
                self.__bestLoop = self.__loop


    #------------------------- Validates Iteration ----------------------
    def validateEvolution(self):

        self.printPopulation((self.__loop == 0 and "Init Population" or "Generation " + str(self.__loop)), True)

        self.__timer.stop()
        if(self.__generations > 0 and self.__loop >= self.__generations):
            self.generateFinalLog()
            return False
        elif(self.__timeStop > 0 and self.__timeStop <= self.__timer.getTimeHours()):
            self.generateFinalLog()
            return False
        else:
            self.__loop += 1
            self.reloadPopulation()
            return True


    #---------------------- Recalculate Population ----------------------
    def initPopulation(self):
        for i in xrange(self.__treeSize):
            cell = self.randomCell()
            cell = self.calculates(cell)
            self.addToCurrent(i,cell)
        self.updateTree(True)


    #---------------------- Recalculate Population ----------------------
    def generation(self):
        rootSize = int(self.__treeSize/self.__treeNodes)
        for i in xrange(rootSize):
            for j in xrange(self.__treeNodes):
                auxNodes = (i * self.__treeNodes) + (j + 1)
                p1 = self.randomSelection(i)
                p2 = self.randomSelection(auxNodes)
                np = self.crossover(p1,p2)
                np = self.mutation(np)
                np = self.calculates(np)
                self.addToCurrent(auxNodes,np)
        self.updateTree()


    #------------------------ Using Local Seach ------------------------
    def calculates(self, np):
        if(self.__isLocalSearch):
            np = self.__LocalSearch.initLocalProcess(np, self.__Ligand)
            auxLigand = copy.deepcopy(self.__Ligand)
            auxLigand.traslateToPoint([(self.__newSpaceCenter[0] + np.x),
                                       (self.__newSpaceCenter[1] + np.y),
                                       (self.__newSpaceCenter[2] + np.z)])
            auxLigand.rotateByVector(Methods.spherePoint(1, np.sph_theta, np.sph_phi), np.theta)
            for r in xrange(len(self.__rotateAtoms)):
                auxLigand.rotateSegment(r, np.rotateAtoms[r])
            self.__ScoreFx.traceScoring(auxLigand,np.score)
        else:
            auxLigand = copy.deepcopy(self.__Ligand)
            auxLigand.traslateToPoint([(self.__newSpaceCenter[0] + np.x),
                                       (self.__newSpaceCenter[1] + np.y),
                                       (self.__newSpaceCenter[2] + np.z)])
            auxLigand.rotateByVector(Methods.spherePoint(1, np.sph_theta, np.sph_phi), np.theta)
            for r in xrange(len(self.__rotateAtoms)):
                auxLigand.rotateSegment(r, np.rotateAtoms[r])

            'Score'
            auxPose = rosetta.core.pose.Pose()
            auxPose.assign(self.__pose_ligand)
            for atom in xrange(auxPose.residue(auxPose.total_residue()).natoms()):
                atomV = rosetta.numeric.xyzVector_Real()
                atomV.x = round(auxLigand.x[atom],3)
                atomV.y = round(auxLigand.y[atom],3)
                atomV.z = round(auxLigand.z[atom],3)
                auxPose.residue(auxPose.total_residue()).set_xyz(auxLigand.atom[atom], atomV)
            np.score = self.__ScoreFx.generateScoringByPose(auxPose, auxLigand)
        return np


    #---------------------- Recalculate Population ----------------------
    def updateTree(self, initGeneration = False ):
        rootSize = int(self.__treeSize/self.__treeNodes)

        'Move currents to pocket'
        for a in xrange(len(self.__treeCurrent)):
            if(self.__treeCurrent[a] != None):
                self.currentToPocket(a)

        '''
        'Move best solution to root'
        if(not initGeneration):
            for i in xrange(rootSize - 1, -1, -1):
                for j in xrange(self.__treeNodes):
                    auxNode = (i * self.__treeNodes) + (j+1)
                    self.replaceBestNodes(i, auxNode)
        '''

        'Move best solution to root'
        if(not initGeneration):
            for i in xrange(rootSize - 1, -1, -1):
                node_ids = [(i * self.__treeNodes) + ( j + 1) for j in xrange(self.__treeNodes)]
                self.replaceFromLeafs(i, node_ids)


        'Recalculates Pockets'
        for j in xrange(self.__treeSize):
            auxPop = self.getPocketsFromNode(j)
            auxPop = self.recalculatesPocket(auxPop)
            self.setPocketsInNode(j,auxPop)


    #--------------------------- Add to Pocket --------------------------
    def currentToPocket(self, nodeId):
        newCell = copy.deepcopy(self.__treeCurrent[nodeId])
        self.__treeCurrent[nodeId] = None
        self.addToPocket(nodeId,newCell)


    #-------------------------- Add to Pocket ---------------------------
    def addToPocket(self, nodeId, newCell):
        population = []
        populationAux = self.getPocketsFromNode(nodeId)

        'Discart all white spaces'
        for x in xrange(len(populationAux)):
            if(populationAux[x]!=None):
                population.append(populationAux[x])

        'Replace '
        if(len(population) == self.__pocketSize ):
            population.sort(key=lambda x: x.score, reverse=False)
            for l in xrange(len(population)-1, -1, -1):
                if(population[l].score > newCell.score):
                    population[l] = newCell
                    break
        else:
            population.append(newCell)

        'Order Pocket'
        population.sort(key=lambda x: x.score, reverse=False)

        'Set Pockets in population'
        self.setPocketsInNode(nodeId, population)


    #----------------------- Recalculates Pockets -----------------------
    def recalculatesPocket(self, population):
        popAux = []

        for q in xrange(len(population)):
            if(population[q] != None):
                popAux.append(copy.deepcopy(population[q]))

        popAux.sort(key=lambda x: x.score, reverse=False)

        return popAux[:]


    #------------------------- Replace Best Node ------------------------
    def replaceBestNodes(self, father, son):
        bestSonCell = son*self.__pocketSize
        bestFatherCell = father*self.__pocketSize

        if(self.__treePopulation[bestSonCell].score < self.__treePopulation[bestFatherCell].score):
            aux = copy.deepcopy(self.__treePopulation[bestSonCell])
            self.__treePopulation[bestSonCell] = self.__treePopulation[bestFatherCell]
            if(father == 0):
                self.addToPocket(father,aux)
            else:
                self.__treePopulation[bestFatherCell] = aux


    #------------------------ Replace from Leafs ------------------------
    def replaceFromLeafs(self, root_id, node_ids):

        bestFatherCell = root_id * self.__pocketSize
        bestSonCell = -1

        for i in xrange(len(node_ids)):
            auxNodeId = node_ids[i]*self.__pocketSize
            if(i == 0 or self.__treePopulation[auxNodeId].score < self.__treePopulation[bestSonCell].score):
                bestSonCell = auxNodeId

        if(bestSonCell != -1 and self.__treePopulation[bestSonCell].score < self.__treePopulation[bestFatherCell].score):
            aux = copy.deepcopy(self.__treePopulation[bestSonCell])
            self.__treePopulation[bestSonCell] = self.__treePopulation[bestFatherCell]
            if(root_id == 0):
                self.addToPocket(root_id, aux)
            else:
                self.__treePopulation[bestFatherCell] = aux


    #--------------------------- Add to Pocket --------------------------
    def addToCurrent(self, nodeId, newCell):
        self.__treeCurrent[nodeId] = copy.deepcopy(newCell)


    #------------------------- Gets Node Pockets ------------------------
    def getPocketsFromNode(self, nodeId):
        return self.__treePopulation[(nodeId*self.__pocketSize) : ((nodeId*self.__pocketSize)+self.__pocketSize)]


    #------------------------- Set Node Pockets -------------------------
    def setPocketsInNode(self, nodeId, pocket):
        for y in xrange(len(pocket)):
            self.__treePopulation[(nodeId*self.__pocketSize)+y] = copy.deepcopy(pocket[y])


    #-------------------------- Calculates Tree -------------------------
    def calculatesTree(self):
        totalNodes = 0
        for i in xrange(self.__treeLevels):
            totalNodes += self.__treeNodes**i
        return totalNodes


    # ------------------------- Apply Changes ---------------------------
    def applyChanges(self, solution):
        auxLigand = copy.deepcopy(self.__Ligand)
        auxLigand.traslateToPoint([(self.__newSpaceCenter[0] + solution.x),
                                   (self.__newSpaceCenter[1] + solution.y),
                                   (self.__newSpaceCenter[2] + solution.z)])
        auxLigand.rotateByVector(Methods.spherePoint(1, solution.sph_theta, solution.sph_phi), solution.theta)
        for r in xrange(len(self.__Ligand.rotationsPoints)):
            auxLigand.rotateSegment(r, solution.rotateAtoms[r])

        return auxLigand


    #------------------------- Random Selection --------------------------
    def randomSelection(self, nodeId):
        population = self.getPocketsFromNode(nodeId)
        popAux = []

        for q in xrange(len(population)):
            if(population[q] != None):
                popAux.append(copy.deepcopy(population[q]))

        randomId = random.randint(0, (len(popAux) - 1))
        selectedPop = copy.deepcopy(popAux[randomId])

        return selectedPop


    #---------------------- Recalculate Population ----------------------
    def randomCell(self):
        gen = Gene()
        gen.x = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        gen.y = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        gen.z = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        gen.sph_theta = math.pi * random.randint(0, 200) / 100.0
        gen.sph_phi = math.pi * random.randint(0, 100) / 100.0
        gen.theta = math.pi * random.randint(0, 200) / 100.0

        for r in xrange(len(self.__rotateAtoms)):
            gen.rotateAtoms.append(math.pi * random.randint(0, 200) / 100.0)

        return copy.deepcopy(gen)


    #---------------------------- CrossOver -----------------------------
    def crossover(self, selectedPop1, selectedPop2 ):

        'Uniform Random Position'
        if(random.randint(0,1)==1):
            aux_x = selectedPop1.x
            selectedPop1.x = selectedPop2.x
            selectedPop2.x = aux_x

        if(random.randint(0,1)==1):
            aux_y = selectedPop1.y
            selectedPop1.y = selectedPop2.y
            selectedPop2.y = aux_y

        if(random.randint(0,1)==1):
            aux_z = selectedPop1.z
            selectedPop1.z = selectedPop2.z
            selectedPop2.z = aux_z

        'Uniform Random Vector'
        if(random.randint(0,1)==1):
            auxSphPhi = selectedPop1.sph_phi
            selectedPop1.sph_phi = selectedPop2.sph_phi
            selectedPop2.sph_phi = auxSphPhi

        if(random.randint(0,1)==1):
            auxSphTheta = selectedPop1.sph_theta
            selectedPop1.sph_theta = selectedPop2.sph_theta
            selectedPop2.sph_theta= auxSphTheta

        'Switch Theta'
        if(random.randint(0,1)==1):
            aux_theta = selectedPop1.theta
            selectedPop1.theta = selectedPop2.theta
            selectedPop2.theta = aux_theta

        'Switch Bonds'
        rot_bond0 = len(selectedPop1.rotateAtoms)
        rot_bond1 = len(selectedPop2.rotateAtoms)
        if( rot_bond0 > 0 and rot_bond1 > 0 and rot_bond0 == rot_bond1):
            for i in xrange(rot_bond0):
                if(random.randint(0,1)==1):
                    aux_bond = selectedPop1.rotateAtoms[i]
                    selectedPop1.rotateAtoms[i] = selectedPop2.rotateAtoms[i]
                    selectedPop2.rotateAtoms[i] = aux_bond

        return selectedPop1


    #----------------------------- Mutation -----------------------------
    def mutation(self, selectedPop):

        muationProb = 0.1

        if(random.uniform(0,1) <= muationProb):
            selectedPop.x = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        if(random.uniform(0,1) <= muationProb):
            selectedPop.y = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        if(random.uniform(0,1) <= muationProb):
            selectedPop.z = random.uniform(-self.__searchSpaceSize, self.__searchSpaceSize)
        if(random.uniform(0,1) <= muationProb):
            selectedPop.sph_theta = math.pi * random.randint(0, 200) / 100.0
        if(random.uniform(0,1) <= muationProb):
            selectedPop.sph_phi = math.pi * random.randint(0, 100) / 100.0
        if(random.uniform(0,1) <= muationProb):
            selectedPop.theta = math.pi * random.randint(0, 200) / 100.0

        for i in xrange(len(selectedPop.rotateAtoms)):
            if(random.uniform(0,1) <= muationProb):
                selectedPop.rotateAtoms[i] = math.pi * random.randint(0, 200) / 100.0

        return selectedPop


    #------------------------- Print Population -------------------------
    def printPopulation(self, title, onLog = False):
        textList = []
        dataList = []

        textList.append("\n* ")
        textList.append(title)
        textList.append("\n")
        textList.append("**************************************************************************")
        for nodeId in xrange(self.__treeSize):
            population = self.getPocketsFromNode(nodeId)
            textList.append("\n")
            textList.append("Nodo ")
            textList.append(str(nodeId))
            textList.append("\n")

            for i in xrange(len(population)):
                textList.append(str(i+1))
                textList.append("-\t")
                textList.append("[\t")

                if(population[i] == None):
                    textList.append("None")
                    textList.append(" ]\n")
                    dataList.append("0")
                else:
                    textList.append("{0:.3f}".format(self.__newSpaceCenter[0] + population[i].x))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(self.__newSpaceCenter[1] + population[i].y))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(self.__newSpaceCenter[2] + population[i].z))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(population[i].sph_theta))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(population[i].sph_phi))
                    textList.append("\t")
                    textList.append("{0:.3f}".format(population[i].theta))
                    textList.append("\t")

                    for j in xrange(len(self.__rotateAtoms)):
                        textList.append("{0:.3f}".format(population[i].rotateAtoms[j]))
                        textList.append("\t")

                    textList.append(" ]\t")
                    textList.append("score: ")
                    textList.append("{0:.3f}".format(population[i].score))
                    textList.append("\t \t")
                    textList.append("\n")
                    dataList.append("{0:.3f}".format(population[i].score))

                if(nodeId != (self.__treeSize - 1) or i < (len(population) - 1)):
                    dataList.append(",")

        dataList.append("\n")


        if (onLog == False):
            print "\n* " + title + "\n"
            print ''.join(textList)
        else:
            self.__log += ''.join(textList)
            self.__dataLog += ''.join(dataList)


    #------------------------- Print Final Log --------------------------
    def generateFinalLog(self):
        self.__timer.stop()

        'Writes Files'
        ligand_PDB_name = "Ligand_score_" + "{0:.3f}".format(self.__ScoreFx.getBestScore()) + ".pdb"
        ligandProtein_PDB_name = "LigandProtein_score_" + "{0:.3f}".format(self.__ScoreFx.getBestScore()) + ".pdb"
        IO.writePDB(ligand_PDB_name, self.__ScoreFx.getBestLigand(), self.__fileDir)
        IO.writeAllPDB(ligandProtein_PDB_name, self.__Protein, self.__ScoreFx.getBestLigand(), self.__fileDir)

        initLog = "Molecule Name: " + self.__moleculeName + "\n"
        initLog += "Test Pdb: " + self.__testPdb + "\n"
        initLog += "Pocket: " + str(self.__pocketSize) + "\n"
        initLog += "Generations: " + str(self.__generations) + "\n"
        initLog += "Time to Stop: " + str(self.__timeStop) + "\n"
        initLog += "Space Size: " + str(self.__searchSpaceSize) + "\n"
        initLog += "Point: ( " + str(self.__spaceCenter[0]) + " , " + str(self.__spaceCenter[1]) + " , " + str(self.__spaceCenter[2]) + " )\n"
        initLog += "Local Search: " + str(self.__isLocalSearch) + "\n"
        initLog += "Rotate Atoms: " + str(self.__rotateAtoms) + "\n"
        initLog += "\n- Best Score: " + str(self.__ScoreFx.getBestScore()) + "\n"
        initLog += "- Reach Best Score: " + str(self.__bestLoop) + "\n"
        auxRMSD = str(self.__ScoreFx.getRMSD("original_" + self.__moleculeName + ".pdb", ligandProtein_PDB_name, self.__tmpDir, self.__fileDir))
        initLog += "- RMSD: " + auxRMSD + "\n"
        initLog += "- Timer: " + str(self.__timer.getTimeChronometer()) + "\n" + "\n"
        self.__log = initLog + self.__log

        'Data Log'
        IO.writeLog("iterations.log", self.__log, self.__fileDir)
        IO.writeLog("data.log", self.__dataLog, self.__fileDir)
        IO.writeHistoricLog("historic.log", str(self.__ScoreFx.getBestScore()) + ";" + auxRMSD + ";" + str(self.__bestLoop) + ";" + str(self.__loop) , self.__fileDir)

        'Delete files'
        Methods.deleteTempFiles(self.__ScoreFx.getLigandName(), self.__moleculeName)


    #----------------------- Reload Population --------------------------
    def reloadPopulation(self):
        if(self.__loopReset != 0):
            if(self.__lastReset >= self.__loopReset):
                self.__lastReset = 0
                self.__loopReset = self.__resetCount
                self.__resetCount = 0
                for i in xrange(1,len(self.__treePopulation)):
                    if( i % self.__pocketSize != 0):
                        self.__treePopulation[i] = None
            elif(self.__bestLoop == self.__loop):
                self.__lastReset = 0
                self.__resetCount += 1
            else:
                self.__lastReset += 1
                self.__resetCount += 1