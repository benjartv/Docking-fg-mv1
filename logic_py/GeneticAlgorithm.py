__author__ = 'Felipe'

from objects_py.Gene import *
from objects_py.TimeManager import *
from logic_py import LocalSearch, ROOTPATH
from utils_py import IO, Methods
import rosetta
import random
import copy
import math
import os


# ----------------------- Evolution Function ------------------------
class GeneticAlgorithm(object):


    # ------------------------------ Init -------------------------------
    def __init__(self, params, scoreFx, protein, ligand):
        self.population = []
        self.__loop = 0
        self.__bestLoop = 0
        self.__log = ""
        self.__dataLog = ""
        self.__rmsdLog = ""
        self.__localSearch = None

        'Params'
        self.__spaceLen = params.searchSpaceSize
        self.__spaceCenter = params.searchCenterPoint
        self.__newSpaceCenter = params.newCenterPoint
        self.__rotateAtoms = params.rotatesAtoms
        self.__localS = params.isLocalSearch
        self.__moleculeName = params.moleculeName
        self.__testPdb = params.testPdb
        self.__ligand = ligand
        self.__protein = protein
        self.__scoreFx = scoreFx
        self.__fileDir = params.folderName
        self.__generations = params.generationNumber
        self.__populationCast = params.populationCast
        self.__tmpDir = os.path.join(ROOTPATH,'temp')
        self.__pose_ligand = None
        self.__timeStop = params.timeStop

        'Timer start'
        self.__timer = TimeManager()

        'Generates Population'
        self.population.extend(self.randomPopulation(params.populationSize))

        'Add Rotate Atoms'
        self.__ligand.addRotationsAtoms(self.__rotateAtoms)

        'Pose Rosetta'
        self.__pose_ligand = self.__scoreFx.getPoseFromPdbStr(self.__protein, self.__ligand)

        'Init Local Search'
        if(self.__localS):
            self.__localSearch = LocalSearch.LocalSeach(params.initTemperature,
                                                        params.minTemperature,
                                                        params.alphaTemperature,
                                                        params.innerGenerations,
                                                        copy.copy(self.__scoreFx),
                                                        self.__protein,
                                                        self.__spaceLen,
                                                        self.__newSpaceCenter,
                                                        self.__pose_ligand)

        'Recalculates'
        self.population = self.recalculate(self.population)


    #--------------------------- Init Process ---------------------------
    def initProcess(self):
        self.__timer.start()

        sizePop = len(self.population)
        sizeA = abs((sizePop*self.__populationCast[0])/100)
        sizeB = abs((sizePop*self.__populationCast[1])/100)
        sizeC = sizePop - (sizeA + sizeB)

        while(self.validateEvolution()):
            newPopulation = []
            newPopulation.extend(self.population[0:sizeA])
            newPopulation.extend(self.mutationAndCross(self.population,sizeB))
            newPopulation.extend(self.randomPopulation(sizeC))

            newPopulation = self.recalculate(newPopulation,True)
            self.population = copy.deepcopy(newPopulation)

            if(self.__scoreFx.findBest()):
                self.__bestLoop = self.__loop


    #------------------------- Random Population ------------------------
    def randomPopulation(self, sizePop):
        population = []
        for i in xrange(sizePop):
            gen = Gene()
            gen.id = i
            gen.x = random.uniform(-self.__spaceLen, self.__spaceLen)
            gen.y = random.uniform(-self.__spaceLen, self.__spaceLen)
            gen.z = random.uniform(-self.__spaceLen, self.__spaceLen)
            gen.sph_theta = math.pi * random.randint(0, 200) / 100.0
            gen.sph_phi = math.pi * random.randint(0, 100) / 100.0
            gen.theta = math.pi * random.randint(0, 200) / 100.0

            for r in xrange(len(self.__rotateAtoms)):
                gen.rotateAtoms.append(math.pi * random.randint(0, 200) / 100.0)

            population.append(gen)
        return population[:]


    #---------------------- Recalculate Population ----------------------
    def recalculate(self, selectedPop, changeId=False):
        if(len(selectedPop)>0):

            for j in xrange(len(selectedPop)):
                auxLigand = copy.deepcopy(self.__ligand)
                auxLigand.traslateToPoint([(self.__newSpaceCenter[0] + selectedPop[j].x),
                                           (self.__newSpaceCenter[1] + selectedPop[j].y),
                                           (self.__newSpaceCenter[2] + selectedPop[j].z)])

                auxLigand.rotateByVector(Methods.spherePoint(1, selectedPop[j].sph_theta, selectedPop[j].sph_phi), selectedPop[j].theta)
                for r in xrange(len(self.__rotateAtoms)):
                    auxLigand.rotateSegment(r, selectedPop[j].rotateAtoms[r])

                'Score'
                auxPose = rosetta.core.pose.Pose()
                auxPose.assign(self.__pose_ligand)
                for atom in xrange(auxPose.residue(auxPose.total_residue()).natoms()):
                    atomV = rosetta.numeric.xyzVector_Real()
                    atomV.x = round(auxLigand.x[atom],3)
                    atomV.y = round(auxLigand.y[atom],3)
                    atomV.z = round(auxLigand.z[atom],3)
                    auxPose.residue(auxPose.total_residue()).set_xyz(auxLigand.atom[atom], atomV)
                selectedPop[j].score = self.__scoreFx.generateScoringByPose(auxPose, auxLigand)

            selectedPop.sort(key=lambda x: x.score, reverse=False)

            if(changeId):
                for k in xrange(len(selectedPop)):
                    selectedPop[k].id = k

        return selectedPop


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
            return True


    #----------------------- Generate New Population --------------------
    def mutationAndCross(self, population, popSize):
        newPopulation = []
        while(len(newPopulation)<popSize):
            selectedPop = self.rouletteWheel_selection(population, newPopulation)
            selectedPop = self.crossover(selectedPop)
            selectedPop = self.mutation(selectedPop)

            'localSearch'
            if(self.__localS):
                selectedPop[0] = self.__localSearch.initLocalProcess(selectedPop[0], self.__ligand)

            newPopulation.extend(selectedPop[:])

        return newPopulation[:]


    #---------------------- Roulette Wheel Selection --------------------
    def rouletteWheel_selection(self, population, newPopulation):
        self.calculateProbabilities(population)

        selectedPop = []
        lenPop = len(population)
        lenNewPop = len(newPopulation)

        if((lenNewPop+1) == lenPop):
            for k in xrange(lenPop):
                if(not Methods.containsId(population[k],newPopulation)):
                    selectedPop.append(copy.deepcopy(population[k]))
        else:
            while(1):
                roulettePoint = random.uniform(0, 1)
                for k in xrange(lenPop):
                    if(not Methods.containsId(population[k],selectedPop) and
                       not Methods.containsId(population[k],newPopulation) and
                       roulettePoint <= population[k].probScore):
                        selectedPop.append(copy.deepcopy(population[k]))
                        break
                if(len(selectedPop)==2):
                    break
        return selectedPop


    #---------------------- Calculate Probabilities ---------------------
    def calculateProbabilities(self, selectedPop):
        scoreList = []
        sumAux = 0.0
        sumScore = 0.0
        deltha = abs(selectedPop[-1].score) + abs(selectedPop[1].score)

        for i in xrange(len(selectedPop)):
            sumAux += (deltha - selectedPop[i].score)

        avarage = sumAux/len(selectedPop)

        for i in xrange(len(selectedPop)):
            fitness = math.exp((deltha - selectedPop[i].score)/avarage)
            scoreList.append(fitness)
            sumScore += fitness

        prob = 0.0
        for k in xrange(len(selectedPop)):
            if(sumScore==0):
                prob += 1
            else:
                prob += (scoreList[k] / sumScore)
            selectedPop[k].probScore = prob


    #---------------------------- CrossOver -----------------------------
    def crossover(self,selectedPop):
        if(len(selectedPop)==2):

            'Uniform Random Position'
            if(random.randint(0,1)==1):
                aux_x = selectedPop[0].x
                selectedPop[0].x = selectedPop[1].x
                selectedPop[1].x = aux_x

            if(random.randint(0,1)==1):
                aux_y = selectedPop[0].y
                selectedPop[0].y = selectedPop[1].y
                selectedPop[1].y = aux_y

            if(random.randint(0,1)==1):
                aux_z = selectedPop[0].z
                selectedPop[0].z = selectedPop[1].z
                selectedPop[1].z = aux_z

            'Uniform Random Vector'
            if(random.randint(0,1)==1):
                auxSphPhi = selectedPop[0].sph_phi
                selectedPop[0].sph_phi = selectedPop[1].sph_phi
                selectedPop[1].sph_phi = auxSphPhi

            if(random.randint(0,1)==1):
                auxSphTheta = selectedPop[0].sph_theta
                selectedPop[0].sph_theta = selectedPop[1].sph_theta
                selectedPop[1].sph_theta = auxSphTheta

            'Switch Theta'
            if(random.randint(0,1)==1):
                aux_theta = selectedPop[0].theta
                selectedPop[0].theta = selectedPop[1].theta
                selectedPop[1].theta = aux_theta

            'Switch Bonds'
            rot_bond0 = len(selectedPop[0].rotateAtoms)
            rot_bond1 = len(selectedPop[1].rotateAtoms)
            if( rot_bond0 > 0 and rot_bond1 > 0 and rot_bond0 == rot_bond1):
                for i in xrange(rot_bond0):
                    if(random.randint(0,1)==1):
                        aux_bond = selectedPop[0].rotateAtoms[i]
                        selectedPop[0].rotateAtoms[i] = selectedPop[1].rotateAtoms[i]
                        selectedPop[1].rotateAtoms[i] = aux_bond

        return selectedPop


    #----------------------------- Mutation -----------------------------
    def mutation(self, selectedPop):
        for i in xrange(len(selectedPop)):

            muationProb = 0.1

            if(random.uniform(0,1) <= muationProb):
                selectedPop[i].x = random.uniform(-self.__spaceLen, self.__spaceLen)
            if(random.uniform(0,1) <= muationProb):
                selectedPop[i].y = random.uniform(-self.__spaceLen, self.__spaceLen)
            if(random.uniform(0,1) <= muationProb):
                selectedPop[i].z = random.uniform(-self.__spaceLen, self.__spaceLen)
            if(random.uniform(0,1) <= muationProb):
                selectedPop[i].sph_theta = math.pi * random.randint(0, 200) / 100.0
            if(random.uniform(0,1) <= muationProb):
                selectedPop[i].sph_phi = math.pi * random.randint(0, 100) / 100.0
            if(random.uniform(0,1) <= muationProb):
                selectedPop[i].theta = math.pi * random.randint(0, 200) / 100.0

            for j in xrange(len(selectedPop[i].rotateAtoms)):
                if(random.uniform(0,1) <= muationProb):
                    selectedPop[i].rotateAtoms[j] = math.pi * random.randint(0, 200) / 100.0

        return selectedPop


    #------------------------- Print Population -------------------------
    def printPopulation(self, title, onLog = False):

        textList = []
        dataList = []

        if(len(self.population)>0):
            textList.append("\n* ")
            textList.append(title)
            textList.append("\n")

            for i in xrange(len(self.population)):
                textList.append(str(self.population[i].id))
                textList.append("-\t")
                textList.append("[\t")
                textList.append("{0:.3f}".format(self.__newSpaceCenter[0] + self.population[i].x))
                textList.append("\t")
                textList.append("{0:.3f}".format(self.__newSpaceCenter[1] + self.population[i].y))
                textList.append("\t")
                textList.append("{0:.3f}".format(self.__newSpaceCenter[2] + self.population[i].z))
                textList.append("\t")
                textList.append("{0:.3f}".format(self.population[i].sph_theta))
                textList.append("\t")
                textList.append("{0:.3f}".format(self.population[i].sph_phi))
                textList.append("\t")
                textList.append("{0:.3f}".format(self.population[i].theta))
                textList.append("\t")

                for j in xrange(len(self.__rotateAtoms)):
                    textList.append("{0:.3f}".format(self.population[i].rotateAtoms[j]))
                    textList.append("\t")

                textList.append(" ]\t")
                textList.append("score: ")
                textList.append("{0:.3f}".format(self.population[i].score))
                textList.append("\t \t")

                if (onLog == False):
                    textList.append("P: ")
                    textList.append("{0:.3f}".format(self.population[i].probScore))
                    textList.append("\t \t")

                textList.append("\n")

                dataList.append("{0:.3f}".format(self.population[i].score))

                if(i < (len(self.population)-1)):
                    dataList.append(",")
                else:
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
        ligand_PDB_name = "Ligand_score_" + "{0:.3f}".format(self.__scoreFx.getBestScore()) + ".pdb"
        ligandProtein_PDB_name = "LigandProtein_score_" + "{0:.3f}".format(self.__scoreFx.getBestScore()) + ".pdb"
        IO.writePDB(ligand_PDB_name,self.__scoreFx.getBestLigand(),self.__fileDir)
        IO.writeAllPDB(ligandProtein_PDB_name, self.__protein, self.__scoreFx.getBestLigand(), self.__fileDir)

        initLog = "Molecule Name: " + self.__moleculeName + "\n"
        initLog += "Test Pdb: " + self.__testPdb + "\n"
        initLog += "Population: " + str(len(self.population)) + "\n"
        initLog += "Generations: " + str(self.__generations) + "\n"
        initLog += "Time to Stop: " + str(self.__timeStop) + "\n"
        initLog += "Space Size: " + str(self.__spaceLen) + "\n"
        initLog += "Point: ( " + str(self.__spaceCenter[0]) + " , " + str(self.__spaceCenter[1]) + " , " + str(self.__spaceCenter[2]) + " )\n"
        initLog += "Local Search: " + str(self.__localS) + "\n"
        initLog += "Rotate Atoms: " + str(self.__rotateAtoms) + "\n"
        initLog += "\n- Best Score: " + str(self.__scoreFx.getBestScore()) + "\n"
        initLog += "- Reach Best Score: " + str(self.__bestLoop) + "\n"
        auxRMSD = str(self.__scoreFx.getRMSD("original_" + self.__moleculeName + ".pdb", ligandProtein_PDB_name, self.__tmpDir, self.__fileDir))
        initLog += "- RMSD: " + auxRMSD + "\n"
        initLog += "- Timer: " + str(self.__timer.getTimeChronometer()) + "\n" + "\n"
        self.__log = initLog + self.__log

        'Data Log'
        IO.writeLog("iterations.log", self.__log, self.__fileDir)
        IO.writeLog("data.log", self.__dataLog, self.__fileDir)
        IO.writeHistoricLog("historic.log", str(self.__scoreFx.getBestScore()) + ";" + auxRMSD + ";" + str(self.__bestLoop) + ";" + str(self.__loop) , self.__fileDir)

        'Delete files'
        Methods.deleteTempFiles(self.__scoreFx.getLigandName(), self.__moleculeName)


