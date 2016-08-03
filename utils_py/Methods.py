__author__ = 'Felipe_Gonzalez_Foncea'

from objects_py import ParamsAlgorithms
from utils_py import ROOTPATH,DELETE_COM
import datetime
import random
import math
import os
import copy


# -------------------------- Cross Product --------------------------
def crossProduct(A, B):
    linesA = len(A)
    columnsA = len(A[0])

    linesB = len(B)
    columnsB = len(B[0])

    if columnsA == linesB:
        M = [[sum(A[m][n] * B[n][p] for n in range(columnsA)) \
              for p in range(columnsB)] for m in range(linesA)]
        return M
    else:
        return -1

# ----------------------- Euclidean Distance ------------------------
def euclideanDistance(p1,p2):
    return math.sqrt(math.pow((p2[0] - p1[0]), 2) +
                math.pow((p2[1] - p1[1]), 2) +
                math.pow((p2[2] - p1[2]), 2))


# --------------------------- Visual Mol ----------------------------
def ejecuteVisualMol():
    fileDir = os.path.join(os.path.dirname( __file__ ), 'VisualMol.py' )
    os.system('python ' + fileDir)


# -------------------------- Print Point -----------------------------
def printPoint(point):
    posx = point[0]
    posy = point[1]
    posz = point[2]
    print "Point ( " + str(posx) + " , " + str(posy) + " , " + str(posz) + ")"


# ------------------------- Contain List -----------------------------
def contains(small, big):
    result = False
    for i in range(len(small)):
        if(small[i] in big):
            result = True
            break
    return result


# ----------------------- Generates Matrix ---------------------------
def generatesMatrix(sizeC, sizeR, value):
    matrix = []
    for a in range(sizeC):
        matrix.append([])
        for b in range(sizeR):
            matrix[a].append(value)

    return matrix


# -------------------- Creates Params Object -------------------------
def generatesParams(searchSpaceSize,
                 originalCenterPoint,
                 newCenterPoint,
                 moleculeName,
                 testPdb,
                 algorithmType,
                 treeNodes,
                 treeLevels,
                 castPercent,
                 rotatesAtoms,
                 populationSize,
                 generationNumber,
                 localSearch,
                 initTemperature,
                 minTemperature,
                 alphaTemperature,
                 innerGenerations,
                 timeStop,
                 loopReset):

    paramsAlgorithms = ParamsAlgorithms.ParamsAlgorithms()
    paramsAlgorithms.searchSpaceSize = searchSpaceSize
    paramsAlgorithms.searchCenterPoint = originalCenterPoint
    paramsAlgorithms.newCenterPoint = newCenterPoint
    paramsAlgorithms.rotatesAtoms = rotatesAtoms
    paramsAlgorithms.populationSize = populationSize
    paramsAlgorithms.generationNumber = generationNumber
    paramsAlgorithms.initTemperature = initTemperature
    paramsAlgorithms.minTemperature = minTemperature
    paramsAlgorithms.alphaTemperature = alphaTemperature
    paramsAlgorithms.innerGenerations = innerGenerations
    paramsAlgorithms.moleculeName = moleculeName
    paramsAlgorithms.testPdb = testPdb
    paramsAlgorithms.algorithmType = algorithmType
    paramsAlgorithms.isLocalSearch = localSearch
    paramsAlgorithms.folderName = generatesFolderName(moleculeName)
    paramsAlgorithms.timeStop = timeStop
    paramsAlgorithms.loopReset = loopReset

    'Optional Params'
    if(algorithmType == 1):
        paramsAlgorithms.populationCast = castPercent
    elif(algorithmType == 2):
        paramsAlgorithms.populationCast = [treeNodes, treeLevels]

    return paramsAlgorithms


# -------------------- Generates Folder Name -------------------------
def generatesFolderName(moleculeName):
    timeNow = datetime.datetime.now()
    folderName = 'Result'
    folderName += '_'
    folderName += moleculeName
    folderName += '_'
    folderName += '{:02d}'.format(timeNow.year)
    folderName += '{:02d}'.format(timeNow.month)
    folderName += '{:02d}'.format(timeNow.day)
    folderName += '_'
    folderName += '{:02d}'.format(timeNow.hour)
    folderName += '{:02d}'.format(timeNow.minute)
    folderName += '{:02d}'.format(timeNow.second)
    folder = os.path.join(ROOTPATH, 'results')
    folder = os.path.join(folder, folderName)
    return folder


# ----------------------- Contain ID List ----------------------------
def containsId(itemGene, listGene):
    result = False
    for i in range(len(listGene)):
        if(listGene[i].id == itemGene.id):
            result = True
            break
    return result


# -------------------- Random Sphere Point ---------------------------
def randomSpherePoint(radius):
    theta = math.pi * random.randint(0, 200) / 100.0
    phi = math.pi * random.randint(0, 100) / 100.0
    xAux = radius * math.cos(theta) * math.sin(phi)
    yAux = radius * math.sin(theta) * math.sin(phi)
    zAux = radius * math.cos(phi)
    return [xAux, yAux, zAux]


# ------------------------ Sphere Point ------------------------------
def spherePoint(radius, theta, phi):
    xAux = radius * math.cos(theta) * math.sin(phi)
    yAux = radius * math.sin(theta) * math.sin(phi)
    zAux = radius * math.cos(phi)
    return [xAux, yAux, zAux]


#---------------------------- Print Gen ------------------------------
def printGen(gen):
    cellStr = ""
    formatDecimal = "{0:.3f}"
    cellStr += "[\t"
    cellStr += formatDecimal.format(gen.x) + "\t"
    cellStr += formatDecimal.format(gen.y) + "\t"
    cellStr += formatDecimal.format(gen.z) + "\t"
    cellStr += formatDecimal.format(gen.sph_theta) + "\t"
    cellStr += formatDecimal.format(gen.sph_phi) + "\t"
    cellStr += formatDecimal.format(gen.theta) + "\t"
    cellStr += " ]\t"
    cellStr += "score: " + formatDecimal.format(gen.score) + "\t \t"
    cellStr += "\n"
    print cellStr


#------------------------- Print Population -------------------------
def printPopulation(population, withProb = False):
    cellStr = ""
    formatDecimal = "{0:.3f}"
    if(len(population)>0):
        for i in range(len(population)):
            cellStr += str(population[i].id) + "-\t"
            cellStr += "[\t"
            cellStr += formatDecimal.format(population[i].x) + "\t"
            cellStr += formatDecimal.format(population[i].y) + "\t"
            cellStr += formatDecimal.format(population[i].z) + "\t"
            cellStr += formatDecimal.format(population[i].sph_theta) + "\t"
            cellStr += formatDecimal.format(population[i].sph_phi) + "\t"
            cellStr += formatDecimal.format(population[i].theta) + "\t"
            cellStr += " ]\t"
            cellStr += "score: " + formatDecimal.format(population[i].score) + "\t"
            if(withProb):
                cellStr += "P: " + formatDecimal.format(population[i].probScore) + "\t "
            cellStr += "\n"
    else:
        cellStr += "0-\t[ ]"
        cellStr += "\n"

    print cellStr


#---------------------------- Print Tree -----------------------------
def printTree(populatio):
    cellStr = ""
    formatDecimal = "{0:.3f}"
    for i in range(len(populatio)):
        cellStr += str(i) + "-\t"
        cellStr += "[\t"
        if(populatio[i]==None):
            cellStr +="None"
            cellStr += "\t]\n"
        else:
            cellStr += formatDecimal.format(populatio[i].x) + "\t"
            cellStr += formatDecimal.format(populatio[i].y) + "\t"
            cellStr += formatDecimal.format(populatio[i].z) + "\t"
            cellStr += formatDecimal.format(populatio[i].sph_theta) + "\t"
            cellStr += formatDecimal.format(populatio[i].sph_phi) + "\t"
            cellStr += formatDecimal.format(populatio[i].theta) + "\t"
            cellStr += " ]\t"
            cellStr += "score: " + formatDecimal.format(populatio[i].score) + "\t \t"
            cellStr += "P: " + formatDecimal.format(populatio[i].probScore) + "\t \t"
            cellStr += "\n"
    print cellStr


# ----------------------- Delete All Files --------------------------
def deleteTempFiles(ligandName, moleculeName):
    tmpDir = os.path.join(ROOTPATH,'temp')
    os.system(DELETE_COM + os.path.join(tmpDir, 'ligand.pdb'))
    os.system(DELETE_COM + os.path.join(tmpDir, 'protein.pdb'))
    os.system(DELETE_COM + os.path.join(tmpDir, ligandName + '.params'))
    os.system(DELETE_COM + os.path.join(tmpDir, ligandName + '_0001.pdb'))
    os.system(DELETE_COM + os.path.join(tmpDir, 'Struct.mdl'))
    os.system(DELETE_COM + os.path.join(tmpDir, 'original_' + moleculeName + '.pdb'))


#-------------------- Replace Postions of Ligand --------------------
def replacePositionLigand(newLigand, ligand):
    ligandAux = copy.deepcopy(ligand)
    ligandAux.x = newLigand.x[:]
    ligandAux.y = newLigand.y[:]
    ligandAux.z = newLigand.z[:]
    return ligandAux