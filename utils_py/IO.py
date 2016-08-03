__author__ = 'Felipe'

from utils_py import ROOTPATH
from objects_py.Molecule import *
from decimal import *
import copy
import os


# ------------------- Read All PDB file separate --------------------
def readAllPDB(nameFile, residueName="", pathfile=""):
    if(pathfile==""):
        fileDir = os.path.join(ROOTPATH, 'testPdb')
        fileDir = os.path.join(fileDir, nameFile)
    else:
        fileDir = os.path.join(pathfile, nameFile)

    with open(fileDir, "r") as filePDF:
        protein = Molecule("ATOM")
        ligand = Molecule("HETATM")
        connect = []
        serialMap = {}
        count = 0

        while True:
            fileLine = filePDF.readline()

            if(fileLine[0:6] == "HETATM" and (fileLine[17:20] == residueName or residueName == "")):
                serialMap[fileLine[6:11].strip()] = count
                ligand.atom.append(fileLine[12:16])
                ligand.alterLocation.append(fileLine[16])
                ligand.residueName.append(fileLine[17:20])
                ligand.chainId.append(fileLine[21])
                ligand.residueSeq.append(fileLine[22:26])
                ligand.insertCode.append(fileLine[26])
                ligand.x.append(float(fileLine[30:38]))
                ligand.y.append(float(fileLine[38:46]))
                ligand.z.append(float(fileLine[46:54]))
                ligand.occupancy.append(fileLine[54:60])
                ligand.tempFactor.append(fileLine[60:66])
                ligand.segmentId.append(fileLine[72:76])
                ligand.symbol.append(fileLine[76:78])
                ligand.charge.append(fileLine[78:80])
                count = count + 1

            elif(fileLine[0:6] == "ATOM  " and (fileLine[17:20] != residueName or residueName == "")):
                protein.atom.append(fileLine[12:16])
                protein.alterLocation.append(fileLine[16])
                protein.residueName.append(fileLine[17:20])
                protein.chainId.append(fileLine[21])
                protein.residueSeq.append(fileLine[22:26])
                protein.insertCode.append(fileLine[26])
                protein.x.append(float(fileLine[30:38]))
                protein.y.append(float(fileLine[38:46]))
                protein.z.append(float(fileLine[46:54]))
                protein.occupancy.append(fileLine[54:60])
                protein.tempFactor.append(fileLine[60:66])
                protein.segmentId.append(fileLine[72:76])
                protein.symbol.append(fileLine[76:78])
                protein.charge.append(fileLine[78:80])

            elif(fileLine[0:6] == "CONECT"):
                connect.append(fileLine[6:].strip().replace("  "," ").split(" "))

            elif(fileLine == "END") or (not fileLine):
                filePDF.close()
                break

        if(len(connect)>0):
            connectMatrix = generatesMatrix(count, count, 0)
            for con in connect:
                for j in con:
                    if(con[0].strip() != j.strip()):
                        connectMatrix[serialMap[con[0].strip()]][serialMap[j.strip()]] = 1

            ligand.connectMatrix = connectMatrix

    return [protein,ligand]


# --------------------- Read PDB file separate ----------------------
def readPDB(nameFile, recordName, pathfile=""):
    if(pathfile==""):
        fileDir = os.path.join(ROOTPATH, 'testPdb')
        fileDir = os.path.join(fileDir, nameFile)
    else:
        fileDir = os.path.join(pathfile, nameFile)

    with open(fileDir, "r") as filePDF:
        molecule = Molecule(recordName)
        connect = []
        serialMap = {}
        count = 0

        while True:
            fileLine = filePDF.readline()

            if(fileLine[0:6] == recordName ):
                serialMap[fileLine[6:11].strip()] = count
                molecule.atom.append(fileLine[12:16])
                molecule.alterLocation.append(fileLine[16])
                molecule.residueName.append(fileLine[17:20])
                molecule.chainId.append(fileLine[21])
                molecule.residueSeq.append(fileLine[22:26])
                molecule.insertCode.append(fileLine[26])
                molecule.x.append(float(fileLine[30:38]))
                molecule.y.append(float(fileLine[38:46]))
                molecule.z.append(float(fileLine[46:54]))
                molecule.occupancy.append(fileLine[54:60])
                molecule.tempFactor.append(fileLine[60:66])
                molecule.segmentId.append(fileLine[72:76])
                molecule.symbol.append(fileLine[76:78])
                molecule.charge.append(fileLine[78:80])
                count = count + 1

            elif(fileLine[0:6] == "CONECT" and recordName == "HETATM"):
                connect.append(fileLine[6:].strip().split(" "))

            elif(fileLine == "END") or (not fileLine):
                filePDF.close()
                break

        if(len(connect)>0):
            connectMatrix = generatesMatrix(count,count,0)
            for con in connect:
                for j in con:
                    if(con[0].strip() != j.strip()):
                        connectMatrix[serialMap[con[0].strip()]][serialMap[j.strip()]] = 1

            molecule.connectMatrix = connectMatrix

    return molecule


# ---------------------- Write All In PDB file ----------------------
def writeAllPDB(nameFile, protein, ligand, pathfile="", withConnect=False):
    if(pathfile==""):
        fileDir = os.path.join(ROOTPATH, 'results')
        fileDir = os.path.join(fileDir, nameFile)
    else:
        fileDir = os.path.join(pathfile, nameFile)

    if not os.path.exists(os.path.dirname(fileDir)):
        os.makedirs(os.path.dirname(fileDir))

    with open(fileDir, "w") as filePDF:

        for i in xrange(0,len(protein.x)):
            textLine = []
            textLine.append("{:<6}".format(protein.recordName))
            textLine.append("{:>5}".format(str(i+1)))
            textLine.append("{:>5}".format(protein.atom[i]))
            textLine.append("{:>1}".format(protein.alterLocation[i]))
            textLine.append("{:>3}".format(protein.residueName[i]))
            textLine.append("{:>2}".format(protein.chainId[i]))
            textLine.append("{:>4}".format(protein.residueSeq[i]))
            textLine.append("{:<4}".format(protein.insertCode[i]))
            textLine.append("{:>8}".format(str(Decimal(protein.x[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>8}".format(str(Decimal(protein.y[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>8}".format(str(Decimal(protein.z[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>6}".format(protein.occupancy[i]))
            textLine.append("{:>6}".format(protein.tempFactor[i]))
            textLine.append("{:>10}".format(protein.segmentId[i]))
            textLine.append("{:>2}".format(protein.symbol[i]))
            textLine.append("{:>2}".format(protein.charge[i]))
            textLine.append("\n")
            filePDF.write(''.join(textLine))

        textLine = []
        textLine.append("{:<6}".format("TER"))
        textLine.append("\n")
        filePDF.write(''.join(textLine))

        for i in xrange(0,len(ligand.x)):
            textLine = []
            textLine.append("{:<6}".format(ligand.recordName))
            textLine.append("{:>5}".format(str(i+1+len(protein.x))))
            textLine.append("{:>5}".format(ligand.atom[i]))
            textLine.append("{:>1}".format(ligand.alterLocation[i]))
            textLine.append("{:>3}".format(ligand.residueName[i]))
            textLine.append("{:>2}".format(ligand.chainId[i]))
            textLine.append("{:>4}".format(ligand.residueSeq[i]))
            textLine.append("{:<4}".format(ligand.insertCode[i]))
            textLine.append("{:>8}".format(str(Decimal(ligand.x[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>8}".format(str(Decimal(ligand.y[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>8}".format(str(Decimal(ligand.z[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>6}".format(ligand.occupancy[i]))
            textLine.append("{:>6}".format(ligand.tempFactor[i]))
            textLine.append("{:>10}".format(ligand.segmentId[i]))
            textLine.append("{:>2}".format(ligand.symbol[i]))
            textLine.append("{:>2}".format(ligand.charge[i]))
            textLine.append("\n")
            filePDF.write(''.join(textLine))

        if(withConnect):
            for i in xrange(len(ligand.connectMatrix)):
                aux = " " + str(i + 1+len(protein.x)) + " "
                for j in xrange(len(ligand.connectMatrix[0])):
                    if(ligand.connectMatrix[i][j]==1):
                        aux += str(j + 1+len(protein.x)) + " "

                if(len(ligand.connectMatrix[0]) > 0):
                    textLine = []
                    textLine.append("{:<6}".format("CONECT"))
                    textLine.append("{:<74}".format(aux))
                    textLine.append("\n")
                    filePDF.write(''.join(textLine))

        filePDF.write("END")
        filePDF.close()


# ------------------------- Write PDB file --------------------------
def writePDB(nameFile, molecule, pathfile="", withConnect=False):
    if(pathfile==""):
        fileDir = os.path.join(ROOTPATH, 'results')
        fileDir = os.path.join(fileDir, nameFile)
    else:
        fileDir = os.path.join(pathfile, nameFile)

    if not os.path.exists(os.path.dirname(fileDir)):
        os.makedirs(os.path.dirname(fileDir))

    with open(fileDir, "w") as filePDF:

        for i in xrange(0,len(molecule.x)):
            textLine = []
            textLine.append("{:<6}".format(molecule.recordName))
            textLine.append("{:>5}".format(str(i+1)))
            textLine.append("{:>5}".format(molecule.atom[i]))
            textLine.append("{:>1}".format(molecule.alterLocation[i]))
            textLine.append("{:>3}".format(molecule.residueName[i]))
            textLine.append("{:>2}".format(molecule.chainId[i]))
            textLine.append("{:>4}".format(molecule.residueSeq[i]))
            textLine.append("{:<4}".format(molecule.insertCode[i]))
            textLine.append("{:>8}".format(str(Decimal(molecule.x[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>8}".format(str(Decimal(molecule.y[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>8}".format(str(Decimal(molecule.z[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>6}".format(molecule.occupancy[i]))
            textLine.append("{:>6}".format(molecule.tempFactor[i]))
            textLine.append("{:>10}".format(molecule.segmentId[i]))
            textLine.append("{:>2}".format(molecule.symbol[i]))
            textLine.append("{:>2}".format(molecule.charge[i]))
            textLine.append("\n")
            filePDF.write(''.join(textLine))

        if(withConnect):
            for i in xrange(len(molecule.connectMatrix)):
                aux = " " + str(i + 1+len(molecule.x)) + " "
                for j in xrange(len(molecule.connectMatrix[0])):
                    if(molecule.connectMatrix[i][j]==1):
                        aux += str(j + 1+len(molecule.x)) + " "

                if(len(molecule.connectMatrix[0]) > 0):
                    textLine = []
                    textLine.append("{:<6}".format("CONECT"))
                    textLine.append("{:<74}".format(aux))
                    textLine.append("\n")
                    filePDF.write(''.join(textLine))

        filePDF.write("END")
        filePDF.close()


# -------------------- Create Original Molecule ---------------------
def createsOriginalMolecule(moleculeName, ligandName, connectMatrix):

    new_ligand = readPDB("ligand.pdb", 'HETATM', os.path.join(ROOTPATH, 'temp'))
    moleculeList = readAllPDB(moleculeName + ".pdb", ligandName, os.path.join(ROOTPATH, 'originalPdb'))
    protein = copy.deepcopy(moleculeList[0])
    ligand = copy.deepcopy(moleculeList[1])

    originPoint = protein.findCenter()
    protein.translate([originPoint[0] * -1, originPoint[1] * -1, originPoint[2] * -1])
    ligand.translate([originPoint[0] * -1, originPoint[1] * -1, originPoint[2] * -1])

    new_ligand.x = ligand.x[:]
    new_ligand.y = ligand.y[:]
    new_ligand.z = ligand.z[:]
    new_ligand.connectMatrix = connectMatrix

    writeAllPDB("original_" + moleculeName + ".pdb", protein, new_ligand, os.path.join(ROOTPATH,'temp'))


# ---------------------------- Write Log ----------------------------
def writeLog(nameFile, data,  pathfile=""):
    if(pathfile==""):
        fileDir = os.path.join(ROOTPATH, 'temp')
    else:
        fileDir = os.path.join(pathfile, nameFile)

    if not os.path.exists(os.path.dirname(fileDir)):
        os.makedirs(os.path.dirname(fileDir))

    with open(fileDir, "w") as filePDF:
        filePDF.write(data)
        filePDF.close()


# ---------------------------- Write Log ----------------------------
def writeHistoricLog(nameFile, data,  pathfile=""):
    if(pathfile==""):
        fileDir = os.path.join(ROOTPATH, 'temp')
    else:
        fileDir = os.path.join(pathfile, nameFile)

    if not os.path.exists(os.path.dirname(fileDir)):
        os.makedirs(os.path.dirname(fileDir))

    with open(fileDir, "a+") as filePDF:
        filePDF.write(data)
        filePDF.close()


# ------------------------ Get PDB as String ------------------------
def getPbdStr(protein, ligand, proteinStr):

    resultFinal = ''
    proteinLen = len(protein.x)
    ligandLen = len(ligand.x)

    if(len(proteinStr)==0):
        textLine = []

        for i in xrange(proteinLen):
            textLine.append("{:<6}".format(protein.recordName))
            textLine.append("{:>5}".format(str(i+1)))
            textLine.append("{:>5}".format(protein.atom[i]))
            textLine.append("{:>1}".format(protein.alterLocation[i]))
            textLine.append("{:>3}".format(protein.residueName[i]))
            textLine.append("{:>2}".format(protein.chainId[i]))
            textLine.append("{:>4}".format(protein.residueSeq[i]))
            textLine.append("{:<4}".format(protein.insertCode[i]))
            textLine.append("{:>8}".format(str(Decimal(protein.x[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>8}".format(str(Decimal(protein.y[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>8}".format(str(Decimal(protein.z[i]).quantize(Decimal('1.000')))))
            textLine.append("{:>6}".format(protein.occupancy[i]))
            textLine.append("{:>6}".format(protein.tempFactor[i]))
            textLine.append("{:>10}".format(protein.segmentId[i]))
            textLine.append("{:>2}".format(protein.symbol[i]))
            textLine.append("{:>2}".format(protein.charge[i]))
            textLine.append("\n")

        textLine.append("{:<6}".format("TER"))
        textLine.append("\n")
        resultFinal += ''.join(textLine)
        proteinStr.append(resultFinal)
    else:
        resultFinal = proteinStr[0]

    textLine = []
    proteinLen += 1
    for i in xrange(ligandLen):
        textLine.append("{:<6}".format(ligand.recordName))
        textLine.append("{:>5}".format(str(i + proteinLen)))
        textLine.append("{:>5}".format(ligand.atom[i]))
        textLine.append("{:>1}".format(ligand.alterLocation[i]))
        textLine.append("{:>3}".format(ligand.residueName[i]))
        textLine.append("{:>2}".format(ligand.chainId[i]))
        textLine.append("{:>4}".format(ligand.residueSeq[i]))
        textLine.append("{:<4}".format(ligand.insertCode[i]))
        textLine.append("{:>8}".format(str(Decimal(ligand.x[i]).quantize(Decimal('1.000')))))
        textLine.append("{:>8}".format(str(Decimal(ligand.y[i]).quantize(Decimal('1.000')))))
        textLine.append("{:>8}".format(str(Decimal(ligand.z[i]).quantize(Decimal('1.000')))))
        textLine.append("{:>6}".format(ligand.occupancy[i]))
        textLine.append("{:>6}".format(ligand.tempFactor[i]))
        textLine.append("{:>10}".format(ligand.segmentId[i]))
        textLine.append("{:>2}".format(ligand.symbol[i]))
        textLine.append("{:>2}".format(ligand.charge[i]))
        textLine.append("\n")
    resultFinal += ''.join(textLine)

    return resultFinal

# ------------------------ knowledge base ------------------------

def readAngles(pathfile=""):
    if(pathfile==""):
        fileDir = os.path.join(ROOTPATH, 'anglesBase')
        fileDir = os.path.join(fileDir, "Dihedral-ANGLE_O5B-C5B__C5B-C4B.txt")
    else:
        fileDir = os.path.join(pathfile, "Dihedral-ANGLE_O5B-C5B__C5B-C4B.txt")
    afile = open(fileDir, "r")
    #afile = open("Dihedral-ANGLE_O5B-C5B__C5B-C4B.txt", "r")
    angles = list()
    for line in afile:
        number = line.split("\t")
        number[1] = number[1].split("\n")[0]
        angles.append(number)
    afile.close()

    return angles





