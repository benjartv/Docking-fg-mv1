__author__ = 'Felipe'

from objects_py import Gene, ROOTPATH, DELETE_COM, ScoringFunction
from utils_py import IO, Methods
import random
import copy
import math
import os


'************************** PDB FILES ***************************'

'''
__MoleculeName = '1ENY'
__LigandName = 'NAD'
__SearchCenterPoint = [-4.200, 31.500, 12.000]
__rotateAtoms = [[36, 35], [5, 6]]
#[[12, 13], [5, 6]]
#[[23, 24], [36, 35]]
#[[23, 24], [12, 13]]
#[[23, 24], [5, 6]]
#[[36, 35], [12, 13]]
#[[36, 35], [5, 6]]
# [[23, 24], [36, 35], [12, 13], [5, 6]]


__MoleculeName = '1AJV'
__LigandName = 'NMB'
__SearchCenterPoint = [15.550, 21.000, 4.600]
__rotateAtoms = [[34, 33], [3, 23]]
#[[8, 9], [3, 23]]
#[[34, 33], [16, 41]]
#[[34, 33], [8, 9]]
#[[34, 33], [3, 23]]
#[[16, 41], [3, 23]]
#[[16, 41], [8, 9]]
#[[34, 33], [16, 41], [8, 9], [3, 23]]


__MoleculeName = '1HPX'
__LigandName = 'KNI'
__SearchCenterPoint = [4.000, 2.000, 10.500]
__rotateAtoms = [[26,27], [37,40]]
#[[37,40], [42,43]]
#[[16,17], [26,27]]
#[[16,17], [37,40]]
#[[16,17], [42,43]]
#[[26,27], [37,40]]
#[[26,27], [42,43]]
#[[16,17], [26,27], [37,40], [42,43]]


__MoleculeName = '2UPJ'
__LigandName = 'U02'
__SearchCenterPoint = [13.505, 23.500, 7.510]
__rotateAtoms = []
#[[14,13], [20,24]]
#[[14,13], [28,32]]
#[[14,13], [36,35]]
#[[20,24], [28,32]]
#[[20,24], [36,35]]
#[[28,32], [36,35]]
#[[14,13], [20,24], [28,32], [36,35]]

'''
__MoleculeName = '1BV9'
__LigandName = 'XV6'
__SearchCenterPoint = [-8.000,  16.600,  24.600]
__rotateAtoms = []
#[[7,27], [47,43]]
#[[7,27], [41,12]]
#[[7,27], [13,14]]
#[[13,14], [41,12]]
#[[13,14], [47,43]]
#[[41,12], [47,43]]
#[[7,27], [13,14], [41,12], [47,43]]



'************************** Parameters **************************'

__numberResult = 1000
__SearchSpaceSize = 5


'*********************+*** Other Params *************************'

__onlyLigand = False
__dataLog = True

'****************************************************************'


'Constants'
resultPath = os.path.join(ROOTPATH,'results')
tempPath = os.path.join(ROOTPATH,'temp')
pdbPath = os.path.join(ROOTPATH,'originalPdb')

'Reads data'
moleculeList = IO.readAllPDB(__MoleculeName + '.pdb', __LigandName, pdbPath)
protein = copy.deepcopy(moleculeList[0])
ligand = copy.deepcopy(moleculeList[1])
ligandOri = copy.deepcopy(moleculeList[1])

'Center Molecule'
originPoint = protein.findCenter()
ligandPoint = ligand.findCenter()
#protein.translate([originPoint[0] * -1, originPoint[1] * -1, originPoint[2] * -1])
#ligand.translate([originPoint[0] * -1, originPoint[1] * -1, originPoint[2] * -1])
#ligandOri.translate([originPoint[0] * -1, originPoint[1] * -1, originPoint[2] * -1])

'Traslate search center point'
alterCenterPoint = [__SearchCenterPoint[0] + (originPoint[0]*-1),
                    __SearchCenterPoint[1] + (originPoint[1]*-1),
                    __SearchCenterPoint[2] + (originPoint[2]*-1)]

'Traslate search center point'
alterLigandPoint = [ligandPoint[0] + (originPoint[0]*-1),
                    ligandPoint[1] + (originPoint[1]*-1),
                    ligandPoint[2] + (originPoint[2]*-1)]


print "\n------------------ Values ------------------"
print "-Search Point"
Methods.printPoint(__SearchCenterPoint)
Methods.printPoint(alterCenterPoint)
print "-Ligand Point"
Methods.printPoint(ligandPoint)
Methods.printPoint(alterLigandPoint)
print "-Distance"
print str(Methods.euclideanDistance(ligandPoint,__SearchCenterPoint))
print "--------------------------------------------\n"


'Writes a new clean structure'
IO.writePDB("ligand.pdb", ligand, tempPath)
IO.writePDB("protein.pdb", protein, tempPath)

'Gets Connection Matrix'
connectMatrix = ligand.connectMatrix

'Init Scoring Function'
scoringFx = ScoringFunction.ScoringFunction()
scoringFx.defineLigandStructure("ligand.pdb", __LigandName, tempPath)
ligand = IO.readPDB("ligand.pdb", 'HETATM', tempPath)

'Recover Connection Matrix'
ligand.connectMatrix = connectMatrix
ligand.addRotationsAtoms(__rotateAtoms)
IO.writeAllPDB("original_" + __MoleculeName + ".pdb", protein, ligand, tempPath)

'Random Change'
listaGenes =[]
for r in range(__numberResult):
    gen = Gene.Gene()
    gen.id = r
    #gen.x = random.uniform(-__SearchSpaceSize, __SearchSpaceSize)
    #gen.y = random.uniform(-__SearchSpaceSize, __SearchSpaceSize)
    #gen.z = random.uniform(-__SearchSpaceSize, __SearchSpaceSize)
    gen.sph_theta = math.pi * random.randint(0, 200) / 100.0
    gen.sph_phi = math.pi * random.randint(0, 100) / 100.0
    gen.theta = math.pi * random.randint(0, 200) / 100.0
    for j in range(len(ligand.rotationsPoints)):
        gen.rotateAtoms.append(math.pi * random.randint(0, 200) / 100.0)
    listaGenes.append(gen)

'Molecule Name'
score = scoringFx.generateScoring(protein, ligand)
ligandName = __MoleculeName + "_Score_" + "{0:.3f}".format(score) + "_original.pdb"
ligandOri.x = ligand.x[:]
ligandOri.y = ligand.y[:]
ligandOri.z = ligand.z[:]

'Original Molecule'
if(__onlyLigand):
    IO.writePDB("Ligand_"+ligandName +".pdb",ligandOri,resultPath)
else:
    IO.writeAllPDB(ligandName, protein, ligandOri, resultPath, True)

'Log'
log = ""

'Generate Files'
for i in range(len(listaGenes)):
    auxLigand = None
    auxLigand = copy.deepcopy(ligand)
    #auxLigand.traslateToPoint([listaGenes[i].x, listaGenes[i].y, listaGenes[i].z])
    auxLigand.rotateByVector(Methods.spherePoint(1, listaGenes[i].sph_theta, listaGenes[i].sph_phi), listaGenes[i].theta)
    for l in range(len(ligand.rotationsPoints)):
        auxLigand.rotateSegment(l, listaGenes[i].rotateAtoms[l])

    #if( len(ligand.rotationsPoints) > 0):
    #    auxLigand.traslateToPoint(ligandPoint)


    'RMSD'
    IO.writeAllPDB("result_aux.pdb", protein, auxLigand, resultPath)
    rmsd = scoringFx.getRMSD("original_" + __MoleculeName + ".pdb", "result_aux.pdb", tempPath, resultPath)

    'Create Result'
    ligandOri.x = auxLigand.x[:]
    ligandOri.y = auxLigand.y[:]
    ligandOri.z = auxLigand.z[:]

    'Create File'
    if( len(ligand.rotationsPoints) > 0):
        if(__onlyLigand):
            IO.writePDB("Ligand_"+__MoleculeName + "_Flexible_RMSD_" + "{0:.3f}".format(rmsd) + ".pdb",ligandOri,resultPath)
        else:
            IO.writeAllPDB(__MoleculeName + "_Flexible_RMSD_" + "{0:.3f}".format(rmsd) + ".pdb", protein, ligandOri, resultPath, True)
    else:
        if(__onlyLigand):
            IO.writePDB("Ligand_"+__MoleculeName + "_Rigid_RMSD_" + "{0:.3f}".format(rmsd) + ".pdb",ligandOri,resultPath)
        else:
            IO.writeAllPDB(__MoleculeName + "_Rigid_RMSD_" + "{0:.3f}".format(rmsd) + ".pdb", protein, ligandOri, resultPath, True)

    'Log'
    log = log + str(rmsd) + ","

    'Delete Aux File'
    os.system(DELETE_COM + os.path.join(os.path.join(ROOTPATH,'results'), "result_aux.pdb"))

'Write log file'
if(__dataLog):
    IO.writeLog("log_"+__MoleculeName+".log",log + "\n",resultPath)

'Delete files'
Methods.deleteTempFiles(__LigandName, __MoleculeName)