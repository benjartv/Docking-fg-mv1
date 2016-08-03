__author__ = 'Felipe'

from objects_py import ROOTPATH, ScoringFunction
from utils_py import IO, Methods
import os
import copy

'************************** PDB FILES ***************************'

__MoleculeName = '1BV9'
__LigandName = 'XV6'
__pdbLigand = '1BV9_Flexible_2_RMSD_1.400_ligand.pdb'

'****************************************************************'


'Constants'
tempPath = os.path.join(ROOTPATH,'temp')
pdbPath = os.path.join(ROOTPATH,'originalPdb')

'Reads data'
newLigand = IO.readPDB(__pdbLigand, 'HETATM', pdbPath)
moleculeList = IO.readAllPDB(__MoleculeName + '.pdb', __LigandName, pdbPath)
protein = copy.deepcopy(moleculeList[0])
ligand = copy.deepcopy(moleculeList[1])

'Writes a new clean structure'
IO.writePDB("ligand.pdb", ligand, tempPath)
IO.writePDB("protein.pdb", protein, tempPath)

'Gets Connection Matrix'
connectMatrix = ligand.connectMatrix

'Init Scoring Function'
scoringFx = ScoringFunction.ScoringFunction()
scoringFx.defineLigandStructure("ligand.pdb", __LigandName, tempPath)

'Recover Connection Matrix'
ligand = IO.readPDB("ligand.pdb", 'HETATM', tempPath)
ligand.connectMatrix = connectMatrix
newLigand = Methods.replacePositionLigand(newLigand, ligand)
IO.writeAllPDB("original.pdb", protein, ligand, tempPath)
IO.writeAllPDB("new.pdb", protein, newLigand, tempPath)

'Score'
minScore = scoringFx.generateScoring(protein, ligand, False)
newScore = scoringFx.generateScoring(protein, newLigand, False)
rmsdScore = scoringFx.getRMSD("original.pdb", "new.pdb", tempPath, tempPath)

print "\n************************** Results **************************"
print "- Ligand PDB: " + __pdbLigand
print "- Best Score: " + "{0:.3f}".format(minScore)
print "- Score: " + "{0:.3f}".format(newScore)
print "- RMSD: " + "{0:.3f}".format(rmsdScore)