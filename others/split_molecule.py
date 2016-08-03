__author__ = 'Felipe'

from objects_py import ROOTPATH
from utils_py import IO
import copy
import os


'************************** PDB FILES ***************************'

#__MoleculeName = '1ENY.pdb'
#__LigandName = 'NAD'

#__MoleculeName = '1AJV.pdb'
#__LigandName = 'NMB'

#__MoleculeName = '1HPX.pdb'
#__LigandName = 'KNI'

__MoleculeName = '2UPJ.pdb'
__LigandName = 'U02'

#__MoleculeName = '1BV9.pdb'
#__LigandName = 'XV6'

'****************************************************************'

'Constants'
moleculeName = __MoleculeName.replace(".pdb", "")
resultPath = os.path.join(ROOTPATH,'results')
pdbPath = os.path.join(ROOTPATH,'testPdb')
originalPath = os.path.join(ROOTPATH,'originalPdb')

'Reads data'
moleculeList = IO.readAllPDB(__MoleculeName, __LigandName, originalPath)
protein = copy.deepcopy(moleculeList[0])
ligand = copy.deepcopy(moleculeList[1])

'Writes a new clean structure'
IO.writePDB(moleculeName + "_ligand.pdb", ligand, resultPath)
IO.writePDB(moleculeName + "_protein.pdb", protein, resultPath)
IO.writeAllPDB(moleculeName + "_clean.pdb", protein, ligand, resultPath, True)

'Print Center'
proteinCenter = protein.findCenter()
ligandCenter = ligand.findCenter()
print "Protein Center [ " + str(proteinCenter[0]) + ", " + str(proteinCenter[1]) + ", " + str(proteinCenter[2]) + " ]"
print "Ligand Center  [ " + str(ligandCenter[0]) + ", " + str(ligandCenter[1]) + ", " + str(ligandCenter[2]) + " ]"