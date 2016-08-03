__author__ = 'Felipe'

from objects_py import BABEL_COM, DELETE_COM, PYTHON_COM, ROOTPATH
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

#__MoleculeName = '2UPJ.pdb'
#__LigandName = 'U02'

#__MoleculeName = '1BV9.pdb'
#__LigandName = 'XV6'

__MoleculeName = '2UPJ_Flexible_4_RMSD_1.300.pdb'
__LigandName = 'U02'

'****************************************************************'

ligandFile = 'ligand_' + __LigandName + '.pdb'
proteinFile = 'protein_' + __LigandName + '.pdb'
molFile = 'ligand_' + __LigandName + '.mol'

'Constants'
resultPath = os.path.join(ROOTPATH,'results')
pdbPath = os.path.join(ROOTPATH,'testPdb')
tempPath = os.path.join(ROOTPATH,'temp')
originalPath = os.path.join(ROOTPATH,'originalPdb')
molPath = os.path.join(ROOTPATH,'others')

moleculeList = IO.readAllPDB(__MoleculeName, __LigandName, originalPath)
protein = copy.deepcopy(moleculeList[0])
ligand = copy.deepcopy(moleculeList[1])

IO.writePDB(ligandFile, ligand, tempPath)
IO.writePDB(proteinFile, protein, tempPath)


'Creates mdl file'
os.system(BABEL_COM + ' -ipdb ' + os.path.join(tempPath, ligandFile) + ' -omdl ' + os.path.join(tempPath, molFile))

os.system(DELETE_COM + os.path.join(tempPath, ligandFile))
os.system(DELETE_COM + os.path.join(tempPath, proteinFile))

