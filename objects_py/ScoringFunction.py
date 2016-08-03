__author__ = 'Felipe_Gonzalez_Foncea'

from objects_py import BABEL_COM, DELETE_COM, PYTHON_COM, ROOTPATH
from utils_py import IO
import rosetta
import copy
import os


# -------------------- Scoring Functions Object ---------------------
class ScoringFunction(object):

    # ------------------------------ Init -------------------------------
    def __init__(self):
        self.__residueSet = None
        self.__scoreRosetta = 0
        self.__ligandFile = None
        self.__ligandName = "NONAME"
        self.__bestLigand = None
        self.__bestScore = 999999.9
        self.__isBest = False
        self.__pdbStr = []
        self.__scorefxn = None

        'Files Paths'
        self.__resDir = os.path.join(ROOTPATH,'results')
        self.__molDir = os.path.join(ROOTPATH,'others')
        self.__tmpDir = os.path.join(ROOTPATH,'temp')

        'Init Rosetta'
        rosetta.init()


    # -------------------- Define Ligand Structure ----------------------
    def defineLigandStructure(self, ligandFile, nameLigand, pathfile=""):
        self.__ligandName = nameLigand

        'Validates path parameters'
        if(pathfile == ""):
            pathfile = self.__resDir

        'Delete Files'
        self.deleteRosettaFiles()

        'Creates mdl file'
        os.system(BABEL_COM + ' -ipdb ' + os.path.join(pathfile, ligandFile) + ' -omdl ' + os.path.join(self.__tmpDir, 'Struct.mdl'))

        'Generates param file'
        os.system(PYTHON_COM + os.path.join(self.__molDir, 'molfile_to_params.py') + ' ' + os.path.join(self.__tmpDir, 'Struct.mdl') + ' -d '+ self.__tmpDir + ' -n ' + self.__ligandName)
        self.__residueSet = rosetta.generate_nonstandard_residue_set([os.path.join(self.__tmpDir, self.__ligandName + '.params')])

        'Generates Ligand'
        self.__ligandFile = IO.readPDB(self.__ligandName + '_0001.pdb', 'HETATM',self.__tmpDir)
        IO.writePDB("ligand.pdb", self.__ligandFile, self.__tmpDir)

        'Init Function'
        self.__scorefxn = rosetta.create_score_function("ligand")


    # ----------------------- Generates Scoring -------------------------
    def generateScoring(self, protein, ligand, traceBest = True):

        ligand_p = rosetta.core.pose.Pose()
        rosetta.core.import_pose.pose_from_pdbstring(ligand_p, IO.getPbdStr(protein, ligand, self.__pdbStr))
        score = self.__scorefxn(ligand_p)

        'Best Values'
        if(traceBest):
            self.traceScoring(ligand, score)
        return score


    # ------------------- Generates Scoring by Pose ---------------------
    def generateScoringByPose(self, pose_ligand, ligand, traceBest = True):

        score = self.__scorefxn(pose_ligand)

        '''
        self.__scoreRosetta = score
        if(self.__scoreRosetta < self.__bestScore):
            self.__bestScore = self.__scoreRosetta
            self.__bestLigand = copy.deepcopy(ligand)
            self.__isBest = True
        '''

        'Best Values'
        if(traceBest):
            self.traceScoring(ligand, score)

        return score

    # -------------------------- Trace Scoring --------------------------
    def traceScoring(self, ligand, score):
        self.__scoreRosetta = score

        'Best Values '
        if(self.__scoreRosetta < self.__bestScore):
            self.__bestScore = self.__scoreRosetta
            self.__bestLigand = copy.deepcopy(ligand)
            self.__isBest = True


    # -------------------- Define Ligand Structure ----------------------
    def getRMSD(self, pdbName1, pdbName2, pathfile1="", pathfile2=""):

        if(pathfile1==""):
            fileDir1 = os.path.join(ROOTPATH, 'results')
            fileDir1 = os.path.join(fileDir1, pdbName1)
        else:
            fileDir1 = os.path.join(pathfile1, pdbName1)

        if(pathfile2==""):
            fileDir2 = os.path.join(ROOTPATH, 'results')
            fileDir2 = os.path.join(fileDir2, pdbName2)
        else:
            fileDir2 = os.path.join(pathfile2, pdbName2)

        pose1 = rosetta.pose_from_pdb(fileDir1)
        pose2 = rosetta.pose_from_pdb(fileDir2)
        scoreRMSD = rosetta.all_atom_rmsd(pose1, pose2)
        return scoreRMSD


    # --------------------- Delete Rosetta Files ------------------------
    def deleteRosettaFiles(self):
        os.system(DELETE_COM + os.path.join(self.__tmpDir, self.__ligandName + '.params'))
        os.system(DELETE_COM + os.path.join(self.__tmpDir, self.__ligandName + '_0001.pdb'))
        os.system(DELETE_COM + os.path.join(self.__tmpDir, 'Struct.mdl'))


    # -------------------- Return Ligand Structure ----------------------
    def getDefinedLigandScoring(self):
        return self.__ligandFile


    # ------------------------- Get Last Score --------------------------
    def getLastScore(self):
        return self.__scoreRosetta


    # ------------------------ Get Best Score ---------------------------
    def getBestScore(self):
        return self.__bestScore


    # ----------------------- Get Ligand Name ---------------------------
    def getLigandName(self):
        return self.__ligandName


    # ------------------- Get Best Ligand Structure ---------------------
    def getBestLigand(self):
        return self.__bestLigand


    # ------------------------ Find Best --------------------------------
    def findBest(self):
        if(self.__isBest):
            self.__isBest = False
            return True
        else:
            return False


    # --------------------- --- Get Pose --------------------------------
    def getPoseFromPdbStr(self, protein, ligand):
        ligand_p = rosetta.core.pose.Pose()
        rosetta.core.import_pose.pose_from_pdbstring(ligand_p, IO.getPbdStr(protein, ligand, self.__pdbStr))
        return ligand_p