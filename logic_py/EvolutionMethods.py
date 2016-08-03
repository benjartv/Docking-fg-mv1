__author__ = 'Felipe'

from logic_py.GeneticAlgorithm import *
from logic_py.MemeticAlgorithm import *


def startAlgorithm(params, scoreFx, protein, ligand ):
    if(params.algorithmType == 1):
        GeneticAlgorithm(params, scoreFx, protein, ligand).initProcess()
    elif(params.algorithmType == 2):
        MemeticAlgorithm(params, scoreFx, protein, ligand).initProcess()