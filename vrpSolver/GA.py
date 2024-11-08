import heapq

from .common import *
from .geometry import *
from .ds import *

def algoGA(
	popSize,
	initSeqFunc=None,
	feasibleCheckFunc=None,
	repairFunc=None,
	fitnessFunc=None,
	localSearchDict=None,
	stopCriteriaDict=None,
	*args
	) -> dict:

	# Create the initial solution bank ========================================
	seqHeap = []
	for i in range(popSize):
		seqBank[i] = {
			'encode': initSeqFunc(args),
			'fitness': fitnessFunc(args)
		}


	# Selection ===============================================================

	# Crossover ===============================================================




def initSeqFunc(length) -> list:

	return seq

def feasibleCheck(seq):
	return True

def repairFunc(seq):
	return True

def fitnessFunc(seq: Ring) -> float:

	return True

