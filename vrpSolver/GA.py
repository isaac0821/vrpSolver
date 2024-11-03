
from .common import *
from .geometry import *
from .ds import *

def algoGA(
	popSize,
	initSeqFunc = initSeqFunc,
	feasibleCheckFunc,
	repairFunc,
	fitnessFunc,
	localSearchDict,
	stopCriteriaDict,
	*args
	) -> dict:

	# Create the initial solution bank ========================================
	seqBank = {}
	for 

def initSeqFunc() -> list:
	return

def feasibleCheck(seq):
	return True

def repairFunc(seq):
	return True

def fitnessFunc(seq) -> float:
	return True

