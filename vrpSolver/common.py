import math
import random

from vrpSolver.const import *

def rndInstance(
	N: "Number of vertices",
	xRange: "A 2-tuple with minimum/maximum range of x" = (0, 100),
	yRange: "A 2-tuple with minimum/maximum range of y" = (0, 100),
	) -> "A set of nodes with id start from 0 to N":
	nodeLoc = {}
	for i in range(N):
		x = random.randrange(xRange[0], xRange[1])
		y = random.randrange(yRange[0], yRange[1])
		nodeLoc[i] = (x, y)

	return nodeLoc

def rndSeq(
	N		: "Integer, Length of sequence",
	s0		: "Integer, Staring index of sequence" = 0,
	closed	: "Boolean, If the sequence is closed, if true, the last element is a duplicate of the first" = False,
	) -> "Randomly generate sequence starting from `start'":
	seq = [i for i in range(s0, N + s0)]
	# Randomly swap
	for i in range(N):
		j = random.randint(0, N - 1)
		t = seq[i]
		seq[i] = seq[j]
		seq[j] = t

	# Return to start?
	if (closed):
		seq.append(seq[0])
	return seq

def iterSeq(seqL, i, direction):
	q = None
	j = None
	if (direction == 'next'):
		if (i < seqL - 1):
			j = i + 1
		else:
			j = 0
	elif (direction == 'prev'):
		if (i > 0):
			j = i - 1
		else:
			j = seqL - 1
	else:
		return None
	return j

def randomPick(coefficients):
	totalSum = sum(coefficients)
	tmpSum = 0
	rnd = random.uniform(0, totalSum)
	index = 0
	for i in range(len(coefficients)):
		tmpSum += coefficients[i]
		if rnd <= tmpSum:
			index = i
			break
	return index

def getTauEuclidean(
	dicNodeLoc	: "Dictionary, {nodeID1: (x, y), ...}",
	lstNodeID	: "List, [nodeID1, nodeID2, ...]"
	) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":

	tau = {}
	for i in lstNodeID:
		for j in lstNodeID:
			if (i != j):
				d = math.sqrt((dicNodeLoc[i][0] - dicNodeLoc[j][0]) ** 2 + (dicNodeLoc[i][1] - dicNodeLoc[j][1]) ** 2)
				tau[i, j] = d
				tau[j, i] = d
			else:
				tau[i, j] = CONST_EPSILON
		tau[None, i] = CONST_EPSILON
		tau[i, None] = CONST_EPSILON

	return tau

def getTauSphereEuclidean(
	dicNodeLoc	: "Dictionary, {nodeID1: (lat, lon), ...}",
	lstNodeID	: "List, [nodeID1, nodeID2, ...]"
	) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":
	
	tau = {}
	for i in lstNodeID:
		for j in lstNodeID:
			if (i != j):
				d = sphereEuclidean2D(dicNodeLoc[i], dicNodeLoc[j])
				tau[i, j] = d
				tau[j, i] = d
			else:
				tau[i, j] = CONST_EPSILON
		tau[None, i] = CONST_EPSILON
		tau[i, None] = CONST_EPSILON

	return tau

def sphereEuclidean2D(
	coord1	: "First coordinate, in (lat, lon)", 
	coord2	: "Second coordinate, in (lat, lon)"
	) -> "Gives a Euclidean distance based on two lat/lon coords, if two coordinates are the same, return a small number":
	if (coord1 != None and coord2 != None):
		R = 3958.8  # Earth radius in miles
		lat1, lon1 = coord1
		lat2, lon2 = coord2
		phi1, phi2 = math.radians(lat1), math.radians(lat2) 
		dphi       = math.radians(lat2 - lat1)
		dlambda    = math.radians(lon2 - lon1)
		a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
		return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))
	else:
		return CONST_EPSILON

