import math
import random

from .const import *

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
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None
	) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":
	tau = {}
	lstNodeID = list(nodes.keys())
	for i in lstNodeID:
		for j in lstNodeID:
			if (i != j):
				d = math.sqrt((nodes[i]['loc'][0] - nodes[j]['loc'][0]) ** 2 + (nodes[i]['loc'][1] - nodes[j]['loc'][1]) ** 2)
				tau[i, j] = d
				tau[j, i] = d
			else:
				tau[i, j] = CONST_EPSILON
		tau[None, i] = CONST_EPSILON
		tau[i, None] = CONST_EPSILON

	return tau

def getTauSphereEuclidean(
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None
	) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":
	tau = {}
	lstNodeID = list(nodes.keys())
	for i in lstNodeID:
		for j in lstNodeID:
			if (i != j):
				d = sphereEuclidean2D(nodes[i]['loc'], nodes[j]['loc'])
				tau[i, j] = d
				tau[j, i] = d
			else:
				tau[i, j] = CONST_EPSILON
		tau[None, i] = CONST_EPSILON
		tau[i, None] = CONST_EPSILON

	return tau

def euclidean2D(
	coord1: "First coordinate, in (x, y)", 
	coord2: "Second coordinate, in (x, y)"
	) -> "Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number":
	if (coord1 != None and coord2 != None):
		return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2)
	else:
		return 0

def sphereEuclidean2D(
	coord1: "First coordinate, in (lat, lon)", 
	coord2: "Second coordinate, in (lat, lon)"
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

