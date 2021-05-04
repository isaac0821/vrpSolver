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

def findComponentsUndirected(
	arcs:		"A list of 2-tuples", 
	) -> "A list of components":
	# Create adj list, each vertex start with an empty list
	adjList = {}
	for e in arcs:
		if (e[0] not in adjList):
			adjList[e[0]] = [e[1]]
		else:
			adjList[e[0]].append(e[1])
		if (e[1] not in adjList):
			adjList[e[1]] = [e[0]]
		else:
			adjList[e[1]].append(e[0])

	# Initialize
	found = {}
	for node in adjList:
		found[node] = 0
	components = []

	# Main algorithm, mark neighbors
	for i in adjList:
		comp = []
		q = []
		if (found[i] == 0):
			found[i] = 1
			comp.append(i)
			q.append(i)
			while (q):
				v = q.pop(0)
				for u in adjList[v]:
					if (found[u] == 0):
						found[u] = 1
						comp.append(u)
						q.append(u)
			components.append(comp)

	return components

def calSeqCostArcs(
	weightArcs:	"A list of 3-tuple (nodeID1, nodeID2, weight)",
	seq:		"List, sequence of visiting node ids"
	) -> "Return the cost on the graph given a list of arcs weights":

	# Accumulate costs ========================================================
	cost = 0
	for i in range(len(seq) - 1):
		c = None
		for j in range(len(weightArcs)):
			if (seq[i] == weightArcs[j][0] and seq[i + 1] == weightArcs[j][1]):
				c = weightArcs[j][2]
				break
			elif (seq[i] == weightArcs[j][1] and seq[i + 1] == weightArcs[j][0]):
				c = weightArcs[j][2]
				break
		if (c == None):
			print("Error: Missing arc (%s, %s) in `weightArcs`" % (seq[i], seq[i + 1]))
			return
		else:
			cost += c

	return cost

def calSeqCostMatrix(
	tau: "Dictionary {(nodeID1, nodeID2): dist, ...}", 
	seq: "List, sequence of visiting node ids"
	) -> "Return the cost on the graph given cost matrix/dictionary tau":

	# Accumulate costs ========================================================
	cost = 0
	for i in range(len(seq) - 1):
		cost += tau[seq[i], seq[i + 1]]

	return cost

def convertArcs2Neighbor(
	arcs:	"1) A list of 3-tuple (nodeID1, nodeID2, weight) or, \
			 2) A list of 2-tuple (nodeID1, nodeID2)"
	) -> "Dictionary of neighbors of each node":

	neighbors = {}
	for i in range(len(arcs)):
		if (arcs[i][0] not in neighbors):
			neighbors[arcs[i][0]] = [arcs[i][1]]
		else:
			neighbors[arcs[i][0]].append(arcs[i][1])
		if (arcs[i][1] not in neighbors):
			neighbors[arcs[i][1]] = [arcs[i][0]]
		else:
			neighbors[arcs[i][1]].append(arcs[i][0])

	return neighbors

