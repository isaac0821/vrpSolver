import math

from vrpSolver.const import *

def getTauEuclidean(
	dicNodeLoc: "Dictionary, {nodeID1: (x, y), ...}",
	lstNodeID: "List, [nodeID1, nodeID2, ...]"
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
	dicNodeLoc: "Dictionary, {nodeID1: (lat, lon), ...}",
	lstNodeID: "List, [nodeID1, nodeID2, ...]"
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

