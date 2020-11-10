###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 11/07/2020                             #
# graphMST() - Minimum Spanning Tree                                          #
# - 11/07/2020 Add Krusal Algorithm                                           #
###############################################################################

import heapq

from vrpSolver.graph.basic import *

def graphMST(
	weightArcs:	"A list of 3-tuples, (ID1, ID2, weight), indexes of vertices must start from 0" = None,
	numVertices:"1) Integer, number of vertices, or \
				 2) (default) None, assuming all vertices are mentioned in `weightArcs`" = None,
	algo:		"1) String, (default) 'Krusal' or, \
				 2) String, (not available) 'Prim_AdjList' or, \
				 3) String, (not available) 'Prim_AdjMat' or, \
				 4) String, (not available) 'Boruvka' or, \
				 5) String, (not available) 'ReverseDelete'" = 'Krusal'
	) -> "A list of weightArcs which forms a minimal spanning tree":

	# Number of vertices ======================================================
	if (numVertices == None):
		vertices = []
		for i in range(len(weightArcs)):
			if (weightArcs[i][0] not in vertices):
				vertices.append(weightArcs[i][0])
			if (weightArcs[i][1] not in vertices):
				vertices.append(weightArcs[i][1])
		numVertices = len(vertices)

	# Call MST ================================================================
	if (algo == 'Krusal'):
		res = mstKrusal(weightArcs, numVertices)
	else:
		print("Error: Incorrect or not available MST option!")
	return res

def mstKrusal(weightArcs, numVertices):
	# Initialize ==============================================================
	mst = []
	val = 0
	compList = []

	# Arc ranking =============================================================
	sortedWeightArcs = []
	for i in range(len(weightArcs)):
		heapq.heappush(sortedWeightArcs, (weightArcs[i][2], weightArcs[i]))
	
	# Krusal algorithm, add weightArcs between components =====================
	while(len(mst) < numVertices - 1 and sortedWeightArcs):	
		# Uninserted arc with minimal weight
		currArc = heapq.heappop(sortedWeightArcs)

		# Mark two nodes
		nodeID1 = currArc[1][0]
		nodeID2 = currArc[1][1]
		weight = currArc[1][2]
		compID1 = None
		compID2 = None
		findNodeFlag1 = False
		findNodeFlag2 = False

		# Find component that nodes belong to
		for i in range(len(compList)):
			if (nodeID1 in compList[i]):
				findNodeFlag1 = True
				compID1 = i
			if (nodeID2 in compList[i]):
				findNodeFlag2 = True
				compID2 = i
			if (findNodeFlag1 and findNodeFlag2):
				break

		# If two nodes are not in the same component, merge components
		if ((not findNodeFlag1) and (not findNodeFlag2)):
			mst.append(currArc[1])
			val += weight
			compList.append([nodeID1, nodeID2])
		elif (findNodeFlag1 and (not findNodeFlag2)):
			mst.append(currArc[1])
			val += weight
			compList[compID1].append(nodeID2)
		elif ((not findNodeFlag1) and findNodeFlag2):
			mst.append(currArc[1])
			val += weight
			compList[compID2].append(nodeID1)
		elif (findNodeFlag1 and findNodeFlag2):
			if (compID1 != compID2):
				mst.append(currArc[1])
				val += weight
				compList[compID1].extend(compList[compID2].copy())
				compList.remove(compList[compID2])

	return {
		'mst': mst,
		'value': val
	}

def mstPrimAdjList(weightArcs):

	return {
		'mst': mst,
		'value': val
	}

def mstPrimAdjMat(weightArcs):

	return {
		'mst': mst,
		'value': val
	}

def mstBoruvka(weightArcs):

	return {
		'mst': mst,
		'value': val
	}

def mstReverse(weightArcs):

	return {
		'mst': mst,
		'value': val
	}