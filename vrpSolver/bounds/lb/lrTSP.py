###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 11/09/2020                             #
# lrTSP - Use Lagrangian Relaxation to Give Held & Karp Bound                 #
# - 11/09/2020 Add Held & Karp Bound                                          #
###############################################################################

import heapq
import math
from gurobipy import *

from vrpSolver.const import *
from vrpSolver.common import *
from vrpSolver.graph.basic import *
from vrpSolver.graph.mst import *

def lrTSP(
	nodeLoc:	"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: (lat, lon), \
						nodeID2: (lat, lon), \
						... \
					}" = None, 
	tau:		"1) String 'Euclidean' or \
				 2) String (default) 'SphereEuclidean' or \
				 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "SphereEuclidean",
	nodeIDs:	"1) String (default) 'All', or \
				 2) A list of node IDs" = 'All',
	subgradM:	"Double" = 1,
	subgradRho:	"Double, (0, 1)" = 0.95,
	stopType:	"1) String, (default) 'Epsilon' (`stopEpsilon` will be used) or \
				 2) String, 'IterationNum' (`stopK` will be used) or \
				 3) String, 'Runtime' (`stopTime` will be used)" = 'Epsilon',
	stopEpsilon:"Double, small number" = 0.01,
	stopK:		"Integer, large number" = 200,
	stopTime:	"Double, in seconds" = 600
	) -> "Returns a Held & Karp lower bound of the TSP using Lagrangian Relaxation":

	# Define nodeIDs ==========================================================
	if (type(nodeIDs) is not list):
		if (nodeIDs == 'All'):
			nodeIDs = []
			for i in nodeLoc:
				nodeIDs.append(i)

	# Define tau ==============================================================
	if (type(tau) is not dict):
		lstNodeID = nodeIDs.copy()
		if (tau == 'Euclidean'):
			tau = getTauEuclidean(nodeLoc, lstNodeID)
		elif (tau == 'SphereEuclidean'):
			tau = getTauSphereEuclidean(nodeLoc, lstNodeID)
		else:
			print("Error: Incorrect type `tau`")
			return None

	# Initialize ==============================================================
	k = 0
	u = [0 for i in range(len(nodeIDs))]
	d = None
	costSum = None
	L = None
	oldL = None

	# Calculate 1 tree ========================================================
	def cal1Tree(weightArcs):
		# Separate first node
		arcsWithVertexOne = []
		arcsWithoutVertexOne = []
		for i in range(len(weightArcs)):
			if (weightArcs[i][0] == 0 or weightArcs[i][1] == 0):
				arcsWithVertexOne.append(weightArcs[i])
			else:
				arcsWithoutVertexOne.append(weightArcs[i])

		# MST for the rest of vertices
		mst = graphMST(arcsWithoutVertexOne)['mst']

		# Find two cheapest arcs to vertex one
		sortedArcswithVertexOne = []
		for i in range(len(arcsWithVertexOne)):
			heapq.heappush(sortedArcswithVertexOne, (arcsWithVertexOne[i][2], arcsWithVertexOne[i]))

		# Build 1-tree
		leastTwo = []
		leastTwo.append(heapq.heappop(sortedArcswithVertexOne))
		leastTwo.append(heapq.heappop(sortedArcswithVertexOne))

		m1t = [i for i in mst]
		m1t.append(leastTwo[0][1])
		m1t.append(leastTwo[1][1])

		# Calculate total cost
		costSum = 0
		for i in range(len(m1t)):
			costSum += m1t[i][2]

		# Arcs to neighbors
		neighbors = convertArcs2Neighbor(m1t)
		d = []
		for i in range(len(nodeIDs)):
			d.append(2 - len(neighbors[i]))

		return {
			'costSum': costSum,
			'm1t': m1t,
			'd': d
		}

	# Main iteration ==========================================================
	continueFlag = True
	while (continueFlag):
		# Update cost of each edge
		weightArcs = []
		for i in range(len(nodeIDs)):
			for j in range(len(nodeIDs)):
				if (i != None and j != None and i < j):
					weightArcs.append((i, j, tau[i, j] - u[i] - u[j]))

		# Calculate 1-tree
		oneTree = cal1Tree(weightArcs)

		# Update L and d
		costSum = oneTree['costSum']
		m1t = oneTree['m1t']
		uSum = sum(u)
		if (L != None):
			oldL = L
		L = costSum + 2 * uSum
		d = oneTree['d']

		# update u
		oldU = [i for i in u]
		u = []
		eff = subgradM * math.pow(subgradRho, k)
		for i in range(len(nodeIDs)):
			u.append(oldU[i] + eff * d[i])

		# Check if continue
		def allZero(d):
			for i in d:
				if (i != 0):
					return False
			return True
		if (k >= stopK):
			continueFlag = False
		elif (oldL != None and abs(oldL - L) < stopEpsilon):
			continueFlag = False
		elif (allZero(d)):
			continueFlag = False
		else:
			k += 1

	return {
		'lrLowerBound': costSum
	}
