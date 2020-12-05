###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 11/09/2020                             #
# consTSP() - Constructive Heuristic for TSP                                  #
# - 11/09/2020 DepthFirst Method                                              #
# - 11/09/2020 Christofides Algorithm                                         #
# - 11/11/2020 Nearest Neighborhood Method                                    #
# - 11/11/2020 Farthest Neighborhood Method                                   #
# - 11/11/2020 Random Sequence                                                #
###############################################################################

from vrpSolver.common import *
from vrpSolver.graph.mst import *
from vrpSolver.graph.basic import *
from vrpSolver.graph.search import *
from vrpSolver.graph.matching import *

def consTSP(
	nodeLoc:"Dictionary, returns the coordinate of given nodeID, \
				{\
					nodeID1: (lat, lon), \
					nodeID2: (lat, lon), \
					... \
				}" = None, 
	tau:	"1) String 'Euclidean' or \
			 2) String (default) 'SphereEuclidean' or \
			 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "SphereEuclidean",
	nodeIDs:"1) String (default) 'All', or \
			 2) A list of node IDs" = 'All',
	algo:	"1) String 'NearestNeighbor' or \
			 2) String 'FarthestNeighbor' or \
			 3) String (not available) 'Insertion' or \
			 4) String (not available) 'Patching' or \
			 5) String (not available) 'Sweep' or \
			 6) String 'DepthFirst' or \
			 7) String (default) 'Christofides' or \
			 8) String 'Random'" = 'Christofides'
	) -> "Heuristic solution for TSP":

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

	# Heuristics that don't need to transform arc representation ==============
	res = None
	if (algo == 'NearestNeighbor'):
		res = _consTSPNearestNeighbor(nodeIDs, tau)
	elif (algo == 'FarthestNeighbor'):
		res = _consTSPFarthestNeighbor(nodeIDs, tau)
	elif (algo == 'Random'):
		res = _consTSPRandomSeq(nodeIDs, tau)
	else:
		pass
	if (res != None):
		return res

	# Create arcs =============================================================
	# FIXME! indexes of nodes starts from 0 here!
	# FIXME! Add mapping between 0 to n and nodeIDs
	weightArcs = []
	for (i, j) in tau:
		if (i != None and j != None and i < j):
			weightArcs.append((i, j, tau[i, j]))

	# Constructive Heuristics for TSP =========================================
	res = None
	if (algo == 'DepthFirst'):
		res = _consTSPDepthFirst(weightArcs)
	elif (algo == 'Christofides'):
		res = _consTSPChristofides(weightArcs)

	return res

def _consTSPRandomSeq(nodeIDs, tau):
	# Get random seq ==========================================================
	seqIndex = rndSeq(len(nodeIDs), closed=True)
	seq = []
	for i in range(len(seqIndex)):
		seq.append(nodeIDs[seqIndex[i]])

	# Calculate Ofv ===========================================================
	ofv = calSeqCostMatrix(tau, seq)
	return {
		'ofv': ofv,
		'seq': seq
	}

def _consTSPNearestNeighbor(nodeIDs, tau):
	# Initialize ==============================================================
	seq = [nodeIDs[0]]
	remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
	ofv = 0

	# Accumulate seq ==========================================================
	while (len(remain) > 0):
		nextLeng = None
		nextID = None
		for node in remain:
			if ((node, seq[-1]) in tau):
				if (nextLeng == None or tau[node, seq[-1]] < nextLeng):
					nextID = node
					nextLeng = tau[node, seq[-1]]
			elif ((seq[-1], node) in tau):
				if (nextLeng == None or tau[seq[-1], node] < nextLeng):
					nextID = node
					nextLeng = tau[seq[-1], node]
		seq.append(nextID)
		remain.remove(nextID)
		ofv += nextLeng
	ofv += tau[seq[0], seq[-1]]
	seq.append(seq[0])

	return {
		'ofv': ofv,
		'seq': seq
	}

def _consTSPFarthestNeighbor(nodeIDs, tau):
	# Initialize ==============================================================
	seq = [nodeIDs[0]]
	remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
	ofv = 0

	# Accumulate seq ==========================================================
	while (len(remain) > 0):
		nextLeng = None
		nextID = None
		for node in remain:
			if ((node, seq[-1]) in tau):
				if (nextLeng == None or tau[node, seq[-1]] > nextLeng):
					nextID = node
					nextLeng = tau[node, seq[-1]]
			elif ((seq[-1], node) in tau):
				if (nextLeng == None or tau[seq[-1], node] > nextLeng):
					nextID = node
					nextLeng = tau[seq[-1], node]
		seq.append(nextID)
		remain.remove(nextID)
		ofv += nextLeng	
	ofv += tau[seq[0], seq[-1]]
	seq.append(seq[0])

	return {
		'ofv': ofv,
		'seq': seq
	}

def _consTSPDepthFirst(weightArcs):
	# Create MST ==============================================================
	mst = graphMST(weightArcs)['mst']

	# Seq of visit is the seq of Depth first search on the MST ================
	seq = traversalGraph(mst)['seq']
	seq.append(seq[0])

	# Calculate ofv ===========================================================
	ofv = calSeqCostArcs(weightArcs, seq)

	return {
		'ofv': ofv,
		'mst': mst,
		'seq': seq
	}

def _consTSPChristofides(weightArcs):
	# Create MST ==============================================================
	mst = graphMST(weightArcs)['mst']

	# Derive subgraph of odd degree vertices ==================================
	neighbors = convertArcs2Neighbor(mst)
	oddDegrees = []
	for node in neighbors:
		if (len(neighbors[node]) % 2 != 0):
			oddDegrees.append(node)
	subGraph = []
	for arc in weightArcs:
		if (arc[0] in oddDegrees and arc[1] in oddDegrees):
			subGraph.append(arc)

	# Find minimum cost matching of the subgraph ==============================
	minMatching = graphMatching(weightArcs=subGraph, mType='Minimum', algo='IP')['matching']

	# Add them back to create a new graph =====================================
	newGraph = []
	for arc in minMatching:
		newGraph.append(arc)
	for arc in mst:
		newGraph.append(arc)

	# Traverse graph and get seq ==============================================
	# Try to find a vertex with degree 1
	rootID = None
	for node in neighbors:
		if (len(neighbors[node]) == 1):
			rootID = node
			break
	seq = traversalGraph(newGraph, rootID=rootID)['seq']
	seq.append(seq[0])

	# Calculate ofv ===========================================================
	ofv = calSeqCostArcs(weightArcs, seq)

	return {
		'ofv': ofv,
		'seq': seq,
		'mst': mst,
		'newGraph': newGraph,
		'matching': minMatching,
	}
	