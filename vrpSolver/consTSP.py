import heapq

from .common import *
from .mst import *
from .matching import *
from .traversal import *

def consTSP(
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	edges:		"1) String (default) 'Euclidean' or \
				 2) String 'SphereEuclidean' or \
				 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
	algo:		"1) String 'NearestNeighbor' or \
				 2) String 'FarthestNeighbor' or \
				 3) String (not available) 'Insertion' or \
				 4) String (not available) 'Patching' or \
				 5) String (not available) 'Sweep' or \
				 6) String 'DepthFirst' or \
				 7) String (default) 'Christofides' or \
				 8) String 'Random'" = 'Christofides'
	) -> "Constructive heuristic solution for TSP":

	# Define edges ============================================================
	if (type(edges) is not dict):
		if (edges == 'Euclidean'):
			edges = getTauEuclidean(nodes)
		elif (edges == 'SphereEuclidean'):
			edges = getTauSphereEuclidean(nodes)
		else:
			print("Error: Incorrect type `edges`")
			return None

	# Heuristics that don't need to transform arc representation ==============
	res = None
	nodeIDs = list(nodes.keys())
	if (algo == 'NearestNeighbor'):
		res = _consTSPNearestNeighbor(nodeIDs, edges)
	elif (algo == 'FarthestNeighbor'):
		res = _consTSPFarthestNeighbor(nodeIDs, edges)
	elif (algo == 'Random'):
		res = _consTSPRandomSeq(nodeIDs, edges)
	else:
		pass
	if (res != None):
		return res

	# Create arcs =============================================================
	# FIXME! indexes of nodes starts from 0 here!
	# FIXME! Add mapping between 0 to n and nodeIDs
	weightArcs = []
	for (i, j) in edges:
		if (i != None and j != None and i < j):
			weightArcs.append((i, j, edges[i, j]))

	# Constructive Heuristics for TSP =========================================
	res = None
	if (algo == 'DepthFirst'):
		res = _consTSPDepthFirst(weightArcs)
	elif (algo == 'Christofides'):
		res = _consTSPChristofides(weightArcs)

	return res

def _consTSPRandomSeq(nodeIDs, edges):
	# Get random seq ==========================================================
	seqIndex = rndSeq(len(nodeIDs), closed=True)
	seq = []
	for i in range(len(seqIndex)):
		seq.append(nodeIDs[seqIndex[i]])

	# Calculate Ofv ===========================================================
	ofv = calSeqCostMatrix(edges, seq)
	return {
		'ofv': ofv,
		'seq': seq
	}

def _consTSPNearestNeighbor(nodeIDs, edges):
	# Initialize ==============================================================
	seq = [nodeIDs[0]]
	remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
	ofv = 0

	# Accumulate seq ==========================================================
	while (len(remain) > 0):
		nextLeng = None
		nextID = None
		for node in remain:
			if ((node, seq[-1]) in edges):
				if (nextLeng == None or edges[node, seq[-1]] < nextLeng):
					nextID = node
					nextLeng = edges[node, seq[-1]]
			elif ((seq[-1], node) in edges):
				if (nextLeng == None or edges[seq[-1], node] < nextLeng):
					nextID = node
					nextLeng = edges[seq[-1], node]
		seq.append(nextID)
		remain.remove(nextID)
		ofv += nextLeng
	ofv += edges[seq[0], seq[-1]]
	seq.append(seq[0])

	return {
		'ofv': ofv,
		'seq': seq
	}

def _consTSPFarthestNeighbor(nodeIDs, edges):
	# Initialize ==============================================================
	seq = [nodeIDs[0]]
	remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
	ofv = 0

	# Accumulate seq ==========================================================
	while (len(remain) > 0):
		nextLeng = None
		nextID = None
		for node in remain:
			if ((node, seq[-1]) in edges):
				if (nextLeng == None or edges[node, seq[-1]] > nextLeng):
					nextID = node
					nextLeng = edges[node, seq[-1]]
			elif ((seq[-1], node) in edges):
				if (nextLeng == None or edges[seq[-1], node] > nextLeng):
					nextID = node
					nextLeng = edges[seq[-1], node]
		seq.append(nextID)
		remain.remove(nextID)
		ofv += nextLeng	
	ofv += edges[seq[0], seq[-1]]
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
		'seq': seq
	}
