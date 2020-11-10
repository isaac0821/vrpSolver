###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 11/09/2020                             #
# consTSP() - Constructive Heuristic for TSP                                  #
# - 11/09/2020 DepthFirst Method                                              #
# - 11/09/2020 Christofides Algorithm                                         #
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
	algo:	"1) String (not available) 'NearestNeighborhood' or \
			 2) String (not available) 'Insertion' or \
			 3) String (not available) 'Patching' or \
			 4) String (not available) 'Sweep' or \
			 5) String 'DepthFirst' or \
			 6) String (default) 'Christofides'" = 'Christofides'
	) -> "Exact solution for TSP":

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