###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 02/01/2021                             #
# consTSP() - Constructive Heuristic for TSP                                  #
# - 11/09/2020 DepthFirst Method                                              #
# - 11/09/2020 Christofides Algorithm                                         #
# - 11/11/2020 Nearest Neighborhood Method                                    #
# - 11/11/2020 Farthest Neighborhood Method                                   #
# - 11/11/2020 Random Sequence                                                #
# consVRP() - Constructive Heuristic for VRP, and typical variants            #
# - 12/03/2020 Clark Wright Saving Algorithm                                  #
###############################################################################

import heapq

from .common import *
from .graph import *

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

def consVRP(
	nodes: 			"Dictionary, returns the detail info of given nodeID, \
						{\
							nodeID1: {'loc': (x, y), 'demand': d}, \
							nodeID2: {'loc': (x, y), 'demand': d}, \
							... \
						}" = None, 
	depotID:		"Node ID for depot" = 0, 
	customerID:		"1) String (default) 'NoFisrt', or \
					 2) A list of node IDs" = 'NoFirst',
	edges:			"1) String (default) 'Euclidean' or \
					 2) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
	vehCap:			"Capacity of each vehicle" = None,
	vehNum:			"Number of vehicles" = None,
	fixVehNum: 		"Boolean, if number of vehicles is fixed" = False,
	algo:			"1) String 'CWSaving' or \
					 2) String (not available) 'Sweep' or \
					 3) String (not available) 'CMT'" = 'CWSaving'
	) -> "Exact solution for VRP":

	# FIXME/TDL ===============================================================
	# 1. Need to further check for infeasibility
	# 2. Check for asymmetric/symmetric VRP

	# Define nodes ============================================================
	if (customerID == 'NoFirst'):
		if (depotID not in nodes):
			print("Error: Cannot find depot in nodes")
			return None
		else:
			customerID = [i for i in nodes if i != depotID]
	else:
		defNodes = {}
		if (depotID not in nodes):
			print("Error: Cannot find depot in nodes")
			return None
		else:
			defNodes[depotID] = nodes[depotID]
			for i in customerID:
				if (i not in nodes):
					print("Error: Cannot find customer %s in nodes" % i)
					return None
				else:
					defNodes[i] = nodes[i]
		nodes = defNodes

	# Define edges ============================================================
	if (type(edges) is not dict):
		if (edges == 'Euclidean'):
			edges = {}
			for i in customerID:
				for j in customerID:
					if (i != j):
						edges[i, j] = euclidean2D(nodes[i]['loc'], nodes[j]['loc'])
			for i in customerID:
				edges[depotID, i] = euclidean2D(nodes[depotID]['loc'], nodes[i]['loc'])
				edges[i, depotID] = euclidean2D(nodes[i]['loc'], nodes[depotID]['loc'])
		else:
			print("Error: Incorrect type `edges`")
			return None
	
	# Solve by different formulations =========================================
	res = None
	if (algo == 'CWSaving'):
		res = _consVRPClarkeWright(nodes, depotID, customerID, edges, vehCap, vehNum, fixVehNum)
	else:
		print("Error: Incorrect or unavailable CVRP formulation option!")

	return res

def _consVRPClarkeWright(nodes, depotID, customerID, edges, vehCap, vehNum, fixVehNum):
	# Initial routes ==========================================================
	routes = {}
	for i in range(len(customerID)):
		routes[i] = {
			'route': [depotID, customerID[i], depotID],
			'demand': nodes[customerID[i]]['demand'],
			'length': 2 * edges[depotID, customerID[i]]
		}

	# Initial saving raking ===================================================
	rankSaving = []
	for i in customerID:
		for j in customerID:
			if (i != j):
				# Calculate saving for each pair
				sav = edges[depotID, i] + edges[depotID, j] - edges[i, j]
				# heapq returns the smallest, so add a negative sign
				heapq.heappush(rankSaving, (-sav, (i, j)))

	# Merge routes subroutine =================================================
	def merge(i, j):
		if (i == j):
			return None
		rI = None
		rJ = None
		iLeft = None
		iRight = None
		jLeft = None
		jRight = None
		for r in routes:
			if (i == routes[r]['route'][1]):
				iLeft = True
				rI = r
			if (i == routes[r]['route'][-2]):
				iRight = True
				rI = r
			if (j == routes[r]['route'][1]):
				jLeft = True
				rJ = r
			if (j == routes[r]['route'][-2]):
				jRight = True
				rJ = r
		newRoute = []
		if (iRight == True and jLeft == True):
			newRoute = [i for i in routes[rI]['route']]
			addRoute = [i for i in routes[rJ]['route']]
			newRoute.extend(addRoute)
		elif (iLeft == True and jRight == True):
			newRoute = [i for i in routes[rJ]['route']]
			addRoute = [i for i in routes[rI]['route']]
			newRoute.extend(addRoute)
		elif (iLeft == True and jLeft == True):
			newRoute = [i for i in routes[rI]['route']]
			newRoute.reverse()
			addRoute = [i for i in routes[rJ]['route']]
			newRoute.extend(addRoute)
		elif (iRight == True and jRight == True):
			newRoute = [i for i in routes[rI]['route']]
			addRoute = [i for i in routes[rJ]['route']]
			addRoute.reverse()
			newRoute.extend(addRoute)

		while (depotID in newRoute):
			newRoute.remove(depotID)
		newRoute.insert(0, depotID)
		newRoute.append(depotID)

		newDemand = routes[rI]['demand'] + routes[rJ]['demand']
		newLength = routes[rI]['length'] + routes[rJ]['length'] + edges[i, j] - edges[depotID, i] - edges[depotID, j]
		routes.pop(rI)
		routes.pop(rJ)
		newRouteIndex = max(list(routes.keys())) + 1
		routes[newRouteIndex] = {
			'route': newRoute,
			'demand': newDemand,
			'length': newLength
		}

	# Merge routes ============================================================
	while (len(rankSaving) > 0 and ((not fixVehNum and -rankSaving[0][0] > 0) or (fixVehNum and len(routes.keys()) > vehNum))):
		# Get the biggest saving
		bestSaving = heapq.heappop(rankSaving)
		# Flip it back
		sav = -bestSaving[0]
		# If there is saving, check which two routes can be merged
		routeI = None
		routeJ = None
		for r in routes:
			if (bestSaving[1][0] == routes[r]['route'][1] or bestSaving[1][0] == routes[r]['route'][-2]):
				routeI = r
			if (bestSaving[1][1] == routes[r]['route'][1] or bestSaving[1][1] == routes[r]['route'][-2]):
				routeJ = r
			if (routeI != None and routeJ != None):
				break
		# Two routes has to be different, and satisfied the capacity
		if (routeI != None and routeJ != None and routeI != routeJ and routes[routeI]['demand'] + routes[routeJ]['demand'] <= vehCap):
			merge(bestSaving[1][0], bestSaving[1][1])

	# Rename the route name ===================================================
	ofv = 0
	route = {}
	acc = 1
	for r in routes:
		ofv += routes[r]['length']
		route[acc] = [i for i in routes[r]['route']]
		acc += 1

	return {
		'ofv': ofv,
		'route': route
	}

def localTSP():

	return

def consTSPPD(
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	edges:		"1) String (default) 'Euclidean' or \
				 2) String 'SphereEuclidean' or \
				 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
	depot:		"Node ID of depot" = None,
	reqs:		"Dictionary, returns the pairs of pickup and delivery node IDs, \
					{\
						reqID1: {'pickup': nodeID, 'delivery': nodeID, 'size': size}, \
						reqID2: {'pickup': nodeID, 'delivery': nodeID, 'size': size}, \
						... \
					}" = None,
	algo:		"1) String 'CheapestFeasibleInsertion' or \
				 2) String 'Mosheiov99' or \
				 3) " = None,
	) -> "":

	return res

# [Constructing]
def _consTSPPDCheapestFeasibleInsertion(nodeIDs, edges, reqs, veh):
	ofv = 0
	seq = []

	# Get all delivery locs ===================================================
	deliveryIDs = []
	for r in reqs:
		deliveryIDs.append(reqs[r]['delivery'])

	# Create a route using all delivered locs =================================
	seq = _consTSPNearestNeighbor(deliveryIDs, edges)

	# Create initial action list with all delivery actions ====================
	# for i in 

	# Subroutine to cal the cost of inserting a pickup location ===============
	def insertPickup(req, actions):

		return {
			'newSeq': newSeq
		}

	# Cheapest insertion, each iteration insert a pickup location =============
	finishedReqs = []
	updatedSeq = [i for i in deliveryIDs]
	while (len(finishedReqs) < len(reqs)):
		tmpCost = None
		tmpUpdatedSeq = None
		# For each request, try to insert the pickup location
		for r in reqs:
			if (r not in finishedReqs):
				pickupNode = reqs[r]['pickup']
				res = insertPickup(r, updatedSeq)
				# If this request is feasible to be inserted and the cost is cheaper
				if (res['feasible'] and (res['cost'] < tmpCost or tmpCost == None)):
					tmpCost = res['cost']
					tmpUpdatedSeq = res['newSeq']
				# If at least one request be added into the seq
				if (tmpCost != None):
					updatedSeq = tmpUpdatedSeq
					finishedReqs.append(r)

	return {
		'ofv': ofv,
		'seq': seq
	}

