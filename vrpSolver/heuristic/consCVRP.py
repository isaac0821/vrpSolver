###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 12/04/2020                             #
# consVRP() - Constructive Heuristic for VRP, and typical variants            #
# - 12/03/2020 Clark Wright Saving Algorithm                                  #
###############################################################################

import heapq

from vrpSolver.common import *
from vrpSolver.graph.mst import *
from vrpSolver.graph.basic import *
from vrpSolver.graph.search import *
from vrpSolver.graph.matching import *

def consCVRP(
	nodes: 			"Dictionary, returns the detail info of given nodeID, \
						{\
							nodeID1: {'loc': (x, y), 'demand': d}, \
							nodeID2: {'loc': (x, y), 'demand': d}, \
							... \
						}" = None, 
	depotID:		"Node ID for depot" = 0, 
	customerID:		"1) String (default) 'NoFisrt', or \
					 2) A list of node IDs" = 'NoFirst',
	edges:			"1) String (default) 'Complete_Euclidean' or \
					 2) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Complete_Euclidean",
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
		if (edges == 'Complete_Euclidean'):
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
		res = _consCVRPClarkeWright(nodes, depotID, customerID, edges, vehCap, vehNum, fixVehNum)
	else:
		print("Error: Incorrect or unavailable CVRP formulation option!")

	return res

def _consCVRPClarkeWright(nodes, depotID, customerID, edges, vehCap, vehNum, fixVehNum):
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