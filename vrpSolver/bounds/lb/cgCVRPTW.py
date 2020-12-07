###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 12/06/2020                             #
# cgCVRPTW - Use column generation to find lower bound of CVRPTW              #
###############################################################################

import math
from gurobipy import *

from vrpSolver.const import *
from vrpSolver.common import *

def cgCVRPTW(
	nodes: 			"Dictionary, returns the detail info of given nodeID, \
						{\
							nodeID1: {'loc': (x, y), 'demand': d, 'tStart': t1, 'tEnd': t2, 'serviceTime': t}, \
							nodeID2: {'loc': (x, y), 'demand': d, 'tStart': t1, 'tEnd': t2, 'serviceTime': t}, \
							... \
						}" = None, 
	depotID:		"Node ID for depot" = 0, 
	customerID:		"1) String (default) 'NoFisrt', or \
					 2) A list of node IDs" = 'NoFirst',
	edges:			"1) String (default) 'Complete_Euclidean' or \
					 2) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Complete_Euclidean",
	vehCap:			"Capacity of each vehicle" = None,
	vehNum:			"Number of vehicles" = None
	) -> "Exact solution for VRP":

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

	# Create adjList for inflow and outflow ===================================
	adjIn = {}
	adjOut = {}
	for e in edges:
		if (e[0] not in adjOut):
			adjOut[e[0]] = [e[1]]
		else:
			adjOut[e[0]].append(e[1])
		if (e[1] not in adjIn):
			adjIn[e[1]] = [e[0]]
		else:
			adjIn[e[1]].append(e[0])

	# Pricing problem =========================================================
	def pricing(pi):
		sub = Model('Pricing')

		# Define decision variables
		w = {}
		for i in nodes:
			for j in nodes:
				if (i != j):
					w[i, j] = sub.addVar(vtype=GRB.BINARY, obj=edges[i, j] - pi[i], name='w_%s_%s' % (i, j))
		s = {}
		for i in nodes:
			s[i] = sub.addVar(vtype=GRB.CONTINUOUS, lb=nodes[i]['tStart'], ub=nodes[i]['tEnd'], name='s_%s' % i)
		dupDepotID = max(list(s.keys())) + 1
		s[dupDepotID] = sub.addVar(vtype=GRB.CONTINUOUS, lb=nodes[depotID]['tStart'], ub=nodes[depotID]['tEnd'], name='s_dupDepot')
		sub.update()

		# Assignment constraints
		for i in customerID:
			sub.addConstr(quicksum(w[i, j] for j in adjOut[i]) == quicksum(w[j, i] for j in adjIn[i]))
		sub.addConstr(quicksum(w[depotID, j] for j in adjOut[depotID]) == 1)
		sub.addConstr(quicksum(w[j, depotID] for j in adjIn[depotID]) == 1)
		sub.update()

		# Capacity constraints
		sub.addConstr(quicksum(w[i, j] * nodes[j]['demand'] for i, j in edges if (i != j)) <= vehCap)

		# Time windows
		# REMEMBER: Yes I have spent a lot of time here! s[depotID] cannot take two values at the same time
		for i in nodes:
			for j in customerID:
				if (i != j):
					sub.addConstr(s[i] + edges[i, j] + nodes[i]['serviceTime'] - nodes[depotID]['tEnd'] * (1 - w[i, j]) <= s[j])
		for i in customerID:
			sub.addConstr(s[i] + edges[i, depotID] + nodes[i]['serviceTime'] - nodes[depotID]['tEnd'] * (1 - w[i, depotID]) <= s[dupDepotID])
		sub.update()

		# Solve
		# sub.write('sub.lp')
		sub.setParam("OutputFlag", 0)
		sub.modelSense = GRB.MINIMIZE
		sub.optimize()

		# Construct solution
		ofv = None
		c = 0
		a = []
		route = []
		if (sub.status == GRB.status.OPTIMAL):
			ofv = sub.getObjective().getValue()
			arcSet = []
			for i, j in w:
				if (w[i, j].x > 0.9):
					c += edges[i, j]
					arcSet.append((i, j))

			route = [depotID]
			while (len(arcSet) > 0):
				cur = None
				for arc in arcSet:					
					if (arc[0] == route[-1]):
						route.append(arc[1])
						arcSet.remove(arc)
						cur = arc[1]
						break
					if (arc[1] == route[-1]):
						route.append(arc[0])
						arcSet.remove(arc)
						cur = arc[0]
						break
				if (cur == depotID):
					break
			a = {}
			for i in customerID:
				if (i in route):
					a[i] = 1
				else:
					a[i] = 0
		else:
			print("WTF???")
		return {
			'ofv': ofv,
			'cr': c,
			'ai': a,
			'route': route
		}

	# Master problem initialization ===========================================
	CVRPTW = Model('CVRPTW')
	# DV, routes
	y = {} # Binary for candidate route, y[routeID], start with n routes
	
	# Parameters
	c = {} # Cost of route, c[routeID]
	a = {} # Whether node in route, a[nodeID, routeID]

	# Route sets
	routes = {}

	# Initialize parameters
	for i in customerID:
		for j in customerID:
			if (i != j):
				a[i, j] = 0
			else:
				a[i, i] = 1
		c[i] = edges[depotID, i] + edges[i, depotID]
		routes[i] = [depotID, i, depotID]
	acc = max(customerID) + 1

	# Initial columns
	for i in customerID:
		y[i] = CVRPTW.addVar(vtype=GRB.CONTINUOUS, obj=c[i])

	# Initial constraints
	cons = {}
	cons[depotID] = CVRPTW.addConstr(quicksum(-y[r] for r in y) >= -vehNum)
	for i in customerID:
		cons[i] = CVRPTW.addConstr(quicksum(a[i, r] * y[r] for r in y) == 1)
	CVRPTW.update()

	# Now solve the Master (Set Partition Formulation) problem ================
	canAddVarFlag = True
	while(canAddVarFlag):
		canAddVarFlag = False

		# Optimize
		CVRPTW.setParam("OutputFlag", 0)
		CVRPTW.modelSense = GRB.MINIMIZE
		CVRPTW.optimize()
		
		if (CVRPTW.status == GRB.status.OPTIMAL):
			# Solve subproblem
			pi = {}
			for constraint in cons:
				pi[constraint] = cons[constraint].Pi
			subproblem = pricing(pi)
			print(subproblem['route'])
			
			# Add route into route set
			if (cons[depotID].Pi - subproblem['ofv'] > CONST_EPSILON):
				canAddVarFlag = True
				newRouteIndex = max(list(y.keys())) + 1
				c[newRouteIndex] = subproblem['cr']
				y[newRouteIndex] = CVRPTW.addVar(vtype=GRB.CONTINUOUS, obj=c[newRouteIndex])
				CVRPTW.update()
				for i in customerID:
					a[i, newRouteIndex] = subproblem['ai'][i]
				routes[newRouteIndex] = subproblem['route']

				# Update columns, add one more
				CVRPTW.chgCoeff(cons[depotID], y[newRouteIndex], -1)
				for i in customerID:
					CVRPTW.chgCoeff(cons[i], y[newRouteIndex], a[i, newRouteIndex])
				CVRPTW.update()
			else:
				canAddVarFlag = False
		else:
			break

	# Interpret solution for lower bound ======================================
	lb = CVRPTW.getObjective().getValue()

	# Early branching heuristic ===============================================
	for i in y:
		y[i].vtype = GRB.BINARY
	CVRPTW.update()
	CVRPTW.optimize()

	# Interpret solution ======================================================
	solRoute = {}
	acc = 1
	for i in y:
		if (y[i].x > 0.9):			
			solRoute[acc] = {
				'route': routes[i],
				'length': c[i]
			}
			acc += 1
	ub = CVRPTW.getObjective().getValue()

	return {
		'lb': lb,
		'ub': ub,
		'routes': solRoute
	}


