###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 12/04/2020                             #
# ipCVRP - Basic Integer Programming Method for Solving CVRP                  #
# - 12/03/2020 Add Golden77                                                   #
# - 12/04/2020 Add Two-Index Flow                                             #
###############################################################################

from gurobipy import *

from vrpSolver.common import *
from vrpSolver.graph.basic import *

def ipCVRP(
	nodes: 			"Dictionary, returns the detail info of given nodeID, \
						{\
							nodeID1: {'loc': (x, y), 'demand': d, 'tStart': t1, 'tEnd': t2}, \
							nodeID2: {'loc': (x, y), 'demand': d, 'tStart': t1, 'tEnd': t2}, \
							... \
						}" = None, 
	depotID:		"Node ID for depot" = 0, 
	customerID:		"1) String (default) 'NoFisrt', or \
					 2) A list of node IDs" = 'NoFirst',
	edges:			"1) String (default) 'Complete_Euclidean' or \
					 2) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Complete_Euclidean",
	vehCap:			"Capacity of each vehicle" = None,
	vehNum:			"Number of vehicles" = None,
	fml:			"1) String 'Golden77' or \
					 2) String 'Two-Index' or \
					 3) String (not available) 'MultiCommodity Flow' or \
					 4) String (not available) 'Set Partitioning'" = 'Golden77',

	cutoffTime: 	"1) Double, in seconds or \
					 2) (default) None, no time limit" = None
	) -> "Exact solution for TSP":

	# Define nodes ==========================================================
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

	# Define edges ==========================================================
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
	if (fml == 'Golden77'):
		res = _ipCVRPGolden77(nodes, depotID, customerID, edges, vehCap, vehNum, cutoffTime)
	if (fml == 'Two-Index'):
		res = _ipCVRPTwoIndex(nodes, depotID, customerID, edges, vehCap, vehNum, cutoffTime)
	else:
		print("Error: Incorrect or unavailable CVRP formulation option!")

	return res

def _ipCVRPTwoIndex(nodes, depotID, customerID, edges, vehCap, vehNum, cutoffTime):
	CVRP = Model('CVRP')

	# Decision variables ======================================================
	x = {}
	for i, j in edges:
		if (i != depotID and j != depotID and i < j):
			x[i, j] = CVRP.addVar(
				vtype = GRB.BINARY,
				obj = edges[i, j], name = 'x_%s_%s' % (i, j))
	for j in customerID:
		x[depotID, j] = CVRP.addVar(
			vtype=GRB.INTEGER,
			ub=2,
			obj=edges[depotID, j], name = 'x_%s_%s' % (depotID, j))

	# CVRP objective function =================================================
	CVRP.modelSense = GRB.MINIMIZE
	CVRP.Params.lazyConstraints = 1
	CVRP.update()

	# Leaving depot ===========================================================
	CVRP.addConstr(quicksum(x[depotID, j] for j in customerID) == 2 * vehNum, name='leaving')

	# Balance constraint ======================================================
	for cus in customerID:
		CVRP.addConstr(quicksum(x[i, k] for (i, k) in edges if (i < k and k == cus)) + quicksum(x[k, j] for (k, j) in edges if (k < j and k == cus)) == 2, name='balance_%s' % cus)

	# Subtour elimination =====================================================
	CVRP._x = x
	def subtourelim(model, where):
		if (where == GRB.Callback.MIPSOL):
			x_sol = model.cbGetSolution(model._x)
			arcs = [(i, j) for i, j in x if (i != depotID and j != depotID and i < j and x_sol[i, j] > 0.9)]
			components = findComponentsUndirected(arcs)
			for comp in components:
				sumQ = 0
				for n in comp:
					sumQ += nodes[n]['demand']
				vS = math.ceil(sumQ / float(vehCap))
				edgesInComp = []
				for (i, j) in x:
					if (i in comp and j in comp):
						edgesInComp.append((i, j))
				model.cbLazy(quicksum(x[i, j] for (i, j) in edgesInComp if i < j) <= len(comp) - vS)

	# CVRP with callback ======================================================
	if (cutoffTime != None):
		CVRP.setParam(GRB.Param.TimeLimit, cutoffTime)
	CVRP.optimize(subtourelim)

	# Reconstruct solution ====================================================
	route = {}
	ofv = None
	gap = None
	lb = None
	ub = None
	runtime = None
	if (CVRP.status == GRB.status.OPTIMAL):		
		ofv = CVRP.getObjective().getValue()
		arcSet = []
		for i, j in x:
			if (x[i, j].X > 0.9 and x[i, j].X < 1.1):
				arcSet.append((i, j))
			elif (x[i, j].X > 1.9):
				arcSet.append((i, j))
				arcSet.append((j, i))
		print(arcSet)
		for k in range(1, vehNum + 1):
			route[k] = [depotID]
			while (len(arcSet) > 0):
				cur = None
				for arc in arcSet:
					if (arc[0] == route[k][-1]):
						route[k].append(arc[1])
						arcSet.remove(arc)
						cur = arc[1]
						break
					if (arc[1] == route[k][-1]):
						route[k].append(arc[0])
						arcSet.remove(arc)
						cur = arc[0]
						break
				if (cur == depotID):
					break
		gap = 0
		lb = ofv
		ub = ofv
		runtime = CVRP.runtime
	elif (CVRP.status == GRB.status.TIME_LIMIT):
		ofv = None
		route = {}
		gap = CVRP.MIPGap
		lb = CVRP.ObjBoundC
		ub = CVRP.ObjVal
		runtime = CVRP.Runtime

	return {
		'ofv': ofv,
		'route': route,
		'gap': gap,
		'lb': lb,
		'ub': ub,
		'runtime': runtime
	}

def _ipCVRPGolden77(nodes, depotID, customerID, edges, vehCap, vehNum, cutoffTime):
	CVRP = Model('CVRP')

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

	# Decision variables  =====================================================
	vehicles = [i for i in range(vehNum)]
	x = {}
	for k in vehicles:
		for i, j in edges:
			x[k, i, j] = CVRP.addVar(
				vtype = GRB.BINARY,
				obj = edges[i, j])

	# CVRP objective function =================================================
	CVRP.modelSense = GRB.MINIMIZE
	CVRP.Params.lazyConstraints = 1
	CVRP.update()

	# Leaving Node ============================================================
	for i in customerID:
		CVRP.addConstr(quicksum(quicksum(x[k, i, j] for j in adjOut[i]) for k in vehicles) == 1)

	# Balance Constraint ======================================================
	for i in nodes:
		for k in vehicles:
			CVRP.addConstr(quicksum(x[k, i, j] for j in adjOut[i]) == quicksum(x[k, j, i] for j in adjIn[i]))

	# Vehicle selecting route =================================================
	for k in vehicles:
		CVRP.addConstr(quicksum(x[k, depotID, i] for i in adjOut[depotID]) <= 1)
		CVRP.addConstr(quicksum(x[k, i, depotID] for i in adjIn[depotID]) <= 1)

	# Capacity Constraint =====================================================
	for k in vehicles:
		CVRP.addConstr(quicksum(nodes[i]['demand'] * quicksum(x[k, i, j] for j in adjOut[i]) for i in customerID) <= vehCap)

	# Sub-tour elimination ====================================================
	CVRP._x = x
	def subtourelim(model, where):
		if (where == GRB.Callback.MIPSOL):
			x_sol = model.cbGetSolution(model._x)
			for k in vehicles:
				arcs = [(i, j) for i, j in edges if (x_sol[k, i, j] > 0.9)]
				components = findComponentsUndirected(arcs)
				if (len(components) != 1):
					for component in components:					
						model.cbLazy(quicksum(x[k, i, j] for i in component for j in component if i != j) <= len(component) - 1)
				if (len(components) == 1 and depotID not in components[0]):
					component = components[0]
					model.cbLazy(quicksum(x[k, i, j] for i in component for j in component if i != j) <= len(component) - 1)
	
	# CVRP with callback ======================================================
	if (cutoffTime != None):
		CVRP.setParam(GRB.Param.TimeLimit, cutoffTime)
	CVRP.optimize(subtourelim)

	# Reconstruct solution ====================================================
	route = {}
	ofv = None
	gap = None
	lb = None
	ub = None
	runtime = None
	if (CVRP.status == GRB.status.OPTIMAL):		
		ofv = CVRP.getObjective().getValue()
		for k in vehicles:
			arcSet = []
			for i, j in edges:
				if (x[k, i, j].X > 0.9):
					# print("x[%s, %s, %s] = %s" % (k, i, j, x[k, i, j].X))
					arcSet.append((i, j))
			route[k] = [depotID]
			while (len(arcSet) > 0):
				for arc in arcSet:
					if (arc[0] == route[k][-1]):
						route[k].append(arc[1])
						arcSet.remove(arc)
						break
					if (arc[1] == route[k][-1]):
						route[k].append(arc[0])
						arcSet.remove(arc)
						break
		gap = 0
		lb = ofv
		ub = ofv
		runtime = CVRP.runtime
	elif (CVRP.status == GRB.status.TIME_LIMIT):
		ofv = None
		route = {}
		gap = CVRP.MIPGap
		lb = CVRP.ObjBoundC
		ub = CVRP.ObjVal
		runtime = CVRP.Runtime

	return {
		'ofv': ofv,
		'route': route,
		'gap': gap,
		'lb': lb,
		'ub': ub,
		'runtime': runtime
	}
