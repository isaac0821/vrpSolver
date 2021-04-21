###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 01/05/2021                             #
# ipTSP - Basic Integer Programming Method for Solving TSP                    #
# - 01/05/2021 Add DFJ_Lazy                                                   #
# - 01/05/2021 Add DFJ_PlainLoop                                              #
# - 01/05/2021 Add MTZ                                                        #
# - 01/05/2021 Add MultiCommodityFlow                                         #
# - 01/05/2021 Add ShortestPath                                               #
# - 01/05/2021 Add QAP                                                        #
# lrTSP - Use Lagrangian Relaxation to Give Held & Karp Bound                 #
# - 11/09/2020 Add Held & Karp Bound                                          #
# ipCVRP - Basic Integer Programming Method for Solving CVRP                  #
# - 12/03/2020 Add Golden77                                                   #
# - 12/04/2020 Add Two-Index Flow                                             #
# cgCVRPTW - Use column generation to find lower bound of CVRPTW              #
###############################################################################

import heapq
import math
from gurobipy import *

from .const import *
from .common import *
from .graph import *

def ipTSP(
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	edges:		"1) String (default) 'Euclidean' or \
				 2) String 'SphereEuclidean' or \
				 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
	fml:		"1) String (default) 'DFJ_Lazy' or \
				 2) String 'DFJ_PlainLoop' or \
				 3) String (not available) 'DFJ_MaxFlow' or \
				 4) String 'MTZ' or \
				 5) String 'MultiCommodityFlow' or \
				 6) String 'ShortestPath' or \
				 7) String 'QAP'" = 'DFJ_Lazy',
	timeLimit: 	"1) Double, in seconds or \
				 2) (default) None, no time limit" = None
	) -> "Exact solution for TSP":

	# Define edges ============================================================
	if (type(edges) is not dict):
		if (edges == 'Euclidean'):
			edges = getTauEuclidean(nodes)
		elif (edges == 'SphereEuclidean'):
			edges = getTauSphereEuclidean(nodes)
		else:
			print("Error: Incorrect type `edges`")
			return None

	# Solve by different formulations =========================================
	res = None
	nodeIDs = list(nodes.keys())
	if (fml == 'DFJ_Lazy'):
		res = _ipTSPLazyCuts(edges, nodeIDs, timeLimit)
	elif (fml == 'DFJ_PlainLoop'):
		res = _ipTSPPlainLoop(edges, nodeIDs, timeLimit)
	elif (fml == 'MTZ'):
		res = _ipTSPMTZ(edges, nodeIDs, timeLimit)
	elif (fml == 'ShortestPath'):
		res = _ipTSPShortestPath(edges, nodeIDs, timeLimit)
	elif (fml == 'MultiCommodityFlow'):
		res = _ipTSPMultiCommodityFlow(edges, nodeIDs, timeLimit)
	elif (fml == 'QAP'):
		res = _ipTSPQAP(edges, nodeIDs, timeLimit)
	else:
		print("Error: Incorrect or not available TSP formulation option!")
		return None
	if (res != None):
		res['fml'] = fml

	return res

def _ipTSPQAP(edges, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if i != j:
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					name = 'x_%s_%s' % (i, j))
	w = {}
	for i in range(n):
		for j in range(n):
			if (i != j):
				for k in range(n):
					w[i, j, k] = TSP.addVar(
						vtype = GRB.CONTINUOUS)

	# TSP objective function ==================================================
	TSP.setObjective(
		quicksum(
			quicksum(
				quicksum(
					edges[nodeIDs[i], nodeIDs[j]] * w[i, j, k] for k in range(n - 1)
				) for j in range(n) if j != i
			) for i in range(n)
		) + 
		quicksum(
			quicksum(
				edges[nodeIDs[i], nodeIDs[j]] * w[i, j, n - 1] for j in range(n) if j != i
			) for i in range(n)
		)
	)

	# Assignment constraints ==================================================
	for i in range(n):
		TSP.addConstr(quicksum(x[i, j] for j in range(n) if j != i) == 1)
	for j in range(n):
		TSP.addConstr(quicksum(x[i, j] for i in range(n) if j != i) == 1)

	# Linearized constraints ==================================================
	for i in range(n):
		for j in range(n):
			if (i != j):
				for k in range(n - 1):
					if (k != i and k + 1 != j):
						TSP.addConstr(w[i, j, k] >= x[i, k] + x[j, k + 1] - 1)
	for i in range(n):
		for j in range(n):
			if (i != j):
				if (i != n - 1 and j != 0):
					TSP.addConstr(w[i, j, n - 1] >= x[i, n - 1] + x[j, 0] - 1)

	# Optimize ================================================================
	if (timeLimit != None):
		TSP.setParam(GRB.Param.TimeLimit, timeLimit)
	TSP.optimize()

	# Reconstruct solution ====================================================
	ofv = None
	seq = []
	gap = None
	lb = None
	ub = None
	runtime = None
	if (TSP.status == GRB.status.OPTIMAL):
		ofv = TSP.getObjective().getValue()
		for j in range(n):
			for i in range(n):
				if (i != j and x[i, j].X > 0.9):
					seq.append(nodeIDs[i])
					break
		seq.append(seq[0])
		gap = 0
		lb = ofv
		ub = ofv
		runtime = TSP.Runtime
	elif (TSP.status == GRB.status.TIME_LIMIT):
		ofv = None
		seq = []
		gap = TSP.MIPGap
		lb = TSP.ObjBoundC
		ub = TSP.ObjVal
		runtime = TSP.Runtime

	return {
		'ofv': ofv,
		'seq': seq,
		'gap': gap,
		'lowerBound': lb,
		'upperBound': ub,
		'runtime': runtime
	}

def _ipTSPMultiCommodityFlow(edges, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if i != j:
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					obj = edges[nodeIDs[i], nodeIDs[j]], 
					name = 'x_%s_%s' % (i, j))
	y = {}
	for i in range(n):
		for j in range(n):
			if (i != j):
				for k in range(1, n):
					y[i, j, k] = TSP.addVar(
						vtype = GRB.CONTINUOUS)

	# TSP objective function ==================================================
	TSP.modelSense = GRB.MINIMIZE
	TSP.update()

	# Degree constraints =====================================================
	for i in range(n):
		TSP.addConstr(quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
		TSP.addConstr(quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

	# MCF =====================================================================
	for i in range(n):
		for j in range(n):
			if (i != j):
				for k in range(1, n):
					TSP.addConstr(y[i, j, k] <= x[i, j])

	for k in range(1, n):
		TSP.addConstr(quicksum(y[0, i, k] for i in range(1, n)) == 1)
		TSP.addConstr(quicksum(y[i, 0, k] for i in range(1, n)) == 0)
		TSP.addConstr(quicksum(y[i, k, k] for i in range(n) if i != k) == 1)
		TSP.addConstr(quicksum(y[k, j, k] for j in range(n) if j != k) == 0)
		for j in range(1, n):
			if (j != k):
				TSP.addConstr(quicksum(y[i, j, k] for i in range(n) if i != j) - quicksum(y[j, i, k] for i in range(n) if i != j) == 0)

	# TSP =====================================================================
	if (timeLimit != None):
		TSP.setParam(GRB.Param.TimeLimit, timeLimit)
	TSP.optimize()

	# Reconstruct solution ====================================================
	ofv = None
	seq = []
	arcs = []
	if (TSP.status == GRB.status.OPTIMAL):
		ofv = TSP.getObjective().getValue()
		for i, j in x:
			if (x[i, j].x > 0.5):
				arcs.append([i, j])
		currentNode = 0
		currentTime = 0
		seq.append(nodeIDs[currentNode])
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(nodeIDs[currentNode])
					arcs.pop(i)
					break
		gap = 0
		lb = ofv
		ub = ofv
		runtime = TSP.Runtime
	elif (TSP.status == GRB.status.TIME_LIMIT):
		ofv = None
		seq = []
		gap = TSP.MIPGap
		lb = TSP.ObjBoundC
		ub = TSP.ObjVal
		runtime = TSP.Runtime

	return {
		'ofv': ofv,
		'seq': seq,
		'gap': gap,
		'lowerBound': lb,
		'upperBound': ub,
		'runtime': runtime
	}

def _ipTSPShortestPath(edges, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if (j != i):
				for t in range(n):
					x[i, j, t] = TSP.addVar(
						obj = edges[nodeIDs[i], nodeIDs[j]],
						vtype = GRB.BINARY)

	# Stage constraints =======================================================
	# Start from depot 
	TSP.addConstr(quicksum(x[0, j, 0] for j in range(1, n)) == 1)
	# First stage
	for i in range(1, n):
		TSP.addConstr(quicksum(x[i, j, 1] for j in range(1, n) if i != j) - x[0, i, 0] == 0)
	# In between
	for i in range(1, n):
		for t in range(2, n - 1):
			TSP.addConstr(quicksum(x[i, j, t] for j in range(1, n) if i != j) - quicksum(x[j, i, t - 1] for j in range(1, n) if i != j) == 0)
	# Last stage
	for i in range(1, n):
		TSP.addConstr(x[i, 0, n - 1] - quicksum(x[j, i, n - 2] for j in range(1, n) if i != j) == 0)
	# Return to depot
	TSP.addConstr(quicksum(x[i, 0, n - 1] for i in range(1, n)) == 1)
	# Consequent
	for i in range(1, n):
		TSP.addConstr(quicksum(quicksum(x[i, j, t] for j in range(1, n) if i != j) for t in range(1, n - 1)) + x[i, 0, n - 1] <= 1)

	# TSP =====================================================================
	if (timeLimit != None):
		TSP.setParam(GRB.Param.TimeLimit, timeLimit)
	TSP.optimize()

	# Reconstruct solution ====================================================
	ofv = None
	seq = []
	arcs = []
	if (TSP.status == GRB.status.OPTIMAL):
		ofv = TSP.getObjective().getValue()
		for i, j, t in x:
			if (x[i, j, t].x > 0.5):
				arcs.append([i, j])
		currentNode = 0
		currentTime = 0
		seq.append(nodeIDs[currentNode])
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(nodeIDs[currentNode])
					arcs.pop(i)
					break
		gap = 0
		lb = ofv
		ub = ofv
		runtime = TSP.Runtime
	elif (TSP.status == GRB.status.TIME_LIMIT):
		ofv = None
		seq = []
		gap = TSP.MIPGap
		lb = TSP.ObjBoundC
		ub = TSP.ObjVal
		runtime = TSP.Runtime

	return {
		'ofv': ofv,
		'seq': seq,
		'gap': gap,
		'lowerBound': lb,
		'upperBound': ub,
		'runtime': runtime
	}

def _ipTSPMTZ(edges, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if i != j:
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					obj = edges[nodeIDs[i], nodeIDs[j]], 
					name = 'x_%s_%s' % (i, j))
	u = {}
	for i in range(n):
		u[i] = TSP.addVar(
			vtype = GRB.CONTINUOUS,
			name = 'u_%s' % (i))

	# TSP objective function ==================================================
	TSP.modelSense = GRB.MINIMIZE
	TSP.update()

	# Degree constraints ======================================================
	for i in range(n):
		TSP.addConstr(quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
		TSP.addConstr(quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

	# Sequence constraints ====================================================
	for i in range(1, n):
		for j in range(1, n):
			if (i != j):
				TSP.addConstr(u[i] - u[j] + (n - 1) * x[i, j] <= n - 2, name = 'seq_%s_%s' % (i, j))
	for i in range(1, n):
		TSP.addConstr(1 <= u[i])
		TSP.addConstr(u[i] <= n - 1)

	# TSP =====================================================================
	if (timeLimit != None):
		TSP.setParam(GRB.Param.TimeLimit, timeLimit)
	TSP.optimize()

	# Reconstruct solution ====================================================
	ofv = None
	gap = None
	seq = []
	arcs = []
	if (TSP.status == GRB.status.OPTIMAL):
		ofv = TSP.getObjective().getValue()
		gap = TSP.Params.MIPGapAbs
		for i, j in x:
			if (x[i, j].x > 0.5):
				arcs.append([i, j])
		currentNode = 0
		currentTime = 0
		seq.append(nodeIDs[currentNode])
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(nodeIDs[currentNode])
					arcs.pop(i)
					break
		gap = 0
		lb = ofv
		ub = ofv
		runtime = TSP.Runtime
	elif (TSP.status == GRB.status.TIME_LIMIT):
		ofv = None
		seq = []
		gap = TSP.MIPGap
		lb = TSP.ObjBoundC
		ub = TSP.ObjVal
		runtime = TSP.Runtime

	return {
		'ofv': ofv,
		'seq': seq,
		'gap': gap,
		'lowerBound': lb,
		'upperBound': ub,
		'runtime': runtime
	}

def _ipTSPPlainLoop(edges, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if (i != j):
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					obj = edges[nodeIDs[i], nodeIDs[j]], 
					name = 'x_%s_%s' % (i, j))
				
	# TSP =====================================================================
	TSP.modelSense = GRB.MINIMIZE
	TSP.Params.lazyConstraints = 1
	TSP.update()

	# Degree constraints ======================================================
	for i in range(n):
		TSP.addConstr(quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
		TSP.addConstr(quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

	# Resolve to optimality and try to find sub-tours =========================
	noSubtourFlag = False
	accRuntime = 0
	while (not noSubtourFlag):
		if (timeLimit != None):
			TSP.setParam(GRB.Param.TimeLimit, timeLimit - accRuntime) # FIXME
		TSP.optimize()
		if (TSP.status == GRB.status.OPTIMAL):
			accRuntime += TSP.Runtime
			arcs = tuplelist((i, j) for i, j in x.keys() if x[i, j].X > 0.9)
			components = findComponentsUndirected(arcs)
			if (len(components) == 1):
				noSubtourFlag = True
				break
			else:
				for comp in components:
					TSP.addConstr(quicksum(x[i, j] for i in comp for j in comp if i != j) <= len(comp) - 1)
		elif (TSP.status == GRB.status.TIME_LIMIT):
			accRuntime += TSP.Runtime
			break

	# Reconstruct solution ====================================================
	ofv = None
	seq = []
	arcs = []
	if (TSP.status == GRB.status.OPTIMAL):
		ofv = TSP.getObjective().getValue()
		for i, j in x:
			if (x[i, j].x > 0.5):
				arcs.append([i, j])
		currentNode = 0
		currentTime = 0
		seq.append(nodeIDs[currentNode])
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(nodeIDs[currentNode])
					arcs.pop(i)
					break
		gap = 0
		lb = ofv
		ub = ofv
		runtime = accRuntime
	elif (TSP.status == GRB.status.TIME_LIMIT):
		ofv = None
		seq = []
		gap = TSP.MIPGap
		lb = TSP.ObjBoundC
		ub = TSP.ObjVal
		runtime = accRuntime

	return {
		'ofv': ofv,
		'seq': seq,
		'gap': gap,
		'lowerBound': lb,
		'upperBound': ub,
		'runtime': runtime
	}

def _ipTSPLazyCuts(edges, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if (i != j):
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					obj = edges[nodeIDs[i], nodeIDs[j]], 
					name = 'x_%s_%s' % (i, j))
				
	# TSP objective function ==================================================
	TSP.modelSense = GRB.MINIMIZE
	TSP.Params.lazyConstraints = 1
	TSP.update()

	# Degree constraints ======================================================
	for i in range(n):
		TSP.addConstr(quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
		TSP.addConstr(quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

	# Sub-tour elimination ====================================================
	TSP._x = x
	def subtourelim(model, where):
		if (where == GRB.Callback.MIPSOL):
			x_sol = model.cbGetSolution(model._x)
			arcs = tuplelist((i, j) for i, j in model._x.keys() if x_sol[i, j] > 0.9)
			components = findComponentsUndirected(arcs)
			for component in components:
				if (len(component) < n):
					model.cbLazy(quicksum(x[i,j] for i in component for j in component if i != j) <= len(component) - 1)

	# TSP with callback =======================================================
	if (timeLimit != None):
		TSP.setParam(GRB.Param.TimeLimit, timeLimit)
	TSP.optimize(subtourelim)

	# Reconstruct solution ====================================================
	ofv = None
	seq = []
	arcs = []
	if (TSP.status == GRB.status.OPTIMAL):
		ofv = TSP.getObjective().getValue()
		for i, j in x:
			if (x[i, j].x > 0.5):
				arcs.append([i, j])
		currentNode = 0
		currentTime = 0
		seq.append(nodeIDs[currentNode])
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(nodeIDs[currentNode])
					arcs.pop(i)
					break
		gap = 0
		lb = ofv
		ub = ofv
		runtime = TSP.Runtime
	elif (TSP.status == GRB.status.TIME_LIMIT):
		ofv = None
		seq = []
		gap = TSP.MIPGap
		lb = TSP.ObjBoundC
		ub = TSP.ObjVal
		runtime = TSP.Runtime

	return {
		'ofv': ofv,
		'seq': seq,
		'gap': gap,
		'lowerBound': lb,
		'upperBound': ub,
		'runtime': runtime
	}


def lrTSP(
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	edges:		"1) String (default) 'Euclidean' or \
				 2) String 'SphereEuclidean' or \
				 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
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
			for i in nodes:
				nodeIDs.append(i)

	# Define edges ==============================================================
	if (type(edges) is not dict):
		lstNodeID = nodeIDs.copy()
		if (edges == 'Euclidean'):
			edges = getTauEuclidean(nodes, lstNodeID)
		elif (edges == 'SphereEuclidean'):
			edges = getTauSphereEuclidean(nodes, lstNodeID)
		else:
			print("Error: Incorrect type `edges`")
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
					weightArcs.append((i, j, edges[i, j] - u[i] - u[j]))

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
	vehNum:			"Number of vehicles" = None,
	cutoffTimeSub:	"Cutoff time (s) for subproblem, if no incumbent sol found within cutoff time, stop adding variables" = 10,
	cutoffTimeMas:	"Accumulated runtime (s) for entire problem" = 300,
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
		sub.setParam(GRB.Param.TimeLimit, cutoffTimeSub)
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
		elif (sub.status == GRB.status.TIME_LIMIT):
			try:
				print("Current best solution")
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
			except:
				print("No incumbent for subproblem")
				ofv = None
				c = None
				a = None
				route = None
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

	# Initial Optimize
	CVRPTW.setParam("OutputFlag", 0)
	CVRPTW.modelSense = GRB.MINIMIZE
	CVRPTW.optimize()

	# Now solve the Master (Set Partition Formulation) problem ================
	accTime = 0
	canAddVarFlag = True
	while(canAddVarFlag):
		startIter = datetime.datetime.now()
		canAddVarFlag = False
		
		if (CVRPTW.status == GRB.status.OPTIMAL):
			# Solve subproblem
			pi = {}
			for constraint in cons:
				pi[constraint] = cons[constraint].Pi
			subproblem = pricing(pi)
			print(subproblem['route'])
			if (subproblem['ofv'] != None):	
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
				canAddVarFlag = False
		else:
			canAddVarFlag = False

		# Re-optimize
		CVRPTW.optimize()
		oneIter = (datetime.datetime.now() - startIter).total_seconds()
		accTime += oneIter
		if (accTime > cutoffTimeMas):
			canAddVarFlag = False

	# Interpret solution for lower bound ======================================
	ColGen = CVRPTW.getObjective().getValue()

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
	EarlyBraching = CVRPTW.getObjective().getValue()

	return {
		'ColGen': ColGen,
		'EarlyBraching': EarlyBraching,
		'routes': solRoute
	}

