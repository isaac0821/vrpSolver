###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 11/09/2020                             #
# ipTSP - Basic Integer Programming Method for Solving TSP                    #
# - 11/01/2020 Add DFJ_Lazy                                                   #
# - 11/01/2020 Add DFJ_PlainLoop                                              #
# - 11/01/2020 Add MTZ                                                        #
# - 11/01/2020 Add MultiCommodityFlow                                         #
# - 11/01/2020 Add ShortestPath                                               #
# - 11/01/2020 Add QAP                                                        #
# lrTSP - Use Lagrangian Relaxation to Give Held & Karp Bound                 #
###############################################################################

import heapq
import math
from gurobipy import *

from vrpSolver.const import *
from vrpSolver.common import *
from vrpSolver.graph.basic import *
from vrpSolver.graph.mst import *

def ipTSP(
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

	# Solve by different formulations =========================================
	res = None
	if (fml == 'DFJ_Lazy'):
		res = _ipTSPLazyCuts(tau, nodeIDs, timeLimit)
	elif (fml == 'DFJ_PlainLoop'):
		res = _ipTSPPlainLoop(tau, nodeIDs, timeLimit)
	elif (fml == 'MTZ'):
		res = _ipTSPMTZ(tau, nodeIDs, timeLimit)
	elif (fml == 'ShortestPath'):
		res = _ipTSPShortestPath(tau, nodeIDs, timeLimit)
	elif (fml == 'MultiCommodityFlow'):
		res = _ipTSPMultiCommodityFlow(tau, nodeIDs, timeLimit)
	elif (fml == 'QAP'):
		res = _ipTSPQAP(tau, nodeIDs, timeLimit)
	else:
		print("Error: Incorrect or not available TSP formulation option!")
		return None
	if (res != None):
		res['fml'] = fml

	return res

def _ipTSPQAP(tau, nodeIDs, timeLimit):
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
	TSP.setObjective(quicksum(quicksum(quicksum(tau[nodeIDs[i], nodeIDs[j]] * w[i, j, k] for k in range(n - 1)) for j in range(n) if j != i) for i in range(n)) + quicksum(quicksum(tau[nodeIDs[i], nodeIDs[j]] * w[i, j, n - 1] for j in range(n) if j != i) for i in range(n)))

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
					seq.append(i + 1)
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

def _ipTSPMultiCommodityFlow(tau, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if i != j:
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					obj = tau[nodeIDs[i], nodeIDs[j]], 
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

	# TSP with callback =======================================================
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
		seq.append(currentNode + 1)
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(currentNode + 1)
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

def _ipTSPShortestPath(tau, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if (j != i):
				for t in range(n):
					x[i, j, t] = TSP.addVar(
						obj = tau[nodeIDs[i], nodeIDs[j]],
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
		seq.append(currentNode + 1)
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(currentNode + 1)
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

def _ipTSPMTZ(tau, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if i != j:
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					obj = tau[nodeIDs[i], nodeIDs[j]], 
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
		seq.append(currentNode + 1)
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(currentNode + 1)
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

def _ipTSPPlainLoop(tau, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if (i != j):
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					obj = tau[nodeIDs[i], nodeIDs[j]], 
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
			components = findComponentsUndirected(arcs, n)
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
		seq.append(currentNode + 1)
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(currentNode + 1)
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

def _ipTSPLazyCuts(tau, nodeIDs, timeLimit):
	n = len(nodeIDs)
	TSP = Model('TSP')

	# Decision variables ======================================================
	x = {}
	for i in range(n):
		for j in range(n):
			if (i != j):
				x[i, j] = TSP.addVar(
					vtype = GRB.BINARY, 
					obj = tau[nodeIDs[i], nodeIDs[j]], 
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
		if where == GRB.Callback.MIPSOL:
			x_sol = model.cbGetSolution(model._x)
			arcs = tuplelist((i, j) for i, j in model._x.keys() if x_sol[i, j] > 0.9)
			components = findComponentsUndirected(arcs, n)
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
		seq.append(currentNode + 1)
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(currentNode + 1)
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
		'costSum': costSum,
		'm1t': m1t,
		'oldL': oldL,
		'L': L,
		'd': d,
		'u': u
	}