from gurobipy import *

from vrpSolver.common import *

def ipTSPPDMTZ(
	nodeLoc : "Dictionary, returns the coordinate of given nodeID, \
	            {\
	                nodeID1: (lat, lon), \
	                nodeID2: (lat, lon), \
	                ... \
	            }" = None, 
	tau     : "1) String 'Euclidean' or \
	           2) String (default) 'SphereEuclidean' or \
	           3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "SphereEuclidean",
	speed   : "Speed of the vehicle" = 50,
	now     : "Current time, design for offsetting 'latestDeliverTime' in reqs" = 0,
	reqs    : "List of paired requests, \
	            {\
	                reqID1: {\
	                    'pickupID': nodeID, \
	                    'deliveryID': nodeID \
	                    }, \
	                ... \
	            }" = None,
	depotID : "Depot location, None as no depot, for Open TSP" = None,
	reqIDs  : "1) String (default) 'All', or \
	           2) A list of request IDs" = 'All',
	vehCap  : "Capacity of vehicle" = None
	) -> "Using callback function to calculate TSPPD, possibly with vehicle capacity, node time window, request deadline":


    # General =================================================================
	TSPPD = Model('TSPPD')

	# reqIDs ==================================================================
	if (reqIDs == 'All'):
		reqIDs = []
		for i in reqs:
			reqIDs.append(i)

	# Sets ====================================================================
	# We need to identify pickup location and delivery locations to distinguish 
	# them by index. In the following sets:
	# 0 - Depot
	# 1 to len(reqIDs) - Pickup locs
	# len(reqIDs) + 1 to 2 * len(reqIDs) - Delivery locs
	# 2 * len(reqIDs) + 1 (or w) - Depot
	w = 2 * len(reqIDs) + 1
	P = [i for i in range(1, len(reqIDs) + 1)]
	D = [i for i in range(len(reqIDs) + 1, w)]
	OUP = list(set([0]).union(P))
	OUPUD = list(set(OUP).union(D))
	PUD = list(set(P).union(D))
	PUDUw = list(set(PUD).union([w]))
	DUw = list(set(D).union([w]))
	OUPUDUw = list(set(OUP).union(DUw))
	# Given a location index, find where it could come from
	def beforeNodes(i):
		delta = []
		if (i in P):
			delta = [i for i in OUPUD]
			delta.remove(i)
			delta.remove(i + len(reqIDs))
		elif (i in D):
			delta = PUD.copy()
			delta.remove(i)
		elif (node == 0):
			delta = []
		elif (node == 2 * len(P) + 1):
			delta = [i for i in D]
		return delta
	# Given a location index, find where it can go afterwards
	def afterNodes(i):
		delta = []
		if (i in P):
			delta = [i for i in PUD]
			delta.remove(i)
		elif (i in D):
			delta = [i for i in PUDUw]
			delta.remove(i)
			delta.remove(i - len(reqIDs))
		elif (i == 0):
			delta = P.copy()
		elif (i == w):
			delta = []
		return delta
	# Given a location index, and the reqID list, return nodeID for querying tau
	def nodeIDByIndex(i):
		nodeID = None
		if (i == 0 or i == w):
			nodeID = depotID
		elif (i <= len(reqIDs)):
			nodeID = reqs[reqIDs[i - 1]]['pickupID']
		elif (i > len(reqIDs)):
			nodeID = reqs[reqIDs[i - len(reqIDs) - 1]]['deliveryID']
		return nodeID

	# tau =====================================================================
	if (depotID == None):
		nodeLoc[None] = None
	if (type(tau) is not dict):
		lstNodeID = []
		lstNodeID.append(depotID)
		for req in reqIDs:
			p = reqs[req]['pickupID']
			d = reqs[req]['deliveryID']
			if (p not in lstNodeID):
				lstNodeID.append(p)
			if (d not in lstNodeID):
				lstNodeID.append(d)
		if (tau == 'Euclidean'):
			tau = getTauEuclidean(nodeLoc, lstNodeID)
		elif (tau == 'SphereEuclidean'):
			tau = getTauSphereEuclidean(nodeLoc, lstNodeID)
		else:
			print("Incorrect type: tau")
			return None

	# Decision variables ======================================================
	# x[i, j] defined as whether the route travel from i to j, i \in OUPUD, j \in afterNodes(i)
	x = {}
	for i in OUPUD:
		afterI = afterNodes(i)
		for j in afterI:
			x[i, j] = TSPPD.addVar(vtype = GRB.BINARY)
	t = {}
	for i in OUPUDUw:
		t[i] = TSPPD.addVar(vtype = GRB.CONTINUOUS)
	q = {}
	for i in OUPUDUw:
		q[i] = TSPPD.addVar(vtype = GRB.CONTINUOUS)

	# get big M ===============================================================
	M = 0
	for i in OUPUD:
		M += tau[nodeIDByIndex(i), nodeIDByIndex(i + 1)]

	# Objective function ======================================================
	TSPPD.modelSense = GRB.MINIMIZE
	TSPPD.setObjective(t[w])
	TSPPD.update()

	# Degree constraints ===================================================
	# First arc of a route must be from depot to a pickup location
	TSPPD.addConstr(quicksum(x[0, j] for j in P) == 1)
	# Last arc of a route must be from a delivery location to depot
	TSPPD.addConstr(quicksum(x[i, w] for i in D) == 1)
	# Each pickups and delivery node must be visited only once, flow balance
	for j in PUD:
		beforeJ = beforeNodes(j)
		TSPPD.addConstr(quicksum(x[i, j] for i in beforeJ) == 1)
		afterJ = afterNodes(j)
		TSPPD.addConstr(quicksum(x[j, k] for k in afterJ) == 1)

	# MTZ constraints =========================================================
	TSPPD.addConstr(t[0] == 0)
	for i in OUPUD:
		afterI = afterNodes(i)
		for j in afterI:
			TSPPD.addConstr(t[j] >= t[i] + tau[nodeIDByIndex(i), nodeIDByIndex(j)] * x[i, j] - M * (1 - x[i, j]))
			TSPPD.addConstr(t[j] <= t[i] + tau[nodeIDByIndex(i), nodeIDByIndex(j)] * x[i, j] + M * (1 - x[i, j]))
	for i in P:
		TSPPD.addConstr(t[i] + tau[nodeIDByIndex(i), nodeIDByIndex(i + len(reqIDs))] <= t[i + len(reqIDs)])

	# Capacity constraints ====================================================
	TSPPD.addConstr(q[0] == 0)
	for j in P:
		beforeJ = beforeNodes(j)
		for i in beforeJ:
			TSPPD.addConstr(q[j] >= q[i] + x[i, j] - vehCap * (1 - x[i, j]))
			TSPPD.addConstr(q[j] <= q[i] + x[i, j] + vehCap * (1 - x[i, j]))
	for j in D:
		beforeJ = beforeNodes(j)
		for i in beforeJ:
			TSPPD.addConstr(q[j] >= q[i] - x[i, j] - vehCap * (1 - x[i, j]))
			TSPPD.addConstr(q[j] <= q[i] - x[i, j] + vehCap * (1 - x[i, j]))

	# SLA constraints =========================================================
	for i in D:
		if ('latestDeliverTime' in reqs[reqIDs[i - len(reqIDs)]]):
			TSPPD.addConstr(speed * t[i] <= reqs[reqIDs[i - len(reqIDs)]]['latestDeliverTime'])

	# TSPPD ===================================================================
	TSPPD.optimize()

	# Reconstruct solution ====================================================
	route = []
	if (TSPPD.status == GRB.status.OPTIMAL):
		segs = []
		for i in OUPUD:
			for j in PUDUw:
				if ((i, j) in x):
					if (x[i, j].X >= 0.8):
						segs.append([i, j])
		currentNode = 0
		while (currentNode <= w):
			if (currentNode == 0):
				route.append((
					'depart',
					None,
					nodeLoc[depotID],
					0,
					0
					))
			elif (currentNode in P):
				route.append((
					'pickup', 
					currentNode,
					nodeIDByIndex(currentNode), 
					round(t[currentNode].X, 3),
					abs(round(q[currentNode].X, 0))
					))
			elif (currentNode in D):
				route.append((
					'deliver', 
					currentNode - len(reqIDs),
					nodeIDByIndex(currentNode), 
					round(t[currentNode].X, 3),
					abs(round(q[currentNode].X, 0))
					))
			elif (currentNode == w):
				route.append((
					'return',
					None,
					nodeLoc[depotID],
					round(t[w].X, 3),
					abs(round(q[w].x, 0))
					))
			if (currentNode != w):
				currentNode = segs[currentNode][1]
			else:
				break
	return route
