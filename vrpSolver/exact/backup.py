
def ipTSPPDLazyCuts(
	nodeLoc : "Dictionary, returns the coordinate of given nodeID, \
	            {\
	                nodeID1: (lat, lon), \
	                nodeID2: (lat, lon), \
	                ... \
	            }" = None, 
	tau     : "1) String 'Euclidean' or \
	           2) String (default) 'SphereEuclidean' or \
	           3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "SphereEuclidean",
	speed   : "Speed of the vehicle" = None,
	now     : "Current time, design for offsetting 'latestDeliverTime' in reqs" = 0,
	reqs    : "List of paired requests, \
	            {\
	                reqID1: {\
	                    'pickupID': nodeID, \
	                    'deliveryID': nodeID , \
	                    'latestDeliverTime': time, \
	                    'load': load,\
	                    }, \
	                ... \
	            }" = None,
	depotID : "Depot location, None as no depot" = None,
	reqIDs  : "1) String (default) 'All', or \
	           2) A list of request IDs" = 'All',
	vehCap  : "Capacity of vehicle" = None,
	fml     : "1) String 'RR', adapted from Ruland and Rodin (1997), Asymmetric TSPPD" = 'RR'
	) -> "Using callback function to calculate TSPPD, possibly with vehicle \
          capacity, node time window, request deadline":

	# General =================================================================
	TSPPD = Model('TSPPD')

	# reqIDs ==================================================================
	if (reqIDs == 'All'):
		reqIDs = []
		for i in reqs:
			reqIDs.append(i)

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
			delta = OUPUD.copy()
			delta.remove(i)
			delta.remove(i + len(reqIDs))
		elif (i in D):
			delta = PUD.copy()
			delta.remove(i)
		elif (node == 0):
			delta = []
		elif (node == 2 * len(P) + 1):
			delta = D.copy()
		return delta
	# Given a location index, find where it can go afterwards
	def afterNodes(i):
		delta = []
		if (i in P):
			delta = PUD.copy()
			delta.remove(i)
		elif (i in D):
			delta = PUDUw.copy()
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

	# Decision variables ======================================================
	# x[i, j] defined as whether the route travel from i to j, i \in OUPUD, j \in afterNodes(i)
	x = {}
	for i in OUPUD:
		DeltaIP = afterNodes(i)
		for j in DeltaIP:
			x[i, j] = TSPPD.addVar(
				vtype = GRB.BINARY, 
				name = 'x_%s_%s' % (i, j))

	# Objective function ======================================================
	TSPPD.modelSense = GRB.MINIMIZE
	TSPPD.Params.lazyConstraints = 1
	TSPPD.update()

	# (Initial) Constraints ===================================================
	# First arc of a route must be from depot to a pickup location
	# \sum_{j \in P} x_{0, j} = 1
	TSPPD.addConstr(quicksum(x[0, j] for j in P) == 1)
	# Last arc of a route must be from a delivery location to depot
	# \sum_{i \in D} x_{i, w} = 1
	TSPPD.addConstr(quicksum(x[i, w] for i in D) == 1)
	# Each pickups and delivery node must be visited only once
	# \sum_{i \in beforeNodes(j)} x_{i, j} = 1 \forall j \in PUD
	for j in PUD:
		beforeJ = beforeNodes(j)
		TSPPD.addConstr(quicksum(x[i, j] for i in beforeJ) == 1)
	# Flow balance, has to go somewhere after visiting node
	# \sum_{k \in afterNodes(j)} x_{j, k} = 1 \forall j \in PUD
	for j in PUD:
		afterJ = afterNodes(j)
		TSPPD.addConstr(quicksum(x[j, k] for k in afterJ) == 1)

	# Lazy Constraints ========================================================
	# x(\delta(S)) \ge 1, \forall S \subset V
	# x(\delta(S)) \ge 4, \forall S \subset V, {(0, p)} \subset S, {(d, w)} \subset V \setminus S
	TSPPD._x = x
	def callback(model, where):
		if where == GRB.callback.MIPSOL:
			x_sol = model.cbGetSolution(model._x)

			subtourFlag = True
			loadFlag = False

			# Check sub-tours
			arcs = tuplelist((i, j) for i, j in model._x.keys() if x_sol[i, j] > 0.9)
			components = findComponentsUndirected(arcs, n)

			# Different cutset will lead to different sub-tour constraints
			if (len(components) == 1):
				# If all connected, no sub-tour
				subtourFlag = False
			else:
				pass

			# Check load capacity constraint
			if (not subtourFlag):
				pass
	return
