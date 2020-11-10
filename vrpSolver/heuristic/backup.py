from .common import *

# [Constructing]
def conBestNearestNeighborhood_TSPPD(
	dicNodeLoc: "Dictionary, {nodeID1: (lat, lon), nodeID2: (lat, lon), ...}",
	tau = 'SphereEuclidean': "1) String 'Euclidean' or 'SphereEuclidean', or 2) Dictionary {(nodeID1, nodeID2): dist, ...}",
	speed: "Speed of the vehicle",
	depotID: "Depot location, None as no depot",
	reqs: "List of paired requests, {reqID1: {'pickupID': nodeID, 'deliveryID': nodeID}, ...}",
	nodeTW: "{nodeID: {'earliestOpen': time, 'latestOpen': time}, ...}",
	reqTW: "{reqID1: {'earliestPickup': time, 'latestDeliver': time}, ...}",
	reqCap: "{reqID1: load, reqID2: load, ...}",
	capacity = None: "Double",
	designatedFirst = None: "None if no designated first pickup, otherwise provide a reqID as start"	
	) -> "Construct route using depth-first search in the Nearest Neighborhood search tree":
	
	'''
	Output
	ofv: Objective function
	actions: List of dictionary, [{'action': 'leave', 'nodeID': depotID, 'accOfv': 0, 'curLoad': 0}, ...]
	'''

	# For output
	ofv = 0
	actions = []

	# States for search
	stateStack = [] # Stack for DFS, for now, use List, if its too big, use collections.deque / queue.LifoQueue
	stateList = []  # List storing final list of states

	# For no depot case
	if (depotID == None):
		dicNodeLoc[None] = None

	# Initial tau
	lstNodeID = []
	if (type(tau) not dict):
		for req in reqs:
			if (reqs[req]['pickupID'] not in lstNodeID):
				lstNodeID.append(reqs[req]['pickupID'])
			if (reqs[req]['deliveryID'] not in lstNodeID):
				lstNodeID.append(reqs[req]['deliveryID'])
		if (tau == 'Euclidean'):
			tau = getTauEuclidean(dicNodeLoc, lstNodeID)
		elif (tau == 'SphereEuclidean'):
			tau = getTauSphereEuclidean(dicNodeLoc, lstNodeID)
		else:
			print("Incorrect type: tau")
			return None

	# Check initial load feasibility - no load should exceed capacity
	if (capacity != None):
		for reqID in reqCap:
			if (reqCap[reqID] > capacity):
				return {
					'ofv': None,
					'route': None
				}

	# Check initial time windows - no delivery can exceed time window
	if (reqTW != None):
		for reqID in reqTW:
			isViolateLatestDeliveryFlag = isViolateLatestDelivery(
				curTime = 0,
				cur2Pickup = tau[depotID, dicNodeLoc[reqs[reqID]['pickupID']]],
				cur2Delivery = tau[depotID, dicNodeLoc[reqs[reqID]['deliveryID']]],
				pickup2Delivery = tau[dicNodeLoc[reqs[reqID]['pickupID']], dicNodeLoc[reqs[reqID]['deliveryID']]],
				latestTime = reqTW[reqID],
				speed = speed,
				state = 'pickup')
			if (not isViolateLatestDeliveryFlag):
				return {
					'ofv': None,
					'route': None
				}

	# First load, for initial state
	reqNotYet = []
	reqOnboard = []
	reqDone = []
	for req in reqs:
		reqNotYet.append(req)
	if (designatedFirst == None or (designatedFirst != None and designatedFirst not in reqs)):
		dist2First = None
		for req in reqs:
			if (distFirst == None or tau[depotID, reqs[req]['pickupID']] < dist2First):
				dist2First = tau[depotID, reqs[req]['pickupID']]
				designatedFirst = req
	reqNotYet.remove(designatedFirst)
	reqOnboard.append(designatedFirst)

	# Initial Search tree
	states = {
		'reqNotYet': reqNotYet,
		'reqOnboard': reqOnboard,
		'reqDone': reqDone,
		'curNodeID': reqs[designatedFirst]['pickupID'],
		'curLoad': reqs[designatedFirst]['load'],
		'accDist': 0,
		'accTime': 0,
		'seq': 0
	}
	stateStack.append(states.copy())

	# Loop until every request has fulfilled pickup and delivery
	while (len(stateList) < 2 * len(reqs)):
		# Pop stack, as current node, if the level is smaller than the last element in stateList, store it
		curState = stateStack.pop()
		if (stateList == []):
			stateList.append(curState.copy())
		elif (stateList[-1]['seq'] < curState['seq']):
			stateList.append(curState.copy())
		elif (stateList[-1]['seq'] == curState['seq']):
			stateList.pop()
			stateList.append(curState.copy())

		# Find feasible children - First find all possible next steps
		nextSteps = []
		for reqID in curState['reqNotYet']:
			nextSteps.append(('pickup', reqID))
		for reqID in curState['reqOnboard']:
			nextSteps.append(('deliver', reqID))

		# Find feasible children - Check capacity, for pickup
		if (capacity != None):
			for cap in reqCap:
				if (('pickup', cap) in nextSteps and curState['curLoad'] + reqCap[cap] > capacity):
					nextSteps.remove(('pickup', cap))

		# Find feasible children - Check time windows, for a child, if it is feasible, 
		#                          it should have at least one child satisfy 
		infNextSteps = []
		isViolateLatestDeliveryFlag = False
		for step in nextSteps:
			# Imagine arrive at the next step, now see it is going to let at least one cannot be delivered
			nextReq = step[1]
			nextAct = step[0]
			newTime = None
			if (nextAct == 'pickup'):
				newTime = curTime + tau[dicNodeID[curState['curNodeID']], dicNodeLoc[reqs[nextReq]['pickupID']]] / speed
				newLocID = reqs[nextReq]['pickupID']
			else:
				newTime = curTime + tau[dicNodeID[curState['curNodeID']], dicNodeLoc[reqs[nextReq]['deliveryID']]] / speed
				newLocID = reqs[nextReq]['deliveryID']

			newNotYet = curState['reqNotYet'].copy()
			newOnboard = curState['reqOnboard'].copy()
			if (nextAct == 'pickup'):
				newOnboard.append(nextReq)
				newNotYet.remove(nextReq)
			else:
				newOnboard.remove(nextReq)

			# If I go there, is it going to be late for any item to be picked up?
			for toPickup in newNotYet:
				isViolateLatestDeliveryFlag = isViolateLatestDelivery(
					curTime = newTime,
					cur2Pickup = tau[dicNodeLoc[newLocID], dicNodeLoc[reqs[toPickup]['pickupID']]],
					cur2Delivery = tau[dicNodeLoc[newLocID], dicNodeLoc[reqs[toPickup]['deliveryID']]],
					pickup2Delivery = tau[dicNodeLoc[reqs[toPickup]['pickupID']], dicNodeLoc[reqs[toPickup]['deliveryID']]],
					latestTime = reqTW[toPickup],
					speed = speed,
					state = 'pickup')
				if (isViolateLatestDeliveryFlag):
					break

			# If I go there, is it going to be late for any item to be delivered?
			if (not isViolateLatestDeliveryFlag):
				for toDeliver in newOnboard:
					isViolateLatestDeliveryFlag = isViolateLatestDelivery(
						curTime = newTime,
						cur2Pickup = tau[dicNodeLoc[newLocID], dicNodeLoc[reqs[toPickup]['pickupID']]],
						cur2Delivery = tau[dicNodeLoc[newLocID], dicNodeLoc[reqs[toPickup]['deliveryID']]],
						pickup2Delivery = tau[dicNodeLoc[reqs[toPickup]['pickupID']], dicNodeLoc[reqs[toPickup]['deliveryID']]],
						latestTime = reqTW[toPickup],
						speed = speed,
						state = 'delivery')
					if (isViolateLatestDeliveryFlag):
						break

			if (isViolateLatestDeliveryFlag):
				infNextSteps.append((step[0], step[1]))

		for step in infNextSteps:
			nextSteps.remove(step)

		# Calculate all feasible children states, order them by how close to the curLoc, push them to stack
		for step in nextSteps:
			nextState = {
				
			}
	
	return {
		'dist': dist,
		'stamps': stamps
	}

def isViolateLatestDelivery(
	curTime: "Current time",
	cur2Pickup: "Distance from current location to pickup", 
	cur2Delivery: "Distance from current location to delivery location",
	pickup2Delivery: "Distance from pickup to delivery",
	latestTime: "Latest delivery time",
	speed: "Speed",
	state = 'pickup': "String 'pickup' not picked up yet, or 'delivery' to be delivered"
	) -> "False if will not violate latest delivery time constraints, True otherwise":

	violateFlag = None

	arrTime = 0
	if (state == 'pickup'):
		arrTime = curTime + (cur2Pickup + pickup2Delivery) / speed
	elif (state == 'delivery'):
		arrTime = curTime + (cur2Delivery) / speed
	else:
		print("Incorrect state")
		return None

	if (arrTime > latestTime):
		violateFlag = False
	else:
		violateFlag = True

	return violateFlag
