import heapq
import math

from .const import *
from .common import *
from .graph import *
from .geometry import *
from .msg import *

def heuVRP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y), 'demand': 1, ...}, \
                        nodeID2: {'loc': (x, y), 'demand': 1, ...}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...} or \
                 4) String 'Grid', will need to add arguments using `edgeArgs`"= "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                    {\
                        'colRow': (numCol, numRow),\
                        'barriers': [(coordX, coordY), ...], \
                    }" = None,
    depotID:    "DepotID, default to be 0" = 0,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime:"Service time spent on each customer (will be added into travel matrix)" = 0,
    objective:  "Objective function\
                 1) String, 'Makespan', or\
                 2) String, 'Cost'" = 'Makespan',
    vehicle:    "Dictionary, describing the vehicle situation, \
                 { \
                    'numTruck': number of truck, if not provided, default to be infinite\
                    'capTruck': capacity of truck, if not provided, default to be infinite\
                 }" = None,
    consAlgo:   "1) String 'CWSaving' or \
                 2) String (not available) 'Sweep' or \
                 3) String (not available) " = 'CWSaving',
    impAlgo:    "1) String '3Opt'" = '3Opt'
    ) -> "Use given heuristic methods to get basic Capacitated VRP solution":

    # Define nodeIDs ==========================================================
    customerID = []
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)
        else:
            msgError(ERROR_INCOR_NODEIDS)
            return
    customerID = [i for i in nodeIDs if i != depotID]

    # Initialize demands ======================================================
    demands = {}
    noDemandList = []
    for n in customerID:        
        if ('demand' not in nodes[n]):
            noDemandList.append(n)
            demands[n] = 1
        else:
            demands[n] = nodes[n]['demand']
    if (len(noDemandList) > 0):
        msgWarning("MESSAGE: 'demands' is missing in nodes %s, set to default value (as 1)" % list2String(noDemandList))

    # Define tau ==============================================================
    tau = {}
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            tau = getTauEuclidean(nodes, nodeIDs)
        elif (edges == 'LatLon'):
            tau = getTauLatLon(nodes, nodeIDs)
        elif (edges == 'Grid'):
            tau = getTauGrid(nodes, nodeIDs, edgeArgs['colRow'], edgeArgs['barriers'])
        else:
            msgError(ERROR_INCOR_TAU)
            return None
    else:
        tau = dict(edges)

    # Service time ============================================================
    if (serviceTime != None and serviceTime > 0):
        for (i, j) in tau:
            if (i != depotID and j != depotID and i != j):
                tau[i, j] += serviceTime
            elif (i == depotID or j == depotID and i != j):
                tau[i, j] += serviceTime / 2 
    else:
        serviceTime = 0

    # Subroutines =============================================================
    def _consVRPClarkeWright(nodes, depotID, customerID, tau, vehCap, vehNum):
        # Initial routes
        routes = {}
        for i in range(len(customerID)):
            routes[i] = {
                'route': [depotID, customerID[i], depotID],
                'demand': demands[customerID[i]],
                'length': tau[depotID, customerID[i]] + tau[customerID[i], depotID]
            }

        # vehCap
        if (vehCap == None):
            vehCap = len(nodeIDs) + 1

        # vehNum
        if (vehNum == None):
            vehNum = len(nodeIDs) + 1

        # Initial saving raking
        rankSaving = []
        for i in customerID:
            for j in customerID:
                if (i != j):
                    # Calculate saving for each pair
                    sav = tau[depotID, i] + tau[depotID, j] - tau[i, j]
                    # heapq returns the smallest, so add a negative sign
                    heapq.heappush(rankSaving, (-sav, (i, j)))

        # Merge routes subroutine
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
            newLength = routes[rI]['length'] + routes[rJ]['length'] + tau[i, j] - tau[depotID, i] - tau[depotID, j]
            routes.pop(rI)
            routes.pop(rJ)
            newRouteIndex = max(list(routes.keys())) + 1
            routes[newRouteIndex] = {
                'route': newRoute,
                'demand': newDemand,
                'length': newLength
            }

        # Merge routes
        while (len(rankSaving) > 0):
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
            if (routeI != None 
                and routeJ != None 
                and routeI != routeJ 
                and routes[routeI]['demand'] + routes[routeJ]['demand'] <= vehCap
                and len(routes) >= vehNum):
                merge(bestSaving[1][0], bestSaving[1][1])

        # Rename the route name
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
    
    # Solve by different formulations =========================================
    res = None
    if (consAlgo == 'CWSaving'):
        res = _consVRPClarkeWright(nodes, depotID, customerID, tau, vehicle['capTruck'], vehicle['numTruck'])
    else:
        print("Error: Incorrect or unavailable CVRP formulation option!")

    return res

def heuVRPGeneral_Temp():
    return

def heuVRPConstructing(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y), 'demand': 1, ...}, \
                        nodeID2: {'loc': (x, y), 'demand': 1, ...}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...} or \
                 4) String 'Grid', will need to add arguments using `edgeArgs`"= "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                    {\
                        'colRow': (numCol, numRow),\
                        'barriers': [(coordX, coordY), ...], \
                    }" = None,
    depotID:    "DepotID, default to be 0" = 0,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime:"Service time spent on each customer (will be added into travel matrix)" = 0,
    objective:  "Objective function\
                 1) String, 'Makespan', or\
                 2) String, 'Cost'" = 'Makespan',
    vehicle:    "Dictionary, describing the vehicle situation, \
                 [{\
                    'vehicleID': vehicleID,\
                    'capacity': capacity of vehicle,\
                    'maxDist': maximum travel distance\
                 }]" = None,
    consAlgo:   "1) String 'CWSaving' or \
                 2) String 'Sweep' or \
                 3) String (not available) 'Petal' " = 'CWSaving',
    impAlgo:    "1) String '2Opt'" = '2Opt'
    ) -> "Use given heuristic methods to get basic Capacitated VRP solution":

    # Define nodeIDs ==========================================================
    customerID = []
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)
        else:
            msgError(ERROR_INCOR_NODEIDS)
            return
    customerID = [i for i in nodeIDs if i != depotID]

    # Initialize demands ======================================================
    demands = {}
    noDemandList = []
    for n in customerID:
        if ('demand' not in nodes[n]):
            noDemandList.append(n)
            demands[n] = 1
        else:
            demands[n] = nodes[n]['demand']
    if (len(noDemandList) > 0):
        msgWarning("MESSAGE: 'demands' is missing in nodes %s, set to default value (as 1)" % list2String(noDemandList))

    # Quick estimate for terminating the function =============================
    capableCus = 0
    for v in vehicle:
        capableCus += v['capacity']
    if (capableCus < len(customerIDs)):
        msgError("ERROR: Insufficient delivery capacity")
        return

    # Define tau ==============================================================
    tau = {}
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            tau = getTauEuclidean(nodes, nodeIDs)
        elif (edges == 'LatLon'):
            tau = getTauLatLon(nodes, nodeIDs)
        elif (edges == 'Grid'):
            tau = getTauGrid(nodes, nodeIDs, edgeArgs['colRow'], edgeArgs['barriers'])
        else:
            msgError(ERROR_INCOR_TAU)
            return
    else:
        tau = dict(edges)

    # Check if it is asymmetric ===============================================
    asymmetricFlag = False
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j and tau[i, j] != tau[j, i]):
                asymmetricFlag = True
                break
        if (asymmetricFlag):
            break

    # Service time ============================================================
    if (serviceTime != None and serviceTime > 0):
        for (i, j) in tau:
            if (i != depotID and j != depotID and i != j):
                tau[i, j] += serviceTime
            elif (i == depotID or j == depotID and i != j):
                tau[i, j] += serviceTime / 2 
    else:
        serviceTime = 0

    # Subroutines =============================================================
    def _consVRPClarkeWright(nodes, depotID, customerID, tau, vehicle, objective):
        # Initial routes
        routes = {}
        for i in range(len(customerID)):
            routes[i] = {
                'route': [depotID, customerID[i], depotID],
                'demand': demands[customerID[i]],
                'length': tau[depotID, customerID[i]] + tau[customerID[i], depotID],
                'reverseLength': tau[depotID, customerID[i]] + tau[customerID[i], depotID]
            }

        # Initial dist saving raking - naked saving, for symmetric VRP is already enough 
        rankSaving = []
        for i in customerID:
            for j in customerID:
                if (i != j):
                    # Calculate saving for connecting depot -> i' -> i -> j -> j' -> depot
                    sav = tau[i, depotID] + tau[depotID, j] - tau[i, j]
                    # heapq returns the smallest, so add a negative sign
                    # NOTE: this step is asymmetric
                    heapq.heappush(rankSaving, (-sav, (i, j)))

        # Initial calculation of the saving - for symmetric VRP will be the same all the time
        saving = []
        while (len(rankSaving) > 0):
            tmp = heapq.heappop(rankSaving)
            saving.append([-tmp[0], tmp[1], tmp[2]]) # Revert the saving

        # Merge routes subroutine
        # i and j are index of customer
        # After merging, i is adjacent to j 
        def merge(i, j):
            # Cannot merge the same route
            if (i == j):
                return {
                    'feasible': False,
                    'saving': None,
                    'route': None,
                    'length': None,
                    'revLength': None
                }

            # Step 1: Find routes that has i and j, and check position
            # i and j has to been either the first visiting customer or the last visiting customer
            rI = None
            rJ = None
            posI = None
            potJ = None
            for r in routes:
                if (rI == None and i in routes[r]['route']):
                    rI = r
                    if (i == routes[r]['route'][1] and i != route[r]['route'][-2]):
                        posI = 'Left'
                    elif (i != routes[r]['route'][1] and i == route[r]['route'][-2]):
                        posI = 'Right'
                    elif (i == routes[r]['route'][1] and i == route[r]['route'][-2]):
                        posI = 'Mid'
                    else:
                        return {
                            'feasible': False,
                            'saving': None,
                            'route': None,
                            'length': None,
                            'revLength': None
                        }
                if (rJ == None and j in routes[r]['route']):
                    rJ = r
                    if (j == routes[r]['route'][1] and j != route[r]['route'][-2]):
                        posJ = 'Left'
                    elif (j != routes[r]['route'][1] and j == route[r]['route'][-2]):
                        posJ = 'Right'
                    elif (j == routes[r]['route'][1] and j == route[r]['route'][-2]):
                        posJ = 'Mid'
                    else:
                        return {
                            'feasible': False,
                            'saving': None,
                            'route': None,
                            'length': None,
                            'revLength': None
                        }
                if (rI != None and rJ != None):
                    break
            if (rI == rJ):
                return {
                    'feasible': False,
                    'saving': None,
                    'route': None,
                    'length': None,
                    'revLength': None
                }

        # Merge routes
        continueFlag = True
        while (continueFlag):
            # 

            # Get the biggest saving
            
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
            if (routeI != None 
                and routeJ != None 
                and routeI != routeJ 
                and routes[routeI]['demand'] + routes[routeJ]['demand'] <= vehCap
                and len(routes) >= vehNum):
                merge(bestSaving[1][0], bestSaving[1][1])

        # Rename the route name
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
    def _consVRPSweep(nodes, depotID, customerID, tau, vehCap, vehNum, objective):
        # NOTE: This method can be applied for all objectives

        sweepSeq = getSweepSeq(
            nodes = nodes,
            nodeIDs = customerID,
            centerLoc = nodes[depotID]['locs'])


        return {
            'ofv': ofv,
            'route': route
        }

    # Solve by different formulations =========================================
    res = None
    if (consAlgo == 'CWSaving'):
        res = _consVRPClarkeWright(nodes, depotID, customerID, tau, vehicle['capTruck'], vehicle['numTruck'])
    else:
        msgError("Error: Incorrect or unavailable CVRP formulation option!")

    # Solve the local improvement =============================================
    def _lImpIntraRoute(route):
        return

    return res

def heuCVRPTW():
    return

def heuCVRPPD():
    return

def heuCVRPPDTW():
    return
