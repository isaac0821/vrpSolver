import heapq
import math

from .const import *
from .common import *
from .graph import *
from .geometry import *
from .msg import *
from .operator import *

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
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
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
    tau = getTau(nodes, edges, edgeArgs, depotID, nodeIDs, serviceTime)

    # Subroutines =============================================================
    def _consVRPClarkeWright(nodes, depotID, customerID, tau, vehicle):
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

def heuVRPMakespan(
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
    constraint: "Dictionary, describing the constraints, \
                 [{\
                    'vehCap': capacity of vehicle, \
                    'numVeh': maximum number of vehicle, \
                    'maxCost': maximum travel distance/time of a vehicle \
                 }]" = None,
    consAlgo:   "1) String 'CWSaving' or \
                 2) String (not available) 'Sweep' or \
                 3) String (not available) 'Petal' " = 'CWSaving',
    consAlgoArgs: "Dictionary" = None,
    impAlgo:    "1) String '2Opt'" = '2Opt',
    impAlgoArgs: "Dictionary" = None
    ) -> "Use given heuristic methods to get basic Capacitated VRP solution":

    # FIXME: for now, this function is only for minimizing makespan

    # Define nodeIDs ==========================================================
    customerID = []
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
        else:
            msgError(ERROR_INCOR_NODEIDS)
            return
    customerID = [i for i in nodeIDs if i != depotID]

    # Initialize demands ======================================================
    demand = {}
    noDemandList = []
    for n in customerID:
        if ('demand' not in nodes[n]):
            noDemandList.append(n)
            demand[n] = 1
        else:
            demand[n] = nodes[n]['demand']
    if (len(noDemandList) > 0):
        msgWarning("MESSAGE: 'demands' is missing in nodes %s, set to default value (as 1)" % list2String(noDemandList))

    # Quick estimate for terminating the function =============================
    capableCus = 0
    if ('vehCap' in constraint and 'numVeh' in constraint):
        capableCus = constraint['vehCap'] * constraint['numVeh']
        if (capableCus < len(customerIDs)):
            msgError("ERROR: Insufficient delivery capacity")
            return

    # Check constraint ========================================================
    if (constraint == None or ('vehCap' not in constraint and 'numVeh' not in constraint and 'maxCost' not in constraint)):
        msgError("ERROR: Missing constraint")
        return

    # Define tau ==============================================================
    tau = getTau(nodes, edges, edgeArgs, depotID, nodeIDs, serviceTime)

    # Check if it is asymmetric ===============================================
    asymFlag = False
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j and tau[i, j] != tau[j, i]):
                asymFlag = True
                break
        if (asymFlag):
            break

    # Subroutines =============================================================
    def _consVRPClarkeWright(nodes, depotID, customerID, tau, constraint):
        # Initial routes
        route = {}
        for i in range(len(customerID)):
            route[i] = {
                'route': [depotID, customerID[i], depotID],
                'demand': demand[customerID[i]],
                'cost': tau[depotID, customerID[i]] + tau[customerID[i], depotID],
                'revCost': tau[depotID, customerID[i]] + tau[customerID[i], depotID]
            }

        # Constraints
        maxDemand = None
        if ('maxDemand' in constraint):
            maxDemand = constraint['maxDemand']
        if ('maxCost' in constraint):
            maxCost = constraint['maxCost']

        # Initial dist saving raking
        mergeSaving = {}
        for i in customerID:
            for j in customerID:
                if (i < j):
                    mergeIJ = merge(
                        routeI = route[i]['route'],
                        routeJ = route[j]['route'],
                        depotID = depotID,
                        tau = tau,
                        costI = route[i]['cost'],
                        revCostI = route[i]['revCost'],
                        costJ = route[j]['cost'],
                        revCostJ = route[j]['cost'],
                        demand = demand,
                        maxDemand = maxDemand,
                        maxCost = maxCost,
                        asymFlag = asymFlag)
                    if (mergeIJ['newSeq'] != None):
                        mergeSaving[i, j] = mergeIJ

        # Merge routes
        if (objective == 'Cost'):
            # Constructing
            print("Constructing")
            return None
        elif (objective == 'Makespan'):
            # NOTE: I don't know if this is the correct implementation of CWSaving...
            print("Constructing")
            return Nond

        return {
            'ofv': ofv,
            'route': route
        }


    def _consVRPRandom(nodes, depotID, customerID, tau):
        return res

    # Solve by different formulations =========================================
    res = None
    if (consAlgo == 'CWSaving'):
        res = _consVRPClarkeWright(nodes, depotID, customerID, tau, vehicle['capTruck'], vehicle['numTruck'])
    else:
        msgError("Error: Incorrect or unavailable CVRP formulation option!")

    # Solve the local improvement =============================================
    def _lImpRoute(nodes, depotID, customerID, tau, route, asymFlag, demand, maxDemand, maxCost):

        # First, initialize the saving of removing any node from its existing positing
        # This dictionary will be updated every time two routes are updated
        removal = {}
        for r in route:
            saving = calRemovalSaving(
                route = route[r]['route'],
                cost = route[r]['cost'],
                revCost = route[r]['revCost'],
                tau = tau,
                asymFlag = asymFlag)
            for n in saving:
                removal[n] = saving[n]
        
        canImproveFlag = True
        while(canImproveFlag):
            canImproveFlag = False

            # Stage 1: Try moving a customer to another route
            for n in nodeIDs:
                # Find which route is the node to be removed comes from
                removalRoute = None
                for r in route:
                    if (n in route[r]['route']):
                        removalRoute = r
                        break

                # Try to find a route that can insert the node
                for r in route:
                    # If the node is not in route r, try insert it
                    if (n not in route[r]['route']):
                        insertion = calInsertionCost(
                            route = route[r]['route'],
                            cost = route[r]['cost'],
                            revCost = route[r]['revCost'],
                            tau = tau,
                            nJ = n,
                            demand = demand,
                            maxDemand = maxDemand,
                            maxCost = maxCost,
                            asymFlag = asymFlag)

                        # If insertion is feasible, try to see if it is improving makespan between these two vehicle
                        if (insertion['newSeq'] != None):
                            oldMakespan = max(route[removalRoute]['cost'], route[r]['cost'])
                            newMakespan = max(removal[n]['newCost'], insertion['newCost'])
                            if (newMakespan < oldMakespan):
                                canImproveFlag = True
                                # Update the route that has node removed
                                route[removalRoute]['route'] = removal[n]['newSeq']
                                route[removalRoute]['cost'] = removal[n]['cost']
                                route[removalRoute]['revCost'] = removal[n]['revCost']
                                # Update the route that has node inserted
                                route[r]['route'] = insertion['newSeq']
                                route[r]['cost'] = insertion['newCost']
                                route[r]['revCost'] = insertion['revNewCost']
                                # Update removal dictionary
                                for updateR in [removalRoute, r]:
                                    saving = calRemovalSaving(
                                        route = route[updateR]['route'],
                                        cost = route[updateR]['cost'],
                                        revCost = route[updateR]['revCost'],
                                        tau = tau,
                                        asymFlag = asymFlag)
                                    for updateN in saving:
                                        removal[updateN] = saving[updateN]

            # Stage 2: Try swapping two customers in the route
            for nI in nodeIDs:
                for nJ in nodeIDs:
                    # FIXME: Skip for now
                    pass

            # Stage 3: Try improving each route using TSP
            for r in routes:
                updateTSP = heuTSP(
                    nodes = nodes,
                    edges = tau,
                    depotID = depotID,
                    nodeIDs = [route[r]['route'][i] for i in range(1, len(route[r]['route']) - 1)],
                    serviceTime = 0, # Do not update tau again
                    consAlgo = 'PreDefine',
                    consAlgoArgs = {'initSeq': route[r]['route']})
                if (updateTSP['ofv'] < route[r]['cost']):
                    rotue[r]['route'] = updateTSP['seq']
                    route[r]['cost'] = updateTSP['ofv']
                    rotue[r]['revCost'] = calSeqCostMatrix(
                        tau, 
                        [updateTSP['seq'][len(updateTSP['seq']) - 1 - i] for i in range(len(updateTSP['seq']))])

    return res

