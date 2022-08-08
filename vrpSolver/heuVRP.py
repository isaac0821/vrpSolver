import heapq
import math
import datetime
import warnings

from .const import *
from .common import *
from .graph import *
from .geometry import *
from .msg import *
from .operator import *
from .calculate import *
from .heuTSP import *

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
                 2) A list of node IDs, all customer IDs, may also include the depot(s)" = 'All',
    serviceTime:"Service time spent on each customer (will be added into travel matrix)" = 0,
    objective:  "Objective function\
                 1) String, 'Min_Max', or\
                 2) String, 'Min_Cost' " = 'Min_Max',
    constraint: "Dictionary, describing the constraints, \
                 {\
                    'vehCap': capacity of vehicle, \
                    'numVeh': maximum number of vehicle, \
                    'maxCost': maximum travel distance/time of a vehicle \
                 }" = None,
    consAlgo:   "If `objective == 'Min_Cost'` \
                 1) String 'CWSaving' or \
                 2) String 'Sweep' or \
                 3) String (not available) 'Petal' or\
                 4) String (not available) 'Cluster-First-Route-Second' or\
                 5) String (not available) 'Route-First-Cluster-Second' or\
                 If `objective == 'Min_Max'` \
                 1) String 'Sweep' " = 'Sweep',
    consAlgoArgs: "Dictionary" = None,
    impAlgo:    "1) String '2Opt'" = '2Opt',
    impAlgoArgs: "Dictionary" = None
    ) -> "Use given heuristic methods to get basic Capacitated VRP solution":

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
        warnings.warn("MESSAGE: 'demands' is missing in nodes %s, set to default value (as 1)" % list2String(noDemandList))

    # Decode and pre-check constraints ========================================
    if (constraint == None or (
        'vehCap' not in constraint 
        and 'numVeh' not in constraint 
        and 'maxCost' not in constraint)):
        msgError("ERROR: Missing constraint")
        return
    # List of constraints 
    vehCap = None
    if ('vehCap' in constraint):
        vehCap = constraint['vehCap']
    vehMaxDist = None
    if ('vehMaxDist' in constraint):
        vehMaxDist = constraint['vehMaxDist']
    numVeh = None
    if ('numVeh' in constraint):
        numVeh = constraint['numVeh']
    # Quick estimate for terminating the function
    # 1. Total capacity
    if (vehCap != None and numVeh != None and vehCap * numVeh < len(customerIDs)):
        msgError("ERROR: Insufficient delivery capacity")
        return
    # 2. Single visit v.s. capacity
    if (vehCap != None and max(demand.values()) > vehCap):
        msgError("ERROR: Insufficient capacity to delivery single demand")
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

    # Solve by different formulations =========================================
    route = None
    if (objective == 'Min_Cost'):
        if (consAlgo == 'CWSaving'):
            route = _consVRPClarkeWright(nodes, depotID, customerID, tau, vehCap, numVeh)
    elif (objective == 'Min_Max'):
        if (consAlgo == 'Sweep'):
            route = _consVRPSweep(nodes, depotID, customerID, tau, numVeh)
    else:
        msgError("Error: Incorrect or unavailable CVRP objective option!")

    if (objective == 'Min_Cost'):
        pass
    elif (objective == 'Min_Max'):
        route = _lImpRouteMinMax(
            nodes = nodes,
            depotID = depotID,
            customerID = customerID,
            tau = tau,
            route = route,
            asymFlag = asymFlag)

    # Stat ====================================================================
    cost = 0
    for r in route:
        cost += route[r]['cost']
    makespan = max([route[r]['cost'] for r in route])

    return {
        'cost': cost,
        'makespan': makespan,
        'route': route,
        'serviceTime': serviceTime
    }

# [Constructing]
def _consVRPInsertion():
    # Initial route
    routeBank = {}
    for i in range(len(customerID)):
        if ((depotID, customerID[i]) in tau and (customerID[i], depotID) in tau):
            routeBank[i] = {
                'route': [depotID, customerID[i], depotID],
                'demand': demand[customerID[i]],
                'cost': tau[depotID, customerID[i]] + tau[customerID[i], depotID],
                'revCost': tau[depotID, customerID[i]] + tau[customerID[i], depotID]
            }
    # Initial dist saving ranking
    mergeSaving = {}
    for i in customerID:
        for j in customerID:
            if (i < j):
                mergeIJ = merge(
                    routeI = routeBank[i]['route'],
                    routeJ = routeBank[j]['route'],
                    depotID = depotID,
                    tau = tau,
                    costI = routeBank[i]['cost'],
                    revCostI = route[i]['revCost'],
                    costJ = routeBank[j]['cost'],
                    revCostJ = route[j]['cost'],
                    demand = demand,
                    maxDemand = vehCap,
                    maxCost = None,
                    asymFlag = asymFlag)
                if (mergeIJ['newSeq'] != None):
                    mergeSaving[i, j] = mergeIJ

    # Merge route
    if (objective == 'Cost'):
        # Constructing
        print("Constructing")
        return None
    elif (objective == 'Makespan'):
        # NOTE: I don't know if this is the correct implementation of CWSaving...
        route = {}
        for i in range(numVeh):
            route[i] =  {
                'route': [depotID, customerID[i], depotID],
                'demand': demand[customerID[i]],
                'cost': tau[depotID, customerID[i]] + tau[customerID[i], depotID],
                'revCost': tau[depotID, customerID[i]] + tau[customerID[i], depotID]
            }

        for i in range(numVeh, len(customerID)):
            for v in range(numVeh):
                pass
        return None

    return {
        'ofv': ofv,
        'route': route
    }

def _consVRPClarkeWright(nodes, depotID, customerID, tau, vehCap, numVeh):
    # Initial route
    routeBank = {}
    for i in range(len(customerID)):
        if ((depotID, customerID[i]) in tau and (customerID[i], depotID) in tau):
            if (vehMaxDist == None or vehMaxDist > tau[depotID, customerID[i]] + tau[customerID[i], depotID]):
                routeBank[i] = {
                    'route': [depotID, customerID[i], depotID],
                    'demand': demand[customerID[i]],
                    'length': tau[depotID, customerID[i]] + tau[customerID[i], depotID]
                }

    # Initial saving raking
    rankSaving = []
    for i in customerID:
        for j in customerID:
            if (i != j):
                # Calculate saving for each pair
                sav = tau[depotID, i] + tau[depotID, j] - tau[i, j]
                # heapq returns the smallest, so add a negative sign
                heapq.heappush(rankSaving, (-sav, (i, j)))

    # Merge route subroutine
    def mergeLocal(i, j):
        if (i == j):
            return None
        rI = None
        rJ = None
        iLeft = None
        iRight = None
        jLeft = None
        jRight = None
        for r in routeBank:
            if (i == routeBank[r]['route'][1]):
                iLeft = True
                rI = r
            if (i == routeBank[r]['route'][-2]):
                iRight = True
                rI = r
            if (j == routeBank[r]['route'][1]):
                jLeft = True
                rJ = r
            if (j == routeBank[r]['route'][-2]):
                jRight = True
                rJ = r
        newRoute = []
        if (iRight == True and jLeft == True):
            newRoute = [i for i in routeBank[rI]['route']]
            addRoute = [i for i in routeBank[rJ]['route']]
            newRoute.extend(addRoute)
        elif (iLeft == True and jRight == True):
            newRoute = [i for i in routeBank[rJ]['route']]
            addRoute = [i for i in routeBank[rI]['route']]
            newRoute.extend(addRoute)
        elif (iLeft == True and jLeft == True):
            newRoute = [i for i in routeBank[rI]['route']]
            newRoute.reverse()
            addRoute = [i for i in routeBank[rJ]['route']]
            newRoute.extend(addRoute)
        elif (iRight == True and jRight == True):
            newRoute = [i for i in routeBank[rI]['route']]
            addRoute = [i for i in routeBank[rJ]['route']]
            addRoute.reverse()
            newRoute.extend(addRoute)

        while (depotID in newRoute):
            newRoute.remove(depotID)
        newRoute.insert(0, depotID)
        newRoute.append(depotID)

        newDemand = routeBank[rI]['demand'] + routeBank[rJ]['demand']
        newLength = routeBank[rI]['length'] + routeBank[rJ]['length'] + tau[i, j] - tau[depotID, i] - tau[depotID, j]
        routeBank.pop(rI)
        routeBank.pop(rJ)
        newRouteIndex = max(list(routeBank.keys())) + 1
        routeBank[newRouteIndex] = {
            'route': newRoute,
            'demand': newDemand,
            'length': newLength
        }

    # Merge route
    while (len(rankSaving) > 0):
        # Get the biggest saving
        bestSaving = heapq.heappop(rankSaving)
        # Flip it back
        sav = -bestSaving[0]
        # If there is saving, check which two route can be merged
        routeI = None
        routeJ = None
        for r in routeBank:
            if (bestSaving[1][0] == routeBank[r]['route'][1] or bestSaving[1][0] == routeBank[r]['route'][-2]):
                routeI = r
            if (bestSaving[1][1] == routeBank[r]['route'][1] or bestSaving[1][1] == routeBank[r]['route'][-2]):
                routeJ = r
            if (routeI != None and routeJ != None):
                break
        # Two route has to be different, and satisfied the capacity
        if (routeI != None 
            and routeJ != None 
            and routeI != routeJ 
            and routeBank[routeI]['demand'] + routeBank[routeJ]['demand'] <= vehCap
            and len(routeBank) > numVeh):
            mergeLocal(bestSaving[1][0], bestSaving[1][1])

    # Rename the route name
    route = {}
    acc = 1
    for r in routeBank:
        route[acc] = {
            'route': [i for i in routeBank[r]['route']],
            'cost': calSeqCostMatrix(tau, [i for i in routeBank[r]['route']]),
            'revCost': calSeqCostMatrix(tau, [routeBank[r]['route'][len(routeBank[r]['route']) - 1 - i] for i in range(len(routeBank[r]['route']))])
        }
        acc += 1

    return route

def _consVRPSweep(nodes, depotID, customerID, tau, numVeh):
    # FIXME: currently I've not yet checked the cases where graphs are incomplete
    # Order the customers in the sequence of sweeping
    sweepSeq = getSweepSeq(
        nodes = nodes,
        nodeIDs = customerID,
        centerLoc = nodes[depotID]['loc'])

    # Initialize the solution by (trying to) evenly partition customers to vehicles
    cusPerVeh = math.ceil(len(customerID) / numVeh)

    route = {}
    # for r in range(numVeh):

    split = [sweepSeq[i:i + cusPerVeh] for i in range(0, len(sweepSeq), cusPerVeh)]

    route = {}
    for r in range(numVeh):
        sub = [depotID]
        sub.extend(split[r])
        sub.append(depotID)
        revSub = [i for i in sub]
        revSub.reverse()
        route[r] = {
            'route': sub,
            'cost': calSeqCostMatrix(tau, sub),
            'revCost': calSeqCostMatrix(tau, revSub)
        }
    return route

def _lImpRouteMinMax(nodes, depotID, customerID, tau, route, asymFlag):
    # First, initialize the saving of removing any node from its existing positing
    # This dictionary will be updated every time two route are updated
    removal = {}

    # Some flags
    # For each route, log if the route has been 
    inrouteOptimizeFlag = {}
    for r in route:
        inrouteOptimizeFlag[r] = False
    
    canImproveFlag = True
    while(canImproveFlag):
        canImproveFlag = False

        # Stage 1: Try moving a customer to another route
        if (not canImproveFlag):
            for r in route:
                for i in range(1, len(route[r]['route']) - 1):
                    saving = calRemovalSaving(
                        route = route[r]['route'],
                        tau = tau,
                        i = i,
                        cost = route[r]['cost'],
                        revCost = route[r]['revCost'],
                        asymFlag = asymFlag)
                    removal[route[r]['route'][i]] = saving

            for n in removal:
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
                            tau = tau,
                            nJ = n,
                            cost = route[r]['cost'],
                            revCost = route[r]['revCost'],
                            asymFlag = asymFlag)

                        # If insertion is feasible, try to see if it is improving makespan between these two vehicle
                        if (insertion != None and insertion['newSeq'] != None):
                            oldMakespan = max(route[removalRoute]['cost'], route[r]['cost'])
                            newMakespan = max(removal[n]['newCost'], insertion['newCost'])
                            if (newMakespan + CONST_EPSILON < oldMakespan):                                    
                                # Update the route that has node removed
                                route[removalRoute]['route'] = removal[n]['newSeq']
                                route[removalRoute]['cost'] = removal[n]['newCost']
                                route[removalRoute]['revCost'] = removal[n]['newRevCost']

                                # Update the route that has node inserted
                                route[r]['route'] = insertion['newSeq']
                                route[r]['cost'] = insertion['newCost']
                                route[r]['revCost'] = insertion['newRevCost']

                                print('Move %s from route %s to %s, oldMakespan: %s, newMakespan: %s' % (n, removalRoute, r, oldMakespan, newMakespan))
                                canImproveFlag = True
                                break
                if (canImproveFlag):
                    break

        # Stage 2: Try swapping two customers in the route
        if (not canImproveFlag):
            for nI in customerID:
                for nJ in customerID:
                    # FIXME: Skip for now
                    pass

        # Stage 3: Try improving each route using TSP
        if (not canImproveFlag):
            for r in route:
                if (not inrouteOptimizeFlag[r]):
                    if (len(route[r]['route']) > 4):
                        updateTSP = heuTSP(
                            nodes = nodes,
                            edges = tau,
                            depotID = depotID,
                            nodeIDs = [route[r]['route'][i] for i in range(len(route[r]['route']) - 1)],
                            serviceTime = 0) # Do not update tau again
                        inrouteOptimizeFlag[r] = True
                        if (updateTSP['ofv'] + CONST_EPSILON < route[r]['cost']):
                            canImproveFlag = True
                            print("Improve Route %s from %s to %s" % (r, route[r]['cost'], updateTSP['ofv']))
                            route[r]['route'] = updateTSP['seq']
                            route[r]['cost'] = updateTSP['ofv']
                            revSeq = [i for i in updateTSP['seq']]
                            revSeq.reverse()
                            route[r]['revCost'] = calSeqCostMatrix(tau, revSeq, True)

    return route


