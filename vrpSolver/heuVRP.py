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
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
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
    vehCap:     "Capacity of each vehicle" = None,
    consAlgo:   "1) String 'CWSaving' or \
                 2) String (not available) 'Sweep'" = 'CWSaving',
    impAlgo:    "1) String '2Opt'" = '2Opt'
    ) -> "Use given heuristic methods to get TSP solution":

    # Define nodeIDs ==========================================================
    customerID = []
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)
        else:
            print(ERROR_INCOR_NODEIDS)
            return
    customerID = [i for i in nodeIDs if i != depotID]

    # Initialize demands ======================================================
    demands = {}
    for n in customerID:
        if ('demand' not in nodes[n]):
            demands[n] = 1
        else:
            demands[n] = nodes[n]['demand']

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
            print(ERROR_INCOR_TAU)
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
    def _consVRPClarkeWright(nodes, depotID, customerID, tau, vehCap):
        # Initial routes
        routes = {}
        for i in range(len(customerID)):
            routes[i] = {
                'route': [depotID, customerID[i], depotID],
                'demand': demands[customerID[i]],
                'length': 2 * tau[depotID, customerID[i]]
            }

        # vehCap
        if (vehCap == None):
            vehCap = len(nodeIDs) + 1

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
                and routes[routeI]['demand'] + routes[routeJ]['demand'] <= vehCap):
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
        res = _consVRPClarkeWright(nodes, depotID, customerID, tau, vehCap)
    else:
        print("Error: Incorrect or unavailable CVRP formulation option!")

    return res
