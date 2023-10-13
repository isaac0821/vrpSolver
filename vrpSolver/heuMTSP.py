import heapq
import math
import datetime
import warnings
import networkx as nx

from .const import *
from .common import *
from .geometry import *
from .msg import *

# History =====================================================================
# 20231001 - Start writing script
# =============================================================================

def heuMTSP(
    nodes: dict,
    locFieldName: str = 'loc',
    numVeh: int = 1,
    vehicles: dict = {},
    objective: str = 'MinMakespan',
    edges: dict = {'method': 'Euclidean', 'ratio': 1},
    method: dict = {'cons': 'Insertion', 'impv': '2Opt'},
    depotID: int | str = 0,
    nodeIDs: list[int|str] | str = 'All',
    serviceTime: float = 0,
    returnRouteObjectFlag = False
    ) -> dict | None:

    """Use heuristic methods to find suboptimal MTSP solution


    Parameters
    ----------

    nodes: dictionary, required, default None
        The coordinates of given nodes, in the following format::
            >>> nodes = {
            ...     nodeID1: {'loc': (x, y)},
            ...     nodeID2: {'loc': (x, y)}, # ...
            ... }
    edges: dictionary, required, default as {'method': "Euclidean", 'ratio': 1}
        The traveling matrix. The options are as follows::
            1) (default) Euclidean space
            >>> edge = {
            ...     'method': 'Euclidean',
            ...     'ratio': 1 # Optional, default to be 1
            ... }
            2) By given pairs of lat/lon
            >>> edge = {
            ...     'method': 'LatLon',
            ...     'unit': 'meters' # Optional, default to be 1
            ... }
            3) ManhattenDistance
            >>> edge = {
            ...     'method': 'Manhatten',
            ...     'ratio': 1 # Optional, default to be 1
            ... }
            4) By a given dictionary
            >>> edge = {
            ...     'method': 'Dictionary',
            ...     'dictionary': dictionary,
            ...     'ratio': 1 # Optional, default to be 1
            ... }
            5) On the grids
            >>> edge = {
            ...     'method': 'Grid',
            ...     'grid': grid
            ... }

    """

    # Sanity check ============================================================
    if (nodes == None or type(nodes) != dict):
        raise MissingParameterError(ERROR_MISSING_NODES)
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
        else:
            for i in nodeIDs:
                if (i not in nodes):
                    raise OutOfRangeError("ERROR: Node %s is not in `nodes`." % i)
    if ((type(nodeIDs) == list and depotID not in nodeIDs)
        or (nodeIDs == 'All' and depotID not in nodes)):
        raise OutOfRangeError("ERROR: Cannot find `depotID` in given `nodes`/`nodeIDs`")
    if ((numVeh == None or numVeh <= 0) and (vehicles == {} or vehicles == None)):
        raise MissingParameterError("ERROR: Missing vehicle definitions.")
    if (vehicles == None or vehicles == {}):
        vehicles = {}
        for i in range(numVeh):
            vehicles[i] = {}

    # Define tau ==============================================================
    tau, path = matrixDist(nodes, edges, depotID, nodeIDs, serviceTime)
    
    # Check symmetric =========================================================
    asymFlag = False
    for (i, j) in tau:
        if (tau[i, j] != tau[j, i]):
            asymFlag = True
            break

    # Construction heuristics =================================================
    # NOTE: Output of this phase should be a Route() object
    if (method == None or ('cons' not in method and 'impv' not in method)):
        raise MissingParameterError(ERROR_MISSING_TSP_ALGO)

    nodeObj = {}
    for n in nodeIDs:
        nodeObj[n] = RouteNode(n, value=nodes[n][locFieldName])

    for veh in vehicles:
        vehicles[veh]['routeObj'] = Route(tau, asymFlag)
        vehicles[veh]['routeObj'].append(nodeObj[depotID].clone())

    if ('cons' not in method):
        raise MissingParameterError("ERROR: Missing 'cons' in `method`.")

    elif (method['cons'] == 'Sweep'):
        vehicles = _consMTSPSweep(nodes, nodeIDs, nodeObj, depotID, vehicles, locFieldName)

    elif (method['cons'] == 'Random'):
        vehicles = _consMTSPRandom(nodeIDs, nodeObj, depotID, vehicles)

    elif (method['cons'] == 'Insertion'):
        vehicles = _consMTSPInsertionMinMakespan(nodeIDs, nodeObj, depotID, vehicles)

    elif (method['cons'] == 'RandomInsertion'):
        vehicles = _consMTSPInsertionMinMakespan(nodeIDs, nodeObj, depotID, vehicles, randomInsertionFlag = True)

    else:
        raise UnsupportedInputError("ERROR: Choose 'cons' from ['Sweep', 'Random']")

    # Local improvement =======================================================
    if ('impv' in method and method['impv'] != None and method['impv'] != []):
        canImpvFlag = True
        while (canImpvFlag):
            canImpvFlag = False

            if (not canImpvFlag and '2Opt' in method['impv']):
                canImpvFlag = _impvMTSP2Opt(vehicles)

    return vehicles

def _consMTSPSweep(nodes, nodeIDs, nodeObj, depotID, vehicles, locFieldName):
    sweep = nodeSeqBySweeping(
        nodes = nodes,
        nodeIDs = [i for i in nodes if i != depotID],
        centerLoc = nodes[depotID][locFieldName])
    sweepPerVeh = splitList(sweep, len(vehicles))
    vehIDs = [v for v in vehicles]
    for i in range(len(vehIDs)):
        for k in range(len(sweepPerVeh[i])):
            vehicles[vehIDs[i]]['routeObj'].cheapestInsert(nodeObj[sweepPerVeh[i][k]])
    return vehicles

def _consMTSPRandom(nodeIDs, nodeObj, depotID, vehicles):
    rndNodeSeq = [i for i in nodeIDs if i != depotID]
    random.shuffle(rndNodeSeq)
    rndPerVeh = splitList(rndNodeSeq, len(vehicles))
    vehIDs = [v for v in vehicles]
    for i in range(len(vehIDs)):
        for k in range(len(rndPerVeh[i])):
            vehicles[vehIDs[i]]['routeObj'].cheapestInsert(nodeObj[rndPerVeh[i][k]])
    return vehicles

def _consMTSPInsertionMinMakespan(nodeIDs, nodeObj, depotID, vehicles, randomInsertionFlag=False):
    unInserted = [i for i in nodeIDs if i != depotID]
    if (randomInsertionFlag):
        random.shuffle(unInserted)
    for n in unInserted:
        cost = {}
        for veh in vehicles:
            oldCost = max(vehicles[veh]['routeObj'].dist for veh in vehicles)
            vehicles[veh]['routeObj'].cheapestInsert(nodeObj[n])
            newCost = max(vehicles[veh]['routeObj'].dist for veh in vehicles)
            cost[veh] = newCost - oldCost
            vehicles[veh]['routeObj'].remove(nodeObj[n])
        cheapestAmongVeh = min(cost, key = lambda k:cost[k])
        vehicles[cheapestAmongVeh]['routeObj'].cheapestInsert(nodeObj[n])
    return vehicles

def _impvMTSP2Opt(vehicles):
    for veh in vehicles:
        vehImpvFlag = vehicles[veh]['routeObj'].impv2Opt()
        if (vehImpvFlag):
            return True
    return False

def _impvMTSP2OptStar(vehicles):
    vehIDs = [i for i in vehicles]
    for i in range(len(vehIDs) - 2):
        for j in range(i + 1, len(vehIDs)):
            vehI = vehIDs[i]
            vehJ = vehIDs[j]

    return False