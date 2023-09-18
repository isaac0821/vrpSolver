import heapq
import math
import datetime
import warnings
import networkx as nx

from .const import *
from .common import *
from .geometry import *
from .msg import *
# from .operator import *
# from .calculate import *

# History =====================================================================
# 20230510 - Cleaning
# 20230822 - Rewrite Christofide using networkx
# =============================================================================

def heuTSP(
    nodes: dict, 
    edges: dict = {'method': "Euclidean", 'ratio': 1}, 
    algo: dict = {'cons': 'Insertion', 'impv': '2opt'}, 
    depotID: int | str = 0, 
    nodeIDs: list[int | str] | str = 'All', 
    serviceTime: float = 0,
    returnRouteObjectFlag = False
    ) -> dict | None:

    """Use heuristic methods to find suboptimal TSP solution

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
    algo: dictionary, required, default as {'cons': 'Insertion', 'impv': '2opt'}
        The algorithm configuration. Includes two phases, use 'cons' to specify constructive heuristic, and 'impv' to specify local improvement heurisitc::
            1) (default) Insertion
            >>> algo = {
            ...     'cons': 'Insertion',
            ...     'initSeq': initSeq, # An initial sequence, defalt [depotID]
            ...     'impv': '2Opt' # Options are: 'Reinsert', '2Opt', can select multiple methods by collecting them into a list, e.g. ['Reinsert', '2Opt']
            ... }
            2) Nearest neighborhood / k-nearest neighborhood
            >>> algo = {
            ...     'cons': 'NearestNeighbor',
            ...     'k': 1, # 1: nearest neighbor, 2 ~ K: k-nearest neighbor, -1: farthest neighbor 
            ... }
            3) Sweep
            >>> algo = {
            ...     'cons': 'Sweep'
            ... }
            4) (not available) Christofides
            >>> algo = {
            ...     'cons': 'Christofides'
            ... }
            5) (not available) Cycle cover, particular for Asymmetric TSP
            >>> algo = {
            ...     'cons': 'CycleCover'
            ... }
            6) Random sequence
            >>> algo = {
            ...     'cons': 'Random'
            ... }
            7) Given sequence for further local improvements
            >>> algo = {
            ...     'cons': None, # or skip this
            ...     'initSeq': initSeq, # An initial sequence, cannot be None in this case
            ... }
    depotID: int or string, required, default as 0
        The ID of depot.
    nodeIDs: string 'All' or a list of node IDs, required, default as 'All'
        The following are two options: 1) 'All', all nodes will be visited, 2) A list of node IDs to be visited.
    serviceTime: float, optional, default as 0
        The service time needed at each location.
    returnRouteObjectFlag: bool, optional, default as False
        If true, in the 'seq' field of output, return the Route() object instead of a list

    Returns
    -------

    dictionary
        A TSP solution in the following format::
            >>> solution = {
            ...     'ofv': ofv,
            ...     'route': route,
            ...     'detail': detail
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
    if (algo == None or ('cons' not in algo and 'impv' not in algo)):
        raise MissingParameterError(ERROR_MISSING_TSP_ALGO)

    nodeObj = {}
    for n in nodeIDs:
        nodeObj[n] = RouteNode(n, value=nodes[n]['loc'])

    seq = Route(tau, asymFlag)
    # An initial solution is given
    if ('cons' not in algo or algo['cons'] == 'Initial' or algo['cons'] == None):
        if ('initSeq' not in algo):
            raise MissingParameterError("ERROR: Need 'initSeq' for local improvement")
        elif (len(algo['initSeq']) != len(nodeIDs) + 1):
            raise UnsupportedInputError("ERROR: Length of 'initSeq' is incorrect, check if the sequence starts and ends with `depotID`")
        else:
            notInNodeIDs = [v for v in algo['initSeq'] if v not in nodeIDs]
            if (len(notInNodeIDs) > 0):
                raise OutOfRangeError("ERROR: The following nodes in 'initSeq' is not in `nodeIDs`: %s" % list2String(notInNodeIDs))
        for i in algo['initSeq'][:-1]:
            seq.append(nodeObj[i])

    # Insertion heuristic
    elif (algo['cons'] == 'Insertion'):
        initSeq = None
        if ('initSeq' not in algo):   
            farthestDist = -1
            farthestID = None
            for n in nodeIDs:
                if ((n, depotID) in tau):
                    if (tau[n, depotID] > farthestDist):
                        farthestID = n
                        farthestDist = tau[n, depotID]
                elif ((depotID, n) in tau):
                    if (tau[depotID, n] > farthestDist):
                        farthestID = n
                        farthestDist = tau[depotID, n]
            initSeq = [depotID, farthestID, depotID]
            seq = _consTSPInsertion(nodeIDs, initSeq, tau, asymFlag)
        else:
            notInNodeIDs = [v for v in algo['initSeq'] if v not in nodeIDs]
            if (len(notInNodeIDs) > 0):
                raise OutOfRangeError("ERROR: The following nodes in 'initSeq' is not in `nodeIDs`: %s" % list2String(notInNodeIDs))
            else:
                seq = _consTSPInsertion(nodeIDs, algo['initSeq'], tau, asymFlag)
    
    # Neighborhood based heuristic, including nearest neighborhood, k-nearest neighborhood, and furthest neighborhood
    elif (algo['cons'] == 'NearestNeighbor'):
        nnSeq = None
        if ('k' not in algo or algo['k'] == 1):
            nnSeq = _consTSPkNearestNeighbor(depotID, nodeIDs, tau, 1)
        elif (algo['k'] == -1):
            nnSeq = _consTSPFarthestNeighbor(depotID, nodeIDs, tau)
        elif (algo['k'] >= 1):
            nnSeq = _consTSPkNearestNeighbor(depotID, nodeIDs, tau, algo['k'])
        for i in nnSeq:
            seq.append(nodeObj[i])

    # Sweep heurisitic
    elif (algo['cons'] == 'Sweep'):
        sweepSeq = _consTSPSweep(nodes, depotID, nodeIDs, tau)
        for i in sweepSeq:
            seq.append(nodeObj[i])
        seq.rehead(depotID)

    # Christofides Algorithm, guaranteed <= 1.5 * optimal
    elif (algo['cons'] == 'Christofides'):
        if (not asymFlag):
            cfSeq = _consTSPChristofides(depotID, tau)
            for i in cfSeq:
                seq.append(nodeObj[i])
        else:
            raise UnsupportedInputError("ERROR: 'Christofides' algorithm is not designed for Asymmetric TSP")

    # Cycle Cover Algorithm, specially designed for Asymmetric TSP
    elif (algo['cons'] == 'CycleCover'):
        raise VrpSolverNotAvailableError("ERROR: 'CycleCover' algorithm is not available yet, please stay tune")
        seq = _consTSPCycleCover(depotID, nodeIDs, tau)

    # Randomly create a sequence
    elif (algo['cons'] == 'Random'):
        rndSeq = _consTSPRandom(depotID, nodeIDs, tau)
        for i in rndSeq:
            seq.append(nodeObj[i])
        seq.rehead(depotID)
    else:
        raise UnsupportedInputError(ERROR_MISSING_TSP_ALGO)

    # Cleaning seq before local improving =================================
    consOfv = seq.dist

    # Local improvment phase ==============================================
    # NOTE: Local improvement phase operates by class methods
    # NOTE: For the local improvement, try every local search operator provided in a greedy way
    if ('impv' in algo and algo['impv'] != None and algo['impv'] != []):
        canImpvFlag = True
        while (canImpvFlag):
            canImpvFlag = False

            # if (not canImpvFlag):
            #     canImpvFlag = seq.impv2Swap()

            # if (not canImpvFlag):
            #     canImpvFlag = seq.impvRemovalReinsert()


            if (not canImpvFlag):
                canImpvFlag = seq.impv2Opt()

    ofv = seq.dist
    nodeSeq = [n.key for n in seq.traverse(closeFlag = True)]

    shapepoints = []
    for i in range(len(nodeSeq) - 1):
        shapepoints.extend(path[nodeSeq[i], nodeSeq[i + 1]])

    return {
        'ofv': ofv,
        'consOfv': consOfv,
        'algo': algo,
        'seq': nodeSeq,
        # 'seqRoute': seq,
        'shapepoints': shapepoints,
        'serviceTime': serviceTime
    }

def _consTSPkNearestNeighbor(depotID, nodeIDs, tau, k = 1):
    # Initialize ----------------------------------------------------------
    seq = [depotID]
    remain = [nodeIDs[i] for i in range(len(nodeIDs)) if nodeIDs[i] != depotID]
    # Accumulate seq ------------------------------------------------------
    while (len(remain) > 0):
        currentNodeID = seq[-1]

        # Sort the distance from current node to the rest of nodes
        sortedSeqHeap = []
        for n in remain:
            if ((currentNodeID, n) in tau):
                dist = tau[currentNodeID, n]
                heapq.heappush(sortedSeqHeap, (dist, n))

        # Get the kth of sorted node and append it to seq
        nextNodeID = None
        for _ in range(k):
            if (len(sortedSeqHeap) > 0):
                nextNodeID = heapq.heappop(sortedSeqHeap)[1]
        seq.append(nextNodeID)

        # Update remain
        remain.remove(nextNodeID)
    return seq

def _consTSPFarthestNeighbor(depotID, nodeIDs, tau):
    # Initialize ----------------------------------------------------------
    seq = [depotID]
    remain = [nodeIDs[i] for i in range(len(nodeIDs)) if nodeIDs[i] != depotID]
    # Accumulate seq ------------------------------------------------------
    while (len(remain) > 0):
        nextLeng = None
        nextID = None
        for n in remain:
            if ((n, seq[-1]) in tau):
                if (nextLeng == None or tau[n, seq[-1]] > nextLeng):
                    nextID = n
                    nextLeng = tau[n, seq[-1]]
            elif ((seq[-1], n) in tau):
                if (nextLeng == None or tau[seq[-1], n] > nextLeng):
                    nextID = n
                    nextLeng = tau[seq[-1], n]
        seq.append(nextID)
        remain.remove(nextID)
    return seq

def _consTSPSweep(nodes, depotID, nodeIDs, tau):
    # Sweep seq -----------------------------------------------------------
    sweep = nodeSeqBySweeping(
        nodes = nodes, 
        nodeIDs = nodeIDs,
        centerLoc = nodes[depotID]['loc'])
    return sweep

def _consTSPRandom(depotID, nodeIDs, tau):
    # Get random seq ------------------------------------------------------
    seq = [i for i in nodeIDs]
    random.shuffle(seq)
    return seq

def _consTSPInsertion(nodeIDs, initSeq, tau, asymFlag):
    # Initialize ----------------------------------------------------------
    routeNodes = {}
    for n in nodeIDs:
        routeNodes[n] = RouteNode(key = n, value = n)
    route = Route(tau, asymFlag)

    # NOTE: initSeq should starts and ends with depotID
    for n in initSeq[:-1]:
        route.append(routeNodes[n])

    unInserted = [i for i in nodeIDs if i not in initSeq[:-1]]
    for n in unInserted:
        route.cheapestInsert(routeNodes[n])

    return route
    
def _consTSPChristofides(depotID, tau):
    G = nx.Graph()
    subG = nx.Graph()
    for (i, j) in tau:
        G.add_edge(i, j, weight=tau[i, j])
    mst = nx.minimum_spanning_tree(G)

    degreeOfEachNodes = {}
    for n in nodes:
        degreeOfEachNodes[n] = 0
    for edge in mst.edges:
        degreeOfEachNodes[edge[0]] += 1
        degreeOfEachNodes[edge[1]] += 1
        subG.add_edge(edge[0], edge[1], weight=tau[edge[0], edge[1]])

    matchG = nx.Graph()
    oddDegrees = []
    for n in degreeOfEachNodes:
        if (degreeOfEachNodes[n] % 2 != 0):
            oddDegrees.append(n)
    for i in oddDegrees:
        for j in oddDegrees:
            if ((i, j) in tau and i < j and (i, j) not in subG.edges and (j, i) not in subG.edges):
                matchG.add_edge(i, j, weight=tau[i, j])
    
    mwm = nx.min_weight_matching(matchG)
    for edge in mwm:
        subG.add_edge(edge[0], edge[1], weight=tau[edge[0], edge[1]])
    s = nx.eulerian_circuit(subG, source=depotID)

    seq = []
    for i in s:
        if (i[0] not in seq):
            seq.append(i[0])
        if (i[1] not in seq):
            seq.append(i[1])

    return seq

def _consTSPCycleCover(depotID, tau):
    raise VrpSolverNotAvailableError("ERROR: vrpSolver has not implement this method yet")
    return seq

def _impTSP2Opts(nodeIDs, tau, initSeq, asymFlag):
    # Initialize ==============================================================
    improvedFlag = False
    impSeq = [i for i in initSeq]

    # Initialize accDist/accRevDist ===========================================
    # FIXME: accDist and accRevDist can be improved, but not necessary right now
    # Accumulated distance from depot
    d = 0
    accDist = []
    for i in range(len(impSeq) - 1):
        accDist.append(d)
        if (d != None and (impSeq[i], impSeq[i + 1]) in tau):
            d += tau[impSeq[i], impSeq[i + 1]]
        else:
            d = None
    accDist.append(d)
    oriOfv = accDist[-1]

    # Accumulated distance to depot (reversed seq)
    revD = 0
    accRevDist = []
    if (asymFlag):
        for i in range(len(impSeq) - 1):
            accRevDist.insert(0, revD)
            if (revD != None and (impSeq[len(impSeq) - i - 1], impSeq[len(impSeq) - i - 2]) in tau):
                revD += tau[impSeq[len(impSeq) - i - 1], impSeq[len(impSeq) - i - 2]]
            else:
                revD = None
    accRevDist.insert(0, revD)
    oriRevOfv = accRevDist[0]

    # Main iteration ==========================================================
    # Needs rewrite, when calculating dist, avoid repeated calculation
    if (len(impSeq) >= 4):
        # Try 2-opt
        can2OptFlag = True
        while (can2OptFlag):
            can2OptFlag = False
            # Try 2Opt
            for i in range(len(impSeq) - 2):
                for j in range(i + 2, len(impSeq) - 1):
                    opt = exchange2Arcs(
                        route = impSeq, 
                        tau = tau, 
                        i = i, 
                        j = j, 
                        accDist = accDist,
                        accRevDist = accRevDist,
                        asymFlag = asymFlag)
                    if (opt != None and opt['deltaCost'] + CONST_EPSILON < 0):
                        can2OptFlag = True
                        improvedFlag = True
                        impSeq = opt['route']
                        oriOfv = opt['newCost']
                        oriRevOfv = opt['newRevCost']

                        d = 0
                        accDist = []
                        for i in range(len(impSeq) - 1):
                            accDist.append(d)
                            if (d != None and (impSeq[i], impSeq[i + 1]) in tau):
                                d += tau[impSeq[i], impSeq[i + 1]]
                            else:
                                d = None
                        accDist.append(d)
                        oriOfv = accDist[-1]

                        revD = 0
                        accRevDist = []
                        if (asymFlag):
                            for i in range(len(impSeq) - 1):
                                accRevDist.insert(0, revD)
                                if (revD != None and (impSeq[len(impSeq) - i - 1], impSeq[len(impSeq) - i - 2]) in tau):
                                    revD += tau[impSeq[len(impSeq) - i - 1], impSeq[len(impSeq) - i - 2]]
                                else:
                                    revD = None
                        accRevDist.insert(0, revD)
                        oriRevOfv = accRevDist[0]
                        break
    return {
        'impSeq': impSeq,
        'improvedFlag': improvedFlag,
        'oriOfv': oriOfv,
        'oriRevOfv': oriRevOfv,
    }

def heuTSPEx(
    nodes: dict, 
    predefinedArcs: list[list[tuple[int|str]]] = [],
    edges: dict = {'method': "Euclidean", 'ratio': 1},
    algo: dict = {'cons': 'NearestNeighbor', 'impv': '2Opt'},
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    serviceTime: float = 0,
    ) -> dict|None:

    
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

    # Define tau ==============================================================
    tau, pathLoc = matrixDist(nodes, edges, depotID, nodeIDs, serviceTime)

    # Check symmetric =========================================================
    asymFlag = False
    for (i, j) in tau:
        if (tau[i, j] != tau[j, i]):
            asymFlag = True
            break

    mustGo = {}
    for arcs in predefinedArcs:
        for arc in arcs:
            if (arc[0] not in mustGo):
                mustGo[arc[0]] = arc[1]
            if (arc[1] not in mustGo):
                mustGo[arc[1]] = arc[0]


    if (algo == None or 'cons' not in algo or 'impv' not in algo):
        raise MissingParameterError("ERROR: Not supported (for now)")

    seq = []
    # Neighborhood based heuristic, including nearest neighborhood, k-nearest neighborhood, and furthest neighborhood
    if (algo['cons'] == 'NearestNeighbor'):
        seq = _consTSPExkNearestNeighbor(depotID, nodeIDs, tau, mustGo, 1)

    # Cleaning seq before local improving =================================
    ofv = calSeqCostMatrix(tau, seq, closeFlag = False)
    revOfv = None if not asymFlag else calSeqCostMatrix(tau, [seq[len(seq) - i - 1] for i in range(len(seq))], closeFlag = False)
    consOfv = ofv

    # Local improvment phase ==============================================
    # NOTE: For the local improvement, try every local search operator provided in a greedy way
    if ('impv' in algo and algo['impv'] != None and algo['impv'] != []):
        canImproveFlag = True
        while (canImproveFlag):
            canImproveFlag = False
            # Try 2Opts
            if (not canImproveFlag and (algo['impv'] == '2Opt' or '2Opt' in algo['impv'])):
                imp = _impTSPEx2Opts(nodeIDs, tau, seq, asymFlag, mustGo)
                if (imp['improvedFlag']):
                    canImproveFlag = True
                    seq = imp['impSeq']
                    ofv = imp['oriOfv']
                    revOfv = imp['oriRevOfv']
    return {
        'ofv': ofv,
        'consOfv': consOfv,
        'algo': algo,
        'seq': seq,
        'serviceTime': serviceTime
    }

def _consTSPExkNearestNeighbor(depotID, nodeIDs, tau, mustGo, k = 1):
    # Initialize ----------------------------------------------------------
    seq = [depotID]
    remain = [nodeIDs[i] for i in range(len(nodeIDs)) if nodeIDs[i] != depotID]
    # Accumulate seq ------------------------------------------------------
    while (len(remain) > 0):
        currentNodeID = seq[-1]

        if (currentNodeID in mustGo and mustGo[currentNodeID] not in seq):
            nextNodeID = mustGo[currentNodeID]
        else:
            # Sort the distance from current node to the rest of nodes
            sortedSeqHeap = []
            for n in remain:
                if ((currentNodeID, n) in tau):
                    dist = tau[currentNodeID, n]
                    heapq.heappush(sortedSeqHeap, (dist, n))

            # Get the kth of sorted node and append it to seq
            nextNodeID = None
            for _ in range(k):
                if (len(sortedSeqHeap) > 0):
                    nextNodeID = heapq.heappop(sortedSeqHeap)[1]
        seq.append(nextNodeID)

        # Update remain
        remain.remove(nextNodeID)
    seq.append(depotID)
    return seq

def _impTSPEx2Opts(nodeIDs, tau, initSeq, asymFlag, mustGo):
    # Initialize ==============================================================
    improvedFlag = False
    impSeq = [i for i in initSeq]

    # Initialize accDist/accRevDist ===========================================
    # FIXME: accDist and accRevDist can be improved, but not necessary right now
    # Accumulated distance from depot
    d = 0
    accDist = []
    for i in range(len(impSeq) - 1):
        accDist.append(d)
        if (d != None and (impSeq[i], impSeq[i + 1]) in tau):
            d += tau[impSeq[i], impSeq[i + 1]]
        else:
            d = None
    accDist.append(d)
    oriOfv = accDist[-1]

    # Accumulated distance to depot (reversed seq)
    revD = 0
    accRevDist = []
    if (asymFlag):
        for i in range(len(impSeq) - 1):
            accRevDist.insert(0, revD)
            if (revD != None and (impSeq[len(impSeq) - i - 1], impSeq[len(impSeq) - i - 2]) in tau):
                revD += tau[impSeq[len(impSeq) - i - 1], impSeq[len(impSeq) - i - 2]]
            else:
                revD = None
    accRevDist.insert(0, revD)
    oriRevOfv = accRevDist[0]

    # Main iteration ==========================================================
    # Needs rewrite, when calculating dist, avoid repeated calculation
    if (len(impSeq) >= 4):
        # Try 2-opt
        can2OptFlag = True
        while (can2OptFlag):
            can2OptFlag = False
            # Try 2Opt
            for i in range(len(impSeq) - 2):
                for j in range(i + 2, len(impSeq) - 1):
                    canTryFlag = ((impSeq[i] not in mustGo and impSeq[j] not in mustGo)
                        or (impSeq[i] in mustGo and mustGo[impSeq[i]] == impSeq[j])
                        or (impSeq[j] in mustGo and mustGo[impSeq[j]] == impSeq[i]))
                    if canTryFlag:
                        opt = exchange2Arcs(
                            route = impSeq, 
                            tau = tau, 
                            i = i, 
                            j = j, 
                            accDist = accDist,
                            accRevDist = accRevDist,
                            asymFlag = asymFlag)
                        if (opt != None and opt['deltaCost'] + CONST_EPSILON < 0):
                            can2OptFlag = True
                            improvedFlag = True
                            impSeq = opt['route']
                            oriOfv = opt['newCost']
                            oriRevOfv = opt['newRevCost']

                            d = 0
                            accDist = []
                            for i in range(len(impSeq) - 1):
                                accDist.append(d)
                                if (d != None and (impSeq[i], impSeq[i + 1]) in tau):
                                    d += tau[impSeq[i], impSeq[i + 1]]
                                else:
                                    d = None
                            accDist.append(d)
                            oriOfv = accDist[-1]

                            revD = 0
                            accRevDist = []
                            if (asymFlag):
                                for i in range(len(impSeq) - 1):
                                    accRevDist.insert(0, revD)
                                    if (revD != None and (impSeq[len(impSeq) - i - 1], impSeq[len(impSeq) - i - 2]) in tau):
                                        revD += tau[impSeq[len(impSeq) - i - 1], impSeq[len(impSeq) - i - 2]]
                                    else:
                                        revD = None
                            accRevDist.insert(0, revD)
                            oriRevOfv = accRevDist[0]
                            break
    return {
        'impSeq': impSeq,
        'improvedFlag': improvedFlag,
        'oriOfv': oriOfv,
        'oriRevOfv': oriRevOfv,
    }