import heapq
import math
import datetime
import warnings

from .const import *
from .common import *
# from .graph import *
from .geometry import *
from .msg import *
# from .operator import *
# from .calculate import *

# History =====================================================================
# 20230510 - Cleaning
# =============================================================================

def bkheuTSP(
    nodes: dict, 
    edges: dict = {'method': "Euclidean", 'ratio': 1}, 
    algo: dict = {'cons': 'Insertion', 'impv': '2opt'}, 
    depotID: int | str = 0, 
    nodeIDs: list[int | str] | str = 'All', 
    serviceTime: float = 0
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
    tau, _ = matrixDist(nodes, edges, depotID, nodeIDs, serviceTime)

    # Check symmetric =========================================================
    asymFlag = False
    for (i, j) in tau:
        if (tau[i, j] != tau[j, i]):
            asymFlag = True
            break

    # Construction heuristics =================================================
    # NOTE: Output format: seq = [depotID, xx, xx, xx, depotID]
    if (algo == None or ('cons' not in algo and 'impv' not in algo)):
        raise MissingParameterError(ERROR_MISSING_TSP_ALGO)

    seq = None
    # Constructive phase
    # An initial solution is given
    if ('cons' not in algo):
        if ('initSeq' not in algo):
            raise MissingParameterError("ERROR: Need 'initSeq' for local improvement")
        elif (len(algo['initSeq']) != len(nodeIDs) + 1):
            raise UnsupportedInputError("ERROR: Length of 'initSeq' is incorrect, check if the sequence starts and ends with `depotID`")
        else:
            notInNodeIDs = [v for v in algo['initSeq'] if v not in nodeIDs]
            if (len(notInNodeIDs) > 0):
                raise OutOfRangeError("ERROR: The following nodes in 'initSeq' is not in `nodeIDs`: %s" % list2String(notInNodeIDs))
        seq = [i for i in algo['initSeq']]

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
            seq = _bkconsTSPInsertion(nodeIDs, initSeq, tau)
        else:
            notInNodeIDs = [v for v in algo['initSeq'] if v not in nodeIDs]
            if (len(notInNodeIDs) > 0):
                raise OutOfRangeError("ERROR: The following nodes in 'initSeq' is not in `nodeIDs`: %s" % list2String(notInNodeIDs))
            else:
                seq = _bkconsTSPInsertion(nodeIDs, algo['initSeq'], tau)
    
     # Randomly create a sequence
    elif (algo['cons'] == 'Random'):
        seq = _bkconsTSPRandomSeq(depotID, nodeIDs, tau)
    
    else:
        raise UnsupportedInputError(ERROR_MISSING_TSP_ALGO)

    # Cleaning seq before local improving =================================
    ofv = bkcalSeqCostMatrix(tau, seq, closeFlag = False)
    revOfv = None if not asymFlag else bkcalSeqCostMatrix(tau, [seq[len(seq) - i - 1] for i in range(len(seq))], closeFlag = False)
    consOfv = ofv

    # Local improvment phase ==============================================
    # NOTE: For the local improvement, try every local search operator provided in a greedy way
    if ('impv' in algo and algo['impv'] != None and algo['impv'] != []):
        canImproveFlag = True
        while (canImproveFlag):
            canImproveFlag = False

            # Try 2Opts
            if (not canImproveFlag and (algo['impv'] == '2Opt' or '2Opt' in algo['impv'])):
                imp = _bkimpTSP2Opts(nodeIDs, tau, seq, asymFlag)
                if (imp['improvedFlag']):
                    canImproveFlag = True
                    seq = imp['impSeq']
                    ofv = imp['oriOfv']
                    revOfv = imp['oriRevOfv']

    return {
        'ofv': ofv,
        'consOfv': consOfv,
        'seq': seq,
        'serviceTime': serviceTime
    }

def _bkconsTSPkNearestNeighbor(depotID, nodeIDs, tau, k = 1):
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
    seq.append(depotID)
    return seq

def _bkconsTSPFarthestNeighbor(depotID, nodeIDs, tau):
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
    seq.append(depotID)
    return seq

def _bkconsTSPSweep(nodes, depotID, nodeIDs, tau):
    # Sweep seq -----------------------------------------------------------
    sweep = getSweepSeq(
        nodes = nodes, 
        nodeIDs = nodeIDs,
        centerLoc = nodes[depotID]['loc'])

    startIndex = 0
    seq = []
    for k in range(len(sweep)):
        if (sweep[k] == depotID):
            startIndex = k
    seq.extend([sweep[k] for k in range(startIndex, len(sweep))])
    seq.extend([sweep[k] for k in range(0, startIndex)])
    seq.append(0)

    return seq

def _bkconsTSPRandomSeq(depotID, nodeIDs, tau):
    # Get random seq ------------------------------------------------------
    seq = [i for i in nodeIDs if i != depotID]
    random.shuffle(seq)
    seq.insert(0, depotID)
    seq.append(depotID)
    return seq

def _bkconsTSPInsertion(nodeIDs, initSeq, tau):
    # Initialize ----------------------------------------------------------
    # NOTE: initSeq should starts and ends with depotID
    seq = [i for i in initSeq]
    insertDict = {}
    if (len(nodeIDs) < 1):
        return {
            'ofv': 0,
            'seq': None
        }
    unInserted = [i for i in nodeIDs if i not in seq]
    # Try insertion one by one --------------------------------------------
    while (len(unInserted) > 0):
        bestCus = None
        bestCost = None
        bestInsertionIndex = 0
        for cus in unInserted:
            for i in range(1, len(seq)):
                if ((seq[i - 1], cus, seq[i]) not in insertDict):
                    insertDict[(seq[i - 1], cus, seq[i])] = (
                        tau[seq[i - 1], cus] 
                        + tau[cus, seq[i]] 
                        - (tau[seq[i - 1], seq[i]] if seq[i - 1] != seq[i] else 0))
                cost = insertDict[(seq[i - 1], cus, seq[i])]
                if (bestCost == None or bestCost > cost):
                    bestCost = cost
                    bestCus = cus
                    bestInsertionIndex = i
        if (bestCost != None):
            seq.insert(bestInsertionIndex, bestCus)
            unInserted.remove(bestCus)
    return seq
    
def _bkconsTSPChristofides(depotID, weightArcs, matchingAlgo):
    # Create MST ----------------------------------------------------------
    mst = graphMST(
        weightArcs = weightArcs,
        exportAs = 'Arcs')['mst']
    mstAsTree = graphArcs2AdjList(mst['mst'])

    # Derive subgraph of odd degree vertices ------------------------------
    oddDegrees = []
    for node in mstAsTree:
        if (len(mstAsTree[node]) % 2 != 0):
            oddDegrees.append(node)
    subGraph = []
    for arc in weightArcs:
        if (arc[0] in oddDegrees and arc[1] in oddDegrees):
            subGraph.append(arc)

    # Find minimum cost matching of the subgraph --------------------------
    minMatching = graphMatching(
        weightArcs = subGraph, 
        algo = matchingAlgo,
        mType = 'Minimize')['matching']

    # Add them back to create a new graph ---------------------------------
    newGraph = []
    for arc in minMatching:
        newGraph.append(arc)
    for arc in mst:
        newGraph.append(arc)

    # Traverse graph and get seq ------------------------------------------
    seq = graphTraversal(
        arcs = newGraph, 
        oID = depotID)['seq']
    seq.append(depotID)

    return seq

def _bkimpTSP2Opts(nodeIDs, tau, initSeq, asymFlag):
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
                    opt = bkexchange2Arcs(
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

def bkcalSeqCostMatrix(
    tau:        "Dictionary {(nodeID1, nodeID2): dist, ...}", 
    seq:        "List, sequence of visiting node ids",
    closeFlag:  "True if the seq is closed" = None,
    reverseFlag: "True if calculates the cost of the reversed seq" = False,
    i:          "Start index" = 0,
    j:          "End index" = None
    ) -> "Return the cost on the graph given cost matrix/dictionary tau":
    # Accumulate costs ========================================================
    if (j == None):
        j = len(seq) - 1

    cost = 0
    for k in range(i, j):
        if (not reverseFlag):
            if ((seq[k], seq[k + 1]) in tau):            
                cost += tau[seq[k], seq[k + 1]]
            else:
                return None
        else:
            if ((seq[len(seq) - k - 1], seq[len(seq) - k - 2]) in tau):
                cost += tau[seq[len(seq) - k - 1], seq[len(seq) - k - 2]]
            else:
                return None

    if (i == 0 and j == len(seq) - 1 and closeFlag):
        if (not reverseFlag):
            if ((seq[-1], seq[0]) in tau):
                cost += tau[seq[-1], seq[0]]
            else:
                return None
        else:
            if ((seq[0], seq[-1]) in tau):
                cost += tau[seq[0], seq[-1]]
            else:
                return None
        
    return cost

def bkexchange2Arcs(
    route:      "A given sequence of vehicle route, start/end with the depot, assuming this route is feasible" = None, 
    tau:        "Traveling cost matrix" = None, 
    completeTauFlag: "True if tau is completed, to bypass calculation" = True,
    i:          "Index in the sequence, 0 <= i <= len(route) - 2" = None, 
    j:          "Index in the sequence, 0 <= i <= len(route) - 2, i.Next() != j" = None,
    accDist:    "Accumulated travel distance from depot" = None,
    accRevDist: "Accumulated travel distance for reveres route" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "2Opt, Reconnect arc between i, i+1 and j, j+1, if it is applicable":

    # Before: ... --> nI -> nINext --> ...abcd... --> nJ     -> nJNext --> ...
    # After:  ... --> nI -> nJ     --> ...dcba... --> nINext -> nJNext --> ...
    N = len(route)
    nI = route[i]
    nINext = route[i + 1]
    nJ = route[j]
    nJNext = route[j + 1]

    # Initialize accDist accRevDist ===========================================
    if (accDist == None):
        accDist = [0]
        for i in range(len(route) - 1):
            accDist.append(tau[route[i], route[i + 1]])
    if (accRevDist == None):
        accRevDist = [0]
        for i in range(len(route) - 1):
            accRevDist.append(tau[route[len(route) - 1 - i], route[len(route) - 2 - i]])

    # Check if can 2-opt
    can2OptFlag = True if (accDist[-1] != None) else False
    canRev2OptFlag = asymFlag and (True if (accRevDist[-1] != None) else False)

    # Early quit
    if (not asymFlag and not can2OptFlag):
        return None
    elif (asymFlag and not can2OptFlag and not canRev2OptFlag):
        return None

    # Calculate deltaCost
    newCost = None
    newRevCost = None
    reverseFlag = False

    cost = accDist[-1]
    revCost = accRevDist[0]

    # If not asymmetric, only one situation possible, so directly calculate deltaCost
    if (not asymFlag):
        newCost = (cost - (tau[nI, nINext] + tau[nJ, nJNext])
                        + (tau[nI, nJ] + tau[nINext, nJNext]))
        deltaCost = newCost - cost
        newRevCost = None

    # If asymmetric, the following cases need to be considered:
    # 1. newSeq feasible + newRevSeq feasible => deltaCost: min(newRevCost, newCost) - cost
    #    - newCost is better
    #    - newRevCost is better, need to take the reversed route as a new route
    # 2. newSeq infeasible + newRevSeq feasible => deltaCost: newRevCost - cost
    #    - need to take the reversed route as a new route
    # 3. newSeq feasible + newRevSeq infeasible => deltaCost: newCost - cost
    #    - same as symmetric
    else:
        costABCD = accDist[j] - accDist[i + 1]
        costDCBA = accRevDist[i + 1] - accRevDist[j]

        if (can2OptFlag):
            newCost = (cost - (tau[nI, nINext] + tau[nJ, nJNext] + costABCD)
                            + (tau[nI, nJ] + tau[nINext, nJNext] + costDCBA))
        if (canRev2OptFlag):
            newRevCost = (revCost - (tau[nINext, nI] + tau[nJNext, nJ] + costDCBA)
                                  + (tau[nJ, nI] + tau[nJNext, nINext] + costABCD))

        # Case 1
        if (can2OptFlag and canRev2OptFlag):
            if (newRevCost >= newCost):
                deltaCost = newCost - cost
            else:
                deltaCost = newRevCost - cost
                reverseFlag = True
                newCost, newRevCost = newRevCost, newCost
        # Case 2
        elif (not can2OptFlag and canRev2OptFlag):
            deltaCost = newRevCost - cost
            reverseFlag = True
            newCost, newRevCost = newRevCost, newCost
        # Case 3
        elif (can2OptFlag and not canRev2OptFlag):
            deltaCost = newCost - cost

    newSeq = []
    if (deltaCost < 0):
        newSeq.extend([route[k] for k in range(i + 1)])
        newSeq.extend([route[j - k] for k in range(j - i)])
        newSeq.extend([route[k] for k in range(j + 1, len(route))])
        if (reverseFlag):
            newSeq.reverse()
        return {
            'route': newSeq,
            'deltaCost': deltaCost,
            'reversed': reverseFlag,
            'newCost': newCost,
            'newRevCost': newRevCost
        }
    else:
        return None