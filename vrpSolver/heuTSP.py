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
from .plot import *

def heuTSP(
    nodes:      "The coordinate of given nodeID" = None, 
    edges:      "The traveling matrix" = None,
    algo:       "Algorithm used in the heuristic" = None,
    depotID:    "DepotID, default to be 0" = 0,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0
    ) -> "Use given heuristic methods to get TSP solution":

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
            ...     'ratio': 1
            ... }
            2) By given pairs of lat/lon
            >>> edge = {
            ...     'method': 'LatLon',
            ...     'unit': 'meters'
            ... }
            3) By a given dictionary
            >>> edge = {
            ...     'method': 'Dictionary',
            ...     'dictionary': dictionary
            ... }
            4) On the grids
            >>> edge = {
            ...     'method': 'Grid',
            ...     'grid': grid
            ... }
    algo: dictionary, required, default as {'cons': 'insert', 'impv': '2opt'}
        The algorithm configuration. Includes two phases, use 'cons' to specify constructive heuristic, and 'impv' to specify local improvement heurisitc::
            1) (default) Insertion
            >>> algo = {
            ...     'cons': 'Insertion',
            ...     'initSeq': initSeq, # An initial sequence, defalt [depotID]
            ...     'impv': '2Opt' # Options are: 'Reinsert', '2Opt'
            ... }
            2) Nearest neighborhood / k-nearest neighborhood
            >>> algo = {
            ...     'cons': 'NearestNeighbor',
            ...     'k': 1, # 1: nearest neighbor, 2 ~ K: k-nearest neighbor, -1: farthest neighbor 
            ... }
            3) Sweep
            >>> algo = {
            ...     'cons': 'Sweep',
            ...     'direction': 'clock-wise', # or 'counter-clock-wise', default 'clock-wise'
            ...     'startDeg': 0, # Starting direction, 0 as north/up, default 0
            ... }
            4) (working) Christofides
            >>> algo = {
            ...     'cons': 'Christofides',
            ... }
            5) (planning) Cycle cover, particular for Asymmetric TSP
            >>> algo = {
            ...     'cons': 'CycleCover',
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
        ... 
        ... }
    """

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
        else:
            msgError(ERROR_INCOR_NODEIDS)
            return

    # Define tau ==============================================================
    tau = getTau(nodes, edges, edgeArgs, depotID, nodeIDs, serviceTime)

    # Check symmetric =========================================================
    asymFlag = False
    for (i, j) in tau:
        if (tau[i, j] != tau[j, i]):
            asymFlag = True
            break

    # Construction heuristics =================================================
    # NOTE: These heuristics don't need to transform into arc representation
    # NOTE: Output format: [depotID, xx, xx, xx, depotID]
    seq = None
    if (consAlgo == None):
        if (consAlgoArgs == None or 'initSeq' not in consAlgoArgs):
            msgError("ERROR: Need initial TSP solution for local improvement")
            if (len(consAlgoArgs['initSeq']) != len(nodeIDs) + 1
                or consAlgoArgs['initSeq'][0] != depotID
                or consAlgoArgs['initSeq'][-1] != depotID):
                msgError("ERROR: Incorrect initial TSP, should start and end with depotID")
        seq = consAlgoArgs['initSeq']
    elif (consAlgo == 'NearestNeighbor'):
        # Nearest neighbor: k = 1
        seq = _consTSPkNearestNeighbor(depotID, nodeIDs, tau, 1)
    elif (consAlgo == 'k-NearestNeighbor'):
        if (consAlgoArgs == None):
            warnings.warn("Warning: Missing parameter k for k-NearestNeighbor, set to be default value as k = 1")
            consAlgoArgs = {'k': 1}
        seq = _consTSPkNearestNeighbor(depotID, nodeIDs, tau, consAlgoArgs['k'])
    elif (consAlgo == 'FarthestNeighbor'):
        # NOTE: "Worst" initial solution
        seq = _consTSPFarthestNeighbor(depotID, nodeIDs, tau)
    elif (consAlgo == 'Insertion'):
        if (consAlgoArgs == None or 'initSeq' not in consAlgoArgs):
            # If missing initial sequence, start with "depot -> farthest -> depot"
            farthestDist = None
            farthestID = None
            for n in nodeIDs:
                if ((n, depotID) in tau):
                    if (farthestDist == None or tau[n, depotID] > farthestDist):
                        farthestID = n
                        farthestDist = tau[n, depotID]
                elif ((depotID, n) in tau):
                    if (farthestDist == None or tau[depotID, n] > farthestDist):
                        farthestID = n
                        farthestDist = tau[depotID, n]
            consAlgoArgs = {'initSeq': [depotID, farthestID, depotID]}
            warnings.warn("WARNING: 'initSeq' was not provided, initial with [depotID, farthest, depotID]")
        seq = _consTSPInsertion(nodeIDs, consAlgoArgs['initSeq'], tau)
    elif (consAlgo == 'Sweep'):
        seq = _consTSPSweep(nodes, depotID, nodeIDs, tau)
    elif (consAlgo == 'Random'):
        seq = _consTSPRandomSeq(depotID, nodeIDs, tau)
    else:
        pass

    weightArcs = []
    # Create arcs =============================================================
    if (seq == None):
        for (i, j) in tau:
            if (i != None and j != None and i < j):
                weightArcs.append((i, j, tau[i, j]))

    # Construction Heuristics for TSP =========================================
    if (seq == None and not asymFlag):
        if (consAlgo == 'DepthFirst'):
            seq = _consTSPDepthFirst(weightArcs)
        elif (consAlgo == 'Christofides'):
            if (consAlgoArgs == None or 'matchingAlgo' not in consAlgoArgs):
                consAlgoArgs = {'matchingAlgo': 'IP'}
            seq = _consTSPChristofides(weightArcs, consAlgoArgs['matchingAlgo'])

    # Confirm construction results ============================================
    if (seq == None):
        msgError("ERROR: Incorrect construction algorithm")
        return

    # Local improvement ======================================================= 
    canImproveFlag = False
    
    # NOTE: seq should starts and ends with depotID
    ofv = calSeqCostMatrix(tau, seq, closeFlag = False)
    revOfv = None if not asymFlag else calSeqCostMatrix(tau, [seq[len(seq) - i - 1] for i in range(len(seq))], closeFlag = False)
    consOfv = ofv

    if (impAlgo != None):
        canImproveFlag = True
        if (type(impAlgo) == str):
            impAlgo = [impAlgo]        

    while (canImproveFlag):
        canImproveFlag = False

        # Try 2Opts
        if (not canImproveFlag and '2Opt' in impAlgo):
            imp = _impTSP2Opts(nodeIDs, tau, seq, asymFlag)
            if (imp['improvedFlag']):
                canImproveFlag = True
                seq = imp['impSeq']
                ofv = imp['oriOfv']
                revOfv = imp['oriRevOfv']

        # Try reinsert
        if (not canImproveFlag and 'Reinsert' in impAlgo):
            imp = _impTSPReinsert(nodeIDs, tau, seq, ofv, revOfv, asymFlag)
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

def _consTSPkNearestNeighbor(depotID, nodeIDs, tau, k = 1):
    # Initialize ----------------------------------------------------------
    seq = [depotID]
    remain = [nodeIDs[i] for i in range(len(nodeIDs)) if nodeIDs[i] != depotID]
    # Accumulate seq ------------------------------------------------------
    while (len(remain) > 0):
        currentNodeID = seq[-1]

        # Sort the distance from current node to the rest of nodes
        sortedSeq = []
        sortedSeqHeap = []
        for n in remain:
            if ((currentNodeID, n) in tau):
                dist = tau[currentNodeID, n]
                heapq.heappush(sortedSeqHeap, (dist, n))

        # Get the kth of sorted node and append it to seq
        nextNodeID = None
        for i in range(k):
            if (len(sortedSeqHeap) > 0):
                nextNodeID = heapq.heappop(sortedSeqHeap)[1]
        seq.append(nextNodeID)

        # Update remain
        remain.remove(nextNodeID)
    seq.append(depotID)
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
    seq.append(depotID)
    return seq

def _consTSPSweep(nodes, depotID, nodeIDs, tau):
    # Sweep seq -----------------------------------------------------------
    sweep = getSweepSeq(
        nodes = nodes, 
        nodeIDs = nodeIDs)

    startIndex = None
    seq = []
    for k in range(len(sweep)):
        if (sweep[k] == depotID):
            startIndex = k
    seq.extend([sweep[k] for k in range(startIndex, len(sweep))])
    seq.extend([sweep[k] for k in range(0, startIndex)])
    seq.append(0)

    return seq

def _consTSPRandomSeq(depotID, nodeIDs, tau):
    # Get random seq ------------------------------------------------------
    seq = [i for i in nodeIDs if i != depotID]
    random.shuffle(seq)
    seq.insert(0, depotID)
    seq.append(depotID)
    return seq

def _consTSPInsertion(nodeIDs, initSeq, tau):
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
        bestInsertionIndex = None
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

def _consTSPDepthFirst(depotID, weightArcs):
    # Create MST ----------------------------------------------------------
    mst = graphMST(
        weightArcs = weightArcs,
        exportAs = 'Arcs')['mst']

    # Seq of visit is the seq of Depth first search on the MST ------------
    seq = graphTraversal(
        arcs = mst,
        oID = depotID)['seq']
    seq.append(depotID)
    return seq

def _consTSPChristofides(depotID, weightArcs, matchingAlgo):
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

def _impTSPReinsert(nodeIDs, tau, initSeq, oriOfv, oriRevOfv, asymFlag):
    # Initialize ==============================================================
    improvedFlag = False
    impSeq = [i for i in initSeq]
    
    # Main iteration ==========================================================
    # Needs rewrite, when calculating dist, avoid repeated calculation
    if (len(impSeq) >= 4):
        # Try reinsert
        canReinsertFlag = True
        while (canReinsertFlag):
            canReinsertFlag = False
            # Does not consider the "reinsert" of the depot
            for i in range(1, len(impSeq) - 1):
                # First remove
                nI = impSeq[i]
                removed = calRemovalSaving(
                    route = impSeq,
                    tau = tau,
                    nI = nI,
                    cost = oriOfv,
                    revCost = oriRevOfv,
                    asymFlag = asymFlag)
                if (removed != None):
                    inserted = calInsertionCost(
                        route = removed['newSeq'],
                        tau = tau,
                        nJ = nI,
                        cost = removed['newCost'],
                        revCost = removed['newRevCost'],
                        asymFlag = asymFlag)
                    if (inserted != None and inserted['newCost'] + CONST_EPSILON < oriOfv):
                        canReinsertFlag = True
                        improvedFlag = True
                        impSeq = inserted['newSeq']
                        oriOfv = inserted['newCost']
                        oriRevOfv = inserted['newRevCost']
    return {
        'impSeq': impSeq,
        'improvedFlag': improvedFlag,
        'oriOfv': oriOfv,
        'oriRevOfv': oriRevOfv
    }
