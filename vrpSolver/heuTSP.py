import heapq
import math

from .const import *
from .node import *
from .common import *
from .graph import *
from .geometry import *
from .msg import *
from .operator import *

def heuTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...} or \
                 4) String 'Grid', will need to add arguments using `edgeArgs`" = "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                    {\
                        'colRow': (numCol, numRow),\
                        'barriers': [(coordX, coordY), ...], \
                    }" = None,
    depotID:    "DepotID, default to be 0" = 0,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    consAlgo:   "1) None, if TSP initial solution is given or \
                 2) String 'NearestNeighbor' (will be regarded as k-NearestNeighbor and k = 1) or \
                 2) String 'k-NearestNeighbor' or \
                 3) String 'FarthestNeighbor' or \
                 4) String (default) 'Insertion' or \
                 5) String 'Sweep' or \
                 6) String 'DepthFirst' or \
                 7) String 'Christofides' or \
                 8) String (not available) 'CycleCover', particular for ATSP, also work for TSP or \
                 9) String 'Random'" = 'Insertion',
    consAlgoArgs: "Dictionary, args for constructive heuristic \
                 1) None for unspecified `algo` options, or \
                 2) for None, which an initial route should be given\
                    {\
                        'initSeq': complete initial TSP route\
                    }\
                 2) for 'k-NearestNeighbor' \
                    {\
                        'k': k-th nearest neighbor\
                    } \
                 3) for 'Christofides' \
                    {\
                        'matchingAlgo': algorithm for finding minimum matching\
                    } \
                 4) for 'Insertion' \
                    {\
                        'initSeq': complete/incomplete initial TSP route\
                    }"= None,
    impAlgo:    "1) List of Strings, options are as follows \
                 2) String (not available) 'LKH' or \
                 3) String (default) '2Opt'" = '2Opt',
    impAlgoArgs: "Dictionary, args for improvement heuristic" = None
    ) -> "Use given heuristic methods to get TSP solution":

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

    # Heuristics that don't need to transform arc representation ==============
    seq = None
    if (consAlgo == None):
        if (consAlgoArgs == None or 'initSeq' not in consAlgoArgs):
            msgError("Need initial TSP solution for local improvement")
        seq = consAlgoArgs['initSeq'][:-1]
    elif (consAlgo == 'NearestNeighbor'):
        seq = _consTSPkNearestNeighbor(nodes, nodeIDs, tau, 1)
    elif (consAlgo == 'k-NearestNeighbor'):
        if (consAlgoArgs == None):
            consAlgoArgs = {'k': 1}
        seq = _consTSPkNearestNeighbor(nodes, nodeIDs, tau, consAlgoArgs['k'])
    elif (consAlgo == 'FarthestNeighbor'):
        seq = _consTSPFarthestNeighbor(nodeIDs, tau)
    elif (consAlgo == 'Insertion'):
        if (consAlgoArgs == None or 'initSeq' not in consAlgoArgs):
            consAlgoArgs = {'initSeq': [nodeIDs[0], nodeIDs[1]]}
        seq = _consTSPInsertion(nodeIDs, consAlgoArgs['initSeq'], tau)
    elif (consAlgo == 'Sweep'):
        seq = _consTSPSweep(nodes, nodeIDs, tau)
    elif (consAlgo == 'Random'):
        seq = _consTSPRandomSeq(nodeIDs, tau)
    else:
        pass

    weightArcs = []
    # Create arcs =============================================================
    if (seq == None):    
        for (i, j) in tau:
            if (i != None and j != None and i < j):
                weightArcs.append((i, j, tau[i, j]))

    # Constructive Heuristics for TSP =========================================
    if (seq == None):
        if (consAlgo == 'DepthFirst'):
            seq = _consTSPDepthFirst(weightArcs)
        elif (consAlgo == 'Christofides'):
            if (consAlgoArgs == None or 'matchingAlgo' not in consAlgoArgs):
                consAlgoArgs = {'matchingAlgo': 'IP'}
            seq = _consTSPChristofides(weightArcs, consAlgoArgs['matchingAlgo'])

    # Heuristics that don't need to transform arc representation ==============
    if (impAlgo == '2Opt'):
        seq = _impTSPOpts(nodeIDs, tau, seq, asymFlag)

    # Fix the sequence to make it start from and end with the depot ===========
    # NOTE: nodeID gets duplicated, if nodeID == 0, the sequence starts and ends with a 0
    startIndex = None
    truckSeq = []
    for k in range(len(seq)):
        if (seq[k] == depotID):
            startIndex = k
    if (startIndex <= len(seq) - 1):
        for k in range(startIndex, len(seq)):
            truckSeq.append(seq[k])
    if (startIndex >= 0):
        for k in range(0, startIndex):
            truckSeq.append(seq[k])
    truckSeq.append(depotID)

    ofv = calSeqCostMatrix(tau, truckSeq)

    return {
        'ofv': ofv,
        'seq': truckSeq,
        'serviceTime': serviceTime
    }

def _consTSPkNearestNeighbor(nodes, nodeIDs, tau, k = 1):
    # Initialize ----------------------------------------------------------
    seq = [nodeIDs[0]]
    remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
    # Accumulate seq ------------------------------------------------------
    while (len(remain) > 0):
        currentNodeID = seq[-1]
        sortedNodes = getSortedNodesByDist(
            nodes = nodes,
            edges = tau,
            nodeIDs = remain,
            refNodeID = currentNodeID)
        if (k > len(remain)):
            seq.append(sortedNodes[-1])
        else:
            seq.append(sortedNodes[k - 1])
        remain = [i for i in nodeIDs if i not in seq]
    return seq

def _consTSPFarthestNeighbor(nodeIDs, tau):
    # Initialize ----------------------------------------------------------
    seq = [nodeIDs[0]]
    remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
    # Accumulate seq ------------------------------------------------------
    while (len(remain) > 0):
        nextLeng = None
        nextID = None
        for node in remain:
            if ((node, seq[-1]) in tau):
                if (nextLeng == None or tau[node, seq[-1]] > nextLeng):
                    nextID = node
                    nextLeng = tau[node, seq[-1]]
            elif ((seq[-1], node) in tau):
                if (nextLeng == None or tau[seq[-1], node] > nextLeng):
                    nextID = node
                    nextLeng = tau[seq[-1], node]
        seq.append(nextID)
        remain.remove(nextID)
    return seq

def _consTSPSweep(nodes, nodeIDs, tau):
    # Sweep seq -----------------------------------------------------------
    sweepSeq = getSweepSeq(
        nodes = nodes, 
        nodeIDs = nodeIDs)
    return sweepSeq

def _consTSPRandomSeq(nodeIDs, tau):
    # Get random seq ------------------------------------------------------
    seq = [i for i in nodeIDs]
    N = len(nodeIDs)
    for i in range(N):
        j = random.randint(0, N - 1)
        seq[i], seq[j] = seq[j], seq[i]
    return seq

def _consTSPInsertion(nodeIDs, initSeq, tau):
    # Initialize ----------------------------------------------------------
    seq = None
    if (initSeq == None):
        seq = [nodeIDs[0], nodeIDs[1]]
    else:
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
                if (list2Tuple([seq[i - 1], cus, seq[i]]) not in insertDict):
                    insertDict[list2Tuple([seq[i - 1], cus, seq[i]])] = (
                        tau[seq[i - 1], cus] 
                        + tau[cus, seq[i]] 
                        - (tau[seq[i - 1], seq[i]] if seq[i - 1] != seq[i] else 0))
                cost = insertDict[list2Tuple([seq[i - 1], cus, seq[i]])]
                if (bestCost == None or bestCost > cost):
                    bestCost = cost
                    bestCus = cus
                    bestInsertionIndex = i
        if (bestCost != None):
            seq.insert(bestInsertionIndex, bestCus)
            unInserted.remove(bestCus)
    return seq

def _consTSPDepthFirst(weightArcs):
    # Create MST ----------------------------------------------------------
    mst = graphMST(weightArcs)['mst']

    # Seq of visit is the seq of Depth first search on the MST ------------
    seq = graphTraversal(mst)['seq']
    seq.append(seq[0])
    return seq

def _consTSPChristofides(weightArcs, matchingAlgo):
    # Create MST ----------------------------------------------------------
    mst = graphMST(weightArcs)['mst']

    # Derive subgraph of odd degree vertices ------------------------------
    neighbors = arcs2AdjList(mst)
    oddDegrees = []
    for node in neighbors:
        if (len(neighbors[node]) % 2 != 0):
            oddDegrees.append(node)
    subGraph = []
    for arc in weightArcs:
        if (arc[0] in oddDegrees and arc[1] in oddDegrees):
            subGraph.append(arc)

    # Find minimum cost matching of the subgraph --------------------------
    minMatching = graphMinMatching(
        weightArcs=subGraph, 
        algo=matchingAlgo)['matching']

    # Add them back to create a new graph ---------------------------------
    newGraph = []
    for arc in minMatching:
        newGraph.append(arc)
    for arc in mst:
        newGraph.append(arc)

    # Traverse graph and get seq ------------------------------------------
    # Try to find a vertex with degree 1
    oID = None
    for node in neighbors:
        if (len(neighbors[node]) == 1):
            oID = node
            break
    seq = graphTraversal(newGraph, oID=oID)['seq']
    seq.append(seq[0])

    return seq

def _impTSPOpts(nodeIDs, tau, initSeq, asymFlag):
    # Initialize ----------------------------------------------------------
    canImproveFlag = True
    impSeq = [i for i in initSeq]
    oriOfv = calSeqCostMatrix(tau, impSeq, closeFlag = True)
    oriRevOfv = None
    if (asymFlag):
        oriRevOfv = calSeqCostMatrix(tau, [impSeq[len(impSeq) - i - 1] for i in range(len(impSeq))], closeFlag = True)
    
    # Main iteration ------------------------------------------------------
    # Needs rewrite, when calculating dist, avoid repeated calculation
    if (len(impSeq) >= 4):
        while (canImproveFlag):
            canImproveFlag = False

            # Try node exchange
            canNodeExchangeFlag = True
            while (canNodeExchangeFlag):           
                canNodeExchangeFlag = False        
                for i in range(len(impSeq) - 2):
                    for j in range(i + 1, len(impSeq) - 1):
                        # Saving
                        opt = exchange2Nodes(
                            seq = impSeq,
                            tau = tau,
                            i = i,
                            j = j,
                            cost = oriOfv,
                            revCost = oriRevOfv,
                            asymFlag = asymFlag)
                        if (opt != None and opt['deltaCost'] + CONST_EPSILON < 0):
                            print("2Nodes: [%s, %s]" % (i, j), opt['deltaCost'], oriOfv, opt['newCost'])
                            canNodeExchangeFlag = True
                            canImproveFlag = True
                            impSeq = opt['seq']
                            oriOfv = opt['newCost']
                            oriRevOfv = opt['newRevCost']
                            break
                    if (canNodeExchangeFlag):
                        break

            # Try 2-opt
            can2OptFlag = True
            while (can2OptFlag):
                can2OptFlag = False                
                for i in range(len(impSeq) - 2):
                    for j in range(i + 2, len(impSeq)):
                        # Saving
                        opt = exchange2Arcs(
                            seq = impSeq, 
                            tau = tau, 
                            i = i, 
                            j = j, 
                            cost = oriOfv, 
                            revCost = oriRevOfv,
                            asymFlag = asymFlag)
                        if (opt != None and opt['deltaCost'] + CONST_EPSILON < 0):
                            print("2Opt: [%s, %s]" % (i, j), opt['deltaCost'], oriOfv, opt['newCost'])
                            can2OptFlag = True
                            canImproveFlag = True
                            impSeq = opt['seq']
                            oriOfv = opt['newCost']
                            oriRevOfv = opt['newRevCost']
                            break
                    if (can2OptFlag):
                        break

            # Try reinsert
            canReinsertFlag = True
            while (canReinsertFlag):    
                canReinsertFlag = False               
                for i in range(1, len(impSeq) - 1):
                    # First remove
                    nI = impSeq[i]
                    removed = calRemovalSaving(
                        route = impSeq,
                        tau = tau,
                        i = i,
                        cost = oriOfv,
                        revCost = oriRevOfv,
                        asymFlag = asymFlag)
                    if (removed != None):
                        removedSeq = removed['newSeq']
                        newRemovedOfv = removed['newCost']
                        newRevRemovedOfv = removed['newRevCost']
                        inserted = calInsertionCost(
                            route = removedSeq,
                            tau = tau,
                            nJ = nI,
                            cost = newRemovedOfv,
                            revCost = newRevRemovedOfv,
                            asymFlag = asymFlag)
                        if (inserted != None and inserted['newCost'] + CONST_EPSILON < oriOfv):
                            print("Reinsert: %s" % nI, oriOfv, inserted['newCost'])
                            canReinsertFlag = True
                            canImproveFlag = True
                            impSeq = inserted['newSeq']
                            oriOfv = inserted['newCost']
                            oriRevOfv = inserted['newRevCost']
                            break

    return impSeq
