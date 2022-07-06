import heapq
import math

from .const import *
from .common import *
from .graph import *
from .geometry import *
from .msg import *

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
                 4) String 'Grid', will need to add arguments using `edgeArgs`"= "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                    {\
                        'colRow': (numCol, numRow),\
                        'barriers': [(coordX, coordY), ...], \
                    }" = None,
    depotID:    "DepotID, default to be 0" = 0,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    consAlgo:   "1) String 'k-NearestNeighbor' or \
                 2) String 'FarthestNeighbor' or \
                 3) String (default) 'Insertion' or \
                 4) String 'Sweep' or \
                 5) String 'DepthFirst' or \
                 6) String 'Christofides' or \
                 7) String (not available) 'CycleCover', particular for ATSP, also work for TSP or \
                 8) String 'Random'" = 'Insertion',
    consAlgoArgs: "Dictionary, args for constructive heuristic \
                 1) None for unspecified `algo` options, or \
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
                        'initSeq': initial TSP route\
                    }"= None,
    locImpFlag: "True if call local improvement to improve solution quality" = True,
    impAlgo:    "1) String (not available) 'LKH' or \
                 2) String (default) '2Opt' = '2Opt'" = '2Opt',
    impAlgoArgs: "Dictionary, args for improvement heuristic" = None,
    ) -> "Use given heuristic methods to get TSP solution":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)
        else:
            print(ERROR_INCOR_NODEIDS)
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

    # Constructive heuristic (with out specifying depot) ======================
    # NOTE: So for, the seq does not necessarily start from depot, and the last digit was not the depot
    seq = _consTSP(
        nodes = nodes,
        tau = tau,
        nodeIDs = nodeIDs,
        algo = consAlgo,
        algoArgs = consAlgoArgs)

    # Local improvement heuristic (with out specifying depot) =================
    # NOTE: The last digit is NOT the duplicate of the depot
    if (locImpFlag):
        seq = _impTSP(
            nodes = nodes, 
            tau = tau, 
            initSeq = seq, 
            algo = impAlgo)

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

def _consTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    tau:      "Dictionary {(nodeID1, nodeID2): dist, ...}" = None,
    nodeIDs:    "A list of node IDs, notice that we don't have duplicated depot ID here" = None,
    algo:       "1) String 'k-NearestNeighbor' or \
                 2) String 'FarthestNeighbor' or \
                 3) String (default) 'Insertion' or \
                 4) String (not available) 'Patching' or \
                 5) String 'Sweep' or \
                 6) String 'DepthFirst' or \
                 7) String 'Christofides' or \
                 8) String (not available) 'CycleCover', particular for ATSP, also work for TSP or \
                 9) String 'Random'" = 'Insertion',
    algoArgs:   "Dictionary, needed for some algorithm options, \
                 1) None for unspecified `algo` options, or \
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
                        'initSeq': initial TSP route\
                    }\
                 5) for 'CycleCover' \
                    {\
                        'VCC': 'LP' or (not available) 'Hungarian'\
                    }"= None
    ) -> "Constructive heuristic solution for TSP, the output sequence does not necessarily start from depot":

    # Subroutines for different constructive heuristics =======================
    def _consTSPkNearestNeighbor(nodes, nodeIDs, tau, k = 1):
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
            t = seq[i]
            seq[i] = seq[j]
            seq[j] = t
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

    # Heuristics that don't need to transform arc representation ==============
    seq = None
    if (algo == 'k-NearestNeighbor'):
        if (algoArgs == None):
            algoArgs = {'k': 1}
        seq = _consTSPkNearestNeighbor(nodes, nodeIDs, tau, algoArgs['k'])
    elif (algo == 'FarthestNeighbor'):
        seq = _consTSPFarthestNeighbor(nodeIDs, tau)
    elif (algo == 'Insertion'):
        if (algoArgs == None):
            algoArgs = {'initSeq': [nodeIDs[0], nodeIDs[1]]}
        seq = _consTSPInsertion(nodeIDs, algoArgs['initSeq'], tau)
    elif (algo == 'Sweep'):
        seq = _consTSPSweep(nodes, nodeIDs, tau)
    elif (algo == 'Random'):
        seq = _consTSPRandomSeq(nodeIDs, tau)
    else:
        pass
    if (seq != None):
        return seq

    # Create arcs =============================================================
    # FIXME! indexes of nodes starts from 0 here!
    # FIXME! Add mapping between 0 to n and nodeIDs
    weightArcs = []
    for (i, j) in tau:
        if (i != None and j != None and i < j):
            weightArcs.append((i, j, tau[i, j]))

    # Subroutines for different constructive heuristics =======================
    def _consTSPDepthFirst(weightArcs):
        # Create MST ----------------------------------------------------------
        mst = graphMST(weightArcs)['mst']

        # Seq of visit is the seq of Depth first search on the MST ------------
        seq = traversalGraph(mst)['seq']
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
        minMatching = graphMatching(
            weightArcs=subGraph, 
            mType='Minimum', 
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
        seq = traversalGraph(newGraph, oID=oID)['seq']
        seq.append(seq[0])

        return seq

    # Constructive Heuristics for TSP =========================================
    seq = None
    if (algo == 'DepthFirst'):
        seq = _consTSPDepthFirst(weightArcs)
    elif (algo == 'Christofides'):
        if (algoArgs == None):
            algoArgs = {'matchingAlgo': 'IP'}
        seq = _consTSPChristofides(weightArcs, algoArgs['matchingAlgo'])

    return seq

def _impTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    tau:      "Dictionary {(nodeID1, nodeID2): dist, ...}" = None,
    initSeq:    "Sequence of visiting nodes, as the initial route, generated by constructive heuristic" = [],
    algo:       "1) String (not available) 'LK' or \
                 2) String (not available) 'k-Opt' or\
                 3) String (default) '2Opt', better for symmetric TSP, not as good for ATSP\
                 4) String 'Or', both symmetric and asymmetric" = '2OptSym',
    algoArgs:   "Dictionary that provides the arguments for the algorithms\
                 1) For `k-Opt`\
                 2) For `Or`,\
                    {\
                       's': 3\
                    }" = None,
    ) -> "Local improvement heuristic solution for TSP":

    # Subroutines for different local improvement heuristics ==================
    def _impTSP2Opt(nodeIDs, tau, initSeq):
        # Check symmetric -----------------------------------------------------
        symFlag = True
        for (i, j) in tau:
            if (tau[i, j] != tau[j, i]):
                symFlag = False
                break

        # Initialize ----------------------------------------------------------
        canImproveFlag = True
        impSeq = [i for i in initSeq]

        # Revert part of the sequences ----------------------------------------
        def revert(i, j, seq):
            rSeq = []
            if (j == i + 1 or j == i + 2):
                rSeq = [i for i in seq]
                rSeq[i], rSeq[j] = rSeq[j], rSeq[i]
            elif (j > i + 2):
                rSeq.extend([seq[k] for k in range(i)])
                rSeq.extend([seq[j - k] for k in range(j - i + 1)])
                rSeq.extend([seq[k] for k in range(j + 1, len(seq))])
            return rSeq

        # Main iteration ------------------------------------------------------
        if (len(impSeq) >= 4):
            while (canImproveFlag):
                canImproveFlag = False

                # Try 2-opt
                # First arc: (i, i + 1)
                # Second arc: (j, j + 1)
                # New arcs: (i, j + 1), (j, i + 1)
                oriOfv = None
                if (not symFlag):
                    oriOfv = calSeqCostMatrix(tau, impSeq, closeFlag = True)

                for i in range(len(impSeq) - 3):
                    for j in range(i + 2, len(impSeq) - 1):
                        # Saving
                        saving = 0
                        if ((impSeq[i], impSeq[j]) in tau and (impSeq[i + 1], impSeq[j + 1]) in tau):
                        # For symmetric TSP, no need to recalculate
                            if (symFlag):
                                saving = (tau[impSeq[i], impSeq[i + 1]] + tau[impSeq[j], impSeq[j + 1]]) - (tau[impSeq[i], impSeq[j]] + tau[impSeq[i + 1], impSeq[j + 1]])
                                if (saving > CONST_EPSILON):
                                    impSeq = revert(i + 1, j, impSeq)
                                    canImproveFlag = True
                            else:
                                newSeq = revert(i + 1, j, impSeq)
                                newOfv = calSeqCostMatrix(tau, newSeq, closeFlag = True)
                                if (newOfv < oriOfv):
                                    impSeq = [i for i in newSeq]
                                    canImproveFlag = True
        return impSeq

    # Heuristics that don't need to transform arc representation ==============
    seq = None
    nodeIDs = list(nodes.keys())
    if (algo == '2Opt'):
        seq = _impTSP2Opt(nodeIDs, tau, initSeq)

    return seq
