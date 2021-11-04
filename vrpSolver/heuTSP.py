import heapq
import math

from .const import *
from .common import *
from .graph import *
from .geometry import *

def lrTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    subgradM:   "Double" = 1,
    subgradRho: "Double, (0, 1)" = 0.95,
    stopType:   "1) String, (default) 'Epsilon' (`stopEpsilon` will be used) or \
                 2) String, 'IterationNum' (`stopK` will be used) or \
                 3) String, 'Runtime' (`stopTime` will be used)" = 'Epsilon',
    stopEpsilon:"Double, small number" = 0.01,
    stopK:      "Integer, large number" = 200,
    stopTime:   "Double, in seconds" = 600
    ) -> "Returns a Held & Karp lower bound of the TSP using Lagrangian Relaxation":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'LatLon'):
            edges = getTauLatLon(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Initialize ==============================================================
    k = 0
    u = [0 for i in range(len(nodeIDs))]
    d = None
    costSum = None
    L = None
    oldL = None

    # Calculate 1 tree ========================================================
    def cal1Tree(weightArcs):
        # Separate first node
        arcsWithVertexOne = []
        arcsWithoutVertexOne = []
        for i in range(len(weightArcs)):
            if (weightArcs[i][0] == 0 or weightArcs[i][1] == 0):
                arcsWithVertexOne.append(weightArcs[i])
            else:
                arcsWithoutVertexOne.append(weightArcs[i])

        # MST for the rest of vertices
        mst = graphMST(arcsWithoutVertexOne)['mst']

        # Find two cheapest arcs to vertex one
        sortedArcswithVertexOne = []
        for i in range(len(arcsWithVertexOne)):
            heapq.heappush(sortedArcswithVertexOne, (arcsWithVertexOne[i][2], arcsWithVertexOne[i]))

        # Build 1-tree
        leastTwo = []
        leastTwo.append(heapq.heappop(sortedArcswithVertexOne))
        leastTwo.append(heapq.heappop(sortedArcswithVertexOne))

        m1t = [i for i in mst]
        m1t.append(leastTwo[0][1])
        m1t.append(leastTwo[1][1])

        # Calculate total cost
        costSum = 0
        for i in range(len(m1t)):
            costSum += m1t[i][2]

        # Arcs to neighbors
        neighbors = arcs2AdjList(m1t)
        d = []
        for i in range(len(nodeIDs)):
            d.append(2 - len(neighbors[i]))

        return {
            'costSum': costSum,
            'm1t': m1t,
            'd': d
        }

    # Main iteration ==========================================================
    continueFlag = True
    while (continueFlag):
        # Update cost of each edge
        weightArcs = []
        for i in range(len(nodeIDs)):
            for j in range(len(nodeIDs)):
                if (i != None and j != None and i < j):
                    weightArcs.append((i, j, edges[i, j] - u[i] - u[j]))

        # Calculate 1-tree
        oneTree = cal1Tree(weightArcs)

        # Update L and d
        costSum = oneTree['costSum']
        m1t = oneTree['m1t']
        uSum = sum(u)
        if (L != None):
            oldL = L
        L = costSum + 2 * uSum
        d = oneTree['d']

        # update u
        oldU = [i for i in u]
        u = []
        eff = subgradM * math.pow(subgradRho, k)
        for i in range(len(nodeIDs)):
            u.append(oldU[i] + eff * d[i])

        # Check if continue
        def allZero(d):
            for i in d:
                if (i != 0):
                    return False
            return True
        if (k >= stopK):
            continueFlag = False
        elif (oldL != None and abs(oldL - L) < stopEpsilon):
            continueFlag = False
        elif (allZero(d)):
            continueFlag = False
        else:
            k += 1

    return {
        'lrLowerBound': costSum
    }

def heuTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    consAlgo:   "1) String 'NearestNeighbor' or \
                 2) String 'FarthestNeighbor' or \
                 3) String (not available) 'Insertion' or \
                 4) String (not available) 'Patching' or \
                 5) String (not available) 'Sweep' or \
                 6) String 'DepthFirst' or \
                 7) String (default) 'Christofides' or \
                 8) String 'Random'" = 'Christofides',
    consAlgoArgs: "Dictionary, args for constructive heuristic" = None,
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

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'LatLon'):
            edges = getTauLatLon(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Constructive heuristic ==================================================
    cons = consTSP(nodes, edges, nodeIDs, consAlgo, consAlgoArgs)

    # Local improvement heuristic =============================================
    tsp = impTSP(nodes, edges, nodeIDs, cons['seq'], impAlgo, impAlgoArgs)

    return tsp

def consTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    algo:       "1) String 'k-NearestNeighbor' or \
                 2) String 'FarthestNeighbor' or \
                 3) String (not available) 'Insertion' or \
                 4) String (not available) 'Patching' or \
                 5) String 'Sweep' or \
                 6) String 'DepthFirst' or \
                 7) String (default) 'Christofides' or \
                 8) String 'Random'" = 'Christofides',
    algoArgs:   "Dictionary, needed for some algorithm options, \
                 1) None for unspecified `algo` options, or \
                 2) for 'k-NearestNeighbor' \
                    {\
                        'k': k-th nearest neighbor\
                    } \
                 3) for 'Christofides' \
                    {\
                        'matchingAlgo': algorithm for finding minimum matching\
                    } " = None
    ) -> "Constructive heuristic solution for TSP":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'LatLon'):
            edges = getTauLatLon(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Subroutines for different constructive heuristics =======================
    def consTSPkNearestNeighbor(nodes, nodeIDs, edges, k = 1):
        seq = [nodeIDs[0]]
        remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
        ofv = 0

        # Accumulate seq ------------------------------------------------------
        while (len(remain) > 0):
            currentNodeID = seq[-1]
            sortedNodes = sortNodesByDist(
                nodes = nodes,
                edges = edges,
                nodeIDs = remain,
                refNodeID = currentNodeID)
            if (k > len(remain)):
                seq.append(sortedNodes[-1])
            else:
                seq.append(sortedNodes[k - 1])
            remain = [i for i in nodeIDs if i not in seq]
        seq.append(seq[0])
        ofv = calSeqCostMatrix(edges, seq)
        return {
            'ofv': ofv,
            'seq': seq
        }
    def consTSPFarthestNeighbor(nodeIDs, edges):
        # Initialize ----------------------------------------------------------
        seq = [nodeIDs[0]]
        remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
        ofv = 0

        # Accumulate seq ------------------------------------------------------
        while (len(remain) > 0):
            nextLeng = None
            nextID = None
            for node in remain:
                if ((node, seq[-1]) in edges):
                    if (nextLeng == None or edges[node, seq[-1]] > nextLeng):
                        nextID = node
                        nextLeng = edges[node, seq[-1]]
                elif ((seq[-1], node) in edges):
                    if (nextLeng == None or edges[seq[-1], node] > nextLeng):
                        nextID = node
                        nextLeng = edges[seq[-1], node]
            seq.append(nextID)
            remain.remove(nextID)
            ofv += nextLeng    
        ofv += edges[seq[0], seq[-1]]
        seq.append(seq[0])

        return {
            'ofv': ofv,
            'seq': seq
        }
    def consTSPSweep(nodes, nodeIDs, edges):
        # Find a center loc ---------------------------------------------------
        centerX = 0
        centerY = 0
        for n in nodes:
            centerX += nodes[n]['loc'][0]
            centerY += nodes[n]['loc'][1]
        centerX /= len(nodes)
        centerY /= len(nodes)

        # Sweep seq -----------------------------------------------------------
        sweepSeq = getSweepSeq(
            nodes = nodes, 
            nodeIDs = nodeIDs,
            centerLoc = [centerX, centerY])
        sweepSeq.append(sweepSeq[0])

        # Calculate ofv -------------------------------------------------------
        ofv = calSeqCostMatrix(edges, sweepSeq)

        return {
            'ofv': ofv,
            'seq': sweepSeq
        }
    def consTSPRandomSeq(nodeIDs, edges):
        # Get random seq ------------------------------------------------------
        seqIndex = rndSeq(len(nodeIDs), closed=True)
        seq = []
        for i in range(len(seqIndex)):
            seq.append(nodeIDs[seqIndex[i]])

        # Calculate Ofv -------------------------------------------------------
        ofv = calSeqCostMatrix(edges, seq)
        return {
            'ofv': ofv,
            'seq': seq
        }

    # Heuristics that don't need to transform arc representation ==============
    tsp = None
    if (algo == 'k-NearestNeighbor'):
        if (algoArgs == None):
            algoArgs = {'k': 1}
        tsp = consTSPkNearestNeighbor(nodes, nodeIDs, edges, algoArgs['k'])
    elif (algo == 'FarthestNeighbor'):
        tsp = consTSPFarthestNeighbor(nodeIDs, edges)
    elif (algo == 'Sweep'):
        tsp = consTSPSweep(nodes, nodeIDs, edges)
    elif (algo == 'Random'):
        tsp = consTSPRandomSeq(nodeIDs, edges)
    else:
        pass
    if (tsp != None):
        return tsp

    # Create arcs =============================================================
    # FIXME! indexes of nodes starts from 0 here!
    # FIXME! Add mapping between 0 to n and nodeIDs
    weightArcs = []
    for (i, j) in edges:
        if (i != None and j != None and i < j):
            weightArcs.append((i, j, edges[i, j]))

    # Subroutines for different constructive heuristics =======================
    def consTSPDepthFirst(weightArcs):
        # Create MST ----------------------------------------------------------
        mst = graphMST(weightArcs)['mst']

        # Seq of visit is the seq of Depth first search on the MST ------------
        seq = traversalGraph(mst)['seq']
        seq.append(seq[0])

        # Calculate ofv -------------------------------------------------------
        ofv = calSeqCostArcs(weightArcs, seq)

        return {
            'ofv': ofv,
            'seq': seq
        }
    def consTSPChristofides(weightArcs, matchingAlgo):
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

        # Calculate ofv -------------------------------------------------------
        ofv = calSeqCostArcs(weightArcs, seq)

        return {
            'ofv': ofv,
            'seq': seq
        }

    # Constructive Heuristics for TSP =========================================
    tsp = None
    if (algo == 'DepthFirst'):
        tsp = consTSPDepthFirst(weightArcs)
    elif (algo == 'Christofides'):
        if (algoArgs == None):
            algoArgs = {'matchingAlgo': 'IP'}
        tsp = consTSPChristofides(weightArcs, algoArgs['matchingAlgo'])

    return tsp

def impTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    initSeq:    "Sequence of visiting nodes, as the initial route, generated by constructive heuristic" = [],
    algo:       "1) String (not available) 'LKH' or \
                 2) String (default) '2Opt' = '2Opt'" = '2Opt'
    ) -> "Local improvement heuristic solution for TSP":

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'LatLon'):
            edges = getTauLatLon(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Subroutines for different local improvement heuristics ==================
    def impTSP2Opt(nodeIDs, edges, initSeq):
        # Initialize ----------------------------------------------------------
        canImproveFlag = True
        impSeq = [i for i in initSeq]

        # Revert part of the sequences ----------------------------------------
        def revert(i, j, seq):
            rSeq = []
            rSeq.extend([seq[k] for k in range(i)])
            rSeq.extend([seq[j - k] for k in range(j - i + 1)])
            rSeq.extend([seq[k] for k in range(j + 1, len(seq))])
            return rSeq

        # Main iteration ------------------------------------------------------
        while (canImproveFlag):
            canImproveFlag = False

            # Try 2-opt -------------------------------------------------------
            # First arc: (i, i + 1)
            # Second arc: (j, j + 1)
            # New arcs: (i, j + 1), (j, i + 1)
            bestSaving = 0
            bestMove = None
            for i in range(len(impSeq) - 3):
                for j in range(i + 2, len(impSeq) - 1):
                    # Saving
                    saving = 0
                    if ((impSeq[i], impSeq[j]) in edges and (impSeq[i + 1], impSeq[j + 1]) in edges):
                        saving = (edges[impSeq[i], impSeq[i + 1]] + edges[impSeq[j], impSeq[j + 1]]) - (edges[impSeq[i], impSeq[j]] + edges[impSeq[i + 1], impSeq[j + 1]])
                    if (saving > bestSaving):
                        bestSaving = saving
                        bestMove = (i, j)

            if (bestSaving > 0):
                (i, j) = bestMove
                impSeq = revert(i + 1, j, impSeq)
                canImproveFlag = True

        # Calculate ofv -------------------------------------------------------
        ofv = calSeqCostMatrix(edges, impSeq)

        return {
            'ofv': ofv,
            'seq': impSeq
        }

    # Heuristics that don't need to transform arc representation ==============
    tsp = None
    nodeIDs = list(nodes.keys())
    if (algo == '2Opt'):
        tsp = impTSP2Opt(nodeIDs, edges, initSeq)

    return tsp



