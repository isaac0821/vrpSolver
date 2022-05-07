import heapq
import datetime

from .common import *
from .const import *
from .msg import *
from .graph import *
from .geometry import *

def lbTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) String 'Grid' or \
                 4) Dictionary {(nodeID1, nodeID2): dist, ...}" = 'Euclidean',
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                    {\
                        'colRow': (numCol, numRow),\
                        'barriers': [(coordX, coordY), ...], \
                    }" = None,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    algo:       "1) String (default) 'HeldKarp'" = 'HeldKarp',
    algoArgs:   "Dictionary, args for the lower bound defining algorithm \
                 1) for 'HeldKarp' \
                    {\
                        'subgradM': subgradM,\
                        'subgradRho': subgradRho, (0, 1),\
                        'stopType': stopType, # Stopping criteria, with the following options: 'Epsilon', 'IterNum', 'Runtime'\
                        'stopEpsilon': epsilon, for stopping criteria 'Epsilon'\
                        'stopK': number of iterations, for stopping criteria 'IterNum'\
                        'stopTime': runtime, for stopping criteria 'Runtime'\
                    }" = None
    ) -> "Returns a Held & Karp lower bound of the TSP using Lagrangian Relaxation":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)
        else:
            msgError(ERROR_INCOR_NODEIDS)
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
            msgError(ERROR_INCOR_TAU)
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
        
    # Default configuration ===================================================
    if (algo == 'HeldKarp' and algoArgs == None):
        algoArgs = {
            'subgradM': 1,
            'subgradRho': 0.95,
            'stopType': 'Epsilon',
            'stopEpsilon': 0.01,
            'stopK': 200,
            'stopRuntime': 600
        }

    # Subroutines for different lower bound algorithms ========================
    def _lbTSPHeldKarp():
        # Initialize ----------------------------------------------------------
        k = 0
        u = [0 for i in range(len(nodeIDs))]
        d = None
        costSum = None
        L = None
        oldL = None
        startTime = datetime.datetime.now()

        # Calculate 1 tree ----------------------------------------------------
        def _cal1Tree(weightArcs):
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

        # Main iteration ------------------------------------------------------
        continueFlag = True
        while (continueFlag):
            # Update cost of each edge
            weightArcs = []
            for i in range(len(nodeIDs)):
                for j in range(len(nodeIDs)):
                    if (i != None and j != None and i < j):
                        weightArcs.append((i, j, tau[i, j] - u[i] - u[j]))

            # Calculate 1-tree
            oneTree = _cal1Tree(weightArcs)

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
            eff = algoArgs['subgradM'] * math.pow(algoArgs['subgradRho'], k)
            for i in range(len(nodeIDs)):
                u.append(oldU[i] + eff * d[i])

            # Check if continue
            def _allZero(d):
                for i in d:
                    if (i != 0):
                        return False
                return True

            # Check stop criteria
            if (algoArgs['stopType'] == 'Epsilon'):
                if (oldL != None and abs(oldL - L) <= algoArgs['stopEpsilon']):
                    continueFlag = False
            elif (algoArgs['stopType'] == 'IterNum'):
                if (k >= algoArgs['stopK']):
                    continueFlag = False
            elif (algoArgs['stopType'] == 'Runtime'):
                if ((datetime.datetime.now() - startTime).total_seconds() >= algoArgs['stopRuntime']):
                    continueFlag = False
            else:
                if (_allZero(d)):
                    continueFlag = False
            k += 1

        return {
            'lowerBound': costSum,
            'm1t': m1t,
            'runtime': (datetime.datetime.now() - startTime).total_seconds()
        }

    # Choose lower bound algorithm ============================================
    lowerBound = None
    if (algo == 'HeldKarp'):
        lowerBound = _lbTSPHeldKarp()
    else:
        print("ERROR: Acceptable options for `algo` is: 'HeldKarp'")
        return

    return lowerBound