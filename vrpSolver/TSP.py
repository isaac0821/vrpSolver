import heapq
import math
from gurobipy import *

from .const import *
from .common import *
from .graph import *
from .geometry import *

def ipTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'SphereEuclidean' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    fml:        "1) String (default) 'DFJ_Lazy' or \
                 2) String 'DFJ_PlainLoop' or \
                 3) String (not available) 'DFJ_MaxFlow' or \
                 4) String 'MTZ' or \
                 5) String 'MultiCommodityFlow' or \
                 6) String 'ShortestPath' or \
                 7) String 'QAP'" = 'DFJ_Lazy',
    timeLimit:  "1) Double, in seconds or \
                 2) (default) None, no time limit" = None
    ) -> "Exact solution for TSP":

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'SphereEuclidean'):
            edges = getTauSphereEuclidean(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Solve by different formulations =========================================
    res = None
    nodeIDs = list(nodes.keys())
    if (fml == 'DFJ_Lazy'):
        res = _ipTSPLazyCuts(edges, nodeIDs, timeLimit)
    elif (fml == 'DFJ_PlainLoop'):
        res = _ipTSPPlainLoop(edges, nodeIDs, timeLimit)
    elif (fml == 'MTZ'):
        res = _ipTSPMTZ(edges, nodeIDs, timeLimit)
    elif (fml == 'ShortestPath'):
        res = _ipTSPShortestPath(edges, nodeIDs, timeLimit)
    elif (fml == 'MultiCommodityFlow'):
        res = _ipTSPMultiCommodityFlow(edges, nodeIDs, timeLimit)
    elif (fml == 'QAP'):
        res = _ipTSPQAP(edges, nodeIDs, timeLimit)
    else:
        print("Error: Incorrect or not available TSP formulation option!")
        return None
    if (res != None):
        res['fml'] = fml

    return res

def _ipTSPQAP(edges, nodeIDs, timeLimit):
    n = len(nodeIDs)
    TSP = Model('TSP')

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = TSP.addVar(
                    vtype = GRB.BINARY, 
                    name = 'x_%s_%s' % (i, j))
    w = {}
    for i in range(n):
        for j in range(n):
            if (i != j):
                for k in range(n):
                    w[i, j, k] = TSP.addVar(
                        vtype = GRB.CONTINUOUS)

    # TSP objective function ==================================================
    TSP.setObjective(
        quicksum(
            quicksum(
                quicksum(
                    edges[nodeIDs[i], nodeIDs[j]] * w[i, j, k] for k in range(n - 1)
                ) for j in range(n) if j != i
            ) for i in range(n)
        ) + 
        quicksum(
            quicksum(
                edges[nodeIDs[i], nodeIDs[j]] * w[i, j, n - 1] for j in range(n) if j != i
            ) for i in range(n)
        )
    )

    # Assignment constraints ==================================================
    for i in range(n):
        TSP.addConstr(quicksum(x[i, j] for j in range(n) if j != i) == 1)
    for j in range(n):
        TSP.addConstr(quicksum(x[i, j] for i in range(n) if j != i) == 1)

    # Linearized constraints ==================================================
    for i in range(n):
        for j in range(n):
            if (i != j):
                for k in range(n - 1):
                    if (k != i and k + 1 != j):
                        TSP.addConstr(w[i, j, k] >= x[i, k] + x[j, k + 1] - 1)
    for i in range(n):
        for j in range(n):
            if (i != j):
                if (i != n - 1 and j != 0):
                    TSP.addConstr(w[i, j, n - 1] >= x[i, n - 1] + x[j, 0] - 1)

    # Optimize ================================================================
    if (timeLimit != None):
        TSP.setParam(GRB.Param.TimeLimit, timeLimit)
    TSP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    gap = None
    lb = None
    ub = None
    runtime = None
    if (TSP.status == GRB.status.OPTIMAL):
        ofv = TSP.getObjective().getValue()
        for j in range(n):
            for i in range(n):
                if (i != j and x[i, j].X > 0.9):
                    seq.append(nodeIDs[i])
                    break
        seq.append(seq[0])
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == GRB.status.TIME_LIMIT):
        ofv = None
        seq = []
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPMultiCommodityFlow(edges, nodeIDs, timeLimit):
    n = len(nodeIDs)
    TSP = Model('TSP')

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = TSP.addVar(
                    vtype = GRB.BINARY, 
                    obj = edges[nodeIDs[i], nodeIDs[j]], 
                    name = 'x_%s_%s' % (i, j))
    y = {}
    for i in range(n):
        for j in range(n):
            if (i != j):
                for k in range(1, n):
                    y[i, j, k] = TSP.addVar(
                        vtype = GRB.CONTINUOUS)

    # TSP objective function ==================================================
    TSP.modelSense = GRB.MINIMIZE
    TSP.update()

    # Degree constraints =====================================================
    for i in range(n):
        TSP.addConstr(quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
        TSP.addConstr(quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

    # MCF =====================================================================
    for i in range(n):
        for j in range(n):
            if (i != j):
                for k in range(1, n):
                    TSP.addConstr(y[i, j, k] <= x[i, j])

    for k in range(1, n):
        TSP.addConstr(quicksum(y[0, i, k] for i in range(1, n)) == 1)
        TSP.addConstr(quicksum(y[i, 0, k] for i in range(1, n)) == 0)
        TSP.addConstr(quicksum(y[i, k, k] for i in range(n) if i != k) == 1)
        TSP.addConstr(quicksum(y[k, j, k] for j in range(n) if j != k) == 0)
        for j in range(1, n):
            if (j != k):
                TSP.addConstr(quicksum(y[i, j, k] for i in range(n) if i != j) - quicksum(y[j, i, k] for i in range(n) if i != j) == 0)

    # TSP =====================================================================
    if (timeLimit != None):
        TSP.setParam(GRB.Param.TimeLimit, timeLimit)
    TSP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    if (TSP.status == GRB.status.OPTIMAL):
        ofv = TSP.getObjective().getValue()
        for i, j in x:
            if (x[i, j].x > 0.5):
                arcs.append([i, j])
        currentNode = 0
        currentTime = 0
        seq.append(nodeIDs[currentNode])
        while (len(arcs) > 0):
            for i in range(len(arcs)):
                if (arcs[i][0] == currentNode):
                    currentNode = arcs[i][1]
                    seq.append(nodeIDs[currentNode])
                    arcs.pop(i)
                    break
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == GRB.status.TIME_LIMIT):
        ofv = None
        seq = []
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPShortestPath(edges, nodeIDs, timeLimit):
    n = len(nodeIDs)
    TSP = Model('TSP')

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if (j != i):
                for t in range(n):
                    x[i, j, t] = TSP.addVar(
                        obj = edges[nodeIDs[i], nodeIDs[j]],
                        vtype = GRB.BINARY)

    # Stage constraints =======================================================
    # Start from depot 
    TSP.addConstr(quicksum(x[0, j, 0] for j in range(1, n)) == 1)
    # First stage
    for i in range(1, n):
        TSP.addConstr(quicksum(x[i, j, 1] for j in range(1, n) if i != j) - x[0, i, 0] == 0)
    # In between
    for i in range(1, n):
        for t in range(2, n - 1):
            TSP.addConstr(quicksum(x[i, j, t] for j in range(1, n) if i != j) - quicksum(x[j, i, t - 1] for j in range(1, n) if i != j) == 0)
    # Last stage
    for i in range(1, n):
        TSP.addConstr(x[i, 0, n - 1] - quicksum(x[j, i, n - 2] for j in range(1, n) if i != j) == 0)
    # Return to depot
    TSP.addConstr(quicksum(x[i, 0, n - 1] for i in range(1, n)) == 1)
    # Consequent
    for i in range(1, n):
        TSP.addConstr(quicksum(quicksum(x[i, j, t] for j in range(1, n) if i != j) for t in range(1, n - 1)) + x[i, 0, n - 1] <= 1)

    # TSP =====================================================================
    if (timeLimit != None):
        TSP.setParam(GRB.Param.TimeLimit, timeLimit)
    TSP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    if (TSP.status == GRB.status.OPTIMAL):
        ofv = TSP.getObjective().getValue()
        for i, j, t in x:
            if (x[i, j, t].x > 0.5):
                arcs.append([i, j])
        currentNode = 0
        currentTime = 0
        seq.append(nodeIDs[currentNode])
        while (len(arcs) > 0):
            for i in range(len(arcs)):
                if (arcs[i][0] == currentNode):
                    currentNode = arcs[i][1]
                    seq.append(nodeIDs[currentNode])
                    arcs.pop(i)
                    break
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == GRB.status.TIME_LIMIT):
        ofv = None
        seq = []
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPMTZ(edges, nodeIDs, timeLimit):
    n = len(nodeIDs)
    TSP = Model('TSP')

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = TSP.addVar(
                    vtype = GRB.BINARY, 
                    obj = edges[nodeIDs[i], nodeIDs[j]], 
                    name = 'x_%s_%s' % (i, j))
    u = {}
    for i in range(n):
        u[i] = TSP.addVar(
            vtype = GRB.CONTINUOUS,
            name = 'u_%s' % (i))

    # TSP objective function ==================================================
    TSP.modelSense = GRB.MINIMIZE
    TSP.update()

    # Degree constraints ======================================================
    for i in range(n):
        TSP.addConstr(quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
        TSP.addConstr(quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

    # Sequence constraints ====================================================
    for i in range(1, n):
        for j in range(1, n):
            if (i != j):
                TSP.addConstr(u[i] - u[j] + (n - 1) * x[i, j] <= n - 2, name = 'seq_%s_%s' % (i, j))
    for i in range(1, n):
        TSP.addConstr(1 <= u[i])
        TSP.addConstr(u[i] <= n - 1)

    # TSP =====================================================================
    if (timeLimit != None):
        TSP.setParam(GRB.Param.TimeLimit, timeLimit)
    TSP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    gap = None
    seq = []
    arcs = []
    if (TSP.status == GRB.status.OPTIMAL):
        ofv = TSP.getObjective().getValue()
        gap = TSP.Params.MIPGapAbs
        for i, j in x:
            if (x[i, j].x > 0.5):
                arcs.append([i, j])
        currentNode = 0
        currentTime = 0
        seq.append(nodeIDs[currentNode])
        while (len(arcs) > 0):
            for i in range(len(arcs)):
                if (arcs[i][0] == currentNode):
                    currentNode = arcs[i][1]
                    seq.append(nodeIDs[currentNode])
                    arcs.pop(i)
                    break
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == GRB.status.TIME_LIMIT):
        ofv = None
        seq = []
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPPlainLoop(edges, nodeIDs, timeLimit):
    n = len(nodeIDs)
    TSP = Model('TSP')

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if (i != j):
                x[i, j] = TSP.addVar(
                    vtype = GRB.BINARY, 
                    obj = edges[nodeIDs[i], nodeIDs[j]], 
                    name = 'x_%s_%s' % (i, j))
                
    # TSP =====================================================================
    TSP.modelSense = GRB.MINIMIZE
    TSP.Params.lazyConstraints = 1
    TSP.update()

    # Degree constraints ======================================================
    for i in range(n):
        TSP.addConstr(quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
        TSP.addConstr(quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

    # Resolve to optimality and try to find sub-tours =========================
    noSubtourFlag = False
    accRuntime = 0
    while (not noSubtourFlag):
        if (timeLimit != None):
            TSP.setParam(GRB.Param.TimeLimit, timeLimit - accRuntime) # FIXME
        TSP.optimize()
        if (TSP.status == GRB.status.OPTIMAL):
            accRuntime += TSP.Runtime
            arcs = tuplelist((i, j) for i, j in x.keys() if x[i, j].X > 0.9)
            components = findComponentsUndirected(arcs)
            if (len(components) == 1):
                noSubtourFlag = True
                break
            else:
                for comp in components:
                    TSP.addConstr(quicksum(x[i, j] for i in comp for j in comp if i != j) <= len(comp) - 1)
        elif (TSP.status == GRB.status.TIME_LIMIT):
            accRuntime += TSP.Runtime
            break

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    if (TSP.status == GRB.status.OPTIMAL):
        ofv = TSP.getObjective().getValue()
        for i, j in x:
            if (x[i, j].x > 0.5):
                arcs.append([i, j])
        currentNode = 0
        currentTime = 0
        seq.append(nodeIDs[currentNode])
        while (len(arcs) > 0):
            for i in range(len(arcs)):
                if (arcs[i][0] == currentNode):
                    currentNode = arcs[i][1]
                    seq.append(nodeIDs[currentNode])
                    arcs.pop(i)
                    break
        gap = 0
        lb = ofv
        ub = ofv
        runtime = accRuntime
    elif (TSP.status == GRB.status.TIME_LIMIT):
        ofv = None
        seq = []
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = accRuntime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPLazyCuts(edges, nodeIDs, timeLimit):
    n = len(nodeIDs)
    TSP = Model('TSP')

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if (i != j):
                x[i, j] = TSP.addVar(
                    vtype = GRB.BINARY, 
                    obj = edges[nodeIDs[i], nodeIDs[j]], 
                    name = 'x_%s_%s' % (i, j))
                
    # TSP objective function ==================================================
    TSP.modelSense = GRB.MINIMIZE
    TSP.Params.lazyConstraints = 1
    TSP.update()

    # Degree constraints ======================================================
    for i in range(n):
        TSP.addConstr(quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
        TSP.addConstr(quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

    # Sub-tour elimination ====================================================
    TSP._x = x
    def subtourelim(model, where):
        if (where == GRB.Callback.MIPSOL):
            x_sol = model.cbGetSolution(model._x)
            arcs = tuplelist((i, j) for i, j in model._x.keys() if x_sol[i, j] > 0.9)
            components = findComponentsUndirected(arcs)
            for component in components:
                if (len(component) < n):
                    model.cbLazy(quicksum(x[i,j] for i in component for j in component if i != j) <= len(component) - 1)

    # TSP with callback =======================================================
    if (timeLimit != None):
        TSP.setParam(GRB.Param.TimeLimit, timeLimit)
    TSP.setParam('OutputFlag', 0)
    TSP.optimize(subtourelim)

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    if (TSP.status == GRB.status.OPTIMAL):
        ofv = TSP.getObjective().getValue()
        for i, j in x:
            if (x[i, j].x > 0.5):
                arcs.append([i, j])
        currentNode = 0
        currentTime = 0
        seq.append(nodeIDs[currentNode])
        while (len(arcs) > 0):
            for i in range(len(arcs)):
                if (arcs[i][0] == currentNode):
                    currentNode = arcs[i][1]
                    seq.append(nodeIDs[currentNode])
                    arcs.pop(i)
                    break
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == GRB.status.TIME_LIMIT):
        ofv = None
        seq = []
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def lrTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'SphereEuclidean' or \
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
        elif (edges == 'SphereEuclidean'):
            edges = getTauSphereEuclidean(nodes)
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
                 2) String 'SphereEuclidean' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    consAlgo:   "1) String 'NearestNeighbor' or \
                 2) String 'FarthestNeighbor' or \
                 3) String (not available) 'Insertion' or \
                 4) String (not available) 'Patching' or \
                 5) String (not available) 'Sweep' or \
                 6) String 'DepthFirst' or \
                 7) String (default) 'Christofides' or \
                 8) String 'Random'" = 'Christofides',
    impAlgo:    "1) String (not available) 'LKH' or \
                 2) String (default) '2Opt' = '2Opt'" = '2Opt'
    ) -> "Use given heuristic methods to get TSP solution":

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'SphereEuclidean'):
            edges = getTauSphereEuclidean(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Constructive heuristic ==================================================
    # Start with a constructive heuristic
    cons = consTSP(nodes, edges, consAlgo)

    # Improve the constructive heuristic result ===============================
    res = impTSP(nodes, edges, cons['seq'], impAlgo)

    return res

def consTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'SphereEuclidean' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    algo:       "1) String 'NearestNeighbor' or \
                 2) String 'FarthestNeighbor' or \
                 3) String (not available) 'Insertion' or \
                 4) String (not available) 'Patching' or \
                 5) String 'Sweep' or \
                 6) String 'DepthFirst' or \
                 7) String (default) 'Christofides' or \
                 8) String 'Random'" = 'Christofides'
    ) -> "Constructive heuristic solution for TSP":

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'SphereEuclidean'):
            edges = getTauSphereEuclidean(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Heuristics that don't need to transform arc representation ==============
    res = None
    nodeIDs = list(nodes.keys())
    if (algo == 'NearestNeighbor'):
        res = _consTSPNearestNeighbor(nodeIDs, edges)
    elif (algo == 'FarthestNeighbor'):
        res = _consTSPFarthestNeighbor(nodeIDs, edges)
    elif (algo == 'Sweep'):
        res = _consTSPSweep(nodes, edges)
    elif (algo == 'Random'):
        res = _consTSPRandomSeq(nodeIDs, edges)
    else:
        pass
    if (res != None):
        return res

    # Create arcs =============================================================
    # FIXME! indexes of nodes starts from 0 here!
    # FIXME! Add mapping between 0 to n and nodeIDs
    weightArcs = []
    for (i, j) in edges:
        if (i != None and j != None and i < j):
            weightArcs.append((i, j, edges[i, j]))

    # Constructive Heuristics for TSP =========================================
    res = None
    if (algo == 'DepthFirst'):
        res = _consTSPDepthFirst(weightArcs)
    elif (algo == 'Christofides'):
        res = _consTSPChristofides(weightArcs)

    return res

def _consTSPRandomSeq(nodeIDs, edges):
    # Get random seq ==========================================================
    seqIndex = rndSeq(len(nodeIDs), closed=True)
    seq = []
    for i in range(len(seqIndex)):
        seq.append(nodeIDs[seqIndex[i]])

    # Calculate Ofv ===========================================================
    ofv = calSeqCostMatrix(edges, seq)
    return {
        'ofv': ofv,
        'seq': seq
    }

def _consTSPNearestNeighbor(nodeIDs, edges):
    # Initialize ==============================================================
    seq = [nodeIDs[0]]
    remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
    ofv = 0

    # Accumulate seq ==========================================================
    while (len(remain) > 0):
        nextLeng = None
        nextID = None
        for node in remain:
            if ((node, seq[-1]) in edges):
                if (nextLeng == None or edges[node, seq[-1]] < nextLeng):
                    nextID = node
                    nextLeng = edges[node, seq[-1]]
            elif ((seq[-1], node) in edges):
                if (nextLeng == None or edges[seq[-1], node] < nextLeng):
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

def _consTSPFarthestNeighbor(nodeIDs, edges):
    # Initialize ==============================================================
    seq = [nodeIDs[0]]
    remain = [nodeIDs[i] for i in range(1, len(nodeIDs))]
    ofv = 0

    # Accumulate seq ==========================================================
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

def _consTSPSweep(nodes, edges):
    # Find a center loc =======================================================
    centerX = 0
    centerY = 0
    for n in nodes:
        centerX += nodes[n]['loc'][0]
        centerY += nodes[n]['loc'][1]
    centerX /= len(nodes)
    centerY /= len(nodes)

    # Sweep seq ===============================================================
    sweepSeq = getSweepSeq(nodes = nodes, centerLoc = [centerX, centerY])
    sweepSeq.append(sweepSeq[0])

    # Calculate ofv ===========================================================
    ofv = calSeqCostMatrix(edges, sweepSeq)

    return {
        'ofv': ofv,
        'seq': sweepSeq
    }

def _consTSPDepthFirst(weightArcs):
    # Create MST ==============================================================
    mst = graphMST(weightArcs)['mst']

    # Seq of visit is the seq of Depth first search on the MST ================
    seq = traversalGraph(mst)['seq']
    seq.append(seq[0])

    # Calculate ofv ===========================================================
    ofv = calSeqCostArcs(weightArcs, seq)

    return {
        'ofv': ofv,
        'seq': seq
    }

def _consTSPChristofides(weightArcs):
    # Create MST ==============================================================
    mst = graphMST(weightArcs)['mst']

    # Derive subgraph of odd degree vertices ==================================
    neighbors = arcs2AdjList(mst)
    oddDegrees = []
    for node in neighbors:
        if (len(neighbors[node]) % 2 != 0):
            oddDegrees.append(node)
    subGraph = []
    for arc in weightArcs:
        if (arc[0] in oddDegrees and arc[1] in oddDegrees):
            subGraph.append(arc)

    # Find minimum cost matching of the subgraph ==============================
    minMatching = graphMatching(weightArcs=subGraph, mType='Minimum', algo='IP')['matching']

    # Add them back to create a new graph =====================================
    newGraph = []
    for arc in minMatching:
        newGraph.append(arc)
    for arc in mst:
        newGraph.append(arc)

    # Traverse graph and get seq ==============================================
    # Try to find a vertex with degree 1
    oID = None
    for node in neighbors:
        if (len(neighbors[node]) == 1):
            oID = node
            break
    seq = traversalGraph(newGraph, oID=oID)['seq']
    seq.append(seq[0])

    # Calculate ofv ===========================================================
    ofv = calSeqCostArcs(weightArcs, seq)

    return {
        'ofv': ofv,
        'seq': seq
    }

def impTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'SphereEuclidean' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    initSeq:    "Sequence of visiting nodes, as the initial route, generated by constructive heuristic" = [],
    algo:       "1) String (not available) 'LKH' or \
                 2) String (default) '2Opt' = '2Opt'" = '2Opt'
    ) -> "Local improvement heuristic solution for TSP":

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'SphereEuclidean'):
            edges = getTauSphereEuclidean(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Heuristics that don't need to transform arc representation ==============
    res = None
    nodeIDs = list(nodes.keys())
    if (algo == '2Opt'):
        res = _impTSP2Opt(nodeIDs, edges, initSeq)
    return res

def _impTSP2Opt(nodeIDs, edges, initSeq):
    # Initialize ==============================================================
    canImproveFlag = True
    impSeq = [i for i in initSeq]

    # Revert part of the sequences ============================================
    def revert(i, j, seq):
        rSeq = []

        rSeq.extend([seq[k] for k in range(i)])
        rSeq.extend([seq[j - k] for k in range(j - i + 1)])
        rSeq.extend([seq[k] for k in range(j + 1, len(seq))])

        return rSeq

    # Main iteration ==========================================================
    while (canImproveFlag):
        canImproveFlag = False

        # Try 2-opt -----------------------------------------------------------
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

    # Calculate ofv ===========================================================
    ofv = calSeqCostMatrix(edges, impSeq)

    return {
        'ofv': ofv,
        'seq': impSeq
    }

