import heapq
import math
import warnings
import networkx as nx

from .common import *
from .geometry import *
from .msg import *

def solveTSP(
    nodes: dict, 
    locFieldName: str = 'loc',
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    serviceTime: float = 0,
    predefinedArcs: list[list[tuple[int|str]]] = [],
    vehicles: dict = {0: {'speed': 1}},
    vehicleID: int|str = 0,
    edges: str = 'Euclidean',
    algo: str = 'IP',
    detailsFlag: bool = False,
    metaFlag: bool = False,
    **kwargs
    ) -> dict:


    """Use MIP/Heuristic/Metaheuristic to find optimal/suboptimal TSP solution

    Parameters
    ----------

    nodes: dict, required
        The coordinates of given nodes.
    locFieldName: str, optional, default as 'loc'
        The key value in `nodes` indicating the location of each node.
    depotID: int|string, optional, default as 0
        The ID of depot.
    nodeIDs: string 'All' or a list of node IDs, optional, default as 'All'
        The following are two options: 1) 'All', all nodes will be visited, 2) A list of node IDs to be visited.
    serviceTime: float, optional, default as 0
        The service time needed at each location.
    edges: string, optional, default as 'Euclidean'
        The methods for the calculation of distances between nodes. Options and required additional information are referred to :func:`~vrpSolver.geometry.matrixDist()`.

    Returns
    -------

    dict
        A TSP solution in the following format::
            >>> solution = {
            ...     'ofv': ofv,
            ...     'seq': seq,
            ...     'gap': gap,
            ...     'solType': solType,
            ...     'lowerBound': lb,
            ...     'upperBound': ub,
            ...     'runtime': runtime
            ... }

    """

    # Sanity check ============================================================
    if (nodes == None or type(nodes) != dict):
        raise MissingParameterError(ERROR_MISSING_NODES)
    for i in nodes:
        if (locFieldName not in nodes[i]):
            raise MissingParameterError("ERROR: missing location information in `nodes`.")
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

    if (detailsFlag == True):
        # For animation propose
        if (vehicles == None):
            raise MissingParameterError("ERROR: Missing required field `vehicles`.")
        if (vehicleID not in vehicles):
            raise MissingParameterError("ERROR: Cannot find `vehicleID` in `vehicles`.")

    if (algo == 'IP'):
        if ('solver' not in kwargs):
            raise MissingParameterError("ERROR: Missing required field `solver`.")
        elif (kwargs['solver'] == 'Gurobi' and kwargs['fml'] not in ['DFJ_Lazy', 'DFJ_Plainloop', 'MTZ', 'ShortestPath', 'MultiCommodityFlow', 'QAP']):
            raise OutOfRangeError("ERROR: Gurobi option s upports 'DFJ_Lazy', 'DFJ_Plainloop', 'MTZ', 'ShortestPath', 'MultiCommodityFlow', and 'QAP' formulations", )
        elif (kwargs['solver'] == 'COPT' and kwargs['fml'] not in ['DFJ_Lazy']):
            raise OutOfRangeError("ERROR: COPT option supports 'DFJ_Lazy' formulations", )
    elif (algo == 'Heuristic'):
        if ('cons' not in kwargs and 'impv' not in kwargs):
            raise MissingParameterError(ERROR_MISSING_TSP_ALGO)

    if (predefinedArcs != None and len(predefinedArcs) > 0):
        if (not (algo == 'IP' and 'fml' in kwargs and kwargs['fml'] in ['DFJ_Lazy'])):
            raise OutOfRangeError("ERROR: TSP with pre-defined arcs is supported by DFJ_Lazy only, for now.")

    # Define tau ==============================================================
    tau = None
    path = None
    if (detailsFlag):
        tau, path = matrixDist(
            nodes = nodes, 
            nodeIDs = nodeIDs,
            edges = edges, 
            locFieldName = locFieldName,
            **kwargs)
    else:
        tau, _ = matrixDist(
            nodes = nodes, 
            nodeIDs = nodeIDs,
            edges = edges, 
            locFieldName = locFieldName,
            **kwargs)

    # Check symmetric =========================================================
    asymFlag = False
    for (i, j) in tau:
        if (tau[i, j] != tau[j, i]):
            asymFlag = True
            break

    # TSP =====================================================================
    tsp = None
    if (algo == 'IP'):
        outputFlag = None if 'outputFlag' not in kwargs else kwargs['outputFlag']
        timeLimit = None if 'timeLimit' not in kwargs else kwargs['timeLimit']
        gapTolerance = None if 'gapTolerance' not in kwargs else kwargs['gapTolerance']
        tsp = _ipTSP(
            nodeIDs = nodeIDs, 
            tau = tau, 
            solver = kwargs['solver'], 
            fml = kwargs['fml'], 
            predefinedArcs = predefinedArcs, 
            outputFlag = outputFlag, 
            timeLimit = timeLimit, 
            gapTolerance = gapTolerance)
        tsp['fml'] = kwargs['fml']
        tsp['solver'] = kwargs['solver']
    elif (algo == 'Heuristic'):
        nodeObj = {}
        for n in nodeIDs:
            nodeObj[n] = RouteNode(n, value=nodes[n][locFieldName])
        tsp = _heuTSP(
            nodeObj = nodeObj, 
            tau = tau, 
            depotID = depotID, 
            nodes = nodeIDs, 
            asymFlag = asymFlag, 
            **kwargs)
    else:
        raise OutOfRangeError("ERROR: Select 'algo' from ['IP', 'Heuristic'].")

    # Fix the sequence to make it start from the depot ========================
    startIndex = 0
    seq = [i for i in tsp['seq']]
    nodeSeq = []
    for k in range(len(seq)):
        if (seq[k] == depotID):
            startIndex = k
    if (startIndex <= len(seq) - 1):
        for k in range(startIndex, len(seq) - 1):
            nodeSeq.append(seq[k])
    if (startIndex >= 0):
        for k in range(0, startIndex):
            nodeSeq.append(seq[k])
    nodeSeq.append(depotID)
    tsp['seq'] = nodeSeq

    # Add service time if provided ============================================
    ofv = tsp['ofv'] + (len(nodeIDs) - 1) * serviceTime

    # Post optimization (for detail information) ==============================
    if (detailsFlag):
        # 返回一个数组，表示路径中的每个点的位置，不包括时间信息
        shapepoints = []        
        for i in range(len(nodeSeq) - 1):
            shapepoints.extend(path[nodeSeq[i], nodeSeq[i + 1]][:-1])
        shapepoints.append(path[nodeSeq[-2], nodeSeq[-1]][-1])

        # 返回一个数组，其中每个元素为二元数组，表示位置+时刻
        curTime = 0
        curLoc = nodes[depotID][locFieldName]
        timedSeq = [(curLoc, curTime)]
        # 对每个leg检索path中的shapepoints，涉及到serviceTime，先不看最后一段leg
        for i in range(1, len(nodeSeq) - 1):
            # 对于Euclidean型的，没有中间节点
            if (edges in ['Euclidean', 'LatLon']):
                curTime += tau[nodeSeq[i - 1], nodeSeq[i]] / vehicles[vehicleID]['speed']
                curLoc = nodes[nodeSeq[i]][locFieldName]
                timedSeq.append((curLoc, curTime))
            else:
                shapepointsInBtw = path[nodeSeq[i - 1], nodeSeq[i]]
                for j in range(1, len(shapepointsInBtw)):
                    curTime += distEuclideanXY(shapepointsInBtw[j - 1], shapepointsInBtw[j])['dist'] / vehicles[vehicleID]['speed']
                    curLoc = shapepointsInBtw[j]
                    timedSeq.append((curLoc, curTime))
            # 如果有service time，则加上一段在原处等待的时间
            if (serviceTime != None and serviceTime > 0):
                curTime += serviceTime
                # curLoc = curLoc
                timedSeq.append((curLoc, curTime))
        # 现在补上最后一段leg
        if (edges in ['Euclidean', 'LatLon']):
            curTime += tau[nodeSeq[-2], nodeSeq[-1]] / vehicles[vehicleID]['speed']
            curLoc = nodes[nodeSeq[-1]][locFieldName]
            timedSeq.append((curLoc, curTime))
        else:
            shapepointsInBtw = path[nodeSeq[-2], nodeSeq[-1]]
            for j in range(1, len(shapepointsInBtw)):
                curTime += distEuclideanXY(shapepointsInBtw[j - 1], shapepointsInBtw[j])['dist'] / vehicles[vehicleID]['speed']
                curLoc = shapepointsInBtw[j]
                timedSeq.append((curLoc, curTime))

        # Add detail information to `vehicles`
        vehicles[vehicleID]['shapepoints'] = shapepoints
        vehicles[vehicleID]['timedSeq'] = timedSeq

    # Add service time info ===================================================
    res = {
        'ofv': ofv,
        'seq': nodeSeq,
    }
    if (algo == 'IP'):
        res['gap'] = tsp['gap']
        res['solType'] = tsp['solType']
        res['lowerBound'] = tsp['lowerBound']
        res['upperBound'] = tsp['upperBound']
        res['runtime'] = tsp['runtime']    
    if (metaFlag):
        res['serviceTime'] = serviceTime
    if (detailsFlag):
        res['vehicles'] = vehicles

    return res

def _ipTSP(nodeIDs, tau, solver, fml, predefinedArcs, outputFlag, timeLimit, gapTolerance) -> dict|None:
    tsp = None
    if (solver == 'Gurobi'):
        if (fml == 'DFJ_Lazy'):
            tsp = _ipTSPGurobiLazyCuts(nodeIDs, tau, predefinedArcs, outputFlag, timeLimit, gapTolerance)
        elif (fml == 'DFJ_Plainloop'):
            tsp = _ipTSPGurobiPlainLoop(nodeIDs, tau, outputFlag, timeLimit, gapTolerance)
        elif (fml == 'MTZ'):
            tsp = _ipTSPGurobiMTZ(nodeIDs, tau, outputFlag, timeLimit, gapTolerance)
        elif (fml == 'ShortestPath'):
            tsp = _ipTSPGurobiShortestPath(nodeIDs, tau, outputFlag, timeLimit, gapTolerance)
        elif (fml == 'MultiCommodityFlow'):
            tsp = _ipTSPGurobiMultiCommodityFlow(nodeIDs, tau, outputFlag, timeLimit, gapTolerance)
        elif (fml == 'QAP'):
            tsp = _ipTSPGurobiQAP(nodeIDs, tau, outputFlag, timeLimit, gapTolerance)
    elif (solver == 'COPT'):
        if (fml == 'DFJ_Lazy'):
            tsp = _ipTSPCOPTLazyCuts(nodeIDs, tau, outputFlag, timeLimit)
    if (tsp == None):
        raise UnsupportedInputError("ERROR: Incorrect or not available TSP formulation option!")
    
    return tsp

def _ipTSPGurobiLazyCuts(nodeIDs, tau, predefinedArcs, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except:
        raise ImportError("ERROR: Cannot find Gurobi")

    # Initialize
    n = len(nodeIDs)
    TSP = grb.Model('TSP')
    if (outputFlag == False):
        TSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                x[i, j] = TSP.addVar(
                    vtype = grb.GRB.BINARY, 
                    obj = tau[i, j], 
                    name = 'x_%s_%s' % (i, j))
                
    # TSP objective function ==================================================
    TSP.modelSense = grb.GRB.MINIMIZE
    TSP.Params.lazyConstraints = 1
    TSP.update()

    # Degree constraints ======================================================
    for i in nodeIDs:
        TSP.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) == 1, name = 'leave_%s' % str(i))
        TSP.addConstr(grb.quicksum(x[j, i] for j in nodeIDs if i != j) == 1, name = 'enter_%s' % str(i))

    # Predefined arcs =========================================================
    if (predefinedArcs != None and len(predefinedArcs) > 0):
        for arcs in predefinedArcs:
            for arc in arcs:
                if ((arc[0], arc[1]) not in x or (arc[1], arc[0]) not in x):
                    raise UnsupportedInputError("ERROR: Cannot find arc[%s, %s] and/or arc[%s, %s]" % (arc[0], arc[1], arc[1], arc[0]))
            TSP.addConstr(grb.quicksum(x[arc[0], arc[1]] for arc in arcs) == 1)

    # Sub-tour elimination ====================================================
    TSP._x = x
    def subtourelim(model, where):
        if (where == grb.GRB.Callback.MIPSOL):
            x_sol = model.cbGetSolution(model._x)
            G = nx.Graph()
            for (i, j) in x.keys():
                if (x_sol[i, j] > 0.9):
                    G.add_edge(i, j, weight = tau[i, j])
            components = [list(c) for c in nx.connected_components(G)]
            for component in components:
                if (len(component) < n):
                    model.cbLazy(grb.quicksum(x[i, j] for i in component for j in component if i != j) <= len(component) - 1)

    # TSP with callback =======================================================
    TSP.optimize(subtourelim)

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None

    ofv = TSP.getObjective().getValue()
    for i, j in x:
        if (x[i, j].x > 0.5):
            arcs.append([i, j])
    currentNode = nodeIDs[0]
    seq.append(currentNode)
    while (len(arcs) > 0):
        for i in range(len(arcs)):
            if (arcs[i][0] == currentNode):
                currentNode = arcs[i][1]
                seq.append(currentNode)
                arcs.pop(i)
                break
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPGurobiPlainLoop(nodeIDs, tau, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
    n = len(nodeIDs)
    TSP = grb.Model('TSP')
    if (outputFlag == False):
        TSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if (i != j):
                x[i, j] = TSP.addVar(
                    vtype = grb.GRB.BINARY, 
                    obj = tau[nodeIDs[i], nodeIDs[j]], 
                    name = 'x_%s_%s' % (i, j))
                
    # TSP =====================================================================
    TSP.modelSense = grb.GRB.MINIMIZE
    TSP.Params.lazyConstraints = 1
    TSP.update()

    # Degree constraints ======================================================
    for i in range(n):
        TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % str(i))
        TSP.addConstr(grb.quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % str(i))

    # Resolve to optimality and try to find sub-tours =========================
    noSubtourFlag = False
    accRuntime = 0
    while (not noSubtourFlag):
        if (timeLimit != None):
            TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit - accRuntime)
        TSP.optimize()
        if (TSP.status == grb.GRB.status.OPTIMAL):
            accRuntime += TSP.Runtime
            G = nx.Graph()
            for i, j in x.keys():
                if (x[i, j].x > 0.9):
                    G.add_edge(i, j, weight = tau[i, j])
            components = [list(c) for c in nx.connected_components(G)]
            if (len(components) == 1):
                noSubtourFlag = True
                break
            else:
                for comp in components:
                    TSP.addConstr(grb.quicksum(x[i, j] for i in comp for j in comp if i != j) <= len(comp) - 1)
        elif (TSP.status == grb.GRB.status.TIME_LIMIT):
            accRuntime += TSP.Runtime
            break

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None

    ofv = TSP.getObjective().getValue()
    for i, j in x:
        if (x[i, j].x > 0.5):
            arcs.append([i, j])
    currentNode = 0
    seq.append(nodeIDs[currentNode])
    while (len(arcs) > 0):
        for i in range(len(arcs)):
            if (arcs[i][0] == currentNode):
                currentNode = arcs[i][1]
                seq.append(nodeIDs[currentNode])
                arcs.pop(i)
                break
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        gap = 0
        lb = ofv
        ub = ofv
        runtime = accRuntime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = accRuntime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPGurobiMTZ(nodeIDs, tau, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
    n = len(nodeIDs)
    TSP = grb.Model('TSP')
    if (outputFlag == False):
        TSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = TSP.addVar(
                    vtype = grb.GRB.BINARY, 
                    obj = tau[nodeIDs[i], nodeIDs[j]], 
                    name = 'x_%s_%s' % (i, j))
    u = {}
    for i in range(n):
        u[i] = TSP.addVar(
            vtype = grb.GRB.CONTINUOUS,
            name = 'u_%s' % (i))

    # TSP objective function ==================================================
    TSP.modelSense = grb.GRB.MINIMIZE
    TSP.update()

    # Degree constraints ======================================================
    for i in range(n):
        TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % str(i))
        TSP.addConstr(grb.quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % str(i))

    # Sequence constraints ====================================================
    for i in range(1, n):
        for j in range(1, n):
            if (i != j):
                TSP.addConstr(u[i] - u[j] + (n - 1) * x[i, j] <= n - 2, name = 'seq_%s_%s' % (i, j))
    for i in range(1, n):
        TSP.addConstr(1 <= u[i])
        TSP.addConstr(u[i] <= n - 1)

    # TSP =====================================================================
    TSP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    gap = None
    seq = []
    arcs = []
    solType = None
    lb = None
    ub = None
    runtime = None

    ofv = TSP.getObjective().getValue()
    gap = TSP.Params.MIPGapAbs
    for i, j in x:
        if (x[i, j].x > 0.5):
            arcs.append([i, j])
    currentNode = 0
    seq.append(nodeIDs[currentNode])
    while (len(arcs) > 0):
        for i in range(len(arcs)):
            if (arcs[i][0] == currentNode):
                currentNode = arcs[i][1]
                seq.append(nodeIDs[currentNode])
                arcs.pop(i)
                break    
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPGurobiShortestPath(nodeIDs, tau, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
    n = len(nodeIDs)
    TSP = grb.Model('TSP')
    if (outputFlag == False):
        TSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if (j != i):
                for t in range(n):
                    x[i, j, t] = TSP.addVar(
                        obj = tau[nodeIDs[i], nodeIDs[j]],
                        vtype = grb.GRB.BINARY)

    # Stage constraints =======================================================
    # Start from depot 
    TSP.addConstr(grb.quicksum(x[0, j, 0] for j in range(1, n)) == 1)
    # First stage
    for i in range(1, n):
        TSP.addConstr(grb.quicksum(x[i, j, 1] for j in range(1, n) if i != j) - x[0, i, 0] == 0)
    # In between
    for i in range(1, n):
        for t in range(2, n - 1):
            TSP.addConstr(grb.quicksum(x[i, j, t] for j in range(1, n) if i != j) - grb.quicksum(x[j, i, t - 1] for j in range(1, n) if i != j) == 0)
    # Last stage
    for i in range(1, n):
        TSP.addConstr(x[i, 0, n - 1] - grb.quicksum(x[j, i, n - 2] for j in range(1, n) if i != j) == 0)
    # Return to depot
    TSP.addConstr(grb.quicksum(x[i, 0, n - 1] for i in range(1, n)) == 1)
    # Consequent
    for i in range(1, n):
        TSP.addConstr(grb.quicksum(grb.quicksum(x[i, j, t] for j in range(1, n) if i != j) for t in range(1, n - 1)) + x[i, 0, n - 1] <= 1)

    # TSP =====================================================================
    TSP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None

    ofv = TSP.getObjective().getValue()
    for i, j, t in x:
        if (x[i, j, t].x > 0.5):
            arcs.append([i, j])
    currentNode = 0
    seq.append(nodeIDs[currentNode])
    while (len(arcs) > 0):
        for i in range(len(arcs)):
            if (arcs[i][0] == currentNode):
                currentNode = arcs[i][1]
                seq.append(nodeIDs[currentNode])
                arcs.pop(i)
                break
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPGurobiMultiCommodityFlow(nodeIDs, tau, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
    n = len(nodeIDs)
    TSP = grb.Model('TSP')
    if (outputFlag == False):
        TSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = TSP.addVar(
                    vtype = grb.GRB.BINARY, 
                    obj = tau[nodeIDs[i], nodeIDs[j]], 
                    name = 'x_%s_%s' % (i, j))
    y = {}
    for i in range(n):
        for j in range(n):
            if (i != j):
                for k in range(1, n):
                    y[i, j, k] = TSP.addVar(
                        vtype = grb.GRB.CONTINUOUS)

    # TSP objective function ==================================================
    TSP.modelSense = grb.GRB.MINIMIZE
    TSP.update()

    # Degree constraints ======================================================
    for i in range(n):
        TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % str(i))
        TSP.addConstr(grb.quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % str(i))

    # MCF ====================================================================
    for i in range(n):
        for j in range(n):
            if (i != j):
                for k in range(1, n):
                    TSP.addConstr(y[i, j, k] <= x[i, j])

    for k in range(1, n):
        TSP.addConstr(grb.quicksum(y[0, i, k] for i in range(1, n)) == 1)
        TSP.addConstr(grb.quicksum(y[i, 0, k] for i in range(1, n)) == 0)
        TSP.addConstr(grb.quicksum(y[i, k, k] for i in range(n) if i != k) == 1)
        TSP.addConstr(grb.quicksum(y[k, j, k] for j in range(n) if j != k) == 0)
        for j in range(1, n):
            if (j != k):
                TSP.addConstr(grb.quicksum(y[i, j, k] for i in range(n) if i != j) 
                    - grb.quicksum(y[j, i, k] for i in range(n) if i != j) == 0)

    # TSP =====================================================================
    TSP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None

    ofv = TSP.getObjective().getValue()
    for i, j in x:
        if (x[i, j].x > 0.5):
            arcs.append([i, j])
    currentNode = 0
    seq.append(nodeIDs[currentNode])
    while (len(arcs) > 0):
        for i in range(len(arcs)):
            if (arcs[i][0] == currentNode):
                currentNode = arcs[i][1]
                seq.append(nodeIDs[currentNode])
                arcs.pop(i)
                break
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPGurobiQAP(nodeIDs, tau, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
    n = len(nodeIDs)
    TSP = grb.Model('TSP')
    if (outputFlag == False):
        TSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = TSP.addVar(
                    vtype = grb.GRB.BINARY, 
                    name = 'x_%s_%s' % (i, j))
    w = {}
    for i in range(n):
        for j in range(n):
            if (i != j):
                for k in range(n):
                    w[i, j, k] = TSP.addVar(
                        vtype = grb.GRB.CONTINUOUS)

    # TSP objective function ==================================================
    TSP.setObjective(
        grb.quicksum(
            grb.quicksum(
                grb.quicksum(
                    tau[nodeIDs[i], nodeIDs[j]] * w[i, j, k] for k in range(n - 1)
                ) for j in range(n) if j != i
            ) for i in range(n)
        ) + 
        grb.quicksum(
            grb.quicksum(
                tau[nodeIDs[i], nodeIDs[j]] * w[i, j, n - 1] for j in range(n) if j != i
            ) for i in range(n)
        )
    )

    # Assignment constraints ==================================================
    for i in range(n):
        TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if j != i) == 1)
    for j in range(n):
        TSP.addConstr(grb.quicksum(x[i, j] for i in range(n) if j != i) == 1)

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
    TSP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    gap = None
    lb = None
    ub = None
    runtime = None
    solType = None

    ofv = TSP.getObjective().getValue()
    for j in range(n):
        for i in range(n):
            if (i != j and x[i, j].X > 0.9):
                seq.append(nodeIDs[i])
                break
    seq.append(seq[0])
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal
        runtime = TSP.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPCOPTLazyCuts(nodeIDs, tau, outputFlag, timeLimit):
    try:
        import coptpy as cp            
    except(ImportError):
        print("ERROR: Cannot find COPT")
        return
    
    env = cp.Envr(envconfig)
    TSP = env.createModel("TSP")

    # Initialize
    n = len(nodeIDs)
    if (outputFlag != None):
        TSP.setParam(cp.COPT.Param.Logging, outputFlag)
        TSP.setParam(cp.COPT.Param.LogToConsole, outputFlag)
    if (timeLimit != None):
        TSP.setParam(cp.COPT.Param.TimeLimit, timeLimit)

    # Decision variables ======================================================
    x = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                x[i, j] = TSP.addVar(
                    vtype = cp.COPT.BINARY, 
                    obj = tau[i, j], 
                    name = 'x_%s_%s' % (i, j))
                
    # TSP objective function ==================================================
    TSP.ObjSense = cp.COPT.MINIMIZE

    # Degree constraints ======================================================
    for i in nodeIDs:
        TSP.addConstr(cp.quicksum(x[i, j] for j in nodeIDs if i != j) == 1, name = 'leave_%s' % str(i))
        TSP.addConstr(cp.quicksum(x[j, i] for j in nodeIDs if i != j) == 1, name = 'enter_%s' % str(i))

    # Callback ================================================================
    class CoptCallback(cp.CallbackBase):
        def __init__(self):
            super().__init__()
        def callback(self):
            if (self.where() == cp.COPT.CBCONTEXT_MIPSOL):
                # print("Called me!")
                x_sol = self.getSolution(x)
                G = nx.Graph()
                for (i, j) in x:
                    if (x_sol[i, j] > 0.9):
                        G.add_edge(i, j, weight = tau[i, j])
                components = [list(c) for c in nx.connected_components(G)]
                # print(components)
                for component in components:
                    if (len(component) < n):
                        self.addLazyConstr(cp.quicksum(x[i, j] for i in component for j in component if i != j) <= len(component) - 1)
    cb = CoptCallback()

    # TSP with callback =======================================================
    TSP.setCallback(cb, cp.COPT.CBCONTEXT_MIPSOL)
    TSP.solve()

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None
    if (TSP.status == cp.COPT.OPTIMAL):
        solType = 'IP_Optimal'
        ofv = TSP.getObjective().getValue()
        for i, j in x:
            if (x[i, j].x > 0.5):
                arcs.append([i, j])
        currentNode = nodeIDs[0]
        seq.append(currentNode)
        # print(arcs)
        while (len(arcs) > 0):
            for i in range(len(arcs)):
                if (arcs[i][0] == currentNode):
                    currentNode = arcs[i][1]
                    seq.append(currentNode)
                    arcs.pop(i)
                    break
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.SolvingTime
    elif (TSP.status == cp.COPT.TIMEOUT):
        solType = 'IP_TimeLimit'
        ofv = None
        seq = []
        gap = TSP.BestGap
        lb = TSP.BestBnd
        ub = TSP.BestObj
        runtime = TSP.SolvingTime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _heuTSP(nodeObj, tau, depotID, nodeIDs, asymFlag, **kwargs) -> dict|None:

    """Use heuristic kwargss to find suboptimal TSP solution
    
    Parameters
    
    ----------

    nodes: dictionary, required, default None
        The coordinates of given nodes, in the following format::
            >>> nodes = {
            ...     nodeID1: {'loc': (x, y)},
            ...     nodeID2: {'loc': (x, y)}, # ...
            ... }
    edges: dictionary, required, default as {'kwargs': "Euclidean", 'ratio': 1}
        The traveling matrix. The options are as follows::
            1) (default) Euclidean space
            >>> edge = {
            ...     'kwargs': 'Euclidean',
            ...     'ratio': 1 # Optional, default to be 1
            ... }
            2) By given pairs of lat/lon
            >>> edge = {
            ...     'kwargs': 'LatLon',
            ...     'unit': 'meters' # Optional, default to be 1
            ... }
            3) ManhattenDistance
            >>> edge = {
            ...     'kwargs': 'Manhatten',
            ...     'ratio': 1 # Optional, default to be 1
            ... }
            4) By a given dictionary
            >>> edge = {
            ...     'kwargs': 'Dictionary',
            ...     'dictionary': dictionary,
            ...     'ratio': 1 # Optional, default to be 1
            ... }
            5) On the grids
            >>> edge = {
            ...     'kwargs': 'Grid',
            ...     'grid': grid
            ... }
    kwargs: dictionary, required, default as {'cons': 'Insertion', 'impv': '2opt'}
        The algorithm configuration. Includes two phases, use 'cons' to specify constructive heuristic, and 'impv' to specify local improvement heurisitc::
            1) (default) Insertion
            >>> kwargs = {
            ...     'cons': 'Insertion',
            ...     'initSeq': initSeq, # An initial sequence, defalt [depotID]
            ...     'impv': '2Opt' # Options are: 'Reinsert', '2Opt', can select multiple kwargss by collecting them into a list, e.g. ['Reinsert', '2Opt']
            ... }
            2) Nearest neighborhood / k-nearest neighborhood
            >>> kwargs = {
            ...     'cons': 'NearestNeighbor',
            ...     'k': 1, # 1: nearest neighbor, 2 ~ K: k-nearest neighbor, -1: farthest neighbor 
            ... }
            3) Sweep
            >>> kwargs = {
            ...     'cons': 'Sweep'
            ... }
            4) (not available) Christofides
            >>> kwargs = {
            ...     'cons': 'Christofides'
            ... }
            5) (not available) Cycle cover, particular for Asymmetric TSP
            >>> kwargs = {
            ...     'cons': 'CycleCover'
            ... }
            6) Random sequence
            >>> kwargs = {
            ...     'cons': 'Random'
            ... }
            7) Given sequence for further local improvements
            >>> kwargs = {
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

    # Construction heuristics =================================================
    # NOTE: Output of this phase should be a Route() object
    seqObj = Route(tau, asymFlag)
    # An initial solution is given
    if ('cons' not in kwargs or kwargs['cons'] == 'Initial' or kwargs['cons'] == None):
        if ('initSeq' not in kwargs):
            raise MissingParameterError("ERROR: Need 'initSeq' for local improvement")
        elif (len(kwargs['initSeq']) != len(nodeIDs) + 1):
            raise UnsupportedInputError("ERROR: Length of 'initSeq' is incorrect, check if the sequence starts and ends with `depotID`")
        else:
            notInNodeIDs = [v for v in kwargs['initSeq'] if v not in nodeIDs]
            if (len(notInNodeIDs) > 0):
                raise OutOfRangeError("ERROR: The following nodes in 'initSeq' is not in `nodeIDs`: %s" % list2String(notInNodeIDs))
        for i in kwargs['initSeq'][:-1]:
            seqObj.append(nodeObj[i])

    # Insertion heuristic
    elif (kwargs['cons'] == 'Insertion' or kwargs['cons'] == 'RandomInsertion'):
        initSeq = None
        
        randomInsertionFlag = False
        if (kwargs['cons'] == 'RandomInsertion'):
            randomInsertionFlag = True
        
        if ('initSeq' not in kwargs):   
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
            seqObj = _consTSPInsertion(nodeIDs, initSeq, nodeObj, tau, asymFlag, randomInsertionFlag)
        else:
            notInNodeIDs = [v for v in kwargs['initSeq'] if v not in nodeIDs]
            if (len(notInNodeIDs) > 0):
                raise OutOfRangeError("ERROR: The following nodes in 'initSeq' is not in `nodeIDs`: %s" % list2String(notInNodeIDs))
            else:
                seqObj = _consTSPInsertion(nodeIDs, kwargs['initSeq'], tau, asymFlag, randomInsertionFlag)
    
    # Neighborhood based heuristic, including nearest neighborhood, k-nearest neighborhood, and furthest neighborhood
    elif (kwargs['cons'] == 'NearestNeighbor'):
        nnSeq = None
        if ('k' not in kwargs or kwargs['k'] == 1):
            nnSeq = _consTSPkNearestNeighbor(depotID, nodeIDs, tau, 1)
        elif (kwargs['k'] == -1):
            nnSeq = _consTSPFarthestNeighbor(depotID, nodeIDs, tau)
        elif (kwargs['k'] >= 1):
            nnSeq = _consTSPkNearestNeighbor(depotID, nodeIDs, tau, kwargs['k'])
        for i in nnSeq:
            seqObj.append(nodeObj[i])

    # Sweep heuristic
    elif (kwargs['cons'] == 'Sweep'):
        sweepSeq = _consTSPSweep(nodes, depotID, nodeIDs, locFieldName)
        for i in sweepSeq:
            seqObj.append(nodeObj[i])
        seqObj.rehead(depotID)

    # Christofides Algorithm, guaranteed <= 1.5 * optimal
    elif (kwargs['cons'] == 'Christofides'):
        if (not asymFlag):
            cfSeq = _consTSPChristofides(depotID, tau)
            for i in cfSeq:
                seqObj.append(nodeObj[i])
        else:
            raise UnsupportedInputError("ERROR: 'Christofides' algorithm is not designed for Asymmetric TSP")

    # Cycle Cover Algorithm, specially designed for Asymmetric TSP
    elif (kwargs['cons'] == 'CycleCover'):
        raise VrpSolverNotAvailableError("ERROR: 'CycleCover' algorithm is not available yet, please stay tune")
        seqObj = _consTSPCycleCover(depotID, nodeIDs, tau)

    # Randomly create a sequence
    elif (kwargs['cons'] == 'Random'):
        rndSeq = _consTSPRandom(depotID, nodeIDs)
        for i in rndSeq:
            seqObj.append(nodeObj[i])
        seqObj.rehead(depotID)

    else:
        raise UnsupportedInputError(ERROR_MISSING_TSP_ALGO)

    # Cleaning seq before local improving =====================================
    consOfv = seqObj.dist

    # Local improvement phase =================================================
    # NOTE: Local improvement phase operates by class kwargss
    # NOTE: For the local improvement, try every local search operator provided in a greedy way
    if ('impv' in kwargs and kwargs['impv'] != None and kwargs['impv'] != []):
        canImpvFlag = True
        while (canImpvFlag):
            canImpvFlag = False

            # 2Opt
            if (not canImpvFlag and '2Opt' in kwargs['impv']):
                canImpvFlag = seqObj.impv2Opt()

    ofv = seqObj.dist
    seq = [n.key for n in seqObj.traverse(closeFlag = True)]

    return {
        'ofv': ofv,
        'seq': seq
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

def _consTSPSweep(nodes, depotID, nodeIDs, locFieldName):
    # Sweep seq -----------------------------------------------------------
    sweep = nodeSeqBySweeping(
        nodes = nodes, 
        nodeIDs = nodeIDs,
        refLoc = nodes[depotID][locFieldName])
    return sweep

def _consTSPRandom(depotID, nodeIDs):
    # Get random seq ------------------------------------------------------
    seq = [i for i in nodeIDs]
    random.shuffle(seq)
    return seq

def _consTSPInsertion(nodeIDs, initSeq, nodeObj, tau, asymFlag, randomInsertionFlag):
    # Initialize ----------------------------------------------------------
    route = Route(tau, asymFlag)

    # NOTE: initSeq should starts and ends with depotID
    for n in initSeq[:-1]:
        route.append(nodeObj[n])

    unInserted = [i for i in nodeIDs if i not in initSeq[:-1]]
    if (randomInsertionFlag):
        random.shuffle(unInserted)
    for n in unInserted:
        route.cheapestInsert(nodeObj[n])

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
    raise VrpSolverNotAvailableError("ERROR: vrpSolver has not implement this kwargs yet")
