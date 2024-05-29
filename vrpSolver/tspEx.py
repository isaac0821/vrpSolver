import heapq
import math
import warnings
import networkx as nx

from .common import *
from .msg import *
from .geometry import *

def ipTSPEx(
    nodes: dict, 
    locFieldName: str = 'loc',
    predefinedArcs: list[list[tuple[int|str]]] = [],
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    serviceTime: float = 0,
    vehicles: dict = {
        0: {'speed': 1}
    },
    vehicleID: int|str = 0,
    edges: dict = {
        'method': "Euclidean", 
        'ratio': 1
    },
    method: dict = {
        'fml': 'DFJ_Lazy',
        'solver': 'Gurobi',
        'timeLimit': None,
        'outputFlag': None,
        'env': None
    },
    detailsFlag: bool = False,
    metaFlag: bool = False
    ) -> dict|None:

    """Use IP formulation to find optimal TSP solution with extra predefined arcs

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
    predefinedArcs: list[list[tuple[int|str]]], optional, default None
        A set of arcs that are predetermined. For example, ``predefinedArcs = [[(1, 2), (2, 1)], [(3, 4)]]'' means arc (1, 2) or (2, 1) is included in the TSP route, and arc (3, 4) is included as well
    method: dictionary, required, default as {'solver': 'Gurobi', 'timeLimit': None, 'gapTolerance': None, 'outputFlag': False}
        The settings for the MILP method, right now Gurobi supports all formulation and COPT supports 'DFJ_Lazy', the format is as following:
            >>> method = {
            ...     'fml': 'DFJ_Lazy',
            ...     'solver': 'Gurobi',
            ...     'timeLimit': timeLimit, # Time limit in seconds, for 'DFJ_Plainloop' is the total time limit
            ...     'gapTolerance': gapTolerance,
            ...     'outputFlag': False # Turn off method log output by default
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
    if (vehicles == None):
        raise MissingParameterError("ERROR: Missing required field `vehicles`.")
    if (vehicleID not in vehicles):
        raise MissingParameterError("ERROR: Cannot find `vehicleID` in `vehicles`.")
    if (method == None or 'solver' not in method):
        raise MissingParameterError("ERROR: Missing required field `method`.")
    elif (method['solver'] == 'Gurobi' and method['fml'] not in ['DFJ_Lazy', 'DFJ_Plainloop', 'MTZ', 'ShortestPath', 'MultiCommodityFlow', 'QAP']):
        raise OutOfRangeError("ERROR: Gurobi option supports 'DFJ_Lazy', 'DFJ_Plainloop', 'MTZ', 'ShortestPath', 'MultiCommodityFlow', and 'QAP' formulations", )
    elif (method['solver'] == 'COPT' and method['fml'] not in ['DFJ_Lazy']):
        raise OutOfRangeError("ERROR: COPT option supports 'DFJ_Lazy' formulations", )

    # Define tau ==============================================================
    tau = None
    path = None
    if (detailsFlag):
        tau, path = matrixDist(nodes, edges, nodeIDs, locFieldName)
    else:
        tau, _ = matrixDist(nodes, edges, nodeIDs, locFieldName)

    # Solve by different formulations =========================================
    tspEx = None
    if (method['solver'] == 'Gurobi'):
        if (method['fml'] == 'DFJ_Lazy'):
            tspEx = _ipTSPExGurobiLazyCuts(
                nodeIDs, 
                predefinedArcs,
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None, 
                method['gapTolerance'] if 'gapTolerance' in method else None) 
    elif (method['solver'] == 'COPT'):
        if (method['fml'] == 'DFJ_Lazy'):
            tspEx = _ipTSPExCOPTLazyCuts(
                nodeIDs, 
                predefinedArcs,
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None,
                method['env'] if 'env' in method else None)
    if (tspEx == None):
        raise UnsupportedInputError("ERROR: Incorrect or not available TSP formulation option!")
    
    tspEx['fml'] = method['fml']

    # Fix the sequence to make it start from the depot ========================
    startIndex = 0
    seq = [i for i in tspEx['seq']]
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
    tspEx['seq'] = nodeSeq

    # Add service time if provided ============================================
    ofv = tspEx['ofv'] + (len(nodeIDs) - 1) * serviceTime

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
            if (edges['method'] in ['Euclidean', 'LatLon']):
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
        if (edges['method'] in ['Euclidean', 'LatLon']):
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
        'gap': tspEx['gap'],
        'solType': tspEx['solType'],
        'lowerBound': tspEx['lowerBound'],
        'upperBound': tspEx['upperBound'],
        'runtime': tspEx['runtime']
    }
    if (metaFlag):
        res['serviceTime'] = serviceTime
        res['method'] = method
    if (detailsFlag):
        res['vehicles'] = vehicles

    return res

def _ipTSPExGurobiLazyCuts(nodeIDs, predefinedArcs, tau, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
    n = len(nodeIDs)
    TSPEx = grb.Model('TSPEx')
    if (outputFlag == False):
        TSPEx.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSPEx.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSPEx.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                x[i, j] = TSPEx.addVar(
                    vtype = grb.GRB.BINARY, 
                    obj = tau[i, j], 
                    name = 'x_%s_%s' % (i, j))
                
    # TSPEx objective function ================================================
    TSPEx.modelSense = grb.GRB.MINIMIZE
    TSPEx.Params.lazyConstraints = 1
    TSPEx.update()

    # Degree constraints ======================================================
    for i in nodeIDs:
        TSPEx.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) == 1, name = 'leave_%s' % str(i))
        TSPEx.addConstr(grb.quicksum(x[j, i] for j in nodeIDs if i != j) == 1, name = 'enter_%s' % str(i))

    # Predefined arcs =========================================================
    for arcs in predefinedArcs:
        for arc in arcs:
            if ((arc[0], arc[1]) not in x or (arc[1], arc[0]) not in x):
                raise UnsupportedInputError("ERROR: Cannot find arc[%s, %s] and/or arc[%s, %s]" % (arc[0], arc[1], arc[1], arc[0]))
        TSPEx.addConstr(grb.quicksum(x[arc[0], arc[1]] for arc in arcs) == 1)

    # Sub-tour elimination ====================================================
    TSPEx._x = x
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
                    model.cbLazy(grb.quicksum(x[i,j] for i in component for j in component if i != j) <= len(component) - 1)

    # TSPEx with callback =======================================================
    TSPEx.optimize(subtourelim)

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None
    if (TSPEx.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        ofv = TSPEx.getObjective().getValue()
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
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSPEx.Runtime
    elif (TSPEx.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        ofv = None
        seq = []
        gap = TSPEx.MIPGap
        lb = TSPEx.ObjBoundC
        ub = TSPEx.ObjVal
        runtime = TSPEx.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPExGurobiMTZ(nodeIDs, predefinedArcs, tau, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
    n = len(nodeIDs)
    TSPEx = grb.Model('TSPEx')
    if (outputFlag == False):
        TSPEx.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSPEx.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSPEx.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = TSPEx.addVar(
                    vtype = grb.GRB.BINARY, 
                    obj = tau[nodeIDs[i], nodeIDs[j]], 
                    name = 'x_%s_%s' % (i, j))
    u = {}
    for i in range(n):
        u[i] = TSPEx.addVar(
            vtype = grb.GRB.CONTINUOUS,
            name = 'u_%s' % (i))

    # TSPEx objective function ================================================
    TSPEx.modelSense = grb.GRB.MINIMIZE
    TSPEx.update()

    # Degree constraints ======================================================
    for i in range(n):
        TSPEx.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % str(i))
        TSPEx.addConstr(grb.quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % str(i))

    # Predefined arcs =========================================================
    for arcs in predefinedArcs:
        for arc in arcs:
            if ((arc[0], arc[1]) not in x or (arc[1], arc[0]) not in x):
                raise UnsupportedInputError("ERROR: Cannot find arc[%s, %s] and/or arc[%s, %s]" % (arc[0], arc[1], arc[1], arc[0]))
        TSPEx.addConstr(grb.quicksum(x[arc[0], arc[1]] for arc in arcs) == 1, name = 'predefind_%s_%s' % list2String(arcs))

    # Sequence constraints ====================================================
    for i in range(1, n):
        for j in range(1, n):
            if (i != j):
                TSPEx.addConstr(u[i] - u[j] + (n - 1) * x[i, j] <= n - 2, name = 'seq_%s_%s' % (i, j))
    for i in range(1, n):
        TSPEx.addConstr(1 <= u[i])
        TSPEx.addConstr(u[i] <= n - 1)

    # TSPEx ===================================================================
    TSPEx.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    gap = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None
    if (TSPEx.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        ofv = TSPEx.getObjective().getValue()
        gap = TSPEx.Params.MIPGapAbs
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
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSPEx.Runtime
    elif (TSPEx.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        ofv = None
        seq = []
        gap = TSPEx.MIPGap
        lb = TSPEx.ObjBoundC
        ub = TSPEx.ObjVal
        runtime = TSPEx.Runtime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPExCOPTLazyCuts(nodeIDs, predefinedArcs, tau, outputFlag, timeLimit, env):
    try:
        import coptpy as cp
        envconfig = cp.EnvrConfig()
        envconfig.set('nobanner', '1')
        AVAIL_SOLVER = 'COPT'
        if (env == None):
            env = cp.Envr(envconfig)
    except(ImportError):
        print("ERROR: Cannot find COPT")
        return

    TSPEx = env.createModel("TSPEx")

    # Initialize
    n = len(nodeIDs)
    if (outputFlag != None):
        TSPEx.setParam(cp.COPT.Param.Logging, outputFlag)
        TSPEx.setParam(cp.COPT.Param.LogToConsole, outputFlag)
    if (timeLimit != None):
        TSPEx.setParam(cp.COPT.Param.TimeLimit, timeLimit)

    # Decision variables ======================================================
    x = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                x[i, j] = TSPEx.addVar(
                    vtype = cp.COPT.BINARY, 
                    obj = tau[i, j], 
                    name = 'x_%s_%s' % (i, j))
                
    # TSPEx objective function ==================================================
    TSPEx.ObjSense = cp.COPT.MINIMIZE

    # Degree constraints ======================================================
    for i in nodeIDs:
        TSPEx.addConstr(cp.quicksum(x[i, j] for j in nodeIDs if i != j) == 1, name = 'leave_%s' % str(i))
        TSPEx.addConstr(cp.quicksum(x[j, i] for j in nodeIDs if i != j) == 1, name = 'enter_%s' % str(i))

    # Predefined arcs =========================================================
    for arcs in predefinedArcs:
        for arc in arcs:
            if ((arc[0], arc[1]) not in x or (arc[1], arc[0]) not in x):
                raise UnsupportedInputError("ERROR: Cannot find arc[%s, %s] and/or arc[%s, %s]" % (arc[0], arc[1], arc[1], arc[0]))
        TSPEx.addConstr(cp.quicksum(x[arc[0], arc[1]] for arc in arcs) == 1, name = 'predefind_%s_%s' % (arc[0], arc[1]))

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

    # TSPEx with callback =======================================================
    TSPEx.setCallback(cb, cp.COPT.CBCONTEXT_MIPSOL)
    TSPEx.solve()

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None
    if (TSPEx.status == cp.COPT.OPTIMAL):
        solType = 'IP_Optimal'
        ofv = TSPEx.getObjective().getValue()
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
        runtime = TSPEx.SolvingTime
    elif (TSPEx.status == cp.COPT.TIMEOUT):
        solType = 'IP_TimeLimit'
        ofv = None
        seq = []
        gap = TSPEx.BestGap
        lb = TSPEx.BestBnd
        ub = TSPEx.BestObj
        runtime = TSPEx.SolvingTime

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def heuTSPEx(
    nodes: dict, 
    predefinedArcs: list[list[tuple[int|str]]] = [],
    edges: dict = {'method': "Euclidean", 'ratio': 1},
    method: dict = {'cons': 'NearestNeighbor', 'impv': '2Opt'},
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    serviceTime: float = 0,
    ) -> dict|None:

    # FIXME: TO BE REWRITEN
    
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


    if (method == None or 'cons' not in method or 'impv' not in method):
        raise MissingParameterError("ERROR: Not supported (for now)")

    seq = []
    # Neighborhood based heuristic, including nearest neighborhood, k-nearest neighborhood, and furthest neighborhood
    if (method['cons'] == 'NearestNeighbor'):
        seq = _consTSPExkNearestNeighbor(depotID, nodeIDs, tau, mustGo, 1)

    # Cleaning seq before local improving =====================================
    ofv = calSeqCostMatrix(tau, seq, closeFlag = False)
    revOfv = None if not asymFlag else calSeqCostMatrix(tau, [seq[len(seq) - i - 1] for i in range(len(seq))], closeFlag = False)
    consOfv = ofv

    # Local improvement phase =================================================
    # NOTE: For the local improvement, try every local search operator provided in a greedy way
    if ('impv' in method and method['impv'] != None and method['impv'] != []):
        canImproveFlag = True
        while (canImproveFlag):
            canImproveFlag = False
            # Try 2Opts
            if (not canImproveFlag and (method['impv'] == '2Opt' or '2Opt' in method['impv'])):
                imp = _impTSPEx2Opts(nodeIDs, tau, seq, asymFlag, mustGo)
                if (imp['improvedFlag']):
                    canImproveFlag = True
                    seq = imp['impSeq']
                    ofv = imp['oriOfv']
                    revOfv = imp['oriRevOfv']
    return {
        'ofv': ofv,
        'consOfv': consOfv,
        'method': method,
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