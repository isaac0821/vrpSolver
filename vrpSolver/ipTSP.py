import heapq
import math
import warnings
import networkx as nx

from .common import *
from .const import *
from .msg import *
from .error import *
from .geometry import *

# History =====================================================================
# 20230510 - use networkx to replace .graph
# 20230517 - Implement PEP 3107
# 20230519 - Now support COPT solver for 'DFJ_Lazy' formulation
# 20231022 - Add detail output information for animation
# =============================================================================

def ipTSP(
    nodes: dict, 
    locFieldName: str = 'loc',
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

    """Use IP formulation to find optimal TSP solution

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
    method: dictionary, required, default as {'solver': 'Gurobi', 'timeLimit': None, 'gapTolerance': None, 'outputFlag': False}
        The settings for the MILP solver, right now Gurobi supports all formulation and COPT supports 'DFJ_Lazy', the format is as following:
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
        tau, path = matrixDist(nodes, edges, depotID, nodeIDs, locFieldName)
    else:
        tau, _ = matrixDist(nodes, edges, depotID, nodeIDs, locFieldName)

    # Solve by different formulations =========================================
    tsp = None
    if (method['solver'] == 'Gurobi'):
        if (method['fml'] == 'DFJ_Lazy'):
            tsp = _ipTSPGurobiLazyCuts(
                nodeIDs, 
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None, 
                method['gapTolerance'] if 'gapTolerance' in method else None)
        elif (method['fml'] == 'DFJ_Plainloop'):
            tsp = _ipTSPGurobiPlainLoop(
                nodeIDs, 
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None, 
                method['gapTolerance'] if 'gapTolerance' in method else None)
        elif (method['fml'] == 'MTZ'):
            tsp = _ipTSPGurobiMTZ(
                nodeIDs, 
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None, 
                method['gapTolerance'] if 'gapTolerance' in method else None)
        elif (method['fml'] == 'ShortestPath'):
            tsp = _ipTSPGurobiShortestPath(
                nodeIDs, 
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None, 
                method['gapTolerance'] if 'gapTolerance' in method else None)
        elif (method['fml'] == 'MultiCommodityFlow'):
            tsp = _ipTSPGurobiMultiCommodityFlow(
                nodeIDs, 
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None, 
                method['gapTolerance'] if 'gapTolerance' in method else None)
        elif (method['fml'] == 'QAP'):
            tsp = _ipTSPGurobiQAP(
                nodeIDs, 
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None, 
                method['gapTolerance'] if 'gapTolerance' in method else None)    
    elif (method['solver'] == 'COPT'):
        if (method['fml'] == 'DFJ_Lazy'):
            tsp = _ipTSPCOPTLazyCuts(
                nodeIDs, 
                tau, 
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None,
                method['env'] if 'env' in method else None)
    if (tsp == None):
        raise UnsupportedInputError("ERROR: Incorrect or not available TSP formulation option!")
    
    tsp['fml'] = method['fml']

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
        'gap': tsp['gap'],
        'solType': tsp['solType'],
        'lowerBound': tsp['lowerBound'],
        'upperBound': tsp['upperBound'],
        'runtime': tsp['runtime']
    }
    if (metaFlag):
        res['serviceTime'] = serviceTime
        res['method'] = method
    if (detailsFlag):
        res['vehicles'] = vehicles

    return res

def _ipTSPGurobiLazyCuts(nodeIDs, tau, outputFlag, timeLimit, gapTolerance):
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
                    model.cbLazy(grb.quicksum(x[i,j] for i in component for j in component if i != j) <= len(component) - 1)

    # TSP with callback =======================================================
    TSP.optimize(subtourelim)
    # TSP.write("T.lp")

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    arcs = []
    solType = None
    gap = None
    lb = None
    ub = None
    runtime = None
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
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
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
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
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
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
        gap = 0
        lb = ofv
        ub = ofv
        runtime = accRuntime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
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
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
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
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
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
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
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
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
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
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
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
        gap = 0
        lb = ofv
        ub = ofv
        runtime = TSP.Runtime
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
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
    if (TSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
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
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
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
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipTSPCOPTLazyCuts(nodeIDs, tau, outputFlag, timeLimit, env):
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
        tau, path = matrixDist(nodes, edges, depotID, nodeIDs, locFieldName)
    else:
        tau, _ = matrixDist(nodes, edges, depotID, nodeIDs, locFieldName)

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
