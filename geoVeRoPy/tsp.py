import heapq
import math
import warnings
import networkx as nx

from .common import *
from .geometry import *
from .msg import *
from .travel import *

def solveTSP(
    nodes: dict, 
    locFieldName: str = 'loc',
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    vehicles: dict = {0: {'speed': 1}},
    vehicleID: int|str = 0,
    serviceTime: float = 0,
    edges: str = 'Euclidean',
    algo: str = 'IP',
    detailFlag: bool = False,
    metaFlag: bool = False,
    **kwargs
    ) -> dict:


    """Use MIP/Heuristic to find optimal/sub-optimal TSP solution

    Parameters
    ----------

    nodes: dict, required
        The `nodes` dictionary, with coordinates of given nodes. See :ref:`nodes`
    locFieldName: str, optional, default as 'loc'
        The key value in `nodes` indicating the location of each node.
    depotID: int|string, optional, default as 0
        The node ID of the depot.
    nodeIDs: string 'All' or a list of node IDs, optional, default as 'All'
        The following are two options: 1) 'All', all nodes will be visited, 2) A list of node IDs to be visited.
    vehicles: dict, optional, default as {0: {'speed': 1}}
        The vehicle information, since the TSP optimizes route for one vehicle, only the first vehicle with `vehicleID` will be considered. The speed of vehicle needs to be specified.
    vehicleID: int|str, optional, default as 0
        The vehicle that takes the TSP route in `vehicles`.
    serviceTime: float, optional, default as 0
        The service time needed at each location. If `serviceTime` is provided in `nodes`, it will be ignored.
    edges: string, optional, default as 'Euclidean'
        The methods for the calculation of distances between nodes. Options and required additional information are referred to :func:`~vrpSolver.geometry.matrixDist()`.
    algo: string, optional, default as 'IP'
        Select the algorithm for calculating TSP. Options and required additional inputs are as follows:
        
        1) (default) 'IP', use (Mixed) Integer Programming to solve TSP. Needs commercial solver such as Gurobi and COPT.
            - solver: str, supports 'Gurobi' and 'COPT'.
            - fml: str, choose formulation of TSP:
                - if use 'Gurobi' as solver, fml supports the following options: ['DFJ_Lazy', 'DFJ_Plainloop', 'MTZ', 'ShortestPath', 'MultiCommodityFlow', 'QAP']. In fact, 'DFJ_Lazy' will be faster than all other formulations, implementation of different formulation is for education purpose.
                - if use 'COPT' as solver, fml supports 'DFJ_Lazy'.
            - timeLimit: int|float, additional stopping criteria
            - gapTolerance: int|float, additional stopping criteria
            - outputFlag: bool, True if turn on the log output from solver. Default to be False
        2) 'Heuristic', use different 2-phase heuristic methods to solve TSP to sub-optimal. 
            - cons: str, choose a construction heuristic to find a feasible TSP solution, options and addition inputs needed are as follows:
                - cons = 'Initial', use a given initial solution, requires the following information:
                    - initSeq: a list of node IDs
                - cons = 'Insertion', in each iteration, insert a node into the existing route, which has the minimal cost.
                - cons = 'RandomInsertion', in each iteration, randomly insert a not that has not yet been inserted to the route.
                - cons = 'NearestNeighbor', in each iteration, find the (k-)nearest node to the end of existing route and add it to the route.
                    - k: the k-th nearest neighbor, default to be 1
                - cons = 'Sweep', sweep all nodes clock-wise (or counter-clock-wise) and add the node to the route one by one.
                - cons = 'Christofides', use Christofides algorithm
                - cons = 'CycleCover', use CycleCover algorithm, particularly design for Asymmetric TSP
                - cons = 'Random', randomly create a feasible route
            - impv: str, choose the local improvement heuristic to improve the existing feasible route
                - impv = '2Opt', use the 2-opt algorithm
        3) 'Metaheuristic', use metaheuristic methods to solve TSP to sub-optimal, different metaheuristic methods requires different construction phase of heuristic
            - cons: str, construction heuristic for metaheuristic, options depends on `meta`
            - meta: str, choose a metaheuristic improvement method
                - meta = 'SimulatedAnnealing', use Simulated Annealing to improve a solution create by 'cons', choice of 'cons' are all construction heuristic available for 'Heuristic'
                - meta = 'GeneticAlgorithm', use Genetic Algorithm to create solutions. Choice of 'cons' includes ['Random', 'RandomInsertion']

    detailFlag: bool, optional, default as False
        If True, an additional field `vehicle` will be added into the solution, which can be used for animation. See :func:`~vrpSolver.plot.aniRouting()`.
    **kwargs: optional
        Provide additional inputs for different `edges` options and `algo` options

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

    if (detailFlag == True):
        # For animation propose
        if (vehicles == None):
            raise MissingParameterError("ERROR: Missing required field `vehicles`.")
        if (vehicleID not in vehicles):
            raise MissingParameterError("ERROR: Cannot find `vehicleID` in `vehicles`.")

    if (algo == 'IP'):
        if ('solver' not in kwargs and 'fml' not in kwargs):
            kwargs['solver'] = 'Gurobi'
            kwargs['fml'] = 'DFJ_Lazy'
            warnings.warn("WARNING: Missing required field `solver`, set to default 'Gurobi' with DFJ + lazy cuts")
        elif (kwargs['solver'] == 'Gurobi' and kwargs['fml'] not in ['DFJ_Lazy', 'DFJ_Plainloop', 'MTZ', 'ShortestPath', 'MultiCommodityFlow', 'QAP']):
            raise OutOfRangeError("ERROR: Gurobi option s upports 'DFJ_Lazy', 'DFJ_Plainloop', 'MTZ', 'ShortestPath', 'MultiCommodityFlow', and 'QAP' formulations", )
        elif (kwargs['solver'] == 'COPT' and kwargs['fml'] not in ['DFJ_Lazy']):
            raise OutOfRangeError("ERROR: COPT option supports 'DFJ_Lazy' formulations", )
    elif (algo == 'Heuristic'):
        if ('cons' not in kwargs):
            kwargs['cons'] = 'NearestNeighbor'
            warnings.warn("WARNING: No construction heuristic are specified.")
        if ('impv' not in kwargs):
            kwargs['impv'] = '2Opt'
            warnings.warn("WARNING: No local improvement heuristic are specified.")
    elif (algo == 'Metaheuristic'):
        if ('cons' not in kwargs):
            kwargs['cons'] = 'NearestNeighbor'
            warnings.warn("WARNING: No construction heuristic are specified.")
        if ('meta' not in kwargs):
            kwargs['meta'] = 'SimulatedAnnealing'
            warnings.warn("WARNING: No metaheuristic are specified.")

    # Define tau ==============================================================
    tau = None
    pathLoc = None
    if (detailFlag):
        res = matrixDist(
            nodes = nodes, 
            nodeIDs = nodeIDs,
            edges = edges, 
            locFieldName = locFieldName,
            detailFlag = True
            **kwargs)
        tau = res['tau']
        pathLoc = res['pathLoc']
    else:
        tau = matrixDist(
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
    startTime = datetime.datetime.now()

    tsp = None
    if (algo == 'IP'):
        tsp = _ipTSP(
            nodeIDs = nodeIDs, 
            tau = tau, 
            solver = kwargs['solver'], 
            fml = kwargs['fml'], 
            outputFlag = False if 'outputFlag' not in kwargs else kwargs['outputFlag'], 
            timeLimit = None if 'timeLimit' not in kwargs else kwargs['timeLimit'], 
            gapTolerance = None if 'gapTolerance' not in kwargs else kwargs['gapTolerance'])
        tsp['fml'] = kwargs['fml']
        tsp['solver'] = kwargs['solver']

    elif (algo == 'Heuristic'):
        nodeObj = {}
        for n in nodeIDs:
            nodeObj[n] = RouteNode(n, value=nodes[n][locFieldName])
        # Two-phase: construction phase
        seqObj = _consTSP(
            nodes = nodes,
            nodeObj = nodeObj, 
            tau = tau, 
            depotID = depotID, 
            nodeIDs = nodeIDs, 
            asymFlag = asymFlag, 
            **kwargs)
        # Two-phase: local improvement phase
        tsp = _impvTSP(
            seqObj = seqObj,
            **kwargs)
        tsp['cons'] = kwargs['cons']
        tsp['impv'] = kwargs['impv']

    elif (algo == 'Metaheuristic'):
        nodeObj = {}
        for n in nodeIDs:
            nodeObj[n] = RouteNode(n, value=nodes[n][locFieldName])

        # Population based
        if (kwargs['meta'] in ['GeneticAlgorithm', 'PSO', 'ACO']):
            # Two-phase: popularize phase
            popObj = _popTSP(
                nodes = nodes,
                nodeObj = nodeObj, 
                tau = tau, 
                depotID = depotID, 
                nodeIDs = nodeIDs, 
                asymFlag = asymFlag, 
                **kwargs)
            # Two-phase: population based search phase
            tsp = _metaPopSearchTSP(
                popObj = popObj,
                nodeObj = nodeObj,
                tau = tau,
                depotID = depotID,
                asymFlag = asymFlag,
                **kwargs)
            tsp['cons'] = kwargs['cons']
            tsp['meta'] = kwargs['meta']
        # Search based
        elif (kwargs['meta'] in ['SimulatedAnnealing', 'ALNS', 'GRASP']):
            # Two-phase: construction phase
            seqObj = _consTSP(
                nodes = nodes,
                nodeObj = nodeObj, 
                tau = tau, 
                depotID = depotID, 
                nodeIDs = nodeIDs, 
                asymFlag = asymFlag, 
                **kwargs)
            # Two-phase: local search phase
            tsp = _metaLocalSearchTSP( 
                seqObj = seqObj, 
                **kwargs)
            tsp['cons'] = kwargs['cons']
            tsp['meta'] = kwargs['meta']

    else:
        raise OutOfRangeError("ERROR: Select 'algo' from ['IP', 'Heuristic', 'Metaheuristic'].")

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
    hasServiceTimeInfoFlag = False
    sumServiceTime = 0
    for n in nodeIDs:
        if ('serviceTime' in nodes[n]):
            sumServiceTime += nodes[n]['serviceTime']
            hasServiceTimeInfoFlag = True
    if (not hasServiceTimeInfoFlag):
        ofv = tsp['ofv'] + (len(nodeIDs) - 1) * serviceTime
    else:
        ofv = tsp['ofv'] + sumServiceTime

    # Post optimization (for detail information) ==============================
    if (detailFlag):
        # 返回一个数组，表示路径中的每个点的位置，不包括时间信息
        shapepoints = []        
        for i in range(len(nodeSeq) - 1):
            shapepoints.extend(pathLoc[nodeSeq[i], nodeSeq[i + 1]][:-1])
        shapepoints.append(pathLoc[nodeSeq[-2], nodeSeq[-1]][-1])

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
                shapepointsInBtw = pathLoc[nodeSeq[i - 1], nodeSeq[i]]
                for j in range(1, len(shapepointsInBtw)):
                    curTime += distEuclideanXY(shapepointsInBtw[j - 1], shapepointsInBtw[j]) / vehicles[vehicleID]['speed']
                    curLoc = shapepointsInBtw[j]
                    timedSeq.append((curLoc, curTime))
            # 如果有service time，则加上一段在原处等待的时间
            if ('serviceTime' in nodes[nodeSeq[i]]):
                curTime += nodes[nodeSeq[i]]['serviceTime']
                timedSeq.append((curLoc, curTime))
            elif (serviceTime != None and serviceTime > 0):
                curTime += serviceTime
                # curLoc = curLoc
                timedSeq.append((curLoc, curTime))
        # 现在补上最后一段leg
        if (edges in ['Euclidean', 'LatLon']):
            curTime += tau[nodeSeq[-2], nodeSeq[-1]] / vehicles[vehicleID]['speed']
            curLoc = nodes[nodeSeq[-1]][locFieldName]
            timedSeq.append((curLoc, curTime))
        else:
            shapepointsInBtw = pathLoc[nodeSeq[-2], nodeSeq[-1]]
            for j in range(1, len(shapepointsInBtw)):
                curTime += distEuclideanXY(shapepointsInBtw[j - 1], shapepointsInBtw[j]) / vehicles[vehicleID]['speed']
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
    if (metaFlag):
        res['algo'] = algo
        res['serviceTime'] = serviceTime
    if (detailFlag):
        res['vehicles'] = vehicles

    res['runtime'] = (datetime.datetime.now() - startTime).total_seconds()

    return res

def _ipTSP(nodeIDs, tau, solver, fml, outputFlag, timeLimit, gapTolerance) -> dict|None:
    tsp = None
    if (solver == 'Gurobi'):
        if (fml == 'DFJ_Lazy'):
            tsp = _ipTSPGurobiLazyCuts(nodeIDs, tau, outputFlag, timeLimit, gapTolerance)
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

def _ipTSPGurobiLazyCuts(nodeIDs, tau, outputFlag, timeLimit, gapTolerance):
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
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub
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
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
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
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub
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
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
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
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
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
    elif (TSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = TSP.MIPGap
        lb = TSP.ObjBoundC
        ub = TSP.ObjVal

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
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
    elif (TSP.status == cp.COPT.TIMEOUT):
        solType = 'IP_TimeLimit'
        ofv = None
        seq = []
        gap = TSP.BestGap
        lb = TSP.BestBnd
        ub = TSP.BestObj

    return {
        'ofv': ofv,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
    }

def _consTSP(nodes, nodeObj, tau, depotID, nodeIDs, asymFlag, **kwargs) -> "Route":
    # Construction heuristics =================================================
    # NOTE: Output of this phase should be a Route() object
    seqObj = Route(tau, asymFlag)
    cons = kwargs['cons']
    # An initial solution is given
    if (cons == 'Initial'):
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
    elif (cons == 'Insertion' or cons == 'RandomInsertion'):
        initSeq = None
        
        randomInsertionFlag = False
        if (cons == 'RandomInsertion'):
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
    elif (cons == 'NearestNeighbor'):
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
    elif (cons == 'Sweep'):
        sweepSeq = _consTSPSweep(nodes, depotID, nodeIDs, locFieldName)
        for i in sweepSeq:
            seqObj.append(nodeObj[i])
        seqObj.rehead(depotID)

    # Christofides Algorithm, guaranteed <= 1.5 * optimal
    elif (cons == 'Christofides'):
        if (not asymFlag):
            cfSeq = _consTSPChristofides(nodes, depotID, tau)
            for i in cfSeq:
                seqObj.append(nodeObj[i])
        else:
            raise UnsupportedInputError("ERROR: 'Christofides' algorithm is not designed for Asymmetric TSP")

   # Cycle Cover Algorithm, specially designed for Asymmetric TSP
    elif (cons == 'CycleCover'):
        raise VrpSolverNotAvailableError("ERROR: 'CycleCover' algorithm is not available yet, please stay tune")
        seqObj = _consTSPCycleCover(depotID, nodeIDs, tau)

    # Randomly create a sequence
    elif (cons == 'Random'):
        rndSeq = _consTSPRandom(depotID, nodeIDs)
        for i in rndSeq:
            seqObj.append(nodeObj[i])
        seqObj.rehead(depotID)

    else:
        raise UnsupportedInputError(ERROR_MISSING_TSP_ALGO)

    # Cleaning seq before local improving =====================================
    seqObj.rehead(depotID)
    return seqObj

def _popTSP(popSize, nodes, nodeObj, tau, depotID, nodeIDs, asymFlag, **kwargs) -> list["Route"]:
    pop = []
    if (kwargs['cons'] in ['Random', 'RandomInsertion']):    
        for i in range(popSize):
            seqObj = _consTSP(nodes, nodeObj, tau, depotID, nodeIDs, asymFlag, **kwargs)
            pop.append(seqObj.clone())
    else:
        raise UnsupportedInputError("ERROR: Population based metaheuristic only supports 'Random' or 'RandomInsertion' as construction heuristic.")
    return pop

def _impvTSP(seqObj, **kwargs):
    # NOTE: For the local improvement, try every local search operator provided in a greedy way
    canImpvFlag = True
    while (canImpvFlag):
        canImpvFlag = False

        # 2Opt
        if (not canImpvFlag and kwargs['impv'] == '2Opt' or '2Opt' in kwargs['impv']):
            canImpvFlag = seqObj.impv2Opt()

    ofv = seqObj.dist
    seq = [n.key for n in seqObj.traverse(closeFlag = True)]

    return {
        'ofv': ofv,
        'seq': seq
    }

def _metaLocalSearchTSP(seqObj, **kwargs):
    if (kwargs['meta'] == 'SimulatedAnnealing'):
        seqObj = _metaTSPSimulatedAnnealing(
            seqObj = seqObj, 
            initTemp = kwargs['initTemp'], 
            lengTemp = kwargs['lengTemp'], 
            neighRatio = kwargs['neighRatio'], 
            coolRate = kwargs['coolRate'], 
            stop = kwargs['stop'])
    else:
        raise UnsupportedInputError("ERROR: Currently not support.")

    ofv = seqObj.dist
    seq = [n.key for n in seqObj.traverse(closeFlag = True)]

    return {
        'ofv': ofv,
        'seq': seq
    }

def _metaPopSearchTSP(popObj, nodeObj, depotID, tau, asymFlag, **kwargs):
    if (kwargs['meta'] == 'GeneticAlgorithm'):
        seqObj = _metaTSPGeneticAlgorithm(
            popObj = popObj, 
            nodeObj = nodeObj, 
            depotID = depotID, 
            tau = tau, 
            asymFlag = asymFlag, 
            neighRatio = kwargs['neighRatio'], 
            stop = kwargs['stop'])
    else:
        raise UnsupportedInputError("ERROR: Currently not support.")

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
    
def _consTSPChristofides(nodes, depotID, tau):
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

def _metaTSPSimulatedAnnealing(seqObj, initTemp, lengTemp, neighRatio, coolRate, stop) -> dict:

    # Subroutines to generate neighborhoods ===================================
    # Swap i and i + 1
    def swap(keyI):
        oldDist = seqObj.dist
        nI = seqObj.query(keyI)
        seqObj.swap(nI)
        newDist = seqObj.dist
        return (newDist - oldDist)

    # Randomly exchange two vertices
    def exchange(keyI, keyJ):
        oldDist = seqObj.dist
        nI = seqObj.query(keyI)
        nJ = seqObj.query(keyJ)
        seqObj.exchange(nI, nJ)
        newDist = seqObj.dist
        return (newDist - oldDist)
        
    # Randomly rotate part of seq
    def rotate(keyI, keyJ):
        oldDist = seqObj.dist
        nI = seqObj.query(keyI)
        nJ = seqObj.query(keyJ)
        seqObj.rotate(nI, nJ)
        newDist = seqObj.dist
        return (newDist - oldDist)

    # Initialize ==============================================================
    # Initial temperature
    T = initTemp
    # Temperature length (maximum temperature iteration)
    L = lengTemp

    # Initial Solution
    curSeq = [i.key for i in seqObj.traverse()]
    ofv = seqObj.dist
    startTime = datetime.datetime.now()

    # Main cooling ============================================================
    contFlag = True

    iterTotal = 0
    iterNoImp = 0
    iterAcc = 0

    while (contFlag):
        # Repeat in the same temperature
        for l in range(L):
            # Increment iterator
            iterTotal += 1

            # Generate a neighbor using different type
            typeOfNeigh = rndPickFromDict(neighRatio)

            deltaC = None
            revAction = {}

            # Randomly swap
            if (typeOfNeigh == 'swap'):
                curSeq = [f.key for f in seqObj.traverse()]
                i = random.randint(0, len(curSeq) - 1)
                keyI = curSeq[i]
                keyINext = seqObj.query(keyI).next.key
                revAction = {
                    'opt': 'swap',
                    'key': keyINext
                }
                deltaC = swap(keyI)

            # Randomly exchange two digits
            elif (typeOfNeigh == 'exchange'):
                curSeq = [f.key for f in seqObj.traverse()]
                i = None
                j = None
                while (i == None 
                        or j == None 
                        or abs(i - j) <= 2 
                        or (i == 0 and j == len(curSeq) - 1) 
                        or (i == len(curSeq) - 1 and j == 0)):
                    i = random.randint(0, len(curSeq) - 1)
                    j = random.randint(0, len(curSeq) - 1)
                keyI = curSeq[i]
                keyJ = curSeq[j]
                revAction = {
                    'opt': 'exchange',
                    'key': (keyJ, keyI)
                }
                nI = seqObj.query(keyI)
                nJ = seqObj.query(keyJ)
                if (nI.next.key == nJ.key or nJ.next.key == nI.key):
                    raise
                deltaC = exchange(keyI, keyJ)

            # Randomly reverse part of path
            elif (typeOfNeigh == 'rotate'):
                curSeq = [f.key for f in seqObj.traverse()]
                i = None
                j = None
                while (i == None 
                        or j == None 
                        or abs(i - j) <= 2 
                        or (i == 0 and j == len(curSeq) - 1) 
                        or (i == len(curSeq) - 1 and j == 0)):
                    i = random.randint(0, len(curSeq) - 1)
                    j = random.randint(0, len(curSeq) - 1)
                keyI = curSeq[i]
                keyJ = curSeq[j]

                revAction = {
                    'opt': 'rotate',
                    'key': (keyJ, keyI)
                }                
                nI = seqObj.query(keyI)
                nJ = seqObj.query(keyJ)
                if (nI.next.key == nJ.key or nJ.next.key == nI.key):
                    raise
                deltaC = rotate(keyI, keyJ)

            # If this new neighbor is good, accept it, 
            #     otherwise accept it with probability
            if (deltaC <= 0): # deltaC = newC - preC, <0 means improve
                # print("Improved: Accept")
                ofv = seqObj.dist
                iterAcc += 1
                iterNoImp = 0
            else:
                sample = random.random()
                if (sample < math.exp(- deltaC / T)):
                    # print("No improve: Accept by chance", sample, "<", math.exp(- deltaC / T))
                    ofv = seqObj.dist
                    iterAcc += 1
                else:
                    # print("No improve: Refused.")
                    beforeSeq = seqObj.traverse()
                    if (revAction['opt'] == 'swap'):
                        _ = swap(revAction['key'])
                    elif (revAction['opt'] == 'exchange'):
                        _ = exchange(revAction['key'][0], revAction['key'][1])
                    elif (revAction['opt'] == 'rotate'):
                        _ = rotate(revAction['key'][0], revAction['key'][1])
                    iterNoImp += 1

            apRate = iterAcc / iterTotal

            # Check stopping criteria
            endCriteria = None
            if ('finalTemp' in stop):
                if (T < stop['finalTemp']):
                    contFlag = False
                    break
            if ('numNoImproveIter' in stop):
                if (iterNoImp > stop['numNoImproveIter']):
                    contFlag = False
                    break
            if ('ptgAcceptedMove' in stop):
                if (iterTotal > 0 and apRate < stop['ptgAcceptedMove']):
                    contFlag = False
                    break
            if ('numIter' in stop):
                if (iterTotal > stop['numIter']):
                    contFlag = False
                    break
            if ('runtime' in stop):
                if ((datetime.datetime.now() - startTime).total_seconds() > stop['runtime']):
                    contFlag = False
                    break
        # Cool down
        T = coolRate * T
    return seqObj

def _metaTSPGeneticAlgorithm(popObj, nodeObj, depotID, tau, asymFlag, neighRatio, stop) -> dict:

    # Subroutines to generate neighborhoods ===================================
    # Swap i and i + 1
    def swap(seqObj, keyI):
        nI = seqObj.query(keyI)
        seqObj.swap(nI)
        return seqObj

    # Randomly exchange two vertices
    def exchange(seqObj, keyI, keyJ):
        nI = seqObj.query(keyI)
        nJ = seqObj.query(keyJ)
        seqObj.exchange(nI, nJ)
        return seqObj
        
    # Randomly rotate part of seq
    def rotate(seqObj, keyI, keyJ):
        nI = seqObj.query(keyI)
        nJ = seqObj.query(keyJ)
        seqObj.rotate(nI, nJ)
        return seqObj

    # Randomly crossover two sequence
    def crossover(seqObj1, seqObj2, idx1, idx2):
        # 原始序列
        seq1 = [i.key for i in seqObj1.traverse()]
        seq2 = [i.key for i in seqObj2.traverse()]
        # 把idx1和idx2排个序换一下，保证idx1在前面
        if (idx1 > idx2):
            idx1, idx2 = idx2, idx1
        # 构造新序列
        newSeq1 = [seq2[i] for i in range(idx1, idx2)]
        for i in range(idx2, len(seq1)):
            if (seq1[i] not in newSeq1):
                newSeq1.append(seq1[i])
        for i in range(idx2):
            if (seq1[i] not in newSeq1):
                newSeq1.append(seq1[i])
        newSeq2 = [seq1[i] for i in range(idx1, idx2)]
        for i in range(idx2, len(seq2)):
            if (seq2[i] not in newSeq2):
                newSeq2.append(seq2[i])
        for i in range(idx2):
            if (seq2[i] not in newSeq2):
                newSeq2.append(seq2[i])
        # 构造新序列实体
        newSeqObj1 = Route(tau, asymFlag)
        for i in newSeq1:
            newSeqObj1.append(nodeObj[i].clone())
        newSeqObj1.rehead(depotID)
        newSeqObj2 = Route(tau, asymFlag)
        for i in newSeq2:
            newSeqObj2.append(nodeObj[i].clone())
        newSeqObj2.rehead(depotID)

        return newSeqObj1, newSeqObj2

    # Initialize ==============================================================
    dashboard = {
        'bestOfv': float('inf'),
        'bestSeq': None
    }
    startTime = datetime.datetime.now()
    for seq in popObj:
        if (seq.dist < dashboard['bestOfv']):
            dashboard['bestOfv'] = seq.dist
            dashboard['bestSeq'] = [i.key for i in seq.traverse()]

    popSize = len(popObj)
    geneLen = len(popObj[0].traverse())

    # Main loop ===============================================================
    contFlag = True

    iterTotal = 0
    iterNoImp = 0   

    while (contFlag):
        # Crossover and create offspring
        while (len(popObj) <= (int)((1 + neighRatio['crossover']) * popSize)):
            # Randomly select two genes, the better gene has better chance to have offspring
            
            rnd1 = None
            rnd2 = None
            while (rnd1 == None or rnd2 == None or rnd1 == rnd2):
                coeff = []
                for i in range(len(popObj)):
                    coeff.append(1 / (popObj[i].dist - 0.8 * dashboard['bestOfv']))
                rnd1 = rndPick(coeff)
                rnd2 = rndPick(coeff)
            # Randomly select a window
            idx1 = None
            idx2 = None
            while (idx1 == None or idx2 == None 
                    or abs(idx1 - idx2) <= 2
                    or (idx1 == 0 and idx2 == geneLen - 1)
                    or (idx1 == geneLen - 1 and idx2 == 0)):
                idx1 = random.randint(0, geneLen - 1)
                idx2 = random.randint(0, geneLen - 1)
        
            newSeq1, newSeq2 = crossover(popObj[rnd1], popObj[rnd2], idx1, idx2)
            popObj.append(newSeq1)
            popObj.append(newSeq2)

        # Mutation
        # NOTE: will always keep the worse outcome
        # swap
        numSwap = (int)(neighRatio['swap'] * popSize)
        for i in range(numSwap):
            rnd = random.randint(0, len(popObj) - 1)
            curSeq = [f.key for f in popObj[rnd].traverse()]
            idx = random.randint(0, geneLen - 1)
            keyI = curSeq[idx]
            popObj[rnd] = swap(popObj[rnd], keyI)

        # exchange
        numExchange = (int)(neighRatio['exchange'] * popSize)
        for i in range(numExchange):
            rnd = random.randint(0, len(popObj) - 1)
            curSeq = [f.key for f in popObj[rnd].traverse()]
            idx1 = None
            idx2 = None
            while (idx1 == None 
                    or idx2 == None 
                    or abs(idx1 - idx2) <= 2 
                    or (idx1 == 0 and idx2 == len(curSeq) - 1) 
                    or (idx1 == len(curSeq) - 1 and idx2 == 0)):
                idx1 = random.randint(0, len(curSeq) - 1)
                idx2 = random.randint(0, len(curSeq) - 1)
            keyI = curSeq[idx1]
            keyJ = curSeq[idx2]
            nI = popObj[rnd].query(keyI)
            nJ = popObj[rnd].query(keyJ)
            if (nI.next.key == nJ.key or nJ.next.key == nI.key):
                raise
            popObj[rnd] = exchange(popObj[rnd], keyI, keyJ)

        # rotate
        numRotate = (int)(neighRatio['rotate'] * popSize)
        for i in range(numRotate):
            rnd = random.randint(0, len(popObj) - 1)
            curSeq = [f.key for f in popObj[rnd].traverse()]
            idx1 = None
            idx2 = None
            while (idx1 == None 
                    or idx2 == None 
                    or abs(idx1 - idx2) <= 2 
                    or (idx1 == 0 and idx2 == len(curSeq) - 1) 
                    or (idx1 == len(curSeq) - 1 and idx2 == 0)):
                idx1 = random.randint(0, len(curSeq) - 1)
                idx2 = random.randint(0, len(curSeq) - 1)
            keyI = curSeq[idx1]
            keyJ = curSeq[idx2]
            nI = popObj[rnd].query(keyI)
            nJ = popObj[rnd].query(keyJ)
            if (nI.next.key == nJ.key or nJ.next.key == nI.key):
                raise
            popObj[rnd] = rotate(popObj[rnd], keyI, keyJ)

        # Tournament
        while (len(popObj) > popSize):
            # Randomly select two genes
            rnd1 = None
            rnd2 = None
            while (rnd1 == None or rnd2 == None or rnd1 == rnd2):
                rnd1 = random.randint(0, len(popObj) - 1)
                rnd2 = random.randint(0, len(popObj) - 1)
            # kill the loser
            if (popObj[rnd1].dist > popObj[rnd2].dist):
                del popObj[rnd1]
            else:
                del popObj[rnd2]

        # Update dashboard
        newOfvFound = False
        for seq in popObj:
            if (seq.dist < dashboard['bestOfv']):
                newOfvFound = True
                dashboard['bestOfv'] = seq.dist
                dashboard['bestSeq'] = [i.key for i in seq.traverse()]
        if (newOfvFound):
            iterNoImp = 0
        else:
            iterNoImp += 1
        iterTotal += 1

        # Check stopping criteria
        if ('numNoImproveIter' in stop):
            if (iterNoImp > stop['numNoImproveIter']):
                contFlag = False
                break
        if ('numIter' in stop):
            if (iterTotal > stop['numIter']):
                contFlag = False
                break
        if ('runtime' in stop):
            if ((datetime.datetime.now() - startTime).total_seconds() > stop['runtime']):
                contFlag = False
                break

    seqObj = Route(tau, asymFlag)
    for i in dashboard['bestSeq']:
        seqObj.append(nodeObj[i].clone())
    seqObj.rehead(depotID)

    return seqObj