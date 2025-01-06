import heapq
import math
import warnings
import networkx as nx

from .common import *
from .geometry import *
from .msg import *
from .travel import *

def solveTSPTW(
    nodes: dict, 
    locFieldName: str = 'loc',
    timeWindowFieldName = 'timeWindows',
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    allowViolation: bool = False,
    vehicles: dict = {0: {'speed': 1}},
    vehicleID: int|str = 0,
    serviceTime: float = 0,
    predefinedArcs: list[list[tuple[int|str]]] = [],
    edges: str = 'Euclidean',
    algo: str = 'IP',
    detailFlag: bool = False,
    metaFlag: bool = False,
    **kwargs
    ) -> dict:

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
        if ('solver' not in kwargs):
            kwargs['solver'] = 'Gurobi'
            kwargs['fml'] = 'MTZ'
        elif (kwargs['solver'] == 'Gurobi' and kwargs['fml'] not in ['MTZ']):
            raise OutOfRangeError("ERROR: Gurobi option supports 'MTZ' only")
    else:
        raise UnsupportedInputError("ERROR: Not supported yet.")

    # Define tau ==============================================================
    tau = None
    path = None
    if (detailFlag):
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

    # Add service time to field ===============================================
    for i in nodeIDs:
        if ('serviceTime' not in nodes[i]):
            nodes[i]['serviceTime'] = serviceTime

    vehSpeed = vehicles[vehicleID]['speed']

    # TSP =====================================================================
    tsptw = None
    unvisited = []
    findUnvisitedRuntime = 0
    if (algo == 'IP'):
        outputFlag = False if 'outputFlag' not in kwargs else kwargs['outputFlag']
        timeLimit = None if 'timeLimit' not in kwargs else kwargs['timeLimit']
        gapTolerance = None if 'gapTolerance' not in kwargs else kwargs['gapTolerance']
        
        if (allowViolation):
            # NOTE: Use trivial way to remove customers that obviously not reachable in time window
            unvisited = []
            checkViolationFlag = True
            while (checkViolationFlag):
                checkViolationFlag = False                
                for i in nodeIDs:
                    if (i != depotID and i not in unvisited):
                        # Latest arrival time for i
                        latestTimeForI = nodes[i][timeWindowFieldName][1]
                        canReachFlag = False
                        for j in nodeIDs:
                            if (j not in unvisited):
                                earliestTimeJ = nodes[j][timeWindowFieldName][0] + nodes[i]['serviceTime'] + tau[i, j]
                                if (earliestTimeJ <= latestTimeForI):
                                    canReachFlag = True
                                    break
                        if (canReachFlag == False):
                            unvisited.append(i)
                            checkViolationFlag = True
            tsptw = _ipTSPTWGurobiMTZMinViolation(
                nodes = nodes,
                depotID = depotID,
                nodeIDs = [i for i in nodeIDs if i not in unvisited], 
                tau = tau, 
                vehSpeed = vehSpeed,
                timeWindowFieldName = timeWindowFieldName,
                outputFlag = outputFlag, 
                timeLimit = timeLimit, 
                gapTolerance = gapTolerance)
            nodeIDs = [i for i in nodeIDs if i not in tsptw['unvisited']]
            unvisited.extend([i for i in tsptw['unvisited']])
            findUnvisitedRuntime = tsptw['runtime']

        tsptw = _ipTSPTWGurobiMTZMinTime(
            nodes = nodes,
            depotID = depotID,
            nodeIDs = nodeIDs, 
            tau = tau, 
            vehSpeed = vehSpeed,
            timeWindowFieldName = timeWindowFieldName,
            outputFlag = outputFlag, 
            timeLimit = timeLimit, 
            gapTolerance = gapTolerance)

        tsptw['fml'] = kwargs['fml']
        tsptw['solver'] = kwargs['solver']
    else:
        raise OutOfRangeError("ERROR: Select 'algo' from ['IP'].")  

    ofv = tsptw['ofv']
    nodeSeq = tsptw['seq']

    # Post optimization (for detail information) ==============================
    # FIXME: Needs rewrite for TSPTW
    if (detailFlag):
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
            shapepointsInBtw = path[nodeSeq[-2], nodeSeq[-1]]
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
        res['gap'] = tsptw['gap']
        res['unvisited'] = unvisited
        res['solType'] = tsptw['solType']
        res['lowerBound'] = tsptw['lowerBound']
        res['upperBound'] = tsptw['upperBound']
        res['runtime'] = tsptw['runtime']    
        res['findUnvisitedRuntime'] = findUnvisitedRuntime
    if (metaFlag):
        res['algo'] = algo
        res['serviceTime'] = serviceTime
    if (detailFlag):
        res['vehicles'] = vehicles

    return res

def _ipTSPTWGurobiMTZMinTime(nodes, depotID, nodeIDs, tau, vehSpeed, timeWindowFieldName, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
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
            if i != j:
                x[i, j] = TSP.addVar(
                    vtype = grb.GRB.BINARY, 
                    name = 'x_%s_%s' % (i, j))
    u = {}
    for i in nodeIDs:
        u[i] = TSP.addVar(
            vtype = grb.GRB.CONTINUOUS,
            name = 'u_%s' % (i))
    u['ret'] = TSP.addVar(vtype = grb.GRB.CONTINUOUS, name = 'u_ret', obj = 1)

    # TSP objective function ==================================================
    TSP.modelSense = grb.GRB.MINIMIZE
    TSP.update()

    # Degree constraints ======================================================
    for i in nodeIDs:
        TSP.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) == 1, name = 'leave_%s' % str(i))
        TSP.addConstr(grb.quicksum(x[j, i] for j in nodeIDs if i != j) == 1, name = 'enter_%s' % str(i))

    # Sequence constraints ====================================================
    M = 0
    for i in nodeIDs:
        maxFromI = 0
        for j in nodeIDs:
            if (tau[i, j] > maxFromI):
                maxFromI = tau[i, j]
        M += maxFromI / vehSpeed

    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j and j != depotID):
                TSP.addConstr(u[i] + nodes[i]['serviceTime'] + tau[i, j] / vehSpeed <= u[j] + M * (1 - x[i, j]))
    for i in nodeIDs:
        if (i != depotID):
            TSP.addConstr(u[i] + nodes[i]['serviceTime'] + tau[i, depotID] / vehSpeed <= u['ret'] + M * (1 - x[i, depotID]))

    for i in nodeIDs:
        TSP.addConstr(nodes[i][timeWindowFieldName][0] <= u[i])
        TSP.addConstr(u[i] <= nodes[i][timeWindowFieldName][1])

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
    currentNode = depotID
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

def _ipTSPTWGurobiMTZMinViolation(nodes, depotID, nodeIDs, tau, vehSpeed, timeWindowFieldName, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except(ImportError):
        print("ERROR: Cannot find Gurobi")
        return

    # Initialize
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
            if i != j:
                x[i, j] = TSP.addVar(
                    vtype = grb.GRB.BINARY, 
                    name = 'x_%s_%s' % (i, j))
    u = {}
    for i in nodeIDs:
        u[i] = TSP.addVar(
            vtype = grb.GRB.CONTINUOUS,
            name = 'u_%s' % (i))
    u['ret'] = TSP.addVar(vtype = grb.GRB.CONTINUOUS, name = 'u_ret')

    visit = {}
    for i in nodeIDs:
        visit[i] = TSP.addVar(
            vtype = grb.GRB.BINARY,
            obj = nodes[i]['value'] if 'value' in nodes[i] else 1,
            name = 'visit_%s' % (i))
    
    # TSP objective function ==================================================
    TSP.modelSense = grb.GRB.MAXIMIZE
    TSP.update()

    # Degree constraints ======================================================
    TSP.addConstr(visit[depotID] == 1)
    for i in nodeIDs:
        TSP.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) == grb.quicksum(x[k, i] for k in nodeIDs if i != k))

    for i in nodeIDs:
        TSP.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) <= 1)
        TSP.addConstr(grb.quicksum(x[k, i] for k in nodeIDs if i != k) <= 1)
        TSP.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) >= visit[i])
        TSP.addConstr(grb.quicksum(x[k, i] for k in nodeIDs if i != k) >= visit[i])

    # Sequence constraints ====================================================
    M = 0
    for i in nodeIDs:
        maxFromI = 0
        for j in nodeIDs:
            if (tau[i, j] > maxFromI):
                maxFromI = tau[i, j]
        M += maxFromI / vehSpeed

    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j and j != depotID):
                TSP.addConstr(u[i] + nodes[i]['serviceTime'] + tau[i, j] / vehSpeed <= u[j] + M * (1 - x[i, j]))
    for i in nodeIDs:
        if (i != depotID):
            TSP.addConstr(u[i] + nodes[i]['serviceTime'] + tau[i, depotID] / vehSpeed <= u['ret'] + M * (1 - x[i, depotID]))

    for i in nodeIDs:
        TSP.addConstr(nodes[i][timeWindowFieldName][0] <= u[i])
        TSP.addConstr(u[i] <= nodes[i][timeWindowFieldName][1])

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
    unvisited = []
    for i in visit:
        if (visit[i].x < 0.5):
            unvisited.append(i)
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
        'unvisited': unvisited,
        'gap': gap,
        'runtime': runtime
    }
