import heapq
import math
import warnings
import networkx as nx

from .common import *
from .geometry import *
from .msg import *
from .travel import *

def solveOP(
    nodes: dict, 
    maxBudget: float,
    locFieldName: str = 'loc',
    priceFieldName: str = 'price',
    startID: int|str = 0,
    endID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    vehicles: dict = {0: {'speed': 1}},
    vehicleID: int|str = 0,
    edges: str = 'Euclidean',
    algo: str = 'IP',
    detailsFlag: bool = False,
    metaFlag: bool = False,
    **kwargs
    ) -> dict|None:

    """Solve the Orienteering Problem"""

    # Sanity check ============================================================
    if (nodes == None or type(nodes) != dict):
        raise MissingParameterError(ERROR_MISSING_NODES)
    for i in nodes:
        if (locFieldName not in nodes[i]):
            raise MissingParameterError("ERROR: Node %s does not specify location information in `nodes`." % i)
        if (priceFieldName not in nodes[i]):
            raise MissingParameterError("ERROR: Node %s does not specify price information in `nodes`." % i)
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
        else:
            for i in nodeIDs:
                if (i not in nodes):
                    raise OutOfRangeError("ERROR: Node %s is not in `nodes`." % i)

    if ((type(nodeIDs) == list and startID not in nodeIDs)
        or (nodeIDs == 'All' and startID not in nodes)):
        raise OutOfRangeError("ERROR: Cannot find `startID` in given `nodes`/`nodeIDs`")
    if ((type(nodeIDs) == list and endID not in nodeIDs)
        or (nodeIDs == 'All' and endID not in nodes)):
        raise OutOfRangeError("ERROR: Cannot find `endID` in given `nodes`/`nodeIDs`")

    if (detailsFlag == True):
        # For animation propose
        if (vehicles == None):
            raise MissingParameterError("ERROR: Missing required field `vehicles`.")
        if (vehicleID not in vehicles):
            raise MissingParameterError("ERROR: Cannot find `vehicleID` in `vehicles`.")

    if (algo == 'IP'):
        if ('solver' not in kwargs):
            raise MissingParameterError("ERROR: Missing required field `solver`.")
        elif (kwargs['solver'] == 'Gurobi' and kwargs['fml'] not in ['MTZ', 'DFJ_Lazy']):
            raise OutOfRangeError("ERROR: Formulation is not supported.", )

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

    price = {}
    for i in nodeIDs:
        price[i] = nodes[i][priceFieldName]

    # Solve by formulation ====================================================
    op = None
    if (algo == 'IP'):
        if (kwargs['fml'] == 'MTZ'):
            op = _ipOPGurobiMTZ(
                nodeIDs, 
                startID,
                endID,
                tau, 
                price,
                maxBudget,
                kwargs['outputFlag'] if 'outputFlag' in kwargs else False, 
                kwargs['timeLimit'] if 'timeLimit' in kwargs else None, 
                kwargs['gapTolerance'] if 'gapTolerance' in kwargs else None)
        elif (kwargs['fml'] == 'DFJ_Lazy'):
            op = _ipOPGurobiLazyCuts(
                nodeIDs, 
                startID,
                endID,
                tau, 
                price,
                maxBudget,
                kwargs['outputFlag'] if 'outputFlag' in kwargs else False, 
                kwargs['timeLimit'] if 'timeLimit' in kwargs else None, 
                kwargs['gapTolerance'] if 'gapTolerance' in kwargs else None)

    if (op == None):
        raise UnsupportedInputError("ERROR: Incorrect or not available OP formulation option!")

    op['fml'] = kwargs['fml']

    # Fix the sequence to make it start from the depot ========================
    startIndex = 0
    seq = [i for i in op['seq']]
    nodeSeq = op['seq']
    
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

    res = {
        'ofv': op['ofv'],
        'dist': op['dist'],
        'seq': op['seq'],
        'gap': op['gap'],
        'solType': op['solType'],
        'lowerBound': op['lowerBound'],
        'upperBound': op['upperBound'],
        'runtime': op['runtime']
    }
    if (metaFlag):
        res['method'] = kwargs['method']
    if (detailsFlag):
        res['vehicles'] = vehicles

    return res

def _ipOP(nodeIDs, startID, endID, tau, price, maxBudget, outputFlag, timeLimit, gapTolerance):
    op = None
    if (solver == 'Gurobi'):
        if (fml == 'DFJ_Lazy'):
            op = _ipOPGurobiLazyCuts(nodeIDs, startID, endID, tau, price, maxBudget, outputFlag, timeLimit, gapTolerance)
        elif (fml == 'MTZ'):
            op = _ipOPGurobiMTZ(nodeIDs, startID, endID, tau, price, maxBudget, outputFlag, timeLimit, gapTolerance)
    if (op == None):
        raise UnsupportedInputError("ERROR: Currently OP can be calculated by Gurobi using MTZ formulation")

    return op

def _ipOPGurobiMTZ(nodeIDs, startID, endID, tau, price, maxBudget, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except:
        raise ImportError("ERROR: Cannot find Gurobi")

    # Initialize
    n = len(nodeIDs)
    OP = grb.Model('OP')
    if (outputFlag == False):
        OP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        OP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        OP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if i != j:
                x[i, j] = OP.addVar(
                    vtype = grb.GRB.BINARY, 
                    name = 'x_%s_%s' % (i, j))
    u = {}
    for i in nodeIDs:
        u[i] = OP.addVar(
            vtype = grb.GRB.CONTINUOUS,
            name = 'u_%s' % (i))

    # OP objective function ===================================================
    OP.setObjective(
        grb.quicksum(price[i] * x[i, j] for i in nodeIDs for j in nodeIDs if i != j))
    OP.modelSense = grb.GRB.MAXIMIZE
    OP.update()

    # a = [1, 13, 17, 6, 18, 9, 14, 3, 2, 1]
    # for i in range(len(a) - 1):
    #     OP.addConstr(x[a[i], a[i + 1]] == 1)

    # Degree constraints ======================================================
    OP.addConstr(grb.quicksum(x[startID, j] for j in nodeIDs if j != startID) == 1, name = 'leave_%s' % str(i))
    OP.addConstr(grb.quicksum(x[j, endID] for j in nodeIDs if j != endID) == 1, name = 'enter_%s' % str(i))

    for i in nodeIDs:
        if (i != startID):
            OP.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) <= 1, name = 'edge_%s' % str(i))
    
    # Flow balancing ==========================================================
    for i in nodeIDs:
        if (i != startID and i != endID):
            OP.addConstr(
                grb.quicksum(x[j, i] for j in nodeIDs if i != j) 
                == grb.quicksum(x[i, j] for j in nodeIDs if i != j), name = 'balance_%s' % str(i))

    # Budget constraints ======================================================
    OP.addConstr(grb.quicksum(tau[i, j] * x[i, j] for i in nodeIDs for j in nodeIDs if i != j) <= maxBudget)

    # Sequence constraints ====================================================
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                OP.addConstr(u[i] - u[j] + (n - 1) * x[i, j] <= n - 2, name = 'seq_%s_%s' % (i, j))
    for i in nodeIDs:
        OP.addConstr(1 <= u[i])
        OP.addConstr(u[i] <= n - 1)

    # OP =====================================================================
    OP.optimize()

    arcs = []
    dist = 0
    for i, j in x:
        if (x[i, j].x > 0.5):
            arcs.append([i, j])
            dist += tau[i, j]

    currentNode = startID
    seq = []
    seq.append(currentNode)
    while (len(arcs) > 0):
        for i in range(len(arcs)):
            if (arcs[i][0] == currentNode):
                currentNode = arcs[i][1]
                seq.append(nodeIDs[currentNode])
                arcs.pop(i)
                break

    if (OP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        ofv = OP.getObjective().getValue()
        gap = 0
        lb = ofv
        ub = ofv
        runtime = OP.Runtime
    elif (OP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        ofv = OP.getObjective().getValue()
        gap = OP.MIPGap
        lb = OP.ObjBoundC
        ub = OP.ObjVal
        runtime = OP.Runtime
    else:
        return None

    return {
        'ofv': ofv,
        'dist': dist,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }

def _ipOPGurobiLazyCuts(nodeIDs, startID, endID, tau, price, maxBudget, outputFlag, timeLimit, gapTolerance):
    try:
        import gurobipy as grb
    except:
        raise ImportError("ERROR: Cannot find Gurobi")

    # Initialize
    n = len(nodeIDs)
    OP = grb.Model('OP')
    if (outputFlag == False):
        OP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        OP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        OP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                x[i, j] = OP.addVar(
                    vtype = grb.GRB.BINARY, 
                    obj = tau[i, j], 
                    name = 'x_%s_%s' % (i, j))
                
    # OP objective function ==================================================
    OP.setObjective(
        grb.quicksum(price[i] * x[i, j] for i in nodeIDs for j in nodeIDs if i != j))
    OP.modelSense = grb.GRB.MAXIMIZE
    OP.Params.lazyConstraints = 1
    OP.update()

    # Degree constraints ======================================================
    OP.addConstr(grb.quicksum(x[startID, j] for j in nodeIDs if j != startID) == 1, name = 'leave_%s' % str(i))
    OP.addConstr(grb.quicksum(x[j, endID] for j in nodeIDs if j != endID) == 1, name = 'enter_%s' % str(i))

    for i in nodeIDs:
        if (i != startID):
            OP.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) <= 1, name = 'edgeS_%s' % str(i))
        if (j != endID):
            OP.addConstr(grb.quicksum(x[i, j] for j in nodeIDs if i != j) <= 1, name = 'edgeE_%s' % str(i))
 
    # Flow balancing ==========================================================
    for i in nodeIDs:
        if (i != startID and i != endID):
            OP.addConstr(
                grb.quicksum(x[j, i] for j in nodeIDs if i != j) 
                == grb.quicksum(x[i, j] for j in nodeIDs if i != j), name = 'balance_%s' % str(i))

    # Budget constraints ======================================================
    OP.addConstr(grb.quicksum(tau[i, j] * x[i, j] for i in nodeIDs for j in nodeIDs if i != j) <= maxBudget)

    # Sub-tour elimination ====================================================
    OP._x = x
    def subtourelim(model, where):
        if (where == grb.GRB.Callback.MIPSOL):
            x_sol = model.cbGetSolution(model._x)
            G = nx.Graph()
            for (i, j) in x.keys():
                if (x_sol[i, j] > 0.9):
                    G.add_edge(i, j, weight = tau[i, j])
            components = [list(c) for c in nx.connected_components(G)]
            for component in components:
                if (startID not in component or endID not in component):
                    model.cbLazy(grb.quicksum(x[i, j] for i in component for j in component if i != j) <= len(component) - 1)
                # 如果成环
                if (startID != endID and (startID in component or endID in component)):
                    model.cbLazy(grb.quicksum(x[i, j] for i in component for j in component if i != j) <= len(component) - 1)

    # OP with callback =======================================================
    OP.optimize(subtourelim)

    arcs = []
    dist = 0
    for i, j in x:
        if (x[i, j].x > 0.5):
            arcs.append([i, j])
            dist += tau[i, j]

    currentNode = startID
    seq = []
    seq.append(currentNode)
    while (len(arcs) > 0):
        for i in range(len(arcs)):
            if (arcs[i][0] == currentNode):
                currentNode = arcs[i][1]
                seq.append(nodeIDs[currentNode])
                arcs.pop(i)
                break

    if (OP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        ofv = OP.getObjective().getValue()
        gap = 0
        lb = ofv
        ub = ofv
        runtime = OP.Runtime
    elif (OP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        ofv = OP.getObjective().getValue()
        gap = OP.MIPGap
        lb = OP.ObjBoundC
        ub = OP.ObjVal
        runtime = OP.Runtime
    else:
        return None

    return {
        'ofv': ofv,
        'dist': dist,
        'seq': seq,
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
    }
