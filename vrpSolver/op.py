import heapq
import math
import warnings
import networkx as nx

from .common import *
from .geometry import *
from .msg import *

def ipOP(
    nodes: dict, 
    maxBudget: float,
    locFieldName: str = 'loc',
    priceFieldName: str = 'price',
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    vehicles: dict = {
        0: {'speed': 1}
    },
    vehicleID: int|str = 0,
    edges: dict = {
        'method': "Euclidean", 
        'ratio': 1
    },
    method: dict = {
        'fml': 'MTZ',
        'solver': 'Gurobi',
        'timeLimit': None,
        'outputFlag': False,
        'env': None
    },
    detailsFlag: bool = False,
    metaFlag: bool = False
    ) -> dict|None:

    """Orienteering Problem"""

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
    if ((type(nodeIDs) == list and depotID not in nodeIDs)
        or (nodeIDs == 'All' and depotID not in nodes)):
        raise OutOfRangeError("ERROR: Cannot find `depotID` in given `nodes`/`nodeIDs`")
    if (detailsFlag == True):
        # For animation propose
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
    price = {}
    for i in nodeIDs:
        price[i] = nodes[i][priceFieldName]

    # Solve by formulation ====================================================
    op = None
    if (method['solver'] == 'Gurobi'):
        if (method['fml'] == 'MTZ'):
            op = _ipOPGurobiMTZ(
                nodeIDs, 
                depotID,
                tau, 
                price,
                maxBudget,
                method['outputFlag'] if 'outputFlag' in method else None, 
                method['timeLimit'] if 'timeLimit' in method else None, 
                method['gapTolerance'] if 'gapTolerance' in method else None)
    if (op == None):
        raise UnsupportedInputError("ERROR: Incorrect or not available OP formulation option!")

    op['fml'] = method['fml']

    # Fix the sequence to make it start from the depot ========================
    startIndex = 0
    seq = [i for i in op['seq']]
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
    op['seq'] = nodeSeq

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

    res = {
        'ofv': op['ofv'],
        'seq': nodeSeq,
        'gap': op['gap'],
        'solType': op['solType'],
        'lowerBound': op['lowerBound'],
        'upperBound': op['upperBound'],
        'runtime': op['runtime']
    }
    if (metaFlag):
        res['method'] = method
    if (detailsFlag):
        res['vehicles'] = vehicles

    return res

def _ipOPGurobiMTZ(nodeIDs, depotID, tau, price, maxBudget, outputFlag, timeLimit, gapTolerance):
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
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = OP.addVar(
                    vtype = grb.GRB.BINARY, 
                    name = 'x_%s_%s' % (i, j))
    u = {}
    for i in range(n):
        u[i] = OP.addVar(
            vtype = grb.GRB.CONTINUOUS,
            name = 'u_%s' % (i))

    # OP objective function ===================================================
    OP.setObjective(
        grb.quicksum(price[i] * x[i, j] for i in nodeIDs for j in nodeIDs if i != j))
    OP.modelSense = grb.GRB.MAXIMIZE
    OP.update()

    # Degree constraints ======================================================
    OP.addConstr(grb.quicksum(x[depotID, j] for j in range(n) if j != depotID) == 1, name = 'leave_%s' % str(i))
    OP.addConstr(grb.quicksum(x[j, depotID] for j in range(n) if j != depotID) == 1, name = 'enter_%s' % str(i))

    for i in range(n):
        if (i != depotID):
            OP.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) <= 1, name = 'leave_%s' % str(i))
            OP.addConstr(
                grb.quicksum(x[j, i] for j in range(n) if i != j) 
                == grb.quicksum(x[i, j] for j in range(n) if i != j), name = 'enter_%s' % str(i))

    # Budget constraints ======================================================
    OP.addConstr(grb.quicksum(tau[i, j] * x[i, j] for i in nodeIDs for j in nodeIDs if i != j) <= maxBudget)

    # Sequence constraints ====================================================
    for i in range(1, n):
        for j in range(1, n):
            if (i != j):
                OP.addConstr(u[i] - u[j] + (n - 1) * x[i, j] <= n - 2, name = 'seq_%s_%s' % (i, j))
    for i in range(1, n):
        OP.addConstr(1 <= u[i])
        OP.addConstr(u[i] <= n - 1)

    # OP =====================================================================
    OP.write("OP.lp")
    OP.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    gap = None
    seq = []
    arcs = []
    solType = None
    lb = None
    ub = None
    runtime = None
    dist = None
    if (OP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        ofv = OP.getObjective().getValue()
        dist = 0
        for i, j in x:
            if (x[i, j].x > 0.5):
                arcs.append([i, j])
                dist += tau[i, j]
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
        runtime = OP.Runtime
    elif (OP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        ofv = OP.getObjective().getValue()
        gap = OP.Params.MIPGapAbs
        dist = 0
        for i, j in x:
            if (x[i, j].x > 0.5):
                arcs.append([i, j])
                dist += tau[i, j]
        currentNode = 0
        seq.append(nodeIDs[currentNode])
        while (len(arcs) > 0):
            for i in range(len(arcs)):
                if (arcs[i][0] == currentNode):
                    currentNode = arcs[i][1]
                    seq.append(nodeIDs[currentNode])
                    arcs.pop(i)
                    break
        gap = OP.MIPGap
        lb = OP.ObjBoundC
        ub = OP.ObjVal
        runtime = OP.Runtime

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
