import heapq
import math
import warnings
import networkx as nx

import gurobipy as grb

from .common import *
from .geometry import *
from .msg import *

@runtime("solveVRPTW")
def solveVRPTW(
    nodes: dict, 
    numVeh: int,
    vehCap: float,
    locFieldName: str = 'loc',
    timeWindowFieldName: str = 'timeWindow',
    demandFieldName: str = 'demand',
    serviceTimeFieldName: str = 'serviceTime',
    depotID: int|str = 0,
    customerIDs: list[int|str]|str = 'AllButDepot',
    edges: str = "Euclidean", 
    method: dict = {
        'algo': 'CG',
        'subproblemAlgo': 'IP',
        'solver': 'Gurobi',
        'stop': {'maxRuntime': 1200},
        'outputFlag': False
    },
    detailsFlag: bool = False,
    metaFlag: bool = False
    ) -> dict|None:

    # Sanity check ============================================================
    if (nodes == None or type(nodes) != dict):
        raise MissingParameterError(ERROR_MISSING_NODES)
    for i in nodes:
        if (locFieldName not in nodes[i]):
            raise MissingParameterError("ERROR: missing location information in `nodes`.")
    if (depotID == None):
        raise MissingParameterError("ERROR: missing `depotID`.")
    if (depotID not in nodes):
        raise OutOfRangeError("ERROR: Cannot find `depotID` in given `nodes`")
    if (type(customerIDs) is not list):
        if (customerIDs == 'AllButDepot'):
            customerIDs = [i for i in nodes if i != depotID]
        else:
            for i in customerIDs:
                if (i not in nodes):
                    raise OutOfRangeError("ERROR: Customer %s is not in `nodes`." % i)

    # Define tau ==============================================================
    nodePlus = [depotID]
    nodePlus.extend([i for i in customerIDs])
    nodeMinus = [i for i in customerIDs]
    dupDepotID = max(list(nodes.keys())) + 1
    nodeMinus.append(dupDepotID)

    tau = None
    path = None
    if (detailsFlag):
        tau, path = matrixDist(
            nodes = nodes, 
            edges = edges, 
            nodeIDs = customerIDs, 
            locFieldName = locFieldName)
        tauStart, tauEnd, pathStart, pathEnd = vectorDist(
            loc = nodes[depotID][locFieldName],
            nodes = nodes,
            edges = edges,
            nodeIDs = customerIDs,
            locFieldName = locFieldName)
        for i in customerIDs:
            tau[depotID, i] = tauStart[i]
            tau[i, dupDepotID] = tauEnd[i]
            path[depotID, i] = pathStart[i]
            path[i, dupDepotID] = pathEnd[i]
        tau[depotID, dupDepotID] = 0
        tau[dupDepotID, depotID] = 0
    else:
        tau, _ = matrixDist(
            nodes = nodes, 
            edges = edges, 
            nodeIDs = customerIDs, 
            locFieldName = locFieldName)
        tauStart, tauEnd, _, _ = vectorDist(
            loc = nodes[depotID][locFieldName],
            nodes = nodes,
            edges = edges,
            nodeIDs = customerIDs,
            locFieldName = locFieldName)
        for i in customerIDs:
            tau[depotID, i] = tauStart[i]
            tau[i, dupDepotID] = tauEnd[i]

    # Pricing problem =========================================================
    @runtime("pricingIP")
    def _pricingIP(pi):
        # print(pi)
        sub = grb.Model('Pricing')

        # Define decision variables
        w = {}
        for i in nodePlus:
            for j in nodeMinus:
                if (i != j and not (i == depotID and j == dupDepotID)):
                    w[i, j] = sub.addVar(
                        vtype = grb.GRB.BINARY, 
                        obj = tau[i, j] - pi[i], 
                        name = 'w_%s_%s' % (i, j))
        t = {}
        for i in customerIDs:
            estT = 0
            latT = float('inf')
            if (timeWindowFieldName in nodes[i]):
                estT = nodes[i][timeWindowFieldName][0]
                latT = nodes[i][timeWindowFieldName][1]
            t[i] = sub.addVar(
                vtype = grb.GRB.CONTINUOUS, 
                lb = estT,
                ub = latT,
                name = 't_%s' % i)

        estT = 0
        latT = float('inf')
        if (timeWindowFieldName in nodes[depotID]):
            estT = nodes[depotID][timeWindowFieldName][0]
            latT = nodes[depotID][timeWindowFieldName][1]
        t[depotID] = sub.addVar(
            vtype = grb.GRB.CONTINUOUS, 
            lb = estT,
            ub = latT,
            name='t_depot')
        t[dupDepotID] = sub.addVar(
            vtype = grb.GRB.CONTINUOUS, 
            lb = estT,
            ub = latT,
            name='t_dupDepot')

        # Assignment constraints
        for i in customerIDs:
            sub.addConstr(grb.quicksum(w[i, j] for j in nodeMinus if i != j and not (i == depotID and j == dupDepotID)) 
                == grb.quicksum(w[j, i] for j in nodePlus if i != j and not (i == depotID and j == dupDepotID)))
        sub.addConstr(grb.quicksum(w[depotID, j] for j in nodeMinus if j != dupDepotID) == 1)
        sub.addConstr(grb.quicksum(w[i, dupDepotID] for i in nodePlus if i != depotID) == 1)
        sub.update()

        # Capacity constraints
        sub.addConstr(
            grb.quicksum(w[i, j] * nodes[i][demandFieldName] for i, j in w if i != j) <= vehCap)

        # Time windows
        # REMEMBER: Yes I have spent a lot of time here! t[depotID] cannot take two values at the same time
        M = max([nodes[i][timeWindowFieldName][1] for i in customerIDs if timeWindowFieldName in nodes[i]])
        M += max(tau.values())
        for i in nodePlus:
            for j in nodeMinus:
                if (i != j and not (i == depotID and j == dupDepotID)):
                    sub.addConstr(t[i] + tau[i, j] + nodes[i][serviceTimeFieldName] - M * (1 - w[i, j]) <= t[j])
        sub.update()

        # Solve
        sub.setParam("OutputFlag", 0)
        if ('stop' in method and 'maxIterRuntime' in method['stop']):
            sub.setParam(grb.GRB.Param.TimeLimit, method['stop']['maxIterRuntime'])
        sub.modelSense = grb.GRB.MINIMIZE
        sub.optimize()

        # Construct solution
        ofv = None
        c = 0
        a = []
        route = []
        if (sub.status == grb.GRB.status.OPTIMAL):
            ofv = sub.getObjective().getValue()
            arcSet = []
            for i, j in w:
                if (w[i, j].x > 0.9):
                    c += tau[i, j]
                    arcSet.append((i, j))
            route = [depotID]
            while (len(arcSet) > 0):
                cur = None
                for arc in arcSet:                    
                    if (arc[0] == route[-1]):
                        route.append(arc[1])
                        arcSet.remove(arc)
                        cur = arc[1]
                        break
                    if (arc[1] == route[-1]):
                        route.append(arc[0])
                        arcSet.remove(arc)
                        cur = arc[0]
                        break
                if (cur == dupDepotID):
                    break
            a = {}
            for i in customerIDs:
                if (i in route):
                    a[i] = 1
                else:
                    a[i] = 0
        elif (sub.status == grb.GRB.status.TIME_LIMIT):
            try:
                ofv = sub.getObjective().getValue()
                arcSet = []
                for i, j in w:
                    if (w[i, j].x > 0.9):
                        c += tau[i, j]
                        arcSet.append((i, j))

                route = [depotID]
                while (len(arcSet) > 0):
                    cur = None
                    for arc in arcSet:                    
                        if (arc[0] == route[-1]):
                            route.append(arc[1])
                            arcSet.remove(arc)
                            cur = arc[1]
                            break
                        if (arc[1] == route[-1]):
                            route.append(arc[0])
                            arcSet.remove(arc)
                            cur = arc[0]
                            break
                    if (cur == dupDepotID):
                        break
                a = {}
                for i in customerIDs:
                    if (i in route):
                        a[i] = 1
                    else:
                        a[i] = 0
            except:
                ofv = None
                c = None
                a = None
                route = None

        return {
            'ofv': ofv,
            'cr': c,
            'ai': a,
            'route': route
        }

    def _pricingLabeling(pi):

        def _createGraph():
            g = nx.DiGraph()
            estT = 0
            latT = float('inf')
            if (timeWindowFieldName in nodes[depotID]):
                estT = nodes[depotID][timeWindowFieldName][0]
                latT = nodes[depotID][timeWindowFieldName][1]
            g.add_node(depotID,
                nodeType = 'Depot',
                timeWindow = [estT, latT],
                weight = 0,
                serviceTime = nodes[depotID][serviceTimeFieldName],
                arrivalTime = None,     # Label - Arrival time
                minDist = 0,            # Label - Minimum distance from depot
                preNode = None)         # Label - Previous node
            g.add_node(dupDepotID,
                nodeType = 'Depot',
                timeWindow = [estT, latT],
                weight = 0,
                serviceTime = nodes[depotID][serviceTimeFieldName],
                arrivalTime = None,     # Label - Arrival time
                minDist = 0,            # Label - Minimum distance from depot
                preNode = None)         # Label - Previous node
            for i in customerIDs:
                estT = 0
                latT = float('inf')
                if (timeWindowFieldName in nodes[i]):
                    estT = nodes[i][timeWindowFieldName][0]
                    latT = nodes[i][timeWindowFieldName][1]
                g.add_node(
                    i,
                    nodeType = 'Customer',
                    timeWindow = [estT, latT],
                    weight = nodes[i][demandFieldName],
                    serviceTime = nodes[i][serviceTimeFieldName],
                    arrivalTime = None,     # Label - Arrival time
                    minDist = 0,            # Label - Minimum distance from depot
                    preNode = None)         # Label - Previous node
            for i in nodePlus:
                for j in nodeMinus:
                    if (i != j and not (i == depotID and j == dupDepotID)):
                        g.add_edge(i, j, 
                            travelTime = tau[i, j], 
                            travelDist = tau[i, j] - pi[i])
            return g

        def _initializeLabel():
            return [{
                'path': [depotID],
                'dist': 0,
                'load': 0,
                'time': 0
            }]

        def _feasibility(curPath, node):
            # Check if node has been covered
            if (node in curPath['path']):
                return False
            # Check time windows
            lastNode = curPath['path'][-1]
            arrTime = curPath['time'] + g.edges[lastNode, node]['travelTime']
            availTW = []
            if (node != dupDepotID):
                availTW = g.nodes[node]['timeWindow']
            else:
                availTW = g.nodes[depotID]['timeWindow']
            if (arrTime < availTW[0] or arrTime > availTW[1]):
                return False
            # Check loads
            accLoad = curPath['load'] + g.nodes[node]['weight']
            if (accLoad > vehCap):
                return False
            return True

        def _extendPath(curPath, node):
            # Make a copy
            extendPath = {
                'path': [i for i in curPath['path']],
                'time': curPath['time'],
                'load': curPath['load'],
                'dist': curPath['dist']
            }                
            extendPath['path'].append(child)
            # Extend the path
            extendPath['dist'] += g.edges[lastNode, child]['travelDist']
            extendPath['load'] += g.nodes[node]['weight']
            extendPath['time'] = max(
                curPath['time'] + g.edges[lastNode, child]['travelTime'],
                g.nodes[node]['timeWindow'][0])
            return extendPath

        def _dominate(label1, label2):
            # Given two labels
            # - return True if label1 is dominating label2, label2 is useless
            # - return False if label2 is not dominated by label1
            if (label1['path'][-1] == label2['path'][-1]
                and label1['dist'] <= label2['dist']
                and label1['time'] <= label2['time']
                and label1['load'] <= label2['load']):
                return True
            else:
                return False

        # Main ================================================================
        g = _createGraph()

        # 初始label, queue, path
        queue = _initializeLabel()
        path = {}
        accSol = 0

        # 主循环
        while (len(queue) > 0):
            # 从队列中取出第一个
            curPath = queue.pop(0)
            
            # 尝试拓展该标签
            lastNode = curPath['path'][-1]
            for child in g.successors(lastNode):
                if (_feasibility(curPath, child)):
                    extendPath = _extendPath(curPath, child)
                    print("Extended: ", extendPath['path'])
                    queue.append(extendPath)
                    if (child == dupDepotID):
                        path[accSol] = extendPath
                        accSol += 1

            # 移除dominated label
            removeQueueIndex = []
            for i in range(len(queue)):
                for j in range(len(queue)):
                    if (i != j 
                        and _dominate(queue[i], queue[j]) == True
                        and j not in removeQueueIndex):
                        removeQueueIndex.append(j)
            print("Remove %s labels." % len(removeQueueIndex))
            queue = [queue[i] for i in range(len(queue)) if i not in removeQueueIndex]

        optPath = None
        minDist = float('inf')
        for i in path:
            if (path[i]['path'][-1] == dupDepotID
                and path[i]['dist'] < minDist):
                minDist = path[i]['dist']
                optPath = path[i]

        ofv = optPath['dist']
        c = 0
        for i in range(len(optPath['path']) - 1):
            c += tau[optPath['path'][i], optPath['path'][i + 1]]
        a = {}
        for i in customerIDs:
            if (i in optPath['path']):
                a[i] = 1
            else:
                a[i] = 0

        print(ofv)

        return {
            'ofv': ofv,
            'cr': c,
            'ai': a,
            'route': optPath['path']
        }

    # Master problem initialization ===========================================
    CVRPTW = grb.Model('CVRPTW')
    CVRPTW.setParam("OutputFlag", 0)

    # Route sets
    routes = {}
    
    # Parameters
    c = {} # Cost of route, c[routeID]
    a = {} # Whether node in route, a[nodeID, routeID]

    # Initialize parameters
    # 初始解，每个customer分配一辆车
    for i in customerIDs:
        for j in customerIDs:
            if (i != j):
                a[i, j] = 0
            else:
                a[i, j] = 1
        c[i] = tau[depotID, i] + tau[i, dupDepotID]
        routes[i] = [depotID, i, dupDepotID]
    acc = max(customerIDs) + 1

    # Binary for candidate route, y[routeID], start with n routes
    y = {}

    # Initial columns
    for i in customerIDs:
        y[i] = CVRPTW.addVar(vtype=grb.GRB.CONTINUOUS, obj=c[i])

    # Set of constraints
    cons = {}

    # Depot处离开的车辆数量约束
    # NOTE: 取符号以保持reduce cost的符号一致性
    cons[depotID] = CVRPTW.addConstr(grb.quicksum(-y[r] for r in y) >= - numVeh)

    # 每个customer都要被访问
    for i in customerIDs:
        cons[i] = CVRPTW.addConstr(grb.quicksum(a[i, r] * y[r] for r in y) == 1)

    # Initial Optimize
    CVRPTW.modelSense = grb.GRB.MINIMIZE
    CVRPTW.optimize()

    # Now solve the Master (Set Partition Formulation) problem ================
    accTime = 0
    canAddVarFlag = True
    while(canAddVarFlag):
        canAddVarFlag = False

        startIter = datetime.datetime.now()
        
        # 主问题接到最优解
        if (CVRPTW.status == grb.GRB.status.OPTIMAL):
            # Solve subproblem
            pi = {}
            for constraint in cons:
                pi[constraint] = cons[constraint].Pi
            # 将约束的对偶转入给子问题
            subproblem = None
            if (method['subproblemAlgo'] == 'IP'):
                subproblem = _pricingIP(pi)
            else:
                subproblem = _pricingLabeling(pi)
            
            # 如果子问题有解，尝试加入
            if (subproblem['ofv'] != None):
                # 检查reduce cost，如果大于0，加入新列
                if (cons[depotID].Pi - subproblem['ofv'] > CONST_EPSILON):
                    # 首先，将加入标识设为可加入
                    canAddVarFlag = True
                    writeLog("Subroute found: " + list2String(subproblem['route']))
                    # 新路径的cost
                    accRID = max(list(y.keys())) + 1
                    c[accRID] = subproblem['cr']
                    y[accRID] = CVRPTW.addVar(vtype = grb.GRB.CONTINUOUS, obj = c[accRID])
                    for i in customerIDs:
                        a[i, accRID] = subproblem['ai'][i]
                    routes[accRID] = subproblem['route']
                    # Update columns, add one more
                    CVRPTW.chgCoeff(cons[depotID], y[accRID], -1)
                    for i in customerIDs:
                        CVRPTW.chgCoeff(cons[i], y[accRID], a[i, accRID])
                    CVRPTW.update()
                else:
                    canAddVarFlag = False
            else:
                canAddVarFlag = False
        else:
            canAddVarFlag = False

        # Re-optimize
        CVRPTW.optimize()
        oneIter = (datetime.datetime.now() - startIter).total_seconds()
        accTime += oneIter
        if ('stop' in method and 'maxRuntime' in method['stop'] and accTime > method['stop']['maxRuntime']):
            canAddVarFlag = False

    # Interpret solution for lower bound ======================================
    lb = CVRPTW.getObjective().getValue()

    # Early branching heuristic ===============================================
    for i in y:
        y[i].vtype = grb.GRB.BINARY
    CVRPTW.update()
    CVRPTW.optimize()

    # Interpret solution ======================================================
    solRoute = {}
    acc = 1
    for i in y:
        if (y[i].x > 0.9):            
            solRoute[acc] = {
                'route': routes[i],
                'length': c[i]
            }
            acc += 1
    ofv = CVRPTW.getObjective().getValue()

    return {
        'ofv': ofv,
        'lb': lb,
        'ub': ofv,
        'gap': (ofv - lb) / lb,
        'routes': solRoute
    }