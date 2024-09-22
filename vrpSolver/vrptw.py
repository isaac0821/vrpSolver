import heapq
import math
import warnings
import networkx as nx

import gurobipy as grb

from .common import *
from .geometry import *
from .msg import *
from .travel import *

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
    algo: str = 'CG_EarlyBranching',
    detailsFlag: bool = False,
    metaFlag: bool = False,
    **kwargs
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
            nodeIDs = customerIDs, 
            edges = edges,             
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

    # Dynamic columns =========================================================
    columns = {}

    # Pricing problem =========================================================
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
                        obj = (tau[i, j] - pi[i]),
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

    @runtime("pricingLabelSetting")
    def _pricingLabelSetting(pi):
        @runtime("createGraph")
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
                label = None)
            g.add_node(dupDepotID,
                nodeType = 'Depot',
                timeWindow = [estT, latT],
                weight = 0,
                serviceTime = nodes[depotID][serviceTimeFieldName],
                label = None)
            for i in customerIDs:
                estT = 0
                latT = float('inf')
                if (timeWindowFieldName in nodes[i]):
                    estT = nodes[i][timeWindowFieldName][0]
                    latT = nodes[i][timeWindowFieldName][1]
                g.add_node(i,
                    nodeType = 'Customer',
                    timeWindow = [estT, latT],
                    weight = nodes[i][demandFieldName],
                    serviceTime = nodes[i][serviceTimeFieldName],
                    label = None)
            for i in nodePlus:
                for j in nodeMinus:
                    if (i != j and not (i == depotID and j == dupDepotID)):
                        g.add_edge(i, j, travelDist = tau[i, j] - pi[i])
            return g
        def _initializeLabel():
            return {
                'path': [depotID],
                'dist': 0,
                'load': 0,
                'time': 0
            }
        @runtime("feasibility")
        def _feasibility(curPath, nextNode):
            # Check if node has been covered
            if (nextNode in curPath['path']):
                return False
            # Check loads
            accLoad = curPath['load'] + g.nodes[nextNode]['weight']
            if (accLoad > vehCap):
                return False
            # Check time windows
            lastNode = curPath['path'][-1]
            arrTime = max(
                curPath['time'] + g.edges[lastNode, nextNode]['travelDist'],
                g.nodes[nextNode]['timeWindow'][0])
            availTW = []
            if (nextNode != dupDepotID):
                availTW = g.nodes[nextNode]['timeWindow']
            else:
                availTW = g.nodes[depotID]['timeWindow']
            if (arrTime < availTW[0] or arrTime > availTW[1]):
                return False
            return True
        @runtime("extendPath")
        def _extendPath(curPath, nextNode):
            # Make a copy
            extendPath = {
                'path': [i for i in curPath['path']],
                'time': curPath['time'],
                'load': curPath['load'],
                'dist': curPath['dist']
            }                
            extendPath['path'].append(nextNode)
            
            # Extend the path
            extendPath['dist'] += g.edges[lastNode, nextNode]['travelDist']
            extendPath['load'] += g.nodes[nextNode]['weight']
            extendPath['time'] = max(
                curPath['time'] + g.edges[lastNode, nextNode]['travelDist'],
                g.nodes[nextNode]['timeWindow'][0]) + g.nodes[nextNode]['serviceTime']
            return extendPath
        @runtime("dominate")
        def _dominate(label1, label2):
            # Given two labels
            # - return True if label1 is dominating label2, label2 is useless
            # - return False if label2 is not dominated by label1
            if (label1['path'][-1] == label2['path'][-1]
                and label1['dist'] <= label2['dist']
                and label1['load'] <= label2['load']
                and label1['time'] <= label2['time']):
                return True
            else:
                return False
        # 将queue进行排序
        @runtime("sortQueue")
        def _sortQueue(queue):
            queue = sorted(queue, key = lambda d: d['dist'])
            return queue

        # Main ================================================================
        g = _createGraph()

        # 初始label, queue, path
        queue = [_initializeLabel()]

        path = []
        accSol = 0

        # 主循环
        while (len(queue) > 0):
            # print(len(queue))
            # 从队列中取出第一个
            curPath = queue.pop()

            # 取出来的路径是不是非支配解？
            curPathNonDominatedFlag = True
            for i in range(len(queue)):
                if (_dominate(queue[i], curPath)):
                    curPathNonDominatedFlag = False
                    break

            # 如果取出来的路径不是支配解，说明有可能可以拓展
            if (curPathNonDominatedFlag):
                # 尝试拓展该标签
                lastNode = curPath['path'][-1]
                # 最后的节点的后续节点
                for child in g.successors(lastNode):
                    if (_feasibility(curPath, child)):
                        extendPath = _extendPath(curPath, child)
                        nonDominatedFlag = True
                        for i in range(len(queue)):
                            if (_dominate(queue[i], extendPath)):
                                nonDominatedFlag = False
                                break
                        if (nonDominatedFlag):
                            queue.append(extendPath)
                if (lastNode == dupDepotID):
                    path.append(curPath)

                # queue = _cleanQueue(queue)
                queue = _sortQueue(queue)


        path = sorted(path, key = lambda d: d['dist'])
        optPath = path[0]
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
    # NOTE: 取负号以保持reduce cost的符号一致性
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
            if (kwargs['subproblemAlgo'] == 'IP'):
                subproblem = _pricingIP(pi)
            elif (kwargs['subproblemAlgo'] == 'LabelSetting'):
                subproblem = _pricingLabelSetting(pi)
            
            # 如果子问题有解，且reduce cost小于0
            if (subproblem['ofv'] != None and subproblem['ofv'] - cons[depotID].Pi < -CONST_EPSILON):
                # 首先，将加入标识设为可加入
                canAddVarFlag = True
                writeLog(hyphenStr(""))
                writeLog("Subroute found: " + list2String(subproblem['route']))
                writeLog("OFV: %s" % subproblem['ofv'])
                writeLog("Cost: %s" % subproblem['cr'])
                writeLog("Reduce cost: %s" % (cons[depotID].Pi - subproblem['ofv']))
                writeLog("Iteration time: " + str(round((datetime.datetime.now() - startIter).total_seconds(), 2)) + "[s]")
                
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

        # Re-optimize
        CVRPTW.optimize()

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
                'route': routes[i][:-1],
                'length': c[i]
            }
            solRoute[acc]['route'].append(depotID)
            acc += 1
    ofv = CVRPTW.getObjective().getValue()

    return {
        'ofv': ofv,
        'lb': lb,
        'ub': ofv,
        'gap': (ofv - lb) / lb,
        'routes': solRoute
    }