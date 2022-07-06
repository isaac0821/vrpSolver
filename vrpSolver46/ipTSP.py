import heapq
import gurobipy as grb
import math

from .common import *
from .const import *
from .geometry import *
from .graph import *

def ipTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) String 'Grid' or \
                 4) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                {\
                    'colRow': (numCol, numRow),\
                    'barriers': [(coordX, coordY), ...], \
                }" = None,
    depotID:    "Depot ID, truck visiting sequence will start from and end with this nodeID" = 0,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    fml:        "1) String (default) 'DFJ_Lazy' or \
                 2) String 'DFJ_PlainLoop' or \
                 3) String 'MTZ' or \
                 4) String 'MultiCommodityFlow' or \
                 5) String 'ShortestPath' or \
                 6) String 'QAP'" = 'DFJ_Lazy',
    timeLimit:  "1) Double, in seconds or \
                 2) (default) None, no time limit" = None,
    gapTolerance: "1) Double, Stopping gap, or \
                 2) (default) None, no gap limit" = None,
    outputFlag: "Boolean, True if export the gurobi logs" = False
    ) -> "Exact solution for TSP":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)
        else:
            print(ERROR_INCOR_NODEIDS)
            return

    # Define tau ==============================================================
    tau = {}
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            tau = getTauEuclidean(nodes, nodeIDs)
        elif (edges == 'LatLon'):
            tau = getTauLatLon(nodes, nodeIDs)
        elif (edges == 'Grid'):
            tau = getTauGrid(nodes, nodeIDs, edgeArgs['colRow'], edgeArgs['barriers'])
        else:
            print(ERROR_INCOR_TAU)
            return None
    else:
        tau = dict(edges)

    # Service time ============================================================
    if (serviceTime != None and serviceTime > 0):
        for (i, j) in tau:
            if (i != depotID and j != depotID and i != j):
                tau[i, j] += serviceTime
            elif (i == depotID or j == depotID and i != j):
                tau[i, j] += serviceTime / 2 
    else:
        serviceTime = 0

    # Gurobi initialize =======================================================
    n = len(nodeIDs)
    TSP = grb.Model('TSP')
    if (outputFlag == False):
        TSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        TSP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Subroutines for different formulations ==================================
    def _ipTSPQAP():
        # Decision variables --------------------------------------------------
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

        # TSP objective function ----------------------------------------------
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

        # Assignment constraints ----------------------------------------------
        for i in range(n):
            TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if j != i) == 1)
        for j in range(n):
            TSP.addConstr(grb.quicksum(x[i, j] for i in range(n) if j != i) == 1)

        # Linearized constraints ----------------------------------------------
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

        # Optimize ------------------------------------------------------------
        TSP.optimize()

        # Reconstruct solution ------------------------------------------------
        ofv = None
        seq = []
        gap = None
        lb = None
        ub = None
        runtime = None
        if (TSP.status == grb.GRB.status.OPTIMAL):
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
    def _ipTSPMultiCommodityFlow():
        # Decision variables --------------------------------------------------
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

        # TSP objective function ----------------------------------------------
        TSP.modelSense = grb.GRB.MINIMIZE
        TSP.update()

        # Degree constraints --------------------------------------------------
        for i in range(n):
            TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
            TSP.addConstr(grb.quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

        # MCF -----------------------------------------------------------------
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

        # TSP -----------------------------------------------------------------
        TSP.optimize()

        # Reconstruct solution ------------------------------------------------
        ofv = None
        seq = []
        arcs = []
        if (TSP.status == grb.GRB.status.OPTIMAL):
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
        elif (TSP.status == grb.GRB.status.TIME_LIMIT):
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
    def _ipTSPShortestPath():
        # Decision variables --------------------------------------------------
        x = {}
        for i in range(n):
            for j in range(n):
                if (j != i):
                    for t in range(n):
                        x[i, j, t] = TSP.addVar(
                            obj = tau[nodeIDs[i], nodeIDs[j]],
                            vtype = grb.GRB.BINARY)

        # Stage constraints ---------------------------------------------------
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

        # TSP -----------------------------------------------------------------
        TSP.optimize()

        # Reconstruct solution ------------------------------------------------
        ofv = None
        seq = []
        arcs = []
        if (TSP.status == grb.GRB.status.OPTIMAL):
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
        elif (TSP.status == grb.GRB.status.TIME_LIMIT):
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
    def _ipTSPMTZ():
        # Decision variables --------------------------------------------------
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

        # TSP objective function ----------------------------------------------
        TSP.modelSense = grb.GRB.MINIMIZE
        TSP.update()

        # Degree constraints --------------------------------------------------
        for i in range(n):
            TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
            TSP.addConstr(grb.quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

        # Sequence constraints ------------------------------------------------
        for i in range(1, n):
            for j in range(1, n):
                if (i != j):
                    TSP.addConstr(u[i] - u[j] + (n - 1) * x[i, j] <= n - 2, name = 'seq_%s_%s' % (i, j))
        for i in range(1, n):
            TSP.addConstr(1 <= u[i])
            TSP.addConstr(u[i] <= n - 1)

        # TSP -----------------------------------------------------------------
        TSP.optimize()

        # Reconstruct solution ------------------------------------------------
        ofv = None
        gap = None
        seq = []
        arcs = []
        if (TSP.status == grb.GRB.status.OPTIMAL):
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
        elif (TSP.status == grb.GRB.status.TIME_LIMIT):
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
    def _ipTSPPlainLoop():
        # Decision variables --------------------------------------------------
        x = {}
        for i in range(n):
            for j in range(n):
                if (i != j):
                    x[i, j] = TSP.addVar(
                        vtype = grb.GRB.BINARY, 
                        obj = tau[nodeIDs[i], nodeIDs[j]], 
                        name = 'x_%s_%s' % (i, j))
                    
        # TSP -----------------------------------------------------------------
        TSP.modelSense = grb.GRB.MINIMIZE
        TSP.Params.lazyConstraints = 1
        TSP.update()

        # Degree constraints --------------------------------------------------
        for i in range(n):
            TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
            TSP.addConstr(grb.quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

        # Resolve to optimality and try to find sub-tours ---------------------
        noSubtourFlag = False
        accRuntime = 0
        while (not noSubtourFlag):
            if (timeLimit != None):
                TSP.setParam(grb.GRB.Param.TimeLimit, timeLimit - accRuntime) # FIXME
            TSP.optimize()
            if (TSP.status == grb.GRB.status.OPTIMAL):
                accRuntime += TSP.Runtime
                arcs = grb.tuplelist((i, j) for i, j in x.keys() if x[i, j].X > 0.9)
                components = getGraphComponents(arcs)
                if (len(components) == 1):
                    noSubtourFlag = True
                    break
                else:
                    for comp in components:
                        TSP.addConstr(grb.quicksum(x[i, j] for i in comp for j in comp if i != j) <= len(comp) - 1)
            elif (TSP.status == grb.GRB.status.TIME_LIMIT):
                accRuntime += TSP.Runtime
                break

        # Reconstruct solution ------------------------------------------------
        ofv = None
        seq = []
        arcs = []
        if (TSP.status == grb.GRB.status.OPTIMAL):
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
        elif (TSP.status == grb.GRB.status.TIME_LIMIT):
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
    def _ipTSPLazyCuts():
        # Decision variables --------------------------------------------------
        x = {}
        for i in range(n):
            for j in range(n):
                if (i != j):
                    x[i, j] = TSP.addVar(
                        vtype = grb.GRB.BINARY, 
                        obj = tau[nodeIDs[i], nodeIDs[j]], 
                        name = 'x_%s_%s' % (i, j))
                    
        # TSP objective function ----------------------------------------------
        TSP.modelSense = grb.GRB.MINIMIZE
        TSP.Params.lazyConstraints = 1
        TSP.update()

        # Degree constraints --------------------------------------------------
        for i in range(n):
            TSP.addConstr(grb.quicksum(x[i, j] for j in range(n) if i != j) == 1, name = 'leave_%s' % i)
            TSP.addConstr(grb.quicksum(x[j, i] for j in range(n) if i != j) == 1, name = 'enter_%s' % i)

        # Sub-tour elimination ------------------------------------------------
        TSP._x = x
        def subtourelim(model, where):
            if (where == grb.GRB.Callback.MIPSOL):
                x_sol = model.cbGetSolution(model._x)
                arcs = grb.tuplelist((i, j) for i, j in model._x.keys() if x_sol[i, j] > 0.9)
                components = getGraphComponents(arcs)
                for component in components:
                    if (len(component) < n):
                        model.cbLazy(grb.quicksum(x[i,j] for i in component for j in component if i != j) <= len(component) - 1)

        # TSP with callback ---------------------------------------------------
        TSP.optimize(subtourelim)

        # Reconstruct solution ------------------------------------------------
        ofv = None
        seq = []
        arcs = []
        if (TSP.status == grb.GRB.status.OPTIMAL):
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
        elif (TSP.status == grb.GRB.status.TIME_LIMIT):
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

    # Solve by different formulations =========================================
    tsp = None
    if (fml == 'DFJ_Lazy'):
        tsp = _ipTSPLazyCuts()
    elif (fml == 'DFJ_PlainLoop'):
        tsp = _ipTSPPlainLoop()
    elif (fml == 'MTZ'):
        tsp = _ipTSPMTZ()
    elif (fml == 'ShortestPath'):
        tsp = _ipTSPShortestPath()
    elif (fml == 'MultiCommodityFlow'):
        tsp = _ipTSPMultiCommodityFlow()
    elif (fml == 'QAP'):
        tsp = _ipTSPQAP()
    else:
        print("Error: Incorrect or not available TSP formulation option!")
        return None
    if (tsp != None):
        tsp['fml'] = fml

    # Fix the sequence to make it start from the depot ========================
    startIndex = None
    seq = [i for i in tsp['seq']]
    truckSeq = []
    for k in range(len(seq)):
        if (seq[k] == depotID):
            startIndex = k
    if (startIndex <= len(seq) - 1):
        for k in range(startIndex, len(seq) - 1):
            truckSeq.append(seq[k])
    if (startIndex >= 0):
        for k in range(0, startIndex):
            truckSeq.append(seq[k])
    truckSeq.append(depotID)
    tsp['seq'] = truckSeq

    # Add service time info ===================================================
    tsp['serviceTime'] = serviceTime

    return tsp
