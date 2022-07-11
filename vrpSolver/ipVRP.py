import heapq
import math

from .common import *
from .const import *
from .geometry import *
from .graph import *

def ipVRP(
    nodes:      "Dictionary, returns the detail info of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y), 'capNeed': d, 'tStart': t1, 'tEnd': t2}, \
                        nodeID2: {'loc': (x, y), 'capNeed': d, 'tStart': t1, 'tEnd': t2}, \
                        ... \
             }" = None, 
    depotID:    "Node ID for depot" = 0, 
    customerID: "1) String (default) 'AllButDepot', all nodes exclude the depotID, or \
                 2) A list of node IDs" = 'AllButDepot',
    edges:      "1) String (default) 'Euclidean' or \
                 2) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    vehicles:   "Dictionary, information of vehicles" = None,
    vehCap:     "Capacity of each vehicle" = None,
    vehNum:     "Number of vehicles" = None,
    fml:        "1) String 'Golden77' or \
                 2) String 'Two-Index' or \
                 3) String (not available) 'MultiCommodityFlow' or \
                 4) String (not available) 'SetPartitioning'" = 'Golden77',
    obj:		"1) String 'Makespan'\
    			 2) String 'EdgeCost'" = 'Makespan',
    timeLimit:  "1) Double, in seconds or \
                 2) (default) None, no time limit" = None,
    gapTolerance: "1) Double, Stopping gap, or \
                 2) (default) None, no gap limit" = None,
    outputFlag: "Boolean, True if export the gurobi logs" = False
    ) -> "Exact solution for VRP with vehicle capacity (CVRP)":

    # Check if Gurobi exists ==================================================
    try:
        import gurobipy as grb
    except(ImportError):
        msgError("ERROR: Cannot find Gurobi")
        return

    # Define nodes ==========================================================
    if (customerID == 'AllButDepot'):
        if (depotID not in nodes):
            print("Error: Cannot find depot in nodes")
            return None
        else:
            customerID = [i for i in nodes if i != depotID]

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'LatLon'):
            edges = getTauLatLon(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Derive capNeed ==========================================================
    capNeed = {}
    for i in customerID:
    	if ('capNeed' in nodes[i]):
    		capNeed[i] = nodes[i]['capNeed']
    	else:
    		capNeed[i] = 1

    # Gurobi initialize =======================================================
    CVRP = Model('CVRP')
    if (outputFlag == False):
        CVRP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        CVRP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        CVRP.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Subroutines for different formulation ===================================
    def _ipVRPTwoIndex(obj, depotID, customerID, edges, vehCap, vehNum):
        # Decision variables --------------------------------------------------
        x = {}
        for i, j in edges:
            if (i != depotID and j != depotID and i < j):
                x[i, j] = CVRP.addVar(
                    vtype = grb.GRB.BINARY, name = 'x_%s_%s' % (i, j))
        for j in customerID:
            x[depotID, j] = CVRP.addVar(
                vtype=grb.GRB.INTEGER,
                ub=vehNum, name = 'x_%s_%s' % (depotID, j))

        # Choose objective function -------------------------------------------
        if (obj == 'Makespan'):
            CVRP.setObjective(grb.quicksum())
        elif (obj == 'TotalCost'):


        # CVRP objective function ---------------------------------------------
        CVRP.modelSense = grb.GRB.MINIMIZE
        CVRP.Params.lazyConstraints = 1
        CVRP.update()

        # Leaving depot -------------------------------------------------------
        CVRP.addConstr(grb.quicksum(x[depotID, j] for j in customerID) == 2 * vehNum, name='leaving')

        # Balance constraint --------------------------------------------------
        for cus in customerID:
            CVRP.addConstr(grb.quicksum(x[i, k] for (i, k) in edges if (i < k and k == cus)) + grb.quicksum(x[k, j] for (k, j) in edges if (k < j and k == cus)) == 2, name='balance_%s' % cus)

        # Subtour elimination -------------------------------------------------
        CVRP._x = x
        def subtourelim(model, where):
            if (where == grb.GRB.Callback.MIPSOL):
                x_sol = model.cbGetSolution(model._x)
                arcs = [(i, j) for i, j in x if (i != depotID and j != depotID and i < j and x_sol[i, j] > 0.9)]
                components = getGraphComponents(arcs)
                for comp in components:
                    sumQ = 0
                    for n in comp:
                        sumQ += capNeed[n]
                    vS = math.ceil(sumQ / float(vehCap))
                    edgesInComp = []
                    for (i, j) in x:
                        if (i in comp and j in comp):
                            edgesInComp.append((i, j))
                    model.cbLazy(grb.quicksum(x[i, j] for (i, j) in edgesInComp if i < j) <= len(comp) - vS)

        # CVRP with callback --------------------------------------------------
        CVRP.optimize(subtourelim)

        # Reconstruct solution ------------------------------------------------
        route = {}
        ofv = None
        gap = None
        lb = None
        ub = None
        runtime = None
        if (CVRP.status == grb.GRB.status.OPTIMAL):        
            ofv = CVRP.getObjective().getValue()
            arcSet = []
            for i, j in x:
                if (x[i, j].X > 0.9 and x[i, j].X < 1.1):
                    arcSet.append((i, j))
                elif (x[i, j].X > 1.9):
                    arcSet.append((i, j))
                    arcSet.append((j, i))
            for k in range(1, vehNum + 1):
                route[k] = [depotID]
                while (len(arcSet) > 0):
                    cur = None
                    for arc in arcSet:
                        if (arc[0] == route[k][-1]):
                            route[k].append(arc[1])
                            arcSet.remove(arc)
                            cur = arc[1]
                            break
                        if (arc[1] == route[k][-1]):
                            route[k].append(arc[0])
                            arcSet.remove(arc)
                            cur = arc[0]
                            break
                    if (cur == depotID):
                        break
            gap = 0
            lb = ofv
            ub = ofv
            runtime = CVRP.runtime
        elif (CVRP.status == grb.GRB.status.TIME_LIMIT):
            ofv = None
            route = {}
            gap = CVRP.MIPGap
            lb = CVRP.ObjBoundC
            ub = CVRP.ObjVal
            runtime = CVRP.Runtime

        return {
            'ofv': ofv,
            'route': route,
            'gap': gap,
            'lb': lb,
            'ub': ub,
            'runtime': runtime
        }

    def _ipVRPGolden77(depotID, customerID, edges, vehCap, vehNum):
        # Create adjList for inflow and outflow -------------------------------
        adjIn = {}
        adjOut = {}
        for e in edges:
            if (e[0] not in adjOut):
                adjOut[e[0]] = [e[1]]
            else:
                adjOut[e[0]].append(e[1])
            if (e[1] not in adjIn):
                adjIn[e[1]] = [e[0]]
            else:
                adjIn[e[1]].append(e[0])

        # Decision variables  -------------------------------------------------
        vehicles = [i for i in range(vehNum)]
        x = {}
        for k in vehicles:
            for i, j in edges:
                x[k, i, j] = CVRP.addVar(
                    vtype = grb.GRB.BINARY,
                    obj = edges[i, j])

        # CVRP objective function ---------------------------------------------
        CVRP.modelSense = grb.GRB.MINIMIZE
        CVRP.Params.lazyConstraints = 1
        CVRP.update()

        # Leaving Node --------------------------------------------------------
        for i in customerID:
            CVRP.addConstr(grb.quicksum(grb.quicksum(x[k, i, j] for j in adjOut[i]) and i != j for k in vehicles) == 1)

        # Balance Constraint --------------------------------------------------
        for i in nodes:
            for k in vehicles:
                CVRP.addConstr(grb.quicksum(x[k, i, j] for j in adjOut[i]) == grb.quicksum(x[k, j, i] for j in adjIn[i]))

        # Vehicle selecting route ---------------------------------------------
        for k in vehicles:
            CVRP.addConstr(grb.quicksum(x[k, depotID, i] for i in adjOut[depotID]) <= 1)
            CVRP.addConstr(grb.quicksum(x[k, i, depotID] for i in adjIn[depotID]) <= 1)

        # Capacity Constraint -------------------------------------------------
        for k in vehicles:
            CVRP.addConstr(grb.quicksum(capNeed[i] * grb.quicksum(x[k, i, j] for j in adjOut[i]) for i in customerID) <= vehCap)

        # Sub-tour elimination ------------------------------------------------
        CVRP._x = x
        def subtourelim(model, where):
            if (where == grb.GRB.Callback.MIPSOL):
                x_sol = model.cbGetSolution(model._x)
                for k in vehicles:
                    arcs = [(i, j) for i, j in edges if (x_sol[k, i, j] > 0.9)]
                    components = getGraphComponents(arcs)
                    if (len(components) != 1):
                        for component in components:                    
                            model.cbLazy(grb.quicksum(x[k, i, j] for i in component for j in component if i != j) <= len(component) - 1)
                    if (len(components) == 1 and depotID not in components[0]):
                        component = components[0]
                        model.cbLazy(grb.quicksum(x[k, i, j] for i in component for j in component if i != j) <= len(component) - 1)
        
        # CVRP with callback --------------------------------------------------
        CVRP.optimize(subtourelim)

        # Reconstruct solution ------------------------------------------------
        route = {}
        ofv = None
        gap = None
        lb = None
        ub = None
        runtime = None
        if (CVRP.status == grb.GRB.status.OPTIMAL):        
            ofv = CVRP.getObjective().getValue()
            for k in vehicles:
                arcSet = []
                for i, j in edges:
                    if (x[k, i, j].X > 0.9):
                        print("x[%s, %s, %s] = %s" % (k, i, j, x[k, i, j].X))
                        arcSet.append((i, j))
                route[k] = [depotID]
                while (len(arcSet) > 0):
                    for arc in arcSet:
                        if (arc[0] == route[k][-1]):
                            route[k].append(arc[1])
                            arcSet.remove(arc)
                            break
                        if (arc[1] == route[k][-1]):
                            route[k].append(arc[0])
                            arcSet.remove(arc)
                            break
            gap = 0
            lb = ofv
            ub = ofv
            runtime = CVRP.runtime
        elif (CVRP.status == grb.GRB.status.TIME_LIMIT):
            ofv = None
            route = {}
            gap = CVRP.MIPGap
            lb = CVRP.ObjBoundC
            ub = CVRP.ObjVal
            runtime = CVRP.Runtime

        return {
            'ofv': ofv,
            'route': route,
            'gap': gap,
            'lb': lb,
            'ub': ub,
            'runtime': runtime
        }

    # Solve by different formulations =========================================
    res = None
    if (fml == 'Golden77'):
        res = _ipVRPGolden77(depotID, customerID, edges, vehCap, vehNum)
    elif (fml == 'Two-Index'):
        res = _ipVRPTwoIndex(depotID, customerID, edges, vehCap, vehNum)
    else:
        print("Error: Incorrect or unavailable CVRP formulation option!")

    return res
