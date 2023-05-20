import shapely
from shapely.geometry import mapping

from .ipTSP import *
from .const import *
from .common import *
from .neighbor import *
from .graph import *
from .geometry import *
from .msg import *

# History =====================================================================
# 20230519 - Initialize
# 20230520 - Naive implementation of CETSP using COPT
# =============================================================================

def heuCETSP(
    nodes: dict,
    solver: dict = {'solver': 'COPT'}
    ) -> dict | None:

    """Use heuristic method to find suboptimal CETSP solution

    Parameters
    ----------

    nodes: dictionary, required, default None
        The coordinates of given nodes, in the following format::
            >>> nodes = {
            ...     'nodeID1': {'loc': (x, y), 'neighbor': poly}, 
            ...     'nodeID2': {'loc': (x, y), 'neighbor': poly}, # ...
            ... }

    Returns
    -------

    dictionary

    """

    if (solver == None or 'solver' not in solver):
        raise MissingParameterError("ERROR: Missing required field `solver`.")    

    # Create convex hull of the nodes, truncate neighborhoods =================
    ch = cutNodesNeighbor(nodes)

    # Create a list of Steiner zones ==========================================
    sz = createSteinerZone(ch)

    # Constructive phase ======================================================
    # Step 1: Create a list of Steiner zones (SZs) ----------------------------
    # NOTE: Start from SZ that has the highest order, greedy set partition
    remain = [i for i in nodes]
    szIDs = []
    for i in range(len(sz)):
        szID = len(sz) - i - 1

        allInRemain = True
        for j in sz[szID]['nodeIDs']:
            if (j not in remain):
                allInRemain = False
        if (allInRemain):
            szIDs.append(szID)
            for k in sz[szID]['nodeIDs']:
                remain.remove(k)

    # Step 2: Get centroids of selected SZs as repPt and do TSP ---------------
    repNodes = {}
    for szID in szIDs:
        repNodes[szID] = {
            'loc': sz[szID]['repPt'],
            'neighbor': sz[szID]['poly'],
            'nodeIDs': sz[szID]['nodeIDs']
        }    

    repTSP = ipTSP(
        nodes = repNodes,
        depotID = min(repNodes),
        edges = {'method': 'Euclidean'},
        fml = 'DFJ_Lazy',
        solver = {'solver': solver['solver'], 'outputFlag': False})

    # Step 3: Transform into shortest path problem ----------------------------
    repSeq = []
    if (repTSP != None):
        repSeq = repTSP['seq']

    # NOTE: For now, naively use all extreme points of neighbor poly as candidate
    ptByStep = []
    for i in range(len(repSeq)):
        ptByStep.append([])
        for p in range(len(repNodes[repSeq[i]]['neighbor'])):
            ptByStep[i].append(repNodes[repSeq[i]]['neighbor'][p])

    repTau = {}
    for i in range(len(ptByStep) - 1):
        for m in range(len(ptByStep[i])):
            for n in range(len(ptByStep[i + 1])):
                repTau[i, m, i + 1, n] = shapely.distance(shapely.Point(ptByStep[i][m]), shapely.Point(ptByStep[i + 1][n]))

    def minCostFlowCOPT(ptByStep, repTau):
        try:
            import coptpy as cp
        except (ImportError):
            print("ERROR: Cannot find COPT")
            return {}

        cost = None

        env = cp.Envr()
        mcf = env.createModel("Minimum Cost Flow")

        # Decision variables
        x = {}
        for i in range(len(ptByStep) - 1):
            for m in range(len(ptByStep[i])):
                for n in range(len(ptByStep[i + 1])):
                    x[i, m, i + 1, n] = mcf.addVar(vtype = cp.COPT.CONTINUOUS, obj = repTau[i, m, i + 1, n], name = "x_%s_%s_%s_%s" % (i, m, i + 1, n))

        # MCF objective
        mcf.ObjSense = cp.COPT.MINIMIZE

        # Degree constraints
        mcf.addConstr(cp.quicksum(x[0, 0, 1, n] for n in range(len(ptByStep[0]))) == 1)
        for i in range(1, len(ptByStep) - 1):
            mcf.addConstr(cp.quicksum(x[i, m, i + 1, n] for m in range(len(ptByStep[i])) for n in range(len(ptByStep[i + 1]))) == 1)
            mcf.addConstr(cp.quicksum(x[i - 1, m, i, n] for m in range(len(ptByStep[i - 1])) for n in range(len(ptByStep[i]))) == 1)
            for m in range(len(ptByStep[i])):
                mcf.addConstr(cp.quicksum(x[i - 1, k, i, m] for k in range(len(ptByStep[i - 1]))) == cp.quicksum(x[i, m, i + 1, n] for n in range(len(ptByStep[i + 1]))))
        mcf.addConstr(cp.quicksum(x[len(ptByStep) - 2, m, len(ptByStep) - 1, 0] for m in range(len(ptByStep[len(ptByStep) - 2]))) == 1)

        # Solve
        mcf.solve()
        seqID = []
        szSeq = []
        if (mcf.status == cp.COPT.OPTIMAL):
            cost = mcf.getObjective().getValue()
            for i in range(len(ptByStep) - 1):
                for m in range(len(ptByStep[i])):
                    for n in range(len(ptByStep[i + 1])):
                        if (x[i, m, i + 1, n].x == 1):
                            seqID.append(m)
                            szSeq.append(repNodes[repSeq[i]]['nodeIDs'])
            seqID.append(0)
            szSeq.append(repNodes[repSeq[0]]['nodeIDs'])

        return {
            'seqID': seqID,
            'szSeq': szSeq,
            'cost': cost
        }
    def minCostFlowGurobi(ptByStep, repTau):
        try:
            import gurobipy as grb
        except (ImportError):
            print("ERROR: Cannot find Gurobi")
            return {}

        cost = None
        mcf = grb.Model("Minimum Cost Flow")

        # Decision variables
        x = {}
        for i in range(len(ptByStep) - 1):
            for m in range(len(ptByStep[i])):
                for n in range(len(ptByStep[i + 1])):
                    x[i, m, i + 1, n] = mcf.addVar(vtype = grb.GRB.CONTINUOUS, obj = repTau[i, m, i + 1, n], name = "x_%s_%s_%s_%s" % (i, m, i + 1, n))

        # MCF objective
        mcf.ObjSense = grb.GRB.MINIMIZE

        # Degree constraints
        mcf.addConstr(grb.quicksum(x[0, 0, 1, n] for n in range(len(ptByStep[0]))) == 1)
        for i in range(1, len(ptByStep) - 1):
            mcf.addConstr(grb.quicksum(x[i, m, i + 1, n] for m in range(len(ptByStep[i])) for n in range(len(ptByStep[i + 1]))) == 1)
            mcf.addConstr(grb.quicksum(x[i - 1, m, i, n] for m in range(len(ptByStep[i - 1])) for n in range(len(ptByStep[i]))) == 1)
            for m in range(len(ptByStep[i])):
                mcf.addConstr(grb.quicksum(x[i - 1, k, i, m] for k in range(len(ptByStep[i - 1]))) == grb.quicksum(x[i, m, i + 1, n] for n in range(len(ptByStep[i + 1]))))
        mcf.addConstr(grb.quicksum(x[len(ptByStep) - 2, m, len(ptByStep) - 1, 0] for m in range(len(ptByStep[len(ptByStep) - 2]))) == 1)

        # Solve
        mcf.optimize()
        seqID = []
        szSeq = []
        if (mcf.status == grb.GRB.OPTIMAL):
            cost = mcf.getObjective().getValue()
            for i in range(len(ptByStep) - 1):
                for m in range(len(ptByStep[i])):
                    for n in range(len(ptByStep[i + 1])):
                        if (x[i, m, i + 1, n].x == 1):
                            seqID.append(m)
                            szSeq.append(repNodes[repSeq[i]]['nodeIDs'])
            seqID.append(0)
            szSeq.append(repNodes[repSeq[0]]['nodeIDs'])

        return {
            'seqID': seqID,
            'szSeq': szSeq,
            'cost': cost
        }    

    mcf = {}
    if (solver['solver'] == 'COPT'):
        mcf = minCostFlowCOPT(ptByStep, repTau)
    elif (solver['solver'] == 'Gurobi'):
        mcf = minCostFlowGurobi(ptByStep, repTau)

    seq = []
    for i in range(len(mcf['seqID'])):
        seq.append(ptByStep[i][mcf['seqID'][i]])

    return {
        'seq': seq,
        'szSeq': mcf['szSeq'],
        'ofv': mcf['cost']
    }