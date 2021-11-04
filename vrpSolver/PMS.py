from .common import *
from .const import *
import gurobipy as grb

def ipPMS(
    numMachines: "Number of parallel machines" = None,
    jobLen:     "List, length of each job" = None,
    timeLimit:  "1) Double, in seconds or \
                 2) (default) None, no time limit" = None,
    gapTolerance: "1) Double, Stopping gap, or \
                 2) (default) None, no gap limit" = None,
    outputFlag: "Boolean, True if export the gurobi logs" = False
    ) -> "Exact solution for PMS":

    # Initialize ==============================================================
    PMS = grb.Model('PMS')
    if (outputFlag == False):
        PMS.setParam('OutputFlag', 0)
    if (timeLimit != None):
        PMS.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        PMS.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Decision variables ======================================================
    x = {}
    for i in range(numMachines):
        for j in range(len(jobLen)):
            x[i, j] = PMS.addVar(vtype = grb.GRB.BINARY, name = 'x_%s_%s' % (i, j))
    t = PMS.addVar(vtype = grb.GRB.CONTINUOUS, obj = 1, name = 't')

    # Constraints =============================================================
    # Time span
    for i in range(numMachines):
        PMS.addConstr(grb.quicksum(jobLen[j] * x[i, j] for j in range(len(jobLen))) <= t)
    # Unique machine
    for j in range(len(jobLen)):
        PMS.addConstr(grb.quicksum(x[i, j] for i in range(numMachines)) == 1)

    # PMS objective function ==================================================
    PMS.modelSense = grb.GRB.MINIMIZE
    PMS.update()

    # Optimize ================================================================
    PMS.optimize()

    # Reconstruct solution ====================================================
    ofv = None
    if (PMS.status == grb.GRB.status.OPTIMAL):
        ofv = PMS.getObjective().getValue()

    return {
        'timespan': ofv
    }

def heuPMS(
    numMachines: "Number of parallel machines" = None,
    jobLen:     "List, length of each job" = None,
    algo:       "Algorithm for the heuristic \
                 1) String, 'LPT', Longest Processing Time first" = "LPT"
    ) -> "Exact solution for PMS":

    # Initialize ==============================================================
    ts = {}
    for m in range(numMachines):
        ts[m] = 0

    # Sort job length =========================================================
    jobLen.sort(reverse=True)
    for l in jobLen:
        # Get the machine with shortest timespan
        m = sorted(ts, key = lambda x: ts[x])[0]
        ts[m] += l

    # Get timespan ============================================================
    timespan = max(ts.values())

    return {
        'timespan': timespan
    }