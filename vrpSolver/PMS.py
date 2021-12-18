from .common import *
from .const import *
import gurobipy as grb

def ipPMSTW(
    numMachines: "Number of parallel machines" = None,
    jobs:       "Dictionary, with job information, in the format of\
                {\
                    'jID': {\
                        'procT': procT\
                        'tw': [earliestAvail, latestAvail]\
                    }\
                }" = None,
    timeLimit:  "1) Double, in seconds or \
                 2) (default) None, no time limit" = None,
    gapTolerance: "1) Double, Stopping gap, or \
                 2) (default) None, no gap limit" = None,
    outputFlag: "Boolean, True if export the gurobi logs" = False
    ) -> "Exact solution for PMSTW":

    # Initialize ==============================================================
    PMSTW = grb.Model('PMSTW')
    if (outputFlag == False):
        PMSTW.setParam('OutputFlag', 0)
    if (timeLimit != None):
        PMSTW.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    if (gapTolerance != None):
        PMSTW.setParam(grb.GRB.Param.MIPGap, gapTolerance)

    # Big-M ===================================================================
    bigM = 0
    for jID in jobs:
        bigM += jobs[jID]['procT']

    # Sets and parameters =====================================================
    # Parameters are saved in `jobs`
    M = [i for i in range(numMachines)] # Set of Machines (mIDs)
    J = [i for i in jobs]              # Set of jobs (jIDs)
    J_zero = [0]
    J_zero.extend(J)
    J_plus = [i for i in jobs]
    J_plus.append(max(J) + 1)
    J_extend = [0]
    J_extend.extend(J_plus)

    # Decision variables ======================================================
    t = {}              # Start time of job jID
    z = {}              # z[m, i, j], machine mID take job i right before job j

    # Define decision variables ===============================================
    # z
    for mID in M:
        for jIDI in J_zero:
            for jIDJ in J_plus:
                z[mID, jIDI, jIDJ] = PMSTW.addVar(vtype = grb.GRB.BINARY, name = "z_%s_%s_%s" % (mID, jIDI, jIDJ))
    # t
    for jID in J_extend:
        t[jID] = PMSTW.addVar(vtype = grb.GRB.CONTINUOUS, name = "t_%s" % jID)

    # Objective function ======================================================
    # Timespan
    PMSTW.setObjective(t[max(J) + 1])

    # Constraints =============================================================
    # Each job should be assigned to one machine
    for jID in J:        
        PMSTW.addConstr(grb.quicksum(z[mID, jIDI, jID] for mID in M for jIDI in J_zero if jIDI != jID) <= 1)
        PMSTW.addConstr(grb.quicksum(z[mID, jID, jIDJ] for mID in M for jIDJ in J_plus if jIDJ != jID) <= 1)

    # z
    for mID in M:
        for jID in J:
            PMSTW.addConstr(z[mID, jID, jID] == 0)
        PMSTW.addConstr(grb.quicksum(z[mID, 0, jID] for jID in J) == 1)
        PMSTW.addConstr(grb.quicksum(z[mID, jID, max(J) + 1] for jID in J) == 1)
        for jIDJ in J:
            PMSTW.addConstr(grb.quicksum(z[mID, jIDI, jIDJ] for jIDI in J_zero) == grb.quicksum(z[mID, jIDJ, jIDK] for jIDK in J_plus))

    for jIDI in J:
        PMSTW.addConstr(grb.quicksum(z[mID, jIDI, jIDJ] for mID in M for jIDJ in J_plus) == 1)

    for jIDI in J:
        for jIDJ in J_plus:
            if (jIDI != jIDJ):
                PMSTW.addConstr(t[jIDI] + jobs[jIDI]['procT']
                    - bigM * (1 - grb.quicksum(z[mID, jIDI, jIDJ] for mID in M)) <= t[jIDJ])
    for jID in J:
        if (jobs[jID]['tw'][0] != None):
            PMSTW.addConstr(t[jID] >= jobs[jID]['tw'][0])
        if (jobs[jID]['tw'][1] != None):
            PMSTW.addConstr(t[jID] <= jobs[jID]['tw'][1])   

    # TDPMS ===================================================================
    PMSTW.modelSense = grb.GRB.MINIMIZE
    PMSTW.optimize()

    # Interpret solution ======================================================
    ofv = None
    assignment = {}
    # cusAssignment = {}
    if (PMSTW.status == grb.GRB.status.OPTIMAL or PMSTW.status == grb.GRB.status.TIME_LIMIT):
        ofv = PMSTW.getObjective().getValue()
        schedule = {}
        for mID in M:
            schedule[mID] = []
            for jIDI in J_zero:
                for jIDJ in J_plus:
                    if (jIDI != jIDJ and jIDI != 0 and z[mID, jIDI, jIDJ].X > 0.5):
                        schedule[mID].append({
                                'jID': jIDI,
                                'ts': t[jIDI].X,
                                'te': t[jIDI].X + jobs[jIDI]['procT']
                            })
        return {
            'timespan': ofv,
            'schedule': schedule,
            'runtime': PMSTW.Runtime,
            'gap': PMSTW.MIPGap,
            'lb': PMSTW.ObjBoundC
        }

    return None

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
        gap = 0
        lb = ofv
        ub = ofv
        runtime = PMS.Runtime
    elif (PMS.status == grb.GRB.status.TIME_LIMIT):
        ofv = None
        gap = PMS.MIPGap
        lb = PMS.ObjBoundC
        ub = PMS.ObjVal
        runtime = PMS.Runtime

    return {
        'ofv': ofv,
        'gap': gap,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime
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
