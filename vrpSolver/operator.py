import datetime

from .common import *
from .graph import *
from .calculate import *

# Definitions =================================================================
# 1. tau v.s. edge: `tau` is the traveling matrix that considered service time
#                   `edge` could be a string indicating the type of traveling matrices, 
#                   and might need additional information in `edgeArgs`
# 2. i, j v.s. nI, nJ: `i` or `j` indicates the index in a route/route, and `nI`
#                      or `nJ` indicates the ID of the node
# =============================================================================

# Note ========================================================================
# The operations here does not involve any constraints checking
# =============================================================================

def exchange2Arcs(
    route:      "A given sequence of vehicle route, start/end with the depot, assuming this route is feasible" = None, 
    tau:        "Traveling cost matrix" = None, 
    completeTauFlag: "True if tau is completed, to bypass calculation" = True,
    i:          "Index in the sequence, 0 <= i <= len(route) - 2" = None, 
    j:          "Index in the sequence, 0 <= i <= len(route) - 2, i.Next() != j" = None,
    accDist:    "Accumulated travel distance from depot" = None,
    accRevDist: "Accumulated travel distance for reveres route" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "2Opt, Reconnect arc between i, i+1 and j, j+1, if it is applicable":

    # Before: ... --> nI -> nINext --> ...abcd... --> nJ     -> nJNext --> ...
    # After:  ... --> nI -> nJ     --> ...dcba... --> nINext -> nJNext --> ...
    N = len(route)
    nI = route[i]
    nINext = route[i + 1]
    nJ = route[j]
    nJNext = route[j + 1]

    # Initialize accDist accRevDist ===========================================
    if (accDist == None):
        accDist = [0]
        for i in range(len(route) - 1):
            accDist.append(tau[route[i], route[i + 1]])
    if (accRevDist == None):
        accRevDist = [0]
        for i in range(len(route) - 1):
            accRevDist.append(tau[route[len(route) - 1 - i], route[len(route) - 2 - i]])

    # Check if can 2-opt
    can2OptFlag = True if (accDist[-1] != None) else False
    canRev2OptFlag = asymFlag and (True if (accRevDist[-1] != None) else False)

    # Early quit
    if (not asymFlag and not can2OptFlag):
        return None
    elif (asymFlag and not can2OptFlag and not canRev2OptFlag):
        return None

    # Calculate deltaCost
    newCost = None
    newRevCost = None
    reverseFlag = False

    cost = accDist[-1]
    revCost = accRevDist[0]

    # If not asymmetric, only one situation possible, so directly calculate deltaCost
    if (not asymFlag):
        newCost = (cost - (tau[nI, nINext] + tau[nJ, nJNext])
                        + (tau[nI, nJ] + tau[nINext, nJNext]))
        deltaCost = newCost - cost
        newRevCost = None

    # If asymmetric, the following cases need to be considered:
    # 1. newSeq feasible + newRevSeq feasible => deltaCost: min(newRevCost, newCost) - cost
    #    - newCost is better
    #    - newRevCost is better, need to take the reversed route as a new route
    # 2. newSeq infeasible + newRevSeq feasible => deltaCost: newRevCost - cost
    #    - need to take the reversed route as a new route
    # 3. newSeq feasible + newRevSeq infeasible => deltaCost: newCost - cost
    #    - same as symmetric
    else:
        costABCD = accDist[j] - accDist[i + 1]
        costDCBA = accRevDist[i + 1] - accRevDist[j]

        if (can2OptFlag):
            newCost = (cost - (tau[nI, nINext] + tau[nJ, nJNext] + costABCD)
                            + (tau[nI, nJ] + tau[nINext, nJNext] + costDCBA))
        if (canRev2OptFlag):
            newRevCost = (revCost - (tau[nINext, nI] + tau[nJNext, nJ] + costDCBA)
                                  + (tau[nJ, nI] + tau[nJNext, nINext] + costABCD))

        # Case 1
        if (can2OptFlag and canRev2OptFlag):
            if (newRevCost >= newCost):
                deltaCost = newCost - cost
            else:
                deltaCost = newRevCost - cost
                reverseFlag = True
                newCost, newRevCost = newRevCost, newCost
        # Case 2
        elif (not can2OptFlag and canRev2OptFlag):
            deltaCost = newRevCost - cost
            reverseFlag = True
            newCost, newRevCost = newRevCost, newCost
        # Case 3
        elif (can2OptFlag and not canRev2OptFlag):
            deltaCost = newCost - cost

    newSeq = []
    if (deltaCost < 0):
        newSeq.extend([route[k] for k in range(i + 1)])
        newSeq.extend([route[j - k] for k in range(j - i)])
        newSeq.extend([route[k] for k in range(j + 1, len(route))])
        if (reverseFlag):
            newSeq.reverse()
        return {
            'route': newSeq,
            'deltaCost': deltaCost,
            'reversed': reverseFlag,
            'newCost': newCost,
            'newRevCost': newRevCost
        }
    else:
        return None

# [Constructing]
def exchange3Arcs(
    route:      "A given sequence of vehicle route, assuming this route is feasible", 
    tau:        "Traveling cost matrix", 
    i:          "Index in the sequence, 0 <= i <= len(route) - 2", 
    j:          "Index in the sequence, 0 <= j <= len(route) - 2, i.Next() != j",
    k:          "Index in the sequence, 0 <= k <= len(route) - 2, j.Next() != k",
    cost:       "Total cost before swapping, if we need to calculate the new cost, input this value, otherwise newCost will be None" = None,
    revCost:    "Reversed cost before swapping, it is possible that after swapping, the revered route is preferable" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "Reconnect arc between i, i+1 and j, j+1, if it is applicable":

    # Before: ... --> nI -> nINext --> ...abcd... --> nJ     -> nJNext --> ...efgh... --> nK     -> nKNext --> ...
    # After1: ... --> nI -> nINext --> ...abcd... --> nJ     -> nK     --> ...hgfe... --> nJNext -> nKNext --> ...
    # After2: ... --> nI -> nJ     --> ...dcba... --> nINext -> nJNext --> ...efgh... --> nK     -> nKNext --> ...
    # After3: ... --> nI -> nJ     --> ...dcba... --> nINext -> nK     --> ...hgfe... --> nJNext -> nKNext --> ...
    # After4: ... --> nI -> nJNext --> ...efgh... --> nK     -> nINext --> ...abcd... --> nJ     -> nKNext --> ...
    # After5: ... --> nI -> nJNext --> ...efgh... --> nK     -> nJ     --> ...dcba... --> nINext -> nKNext --> ...
    # After6: ... --> nI -> nK     --> ...hgfe... --> nJNext -> nINext --> ...abcd... --> nJ     -> nKNext --> ...
    # After7: ... --> nI -> nK     --> ...hgfe... --> nJNext -> nJ     --> ...dcba... --> nINext -> nKNext --> ...
    N = len(route)
    nIPrev = route[i - 1]
    nI = route[i]
    nINext = route[i + 1]
    nJPrev = route[j - 1]
    nJ = route[j]
    nJNext = route[j + 1]
    nKPrev = route[k - 1]
    nK = route[k]
    nKNext = route[k + 1]

    revSeq = [route[len(route) - i - 1] for i in range(len(route))]
    costABCD = calSeqCostMatrix(tau, route, i + 1, j)
    costEFGH = calSeqCostMatrix(tau, route, j + 1, k)
    costDCBA = calSeqCostMatrix(tau, revSeq, N - j, N - i + 1)
    costHGFE = calSeqCostMatrix(tau, revSeq, N - k, N - j + 1)

    # Check feasibility
    availDCBA = True
    availHGFE = True
    if (asymFlag):
        for m in range(N - j, N - i + 1):
            if ((revSeq[m], revSeq[m + 1]) not in tau):
                availDCBA = False
                break
        for m in range(N - k, N - j + 1):
            if ((revSeq[m], revSeq[m + 1] not in tau)):
                availHGFE = False
                break

    # Seq 1
    newSeq1 = []
    newSeq1.extend([route[m] for m in range(i + 1)])
    newSeq1.extend([route[j - m] for m in range(j - i)])
    newSeq1.extend([route[m] for m in range(j + 1, len(route))])

    newSeq1 = []
    newSeq1.extend([route[m] for m in range(j + 1)])


    print(route)
    print(i, j, k)
    print(nI, nJ, nK)
    print(newSeq1)
    print(newSeq2)
    print(newSeq3)
    print(newSeq4)
    print(newSeq5)
    print(newSeq6)
    print(newSeq7)

    return

def conseqArcExchange():
    return

# [Constructing]
def merge2Routes(
    routeI:     "Route I, or abcd, reverse of route I is represented as dcba, must include depot as the first and last element" = None,
    routeJ:     "Route J, or efgh, reverse of route J is represented as hgfe, must include depot as the first and last element" = None,
    tau:        "Traveling cost matrix" = None, 
    costI:      "Length of route I" = None,
    revCostI:   "Length of reverse route I, could be None" = None,
    costJ:      "Length of route J" = None,
    revCostJ:   "Reverse length of route J, could be None" = None,
    demand:     "Dictionary, demand at each node" = None,
    vehCap:     "Maximum demand, if provided, will filter out routes that exceed the demand limit" = None,
    vehMaxDist: "Maximum length, if provided, will filter out routes that exceed the max cost limit" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "Given two routes, return the saving of merging them": 

    # Get the first/last non-depot element
    nA = routeI[1]
    nD = routeI[-2]
    nE = routeJ[1]
    nH = routeJ[-2]
    depotID = routeI[0]

    # Check routes start/end with depot
    if (routeI[-1] != depotID or routeJ[0] != depotID or routeJ[-1] != depotID):
        msgError("ERROR: sequences for merging need to start and end at the depot")
        return

    # Check new demands
    newDemand = 0
    for i in range(1, len(routeI) - 1):
        newDemand += demand[i]
    for i in range(1, len(routeJ) - 1):
        newDemand += demand[i]
    if (vehCap != None and newDemand > vehCap):
        return {
            'newSeq': None,
            'deltaMakespan': None,
            'deltaCost': None,
            'newCost': None,
            'newRevCost': None,
            'demand': None
        }

    # Check feasibility of reversed routeI and reversed routeJ
    dcbaFlag = True
    hgfeFlag = True
    if (asymFlag):
        for k in range(len(routeI) - 1):
            if ((len(routeI) - k - 1, len(routeI) - k - 2) not in tau):
                dcbaFlag = False
                break        
        for k in range(len(routeJ) - 1):
            if ((len(routeJ) - k - 1, len(routeJ) - k - 2) not in tau):
                hgfeFlag = False
                break

    # If reverse of route is feasible, but the cost was not provided, calculate it
    if (asymFlag):
        if (revCostI == None and dcbaFlag):
            revCostI = calSeqCostMatrix(tau, [routeI[len(routeI) - i - 1] for i in range(len(routeI))])
        if (revCostJ == None and hgfeFlag):
            revCostJ = calSeqCostMatrix(tau, [routeJ[len(routeJ) - i - 1] for i in range(len(routeJ))])

    # New options
    opt = {}

    # Case 1:
    # ... abcd ... efgh ...
    if ((nD, nE) in tau):
        # New route
        newSeq = []
        newSeq.extend([i for i in routeI[:-1]])
        newSeq.extend([i for i in routeJ[1:]])
        newCost = costI + costJ + tau[nD, nE] - (tau[nD, depotID] + tau[depotID, nE])
        deltaCost = newCost - costI - costJ
        if (vehMaxDist == None or newCost <= vehMaxDist):
            newRevCost = None
            if (asymFlag and (nE, nD) in tau and dcbaFlag and hgfeFlag):
                newRevCost = revCostI + revCostJ - tau[nE, depotID] - tau[depotID, nD] + tau[nE, nD]
                if (vehMaxDist != None and newRevCost > vehMaxDist):
                    newRevCost = None
            opt['ABCD-EFGH'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'newRevCost': newRevCost,
                'demand': newDemand
            }

    # Case 2:
    # ... abcd ... hgfe ...
    if ((nD, nH) in tau and len(routeJ) > 3 and hgfeFlag):
        newSeq = []
        newSeq.extend([i for i in routeI[:-1]])
        newSeq.extend([routeJ[len(routeJ) - i - 1] for i in range(1, len(routeJ))])
        newCost = costI + revCostJ + tau[nD, nH] - (tau[nD, depotID] + tau[depotID, nH])
        deltaCost = newCost - costI - costJ
        if (vehMaxDist == None or newCost <= vehMaxDist):
            newRevCost = None
            if (asymFlag and (nH, nD) in tau and dcbaFlag):
                newRevCost = revCostI + costJ - tau[nH, depotID] - tau[depotID, nD] + tau[nH, nD]
                if (vehMaxDist != None and newRevCost > vehMaxDist):
                    newRevCost = None
            opt['ABCD-HGFE'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'newRevCost': newRevCost,
                'demand': newDemand
            }

    # Case 3:
    # ... efgh ... abcd ...
    if ((nH, nA) in tau):
        newSeq = []
        newSeq.extend([i for i in routeJ[:-1]])
        newSeq.extend([i for i in routeI[1:]])
        newCost = costJ + costI + tau[nH, nA] - (tau[nH, depotID] + tau[depotID, nA])
        deltaCost = newCost - costI - costJ
        if (vehMaxDist == None or newCost <= vehMaxDist):
            newRevCost = None
            if (asymFlag and (nA, nH) in tau and dcbaFlag and hgfeFlag):
                newRevCost = revCostJ + revCostI - tau[nA, depotID] - tau[depotID, nH] + tau[nA, nH]
                if (vehMaxDist != None and newRevCost > vehMaxDist):
                    newRevCost = None
            opt['EFGH-ABCD'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'newRevCost': newRevCost,
                'demand': newDemand
            }

    # Case 4:
    # ... efgh ... dcba ...
    if ((nH, nD) in tau and len(routeI) > 3 and dcbaFlag):
        newSeq = []
        newSeq.extend([i for i in routeJ[:-1]])
        newSeq.extend([routeI[len(routeI) - i - 1] for i in range(1, len(routeI))])
        newCost = costJ + revCostI + tau[nH, nD] - (tau[nH, depotID] + tau[depotID, nD])
        deltaCost = newCost - costI - costJ
        if (vehMaxDist == None or newCost <= vehMaxDist):
            newRevCost = None
            if (asymFlag and (nD, nH) in tau and hgfeFlag):
                newRevCost = costI + revCostJ - tau[nD, depotID] - tau[depotID, nH] + tau[nD, nH]
                if (vehMaxDist != None and newRevCost > vehMaxDist):
                    newRevCost = None
            opt['EFGH-DCBA'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'newRevCost': newRevCost,
                'demand': newDemand
            }

    # Case 5: (For asymmetric)
    # ... dcba ... efgh ...
    if (asymFlag and (nA, nE) in tau and len(routeI) > 3 and dcbaFlag):
        newSeq = []
        newSeq.extend([routeI[len(routeI) - i - 1] for i in range(len(routeI) - 1)])
        newSeq.extend([i for i in routeJ[1:]])
        newCost = revCostI + costJ + tau[nA, nE] - (tau[nA, depotID] + tau[depotID, nE])
        deltaCost = newCost - costI - costJ
        if (vehMaxDist == None or newCost <= vehMaxDist):
            newRevCost = None
            if (asymFlag and (nE, nA) in tau and hgfeFlag):
                newRevCost = revCostJ + costI - tau[nE, depotID] - tau[depotID, nA] + tau[nE, nA]
                if (vehMaxDist != None and newRevCost > vehMaxDist):
                    newRevCost = None
            opt['DCBA-EFGH'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'newRevCost': newRevCost,
                'demand': newDemand
            }

    # Case 6: (For asymmetric)
    # ... dcba ... hgfe ...
    if (asymFlag and (nA, nH) in tau and len(routeI) > 3 and len(routeJ) > 3 and dcbaFlag and hgfeFlag):
        newSeq = []
        newSeq.extend([routeI[len(routeI) - i - 1] for i in range(len(routeI) - 1)])
        newSeq.extend([routeJ[len(routeJ) - i - 1] for i in range(1, len(routeJ))])
        newCost = revCostI + revCostJ + tau[nA, nH] - (tau[nA, depotID] + tau[depotID, nH])
        deltaCost = newCost - costI - costJ
        if (vehMaxDist == None or newCost <= vehMaxDist):
            newRevCost = None
            if (asymFlag and (nH, nA) in tau):
                newRevCost = costI + costJ - tau[nH, depotID] - tau[depotID, nA] + tau[nH, nA]
                if (vehMaxDist != None and newRevCost > vehMaxDist):
                    newRevCost = None
            opt['DCBA-HGFE'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'newRevCost': newRevCost,
                'demand': newDemand
            }

    # Case 7: (For asymmetric)
    # ... hgfe ... abcd ...
    if (asymFlag and (nE, nA) in tau and len(routeJ) > 3 and hgfeFlag):
        newSeq = []
        newSeq.extend([routeJ[len(routeJ) - i - 1] for i in range(len(routeJ) - 1)])
        newSeq.extend([i for i in routeI[1:]])
        newCost = revCostJ + costI + tau[nE, nA] - (tau[nE, depotID] + tau[depotID, nA])
        deltaCost = newCost - costI - costJ
        if (vehMaxDist == None or newCost <= vehMaxDist):
            newRevCost = None
            if (asymFlag and (nA, nE) in tau and dcbaFlag):
                newRevCost = revCostI + costJ - tau[nA, depotID] - tau[depotID, nE] + tau[nA, nE]
                if (vehMaxDist != None and newRevCost > vehMaxDist):
                    newRevCost = None
            opt['HGFE-ABCD'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'newRevCost': newRevCost,
                'demand': newDemand
            }

    # Case 8: (For asymmetric)
    # ... hgfe ... dcba ...
    if (asymFlag and (nE, nD) in tau and len(routeJ) > 3 and len(routeI) > 3 and hgfeFlag and dcbaFlag):
        newSeq = []
        newSeq.extend([routeJ[len(routeJ) - i - 1] for i in range(len(routeJ) - 1)])
        newSeq.extend([routeI[len(routeI) - i - 1] for i in range(1, len(routeI))])
        newCost = revCostJ + revCostI + tau[nE, nD] - (tau[nE, depotID] + tau[depotID, nD])
        deltaCost = newCost - costI - costJ
        if (vehMaxDist == None or newCost <= vehMaxDist):
            newRevCost = None
            if (asymFlag and (nD, nE) in tau):
                newRevCost = costI + costJ - tau[nD, depotID] - tau[depotID, nE] + tau[nD, nE]
                if (vehMaxDist != None and newRevCost > vehMaxDist):
                    newRevCost = None
            opt['HGFE-DCBA'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'newRevCost': newRevCost,
                'demand': newDemand
            }

    # Return the best saving - could be bad move
    newSeq = None
    deltaCost = None
    newCost = None
    newRevCost = None
    demand = None
    move = None
    if (len(opt) > 0):
        move = min(opt, key = lambda x: opt[x]['newCost'])
        newSeq = opt[move]['newSeq']
        deltaCost = opt[move]['deltaCost']
        newCost = opt[move]['newCost']
        newRevCost = opt[move]['newRevCost']
        demand = opt[move]['demand']
        move = move
    return {
        'newSeq': newSeq,
        'deltaCost': deltaCost,
        'newCost': newCost,
        'newRevCost': newRevCost,
        'demand': demand,
        'move': move
    }
