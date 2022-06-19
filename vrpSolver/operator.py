from .common import *
from .graph import *

# Definitions =================================================================
# 1. seq v.s. route: `seq` is a sequence of visits, may not start/end with the depot, 
#                    `route` is a sequence of visits that starts and ends with the depot
# 2. tau v.s. edge: `tau` is the traveling matrix that considered service time
#                   `edge` could be a string indicating the type of traveling matrices, 
#                   and might need additional information in `edgeArgs`
# 3. i, j v.s. nI, nJ: `i` or `j` indicates the index in a route/seq, and `nI`
#                      or `nJ` indicates the ID of the node
# =============================================================================

# Note ========================================================================
# The operations here does not involve any constraints checking
# =============================================================================

def swapNeighbor(
    seq:        "A given sequence of vehicle route, assuming this route is feasible",
    tau:        "Traveling cost matrix", 
    i:          "Index in the sequence, 0 <= i <= len(seq) - 2",     
    cost:       "Total cost before swapping, if we need to calculate the new cost, input this value, otherwise newCost will be None" = None,
    revCost:    "Reversed cost before swapping, it is possible that after swapping, the revered route is preferable" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "Swap the ith and (i+1)th visit in the sequence, if it is applicable, notice that for asymmetric tau, the result is the same":

    # Before: ... --> nIPrev -> nI -> nJ -> nJNext --> ...
    # After:  ... --> nIPrev -> nJ -> nI -> nJNext --> ...
    N = len(seq)
    nIPrev = seq[iterSeq(N, i, 'prev')]
    nI = seq[i]
    nJ = seq[iterSeq(N, i, 'next')]
    nJNext = seq[iterSeq(N, iterSeq(N, i, 'next'), 'next')]

    # Check if can swap
    # Flag that swapping is feasible
    canSwapFlag = True
    # Flag that swapping in the reversed the seq is feasible and need to consider, will be False if not asymmetric
    canRevSwapFlag = asymFlag
    # If any one of the new links is not available, swapping is not available
    # NOTE: do not need to provide `cost`
    if ((nIPrev, nJ) not in tau or (nJ, nI) not in tau or (nI, nJNext) not in tau):
        canSwapFlag = False
    if (asymFlag):
        # If asymmetric, check the reversed seq after swapping
        if ((nJNext, nI) not in tau or (nI, nJ) not in tau or (nJ, nIPrev) not in tau):
            canRevSwapFlag = False

    # Early quit
    if (not asymFlag and not canSwapFlag):
        return None
    elif (asymFlag and not canSwapFlag and not canRevSwapFlag):
        return None

    # newSeq
    newSeq = [k for k in seq]
    j = iterSeq(N, i, 'next')
    newSeq[i], newSeq[j] = newSeq[j], newSeq[i]
    newRevSeq = []
    if (asymFlag and canRevSwapFlag):
        newRevSeq = [newSeq[len(newSeq) - 1 - i] for i in range(len(newSeq))]

    # Calculate deltaCost
    newCost = None
    newRevCost = None
    reverseFlag = False

    # If not asymmetric, only one situation possible, so directly calculate deltaCost
    if (not asymFlag):
        deltaCost = (tau[nIPrev, nJ] + tau[nJ, nI] + tau[nI, nJNext]
                  - (tau[nIPrev, nI] + tau[nI, nJ] + tau[nJ, nJNext]))
        # For the case that we do not need to calculate the total cost, let them be None
        if (cost == None):
            newCost = None
            newRevCost = None
        else:
            if (cost != None):
                newCost = cost + deltaCost
                newRevCost = None

    # If asymmetric, the following cases need to be considered:
    # 1. newSeq feasible + newRevSeq feasible => deltaCost: min(newRevCost, newCost) - cost
    #    - newCost is better
    #    - newRevCost is better, need to take the reversed seq as a new seq
    # 2. newSeq infeasible + newRevSeq feasible => deltaCost: newRevCost - cost
    #    - need to take the reversed seq as a new seq
    # 3. newSeq feasible + newRevSeq infeasible => deltaCost: newCost - cost
    #    - same as symmetric
    else:
        if (canSwapFlag):
            if (cost != None):
                newCost = (cost + tau[nIPrev, nJ] + tau[nJ, nI] + tau[nI, nJNext]
                        - (tau[nIPrev, nI] + tau[nI, nJ] + tau[nJ, nJNext]))
            else:
                newCost = calSeqCostMatrix(
                    tau = tau,
                    seq = newSeq) 
        if (canRevSwapFlag):
            # NOTE: if revCost != None, it means the reverse seq is feasible, no need to check the links
            if (revCost != None):
                newRevCost = (revCost + tau[nJNext, nI] + tau[nI, nJ] + tau[nJ, nIPrev]
                           - (tau[nJNext, nJ] + tau[nJ, nI] + tau[nI, nIPrev]))
            else:
                newRevCost = calSeqCostMatrix(
                    tau = tau,
                    seq = newRevSeq)
                if (newRevCost == None):
                    canRevSwapFlag = False

        # Case 1
        if (canSwapFlag and canRevSwapFlag):
            deltaCost1 = newCost - cost
            deltaCost2 = newRevCost - cost
            if (deltaCost1 <= deltaCost2):
                deltaCost = deltaCost1
            else:
                deltaCost = deltaCost2
                newSeq = newRevSeq
                reverseFlag = True
                newCost, newRevCost = newRevCost, newCost
        # Case 2
        elif (not canSwapFlag and canRevSwapFlag):
            deltaCost = newRevCost - cost
            newSeq = newRevSeq
            reverseFlag = True
            newCost, newRevCost = newRevCost, newCost
        # Case 3
        elif (canSwapFlag and not canRevSwapFlag):
            deltaCost = newCost - cost

    return {
        'seq': newSeq,
        'deltaCost': deltaCost,
        'reversed': reverseFlag,
        'newCost': newCost,
        'newRevCost': newRevCost
    }

def exchange2Nodes(
    seq:        "A given sequence of vehicle route, assuming this route is feasible", 
    tau:        "Traveling cost matrix", 
    i:          "Index in the sequence, 0 <= i <= len(seq) - 2", 
    j:          "Index in the sequence, 0 <= i <= len(seq) - 2, i != j",     
    cost:       "Total cost before swapping, if we need to calculate the new cost, input this value, otherwise newCost will be None" = None,
    revCost:    "Reversed cost before swapping, it is possible that after swapping, the revered route is preferable" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "Swap the ith and jth visit in the sequence, if it is applicable, notice that for asymmetric tau, the result is the same":

    # Before: ... --> nIPrev -> nI -> nINext --> ... --> nJPrev -> nJ -> nJNext --> ...
    # After:  ... --> nIPrev -> nJ -> nINext --> ... --> nJPrev -> nI -> nJNext --> ...
    N = len(seq)
    nIPrev = seq[iterSeq(N, i, 'prev')]
    nI = seq[i]
    nINext = seq[iterSeq(N, i, 'next')]
    nJPrev = seq[iterSeq(N, j, 'prev')]
    nJ = seq[j]
    nJNext = seq[iterSeq(N, j, 'next')]

    # nINext == nJ
    if (nINext == nJ):
        return swapNeighbor(
            seq = seq,
            tau = tau,
            i = i,
            cost = cost,
            revCost = revCost,
            asymFlag = asymFlag)

    # Check if can exchange
    canExchangeFlag = True
    canRevExchangeFlag = asymFlag
    if ((nIPrev, nJ) not in tau or (nJ, nINext) not in tau or (nJPrev, nI) not in tau or (nI, nJNext) not in tau):
        canExchangeFlag = False
    if (asymFlag):
        if ((nJNext, nI) not in tau or (nI, nJPrev) not in tau or (nINext, nJ) not in tau or (nJ, nIPrev) not in tau):
            canRevExchangeFlag = False

    # Early quit
    if (not asymFlag and not canExchangeFlag):
        return None
    elif (asymFlag and not canExchangeFlag and not canRevExchangeFlag):
        return None

    # New seq
    newSeq = [k for k in seq]
    newSeq[i], newSeq[j] = newSeq[j], newSeq[i]
    newRevSeq = []
    if (asymFlag and canRevExchangeFlag):
        newRevSeq = [newSeq[len(newSeq) - 1 - i] for i in range(len(newSeq))]

    # Calculate deltaCost
    newCost = None
    newRevCost = None
    reverseFlag = False

    # If not asymmetric, only one situation possible, so directly calculate deltaCost
    if (not asymFlag):
        deltaCost = ((tau[nIPrev, nJ] + tau[nJ, nINext] + tau[nJPrev, nI] + tau[nI, nJNext])
                  - (tau[nIPrev, nI] + tau[nI, nINext] + tau[nJPrev, nJ] + tau[nJ, nJNext]))
        # For the case that we do not need to calculate the total cost, let them be None
        if (cost == None):
            newCost = None
            newRevCost = None
        else:
            if (cost != None):
                newCost = cost + deltaCost
                newRevCost = None

    # If asymmetric, the following cases need to be considered:
    # 1. newSeq feasible + newRevSeq feasible => deltaCost: min(newRevCost, newCost) - cost
    #    - newCost is better
    #    - newRevCost is better, need to take the reversed seq as a new seq
    # 2. newSeq infeasible + newRevSeq feasible => deltaCost: newRevCost - cost
    #    - need to take the reversed seq as a new seq
    # 3. newSeq feasible + newRevSeq infeasible => deltaCost: newCost - cost
    #    - same as symmetric
    else:
        if (canExchangeFlag):
            if (cost != None):
                newCost = (cost + tau[nIPrev, nJ] + tau[nJ, nINext] + tau[nJPrev, nI] + tau[nI, nJNext]
                          - (tau[nIPrev, nI] + tau[nI, nINext] + tau[nJPrev, nJ] + tau[nJ, nJNext]))
            else:
                newCost = calSeqCostMatrix(
                    tau = tau,
                    seq = newSeq) 
        if (canRevExchangeFlag):
            # NOTE: if revCost != None, it means the reverse seq is feasible, no need to check the links
            if (revCost != None):
                newRevCost = (revCost + tau[nJNext, nI] + tau[nI, nJPrev] + tau[nINext, nJ] + tau[nJ, nIPrev]
                           - (tau[nJNext, nJ] + tau[nJ, nJPrev] + tau[nINext, nI] + tau[nI, nIPrev]))
            else:
                newRevCost = calSeqCostMatrix(
                    tau = tau,
                    seq = newRevSeq)
                if (newRevCost == None):
                    canRevExchangeFlag = False

        # Case 1
        if (canExchangeFlag and canRevExchangeFlag):
            deltaCost1 = newCost - cost
            deltaCost2 = newRevCost - cost
            if (deltaCost1 <= deltaCost2):
                deltaCost = deltaCost1
            else:
                deltaCost = deltaCost2
                newSeq = newRevSeq
                reverseFlag = True
                newCost, newRevCost = newRevCost, newCost
        # Case 2
        elif (not canExchangeFlag and canRevExchangeFlag):
            deltaCost = newRevCost - cost
            newSeq = newRevSeq
            reverseFlag = True
            newCost, newRevCost = newRevCost, newCost
        # Case 3
        elif (canExchangeFlag and not canRevExchangeFlag):
            deltaCost = newCost - cost

    return {
        'seq': newSeq,
        'deltaCost': deltaCost,
        'reversed': reverseFlag,
        'newCost': newCost,
        'newRevCost': newRevCost
    }

def exchange2Arcs(
    seq:        "A given sequence of vehicle route, assuming this route is feasible", 
    tau:        "Traveling cost matrix", 
    i:          "Index in the sequence, 0 <= i <= len(seq) - 2", 
    j:          "Index in the sequence, 0 <= i <= len(seq) - 2, i.Next() != j",
    cost:       "Total cost before swapping, if we need to calculate the new cost, input this value, otherwise newCost will be None" = None,
    revCost:    "Reversed cost before swapping, it is possible that after swapping, the revered route is preferable" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "2Opt, Reconnect arc between i, i+1 and j, j+1, if it is applicable":

    # Before: ... --> nI -> nINext --> ...abcd... --> nJ     -> nJNext --> ...
    # After:  ... --> nI -> nJ     --> ...dcba... --> nINext -> nJNext --> ...
    N = len(seq)
    nI = seq[i]
    nINext = seq[iterSeq(N, i, 'next')]
    nJ = seq[j]
    nJNext = seq[iterSeq(N, j, 'next')]

    # new seq
    newSeq = []
    newSeq.extend([seq[k] for k in range(i + 1)])
    newSeq.extend([seq[j - k] for k in range(j - i)])
    newSeq.extend([seq[k] for k in range(j + 1, len(seq))])

    # Check if can 2-opt
    can2OptFlag = True
    canRev2OptFlag = asymFlag
    if ((nI, nJ) not in tau or (nINext, nJNext) not in tau):
        can2OptFlag = False
    if (asymFlag):
        # Check dcba
        for k in range(i, min(iterSeq(N, j, 'next'), N - 1)):
            if ((newSeq[k], newSeq[k + 1] not in tau)):
                can2OptFlag = False
                break
        # Check newRevSeq
        for k in range(0, min(i, N - 1)):
            if ((seq[k + 1], seq[k]) not in tau):
                canRev2OptFlag = False
                break
        for k in range(min(iterSeq(N, j, 'next'), N - 1), N - 1):
            if ((seq[k + 1], seq[k]) not in tau):
                canRev2OptFlag = False
                break

    newRevSeq = []
    if (asymFlag and canRev2OptFlag):
        newRevSeq = [newSeq[len(newSeq) - 1 - i] for i in range(len(newSeq))]

    # Early quit
    if (not asymFlag and not can2OptFlag):
        return None
    elif (asymFlag and not can2OptFlag and not canRev2OptFlag):
        return None

    # Calculate deltaCost
    newCost = None
    newRevCost = None
    reverseFlag = False

    # If not asymmetric, only one situation possible, so directly calculate deltaCost
    if (not asymFlag):
        deltaCost = ((tau[nI, nJ] + tau[nINext, nJNext]) 
                - (tau[nI, nINext] + tau[nJ, nJNext]))
        # For the case that we do not need to calculate the total cost, let them be None
        if (cost == None):
            newCost = None
            newRevCost = None
        else:
            if (cost != None):
                newCost = cost + deltaCost
                newRevCost = None

    # If asymmetric, the following cases need to be considered:
    # 1. newSeq feasible + newRevSeq feasible => deltaCost: min(newRevCost, newCost) - cost
    #    - newCost is better
    #    - newRevCost is better, need to take the reversed seq as a new seq
    # 2. newSeq infeasible + newRevSeq feasible => deltaCost: newRevCost - cost
    #    - need to take the reversed seq as a new seq
    # 3. newSeq feasible + newRevSeq infeasible => deltaCost: newCost - cost
    #    - same as symmetric
    else:
        costABCD = calSeqCostMatrix(
            tau = tau,
            seq = seq,
            i = iterSeq(N, i, 'next'),
            j = j)
        costDCBA = calSeqCostMatrix(
            tau = tau,
            seq = newSeq,
            i = iterSeq(N, i, 'next'),
            j = j)

        if (can2OptFlag):
            if (cost != None):
                newCost = (cost + tau[nI, nJ] + tau[nINext, nJNext] + costDCBA
                        - (tau[nI, nINext] + tau[nJ, nJNext] + costABCD))
            else:
                newCost = calSeqCostMatrix(
                    tau = tau,
                    seq = newSeq) 
        if (canRev2OptFlag):
            # NOTE: if revCost != None, it means the reverse seq is feasible, no need to check the links
            if (revCost != None):
                newRevCost = (revCost + tau[nJ, nI] + tau[nJNext, nINext] + costABCD
                           - (tau[nINext, nI] + tau[nJNext, nJ] + costDCBA))
            else:
                newRevCost = calSeqCostMatrix(
                    tau = tau,
                    seq = newRevSeq)
                if (newRevCost == None):
                    canRev2OptFlag = False

        # Case 1
        if (can2OptFlag and canRev2OptFlag):
            deltaCost1 = newCost - cost
            deltaCost2 = newRevCost - cost
            if (deltaCost1 <= deltaCost2):
                deltaCost = deltaCost1
            else:
                deltaCost = deltaCost2
                newSeq = newRevSeq
                reverseFlag = True
                newCost, newRevCost = newRevCost, newCost
        # Case 2
        elif (not can2OptFlag and canRev2OptFlag):
            deltaCost = newRevCost - cost
            newSeq = newRevSeq
            reverseFlag = True
            newCost, newRevCost = newRevCost, newCost
        # Case 3
        elif (can2OptFlag and not canRev2OptFlag):
            deltaCost = newCost - cost

    return {
        'seq': newSeq,
        'deltaCost': deltaCost,
        'reversed': reverseFlag,
        'newCost': newCost,
        'newRevCost': newRevCost
    }

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
        # New seq
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

def calInsertionCost(
    route:      "A given sequence of vehicle route, assuming this route is feasible" = None, 
    tau:        "Traveling cost matrix" = None, 
    nJ:         "Node to be inserted" = None,
    cost:       "Cost of the route" = None,
    revCost:    "Reverse length of the route" = None,
    asymFlag:   "True if asymmetric" = None
    ) -> "Calculate the cost of inserting node nJ into the route":

    # Initialize ==============================================================
    opt = {}

    # FIXME: Potential improvement here - reduce insertion attempts
    # distNI = []
    # for j in range(len(route) - 1):
    #     if ((nJ, nI) in tau):
    #         heapq.heappush(distNI, (tau[nJ, nI], nJ, 'oriDirection'))
    #     if (asymFlag and (nI, nJ) in tau):
    #         heapq.heappush(distNI, (tau[nI, nJ], nJ, 'revDirection'))

    # # Insertion positions
    # insertPosCandi = []
    # if (insertSearchRange == None):
    #     if (not asymFlag):
    #         insertSearchRange = len(route) - 1
    #     else:
    #         insertSearchRange = 2 * len(route) - 2
    # for i in range(min(insertSearchRange, len(distNI))):
    #     insertPosCandi.append(heapq.heappop(distNI))

    # Try to insert between any two existing nodes ============================
    for i in range(0, len(route) - 1):
        # Before: ... --> nI -> xx -> nINext --> ...
        # After:  ... --> nI -> nJ -> nINext --> ...
        nI = route[i]
        nINext = route[i + 1]
        
        # New route
        newSeq = [k for k in route]
        newSeq.insert(i + 1, nJ)

        # Update costs
        newCost = None
        if ((nI, nJ) in tau and (nJ, nINext) in tau):
            # deltaC = newCost - cost
            newCost = cost + tau[nI, nJ] + tau[nJ, nINext] - tau[nI, nINext]
            deltaCost = newCost - cost
        newRevCost = None
        if (asymFlag and (nJ, nI) in tau and (nINext, nJ) in tau):
            newRevCost = revCost + tau[nINext, nJ] + tau[nJ, nI] - tau[nINext, nI]

        # Check which direction is better
        if (not asymFlag):
            if (newCost != None):
                opt['Insert_Btw_%s_%s' % (nI, nINext)] = {
                    'newSeq': newSeq,
                    'deltaCost': deltaCost,
                    'reversed': False,
                    'newCost': newCost,
                    'newRevCost': newRevCost
                }
        else:
            # If only revert route is feasible
            if (newCost == None and newRevCost != None):
                newSeq.reverse()
                deltaCost = newRevCost - cost
                newCost, newRevCost = newRevCost, newCost
                opt['Insert_Btw_%s_%s' % (nINext, nI)] = {
                    'newSeq': newSeq,
                    'deltaCost': deltaCost,
                    'reversed': True,
                    'newCost': newCost,
                    'newRevCost': newRevCost
                }
            elif (newCost != None and newRevCost == None):
                opt['Insert_Btw_%s_%s' % (nI, nINext)] = {
                    'newSeq': newSeq,
                    'deltaCost': deltaCost,
                    'reversed': False,
                    'newCost': newCost,
                    'newRevCost': newRevCost
                }
            elif (newCost != None and newRevCost != None):
                if (newCost < newRevCost):
                    opt['Insert_Btw_%s_%s' % (nI, nINext)] = {
                        'newSeq': newSeq,
                        'deltaCost': deltaCost,
                        'reversed': False,
                        'newCost': newCost,
                        'newRevCost': newRevCost
                    }
                else:
                    newSeq.reverse()
                    deltaCost = newRevCost - cost
                    newCost, newRevCost = newRevCost, newCost 
                    opt['Insert_Btw_%s_%s' % (nINext, nI)] = {
                        'newSeq': newSeq,
                        'deltaCost': deltaCost,
                        'reversed': True,
                        'newCost': newCost,
                        'newRevCost': newRevCost
                    }

    # print(opt)
    if (len(opt) > 0):
        move = min(opt, key = lambda x: opt[x]['deltaCost'])
        newSeq = opt[move]['newSeq']
        deltaCost = opt[move]['deltaCost']
        reverseFlag = opt[move]['reversed']
        newCost = opt[move]['newCost']
        newRevCost = opt[move]['newRevCost']
        return {
            'newSeq': newSeq,
            'deltaCost': deltaCost,
            'reversed': reverseFlag,
            'newCost': newCost,
            'newRevCost': newRevCost,
        }
    else:
        return None

def calRemovalSaving(
    route:      "A given sequence of vehicle route, assuming this route is feasible" = None, 
    tau:        "Traveling cost matrix" = None, 
    i:          "Index in the sequence, 0 <= i <= len(seq) - 2" = None, 
    cost:       "Cost of the route" = None,
    revCost:    "Reverse cost of the route" = None,    
    asymFlag:   "True if asymmetric" = None
    ) -> "Given a route, returns the saving of removing (potentially beneficial) customer(s)":

    # Before: ... --> nIPrev -> nI -> nINext --> ...
    # After:  ... --> nIPrev -> xx -> nINext --> ...
    nIPrev = route[i - 1]
    nI = route[i]
    nINext = route[i + 1]

    # First check saving of removing the customer
    newCost = None
    deltaCost = None
    newSeq = None
    reverseFlag = False
    if ((nIPrev, nINext) in tau):
        # deltaCost = newCost - cost
        newCost = cost + tau[nIPrev, nINext] - (tau[nIPrev, nI] + tau[nI, nINext])
        deltaCost = newCost - cost
        newSeq = [j for j in route if j != nI]

    # Then, check if the reversed route is feasible and can give better saving
    newRevCost = None
    if (asymFlag):
        if ((nINext, nIPrev) in tau):
            if (revCost != None):
                newRevCost = revCost + tau[nINext, nIPrev] - (tau[nINext, nI] + tau[nI, nIPrev])

    # If asymmetric, revert the sequence if needed
    if (asymFlag and newRevCost != None and newCost != None and newRevCost < newCost):
        newSeq.reverse()
        reverseFlag = True
        deltaCost = newRevCost - cost
        newCost, newRevCost = newRevCost, newCost

    return {
        'newSeq': newSeq,
        'deltaCost': deltaCost,
        'reversed': reverseFlag,
        'newCost': newCost,
        'newRevCost': newRevCost,
    }
