from .common import *
from .graph import *

# Note ========================================================================
# 1. seq v.s. route: `seq` is a sequence of visits, may not start/end with the depot, 
#                    `route` is a sequence of visits that starts and ends with the depot
# 2. tau v.s. edge: `tau` is the traveling matrix that considered service time
#                   `edge` could be a string indicating the type of traveling matrices, 
#                   and might need additional information in `edgeArgs`
# 3. i, j v.s. nI, nJ: `i` or `j` indicates the index in a route/seq, and `nI`
#                      or `nJ` indicates the ID of the node
# =============================================================================

def neighborSwap(
    seq:        "A given sequence of vehicle route, assuming this route is feasible",
    tau:        "Traveling cost matrix", 
    i:          "Index in the sequence, 0 <= i <= len(seq) - 2",     
    oldCost:    "Total cost before swapping, if we need to calculate the new cost, input this value, otherwise newCost will be None" = None
    ) -> "Swap the ith and (i+1)th visit in the sequence, if it is applicable, notice that for asymmetric tau, the result is the same":

    # Before: ... --> nIPrev -> nI -> nJ -> nJNext --> ...
    # After:  ... --> nIPrev -> nJ -> nI -> nJNext --> ...
    N = len(seq)
    nIPrev = seq[iterSeq(N, i, 'prev')]
    nI = seq[i]
    nJ = seq[iterSeq(N, i, 'next')]
    nJNext = seq[iterSeq(N, iterSeq(N, i, 'next'), 'next')]

    # newSeq
    newSeq = [k for k in seq]
    j = iterSeq(N, i, 'next')
    newSeq[i], newSeq[j] = newSeq[j], newSeq[i]

    # Check tau
    if ((nIPrev, nJ) not in tau
        or (nJ, nI) not in tau
        or (nI, nJNext) not in tau):
        return {
            'seq': None,
            'deltaCost': None,
            'newCost': None
        }

    # deltaCost = newCost - oldCost
    deltaCost = ((tau[nIPrev, nJ] + tau[nJ, nI] + tau[nI, nJNext])
            - (tau[nIPrev, nI] + tau[nI, nJ] + tau[nJ, nJNext]))
    newCost = None
    if (oldCost != None):
        newCost = oldCost + deltaCost

    return {
        'seq': newSeq,
        'deltaCost': deltaCost,
        'newCost': newCost
    }

def exchange(
    seq:        "A given sequence of vehicle route, assuming this route is feasible", 
    tau:        "Traveling cost matrix", 
    i:          "Index in the sequence, 0 <= i <= len(seq) - 2", 
    j:          "Index in the sequence, 0 <= i <= len(seq) - 2, i != j",     
    oldCost:    "Total cost before swapping, if we need to calculate the new cost, input this value, otherwise newCost will be None" = None
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

    # New seq
    newSeq = [k for k in seq]
    newSeq[i], newSeq[j] = newSeq[j], newSeq[i]

    # Case 1: nINext == nJ
    if (nINext == nJ):
        return neighborSwap(seq, i, tau)

    # Case 2: At least one visit between nI and nJ
    # Check tau
    if ((nIPrev, nJ) not in tau
        or (nJ, nINext) not in tau
        or (nJPrev, nI) not in tau
        or (nI, nJNext) not in tau):
        return {
            'seq': None,
            'deltaCost': None,
            'newCost': None
        }

    # deltaCost = newCost - oldCost
    deltaCost = ((tau[nIPrev, nJ] + tau[nJ, nINext] + tau[nJPrev, nI] + tau[nI, nJNext])
            - (tau[nIPrev, nI] + tau[nI, nINext] + tau[nJPrev, nJ] + tau[nJ, nJNext]))
    newCost = None
    if (oldCost != None):
        newCost = oldCost + deltaCost

    return {
        'seq': newSeq,
        'deltaCost': deltaCost,
        'newCost': newCost
    }

def twoOpt(
    seq:        "A given sequence of vehicle route, assuming this route is feasible", 
    tau:        "Traveling cost matrix", 
    i:          "Index in the sequence, 0 <= i <= len(seq) - 2", 
    j:          "Index in the sequence, 0 <= i <= len(seq) - 2, i.Next() != j",
    oldCost:    "Total cost before swapping, if we need to calculate the new cost, input this value, otherwise newCost will be None" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "Reconnect arc between i, i+1 and j, j+1, if it is applicable":

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

    # Check feasibility
    if (asymFlag):
        for k in range(i, min(iterSeq(N, j, 'next'), N - 1)):
            if ((newSeq[k], newSeq[k + 1] not in tau)):
                return {
                    'seq': None,
                    'deltaCost': None,
                    'newCost': None
                }
    else:
        if ((nI, nJ) not in tau
            or (nINext, nJNext) not in tau):
            return {
                'seq': None,
                'deltaCost': None,
                'newCost': None
            }

    # deltaCost = newCost - oldCost
    # For asymmetric, need to calculate the reverted sub-route
    deltaCost = None
    if (not asymFlag):
        deltaCost = ((tau[nI, nJ] + tau[nINext, nJNext]) 
                - (tau[nI, nINext] + tau[nJ, nJNext]))
    else:
        deltaCost = (calSeqCostMatrix(tau, newSeq, i, iterSeq(N, j, 'next'))
                - calSeqCostMatrix(tau, seq, i, iterSeq(N, j, 'next')))
    newCost = None
    if (oldCost != None):
        newCost = oldCost + deltaCost

    return {
        'seq': newSeq,
        'deltaCost': deltaCost,
        'newCost': newCost
    }

def merge(
    routeI:     "Route I, or abcd, reverse of route I is represented as dcba, must include depot as the first and last element" = None,
    routeJ:     "Route J, or efgh, reverse of route J is represented as hgfe, must include depot as the first and last element" = None,
    depotID:    "Depot ID" = 0,
    tau:        "Traveling cost matrix" = None, 
    costI:      "Length of route I" = None,
    revCostI:   "Length of reverse route I, could be None" = None,
    costJ:      "Length of route J" = None,
    revCostJ:   "Reverse length of route J, could be None" = None,
    demand:     "Dictionary, demand at each node" = None,
    maxDemand:  "Maximum demand, if provided, will filter out routes that exceed the demand limit" = None,
    maxCost:    "Maximum length, if provided, will filter out routes that exceed the max cost limit" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "Given two routes, return the saving of merging them": 

    # Get the first/last non-depot element
    nA = routeI[1]
    nD = routeI[-2]
    nE = routeJ[1]
    nH = routeJ[-2]

    # Check routes start/end with depot
    if (routeI[0] != depotID or routeI[-1] != depotID or routeJ[0] != depotID or routeJ[-1] != depotID):
        msgError("ERROR: sequences for merging need to start and end at the depot")
        return

    # Check new demands
    newDemand = 0
    for i in range(1, len(routeI) - 1):
        newDemand += demand[i]
    for i in range(1, len(routeJ) - 1):
        newDemand += demand[i]
    if (maxDemand != None and newDemand > maxDemand):
        return {
            'newSeq': None,
            'deltaMakespan': None,
            'deltaCost': None,
            'newCost': None,
            'revNewCost': None,
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
        if (maxCost == None or newCost <= maxCost):
            revNewCost = None
            if (asymFlag and (nE, nD) in tau and dcbaFlag and hgfeFlag):
                revNewCost = revCostI + revCostJ - tau[nE, depotID] - tau[depotID, nD] + tau[nE, nD]
                if (maxCost != None and revNewCost > maxCost):
                    revNewCost = None
            opt['ABCD-EFGH'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'revNewCost': revNewCost,
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
        if (maxCost == None or newCost <= maxCost):
            revNewCost = None
            if (asymFlag and (nH, nD) in tau and dcbaFlag):
                revNewCost = revCostI + costJ - tau[nH, depotID] - tau[depotID, nD] + tau[nH, nD]
                if (maxCost != None and revNewCost > maxCost):
                    revNewCost = None
            opt['ABCD-HGFE'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'revNewCost': revNewCost,
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
        if (maxCost == None or newCost <= maxCost):
            revNewCost = None
            if (asymFlag and (nA, nH) in tau and dcbaFlag and hgfeFlag):
                revNewCost = revCostJ + revCostI - tau[nA, depotID] - tau[depotID, nH] + tau[nA, nH]
                if (maxCost != None and revNewCost > maxCost):
                    revNewCost = None
            opt['EFGH-ABCD'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'revNewCost': revNewCost,
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
        if (maxCost == None or newCost <= maxCost):
            revNewCost = None
            if (asymFlag and (nD, nH) in tau and hgfeFlag):
                revNewCost = costI + revCostJ - tau[nD, depotID] - tau[depotID, nH] + tau[nD, nH]
                if (maxCost != None and revNewCost > maxCost):
                    revNewCost = None
            opt['EFGH-DCBA'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'revNewCost': revNewCost,
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
        if (maxCost == None or newCost <= maxCost):
            revNewCost = None
            if (asymFlag and (nE, nA) in tau and hgfeFlag):
                revNewCost = revCostJ + costI - tau[nE, depotID] - tau[depotID, nA] + tau[nE, nA]
                if (maxCost != None and revNewCost > maxCost):
                    revNewCost = None
            opt['DCBA-EFGH'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'revNewCost': revNewCost,
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
        if (maxCost == None or newCost <= maxCost):
            revNewCost = None
            if (asymFlag and (nH, nA) in tau):
                revNewCost = costI + costJ - tau[nH, depotID] - tau[depotID, nA] + tau[nH, nA]
                if (maxCost != None and revNewCost > maxCost):
                    revNewCost = None
            opt['DCBA-HGFE'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'revNewCost': revNewCost,
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
        if (maxCost == None or newCost <= maxCost):
            revNewCost = None
            if (asymFlag and (nA, nE) in tau and dcbaFlag):
                revNewCost = revCostI + costJ - tau[nA, depotID] - tau[depotID, nE] + tau[nA, nE]
                if (maxCost != None and revNewCost > maxCost):
                    revNewCost = None
            opt['HGFE-ABCD'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'revNewCost': revNewCost,
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
        if (maxCost == None or newCost <= maxCost):
            revNewCost = None
            if (asymFlag and (nD, nE) in tau):
                revNewCost = costI + costJ - tau[nD, depotID] - tau[depotID, nE] + tau[nD, nE]
                if (maxCost != None and revNewCost > maxCost):
                    revNewCost = None
            opt['HGFE-DCBA'] = {
                'newSeq': newSeq,
                'deltaCost': deltaCost,
                'newCost': newCost,
                'revNewCost': revNewCost,
                'demand': newDemand
            }

    # Return the best saving - could be bad move
    newSeq = None
    deltaCost = None
    newCost = None
    revNewCost = None
    demand = None
    move = None
    if (len(opt) > 0):
        move = min(opt, key = lambda x: opt[x]['newCost'])
        newSeq = opt[move]['newSeq']
        deltaCost = opt[move]['deltaCost']
        newCost = opt[move]['newCost']
        revNewCost = opt[move]['revNewCost']
        demand = opt[move]['demand']
        move = move
    return {
        'newSeq': newSeq,
        'deltaCost': deltaCost,
        'newCost': newCost,
        'revNewCost': revNewCost,
        'demand': demand,
        'move': move
    }

def calInsertionCost(
    route:      "A given sequence of vehicle route, assuming this route is feasible" = None, 
    cost:       "Cost of the route" = None,
    revCost:    "Reverse length of the route" = None,
    tau:        "Traveling cost matrix" = None, 
    nJ:         "Node to be inserted" = None,
    demand:     "Dictionary, demand at each node" = None,
    maxDemand:  "Maximum demand, if provided, will filter out routes that exceed the demand limit" = None,
    maxCost:    "Maximum length, if provided, will filter out routes that exceed the max cost limit" = None,
    asymFlag:   "True if asymmetric" = None
    ) -> "Insert node i to seq":

    # Initialize ==============================================================
    opt = {}

    # Check new demands
    newDemand = None
    if (maxDemand != None):
        newDemand = 0
        for i in range(1, len(route) - 1):
            newDemand += demand[i]
        newDemand += demand[nJ]
        if (maxDemand != None and newDemand > maxDemand):
            return {
                'newSeq': None,
                'deltaCost': None,
                'newCost': None,
                'revNewCost': None,
                'demand': None
            }

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
    for i in range(1, len(route) - 1):
        nI = route[i]
        nINext = route[i + 1]
        
        # New route
        newSeq = [k for k in route]
        newSeq.insert(i + 1, nJ)

        # Update costs
        newCost = None
        if ((nI, nJ) in tau and (nJ, nINext) in tau):
            # deltaC = newCost - oldCost
            newCost = cost + tau[nI, nJ] + tau[nJ, nINext] - tau[nI, nINext]
            deltaCost = newCost - cost
        revNewCost = None
        if (asymFlag and (nJ, nI) in tau and (nINext, nJ) in tau):
            revNewCost = revCost + tau[nINext, nJ] + tau[nJ, nI] - tau[nINext, nI]

        # Check which direction is better
        if (asymFlag):
            # If only revert route is feasible
            if (newCost == None and revNewCost != None):
                newSeq.reverse()
                deltaCost = revNewCost - cost
                newCost, revNewCost = revNewCost, newCost
                opt['Insert_%s' % i] = {
                    'newSeq': newSeq,
                    'deltaCost': deltaCost,
                    'newCost': newCost,
                    'revNewCost': revNewCost
                }
            elif (newCost != None and revNewCost == None):
                opt['Insert_%s' % i] = {
                    'newSeq': newSeq,
                    'deltaCost': deltaCost,
                    'newCost': newCost,
                    'revNewCost': revNewCost
                }
            elif (newCost != None and revNewCost != None):
                if (newCost < revNewCost):
                    opt['Insert_%s' % i] = {
                        'newSeq': newSeq,
                        'deltaCost': deltaCost,
                        'newCost': newCost,
                        'revNewCost': revNewCost
                    }
                else:
                    newSeq.reverse()
                    deltaCost = revNewCost - cost
                    newCost, revNewCost = revNewCost, newCost 
                    opt['Insert_%s' % i] = {
                        'newSeq': newSeq,
                        'deltaCost': deltaCost,
                        'newCost': newCost,
                        'revNewCost': revNewCost
                    }
        else:
            if (newCost != None):
                opt['Insert_%s' % i] = {
                    'newSeq': newSeq,
                    'deltaCost': deltaCost,
                    'newCost': newCost,
                    'revNewCost': revNewCost
                }

    newSeq = None
    deltaCost = None
    newCost = None
    revNewCost = None
    demand = None
    move = None
    if (len(opt) > 0):
        move = min(opt, key = lambda x: opt[x]['newCost'])
        newSeq = opt[move]['newSeq']
        deltaCost = opt[move]['deltaCost']
        newCost = opt[move]['newCost']
        revNewCost = opt[move]['revNewCost']
        demand = newDemand
        move = move
    return {
        'newSeq': newSeq,
        'deltaCost': deltaCost,
        'newCost': newCost,
        'revNewCost': revNewCost,
        'demand': demand,
        'move': move
    }

def calRemovalSaving(
    route:      "A given sequence of vehicle route, assuming this route is feasible" = None, 
    cost:       "Cost of the route" = None,
    revCost:    "Reverse length of the route" = None,
    tau:        "Traveling cost matrix" = None, 
    asymFlag:   "True if asymmetric" = None
    ) -> "Given a route, returns the saving of removing (potentially beneficial) customer(s)":

    # Initialize ==============================================================
    saving = {}

    # Check saving ============================================================
    for i in range(1, len(route) - 1):
        nIPrev = route[i - 1]
        nI = route[i]
        nINext = route[i + 1]

        # First check saving of removing the customer
        newCost = None
        deltaCost = None
        newSeq = None
        if ((nIPrev, nINext) in tau):
            # deltaCost = newCost - oldCost
            newCost = cost + tau[nIPrev, nINext] - (tau[nIPrev, nI] + tau[nI, nINext])
            deltaCost = newCost - cost
            newSeq = [j for j in route if j != nI]

        # Then, check if the reversed route is feasible and can give better saving
        revNewCost = None
        if (asymFlag):
            if ((nINext, nIPrev) in tau):
                if (revCost != None):
                    revNewCost = revCost + tau[nINext, nIPrev] - (tau[nINext, nI] + tau[nI, nIPrev])

        # Third, try to remove multiple customer that are closed to each other
        pass

        # If asymmetric, revert the sequence if needed
        if (asymFlag and revNewCost != None and newCost != None and revNewCost < newCost):
            newSeq.reverse()
            deltaCost = revNewCost - cost
            newCost, revNewCost = revNewCost, newCost

        # Save the delta
        saving[route[i]] = {
            'newSeq': newSeq,
            'deltaCost': deltaCost,
            'newCost': newCost,
            'revNewCost': revNewCost,
        }

    return saving

