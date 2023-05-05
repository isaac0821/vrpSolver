import math

from .common import *

def calInsertionCost(
    route:      "A given sequence of vehicle route, starts and ends with depotID, assuming this route is feasible" = None, 
    tau:        "Traveling cost matrix" = None, 
    nJ:         "nodeID to be inserted" = None,
    cost:       "Cost of the route" = None,
    revCost:    "Reverse length of the route" = None,
    asymFlag:   "True if asymmetric" = None
    ) -> "Calculate the cheapest cost of inserting node nJ into the route":

    # Initialize ==============================================================
    opt = None
    sofarCheapestCost = max(cost, revCost)
    insertIndex = None

    # Try to insert between any two existing nodes ============================
    for i in range(0, len(route) - 1):
        # Before: ... --> nI -> xx -> nINext --> ...
        # After:  ... --> nI -> nJ -> nINext --> ...
        nI = route[i]
        nINext = route[i + 1]

        # Update costs
        newCost = None
        if ((nI, nJ) in tau and (nJ, nINext) in tau):
            # deltaC = newCost - cost
            newCost = cost + tau[nI, nJ] + tau[nJ, nINext] - tau[nI, nINext]
        newRevCost = None
        if (asymFlag and (nJ, nI) in tau and (nINext, nJ) in tau):
            newRevCost = revCost + tau[nINext, nJ] + tau[nJ, nI] - tau[nINext, nI]

        # Check which direction is better
        if (not asymFlag):
            if (newCost != None):
                deltaCost = newCost - cost
                if (deltaCost < sofarCheapestCost):
                    sofarCheapestCost = deltaCost
                    insertIndex = i
                    opt = {
                        'deltaCost': deltaCost,
                        'reversed': False,
                        'newCost': newCost,
                        'newRevCost': newRevCost
                    }
        else:
            # If only revert route is feasible
            if (newCost == None and newRevCost != None):
                deltaCost = newRevCost - cost
                if (deltaCost < sofarCheapestCost):
                    sofarCheapestCost = deltaCost
                    insertIndex = i
                    opt = {
                        'deltaCost': deltaCost,
                        'reversed': True,
                        'newCost': newRevCost,
                        'newRevCost': newCost
                    }
            elif (newCost != None and newRevCost == None):
                deltaCost = newCost - cost
                if (deltaCost < sofarCheapestCost):
                    sofarCheapestCost = deltaCost
                    insertIndex = i
                    opt = {
                        'deltaCost': deltaCost,
                        'reversed': False,
                        'newCost': newCost,
                        'newRevCost': newRevCost
                    }
            elif (newCost != None and newRevCost != None and newCost < newRevCost):
                deltaCost = newCost - cost 
                if (deltaCost < sofarCheapestCost):
                    sofarCheapestCost = deltaCost
                    insertIndex = i
                    opt = {
                        'deltaCost': deltaCost,
                        'reversed': False,
                        'newCost': newCost,
                        'newRevCost': newRevCost
                    }
            elif (newCost != None and newRevCost != None and newCost > newRevCost):
                deltaCost = newRevCost - cost
                if (deltaCost < sofarCheapestCost):
                    sofarCheapestCost = deltaCost
                    insertIndex = i
                    opt = {
                        'deltaCost': deltaCost,
                        'reversed': True,
                        'newCost': newRevCost,
                        'newRevCost': newCost
                    }

    if (opt != None):
        newSeq = [k for k in route]
        newSeq.insert(insertIndex + 1, nJ)
        if (opt['reversed']):
            newSeq.reverse()
        opt['newSeq'] = newSeq
        return opt
    return None

def calRemovalSaving(
    route:      "A given sequence of vehicle route, assuming this route is feasible" = None, 
    tau:        "Traveling cost matrix" = None, 
    nI:         "nodeID to be removed" = None,
    cost:       "Cost of the route" = None,
    revCost:    "Reverse cost of the route" = None,    
    asymFlag:   "True if asymmetric" = None
    ) -> "Given a route, returns the saving of removing (potentially beneficial) customer(s)":

    # Before: ... --> nIPrev -> nI -> nINext --> ...
    # After:  ... --> nIPrev -> xx -> nINext --> ...
    nIPrev = None
    nINext = None
    if (nI not in route or nI == route[0] or nI == route[-1]):
        return None
    else:
        i = route.index(nI)
        nIPrev = route[i - 1]
        nINext = route[i + 1]
    opt = None

    # Update costs
    newCost = None
    if ((nIPrev, nINext) in tau):
        newCost = cost + tau[nIPrev, nINext] - (tau[nIPrev, nI] + tau[nI, nINext])
    newRevCost = None
    if (asymFlag and (nINext, nIPrev) in tau):
        newRevCost = revCost + tau[nINext, nIPrev] - (tau[nINext, nI] + tau[nI, nIPrev])

    # Check which direction is better
    if (not asymFlag):
        if (newCost != None):
            deltaCost = newCost - cost
            opt = {
                'deltaCost': deltaCost,
                'reversed': False,
                'newCost': newCost,
                'newRevCost': newRevCost
            }
    else:
        # If only revert route is feasible
        if (newCost == None and newRevCost != None):
            deltaCost = newRevCost - cost
            opt = {
                'deltaCost': deltaCost,
                'reversed': True,
                'newCost': newRevCost,
                'newRevCost': newCost
            }
        elif (newCost != None and newRevCost == None):
            deltaCost = newCost - cost
            opt = {
                'deltaCost': deltaCost,
                'reversed': False,
                'newCost': newCost,
                'newRevCost': newRevCost
            }
        elif (newCost != None and newRevCost != None and newCost < newRevCost):
            deltaCost = newCost - cost 
            opt = {
                'deltaCost': deltaCost,
                'reversed': False,
                'newCost': newCost,
                'newRevCost': newRevCost
            }
        elif (newCost != None and newRevCost != None and newCost > newRevCost):
            deltaCost = newRevCost - cost
            opt = {
                'deltaCost': deltaCost,
                'reversed': True,
                'newCost': newRevCost,
                'newRevCost': newCost
            }

    if (opt != None):
        newSeq = [k for k in route]
        newSeq.remove(nI)
        if (opt['reversed']):
            newSeq.reverse()
        opt['newSeq'] = newSeq
        return opt
    return None

def calSeqCostArcs(
    weightArcs: "A list of 3-tuple (nodeID1, nodeID2, weight)",
    seq:        "List, sequence of visiting node ids"
    ) -> "Return the cost on the graph given a list of arcs weights":

    # Accumulate costs ========================================================
    cost = 0
    for i in range(len(seq) - 1):
        c = None
        for j in range(len(weightArcs)):
            if (seq[i] == weightArcs[j][0] and seq[i + 1] == weightArcs[j][1]):
                c = weightArcs[j][2]
                break
            elif (seq[i] == weightArcs[j][1] and seq[i + 1] == weightArcs[j][0]):
                c = weightArcs[j][2]
                break
        if (c == None):
            msgError("Error: Missing arc (%s, %s) in `weightArcs`" % (seq[i], seq[i + 1]))
            return
        else:
            cost += c

    return cost

def calSeqCostMatrix(
    tau:        "Dictionary {(nodeID1, nodeID2): dist, ...}", 
    seq:        "List, sequence of visiting node ids",
    closeFlag:  "True if the seq is closed" = None,
    reverseFlag: "True if calculates the cost of the reversed seq" = False,
    i:          "Start index" = 0,
    j:          "End index" = None
    ) -> "Return the cost on the graph given cost matrix/dictionary tau":
    # Accumulate costs ========================================================
    if (j == None):
        j = len(seq) - 1

    cost = 0
    for k in range(i, j):
        if (not reverseFlag):
            if ((seq[k], seq[k + 1]) in tau):            
                cost += tau[seq[k], seq[k + 1]]
            else:
                return None
        else:
            if ((seq[len(seq) - k - 1], seq[len(seq) - k - 2]) in tau):
                cost += tau[seq[len(seq) - k - 1], seq[len(seq) - k - 2]]
            else:
                return None

    if (i == 0 and j == len(seq) - 1 and closeFlag):
        if (not reverseFlag):
            if ((seq[-1], seq[0]) in tau):
                cost += tau[seq[-1], seq[0]]
            else:
                return None
        else:
            if ((seq[0], seq[-1]) in tau):
                cost += tau[seq[0], seq[-1]]
            else:
                return None
        
    return cost

