import math

from .common import *

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
                newSeq = [k for k in route]
                newSeq.insert(i + 1, nJ)
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
                newSeq = [k for k in route]
                newSeq.insert(i + 1, nJ)
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
                newSeq = [k for k in route]
                newSeq.insert(i + 1, nJ)
                opt['Insert_Btw_%s_%s' % (nI, nINext)] = {
                    'newSeq': newSeq,
                    'deltaCost': deltaCost,
                    'reversed': False,
                    'newCost': newCost,
                    'newRevCost': newRevCost
                }
            elif (newCost != None and newRevCost != None):
                if (newCost < newRevCost):
                    newSeq = [k for k in route]
                    newSeq.insert(i + 1, nJ)
                    opt['Insert_Btw_%s_%s' % (nI, nINext)] = {
                        'newSeq': newSeq,
                        'deltaCost': deltaCost,
                        'reversed': False,
                        'newCost': newCost,
                        'newRevCost': newRevCost
                    }
                else:
                    newSeq = [k for k in route]
                    newSeq.insert(i + 1, nJ)
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
    N = len(route)
    # nIPrev = route[iterSeq(N, i, 'prev')]
    nIPrev = route[i - 1]
    nI = route[i]
    # nINext = route[iterSeq(N, i, 'next')]
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
    i:          "Start index" = 0,
    j:          "End index" = None
    ) -> "Return the cost on the graph given cost matrix/dictionary tau":
    # Accumulate costs ========================================================
    if (j == None):
        j = len(seq) - 1

    cost = 0
    for k in range(i, j):
        if ((seq[k], seq[k + 1]) in tau):
            cost += tau[seq[k], seq[k + 1]]
        else:
            msgError(seq[k], seq[k + 1])
            return None

    if (i == 0 and j == len(seq) - 1 and closeFlag):
        if ((seq[-1], seq[0]) in tau):
            cost += tau[seq[-1], seq[0]]
        else:
            msgError(seq[-1], seq[0])
            return None
        
    return cost

