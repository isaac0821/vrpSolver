import datetime

from .common import *
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
