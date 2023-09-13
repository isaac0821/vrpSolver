import math

from .common import *

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

