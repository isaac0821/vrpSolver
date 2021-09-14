import math
import numpy as np
import random
import pickle

from .const import *

def saveDictionary(
    obj:        "The dictionary to be saved", 
    name:       "Local file name"
    ) -> "Save a dictionary data to a local .pkl file":
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def loadDictionary(
    name:       "Local file name"
    ) -> "Read a dictionary data from local .pkl file":
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def rndSeq(
    N:          "Integer, Length of sequence",
    s0:         "Integer, Staring index of sequence" = 0,
    closed:     "Boolean, If the sequence is closed, if true, the last element is a duplicate of the first" = False,
    ) -> "Randomly generate sequence starting from `start'":

    seq = [i for i in range(s0, N + s0)]
    # Randomly swap
    for i in range(N):
        j = random.randint(0, N - 1)
        t = seq[i]
        seq[i] = seq[j]
        seq[j] = t

    # Return to start?
    if (closed):
        seq.append(seq[0])
    return seq

def rndPick(
    coefficients: "A list of float numbers as the weight"
    ) -> "Given a list of coefficients, randomly return an index according to that coefficient.":
    totalSum = sum(coefficients)
    tmpSum = 0
    rnd = np.random.uniform(0, totalSum)
    idx = 0
    for i in range(len(coefficients)):
        tmpSum += coefficients[i]
        if rnd <= tmpSum:
            idx = i
            break

    return idx

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
            print("Error: Missing arc (%s, %s) in `weightArcs`" % (seq[i], seq[i + 1]))
            return
        else:
            cost += c

    return cost

def calSeqCostMatrix(
    tau:        "Dictionary {(nodeID1, nodeID2): dist, ...}", 
    seq:        "List, sequence of visiting node ids"
    ) -> "Return the cost on the graph given cost matrix/dictionary tau":

    # Accumulate costs ========================================================
    cost = 0
    for i in range(len(seq) - 1):
        cost += tau[seq[i], seq[i + 1]]

    return cost

def arcs2AdjList(
    arcs:       "1) A list of 3-tuple (nodeID1, nodeID2, weight) or, \
                 2) A list of 2-tuple (nodeID1, nodeID2)"
    ) -> "Dictionary of neighbors of each node":

    neighbors = {}
    for i in range(len(arcs)):
        if (arcs[i][0] not in neighbors):
            neighbors[arcs[i][0]] = [arcs[i][1]]
        else:
            neighbors[arcs[i][0]].append(arcs[i][1])
        if (arcs[i][1] not in neighbors):
            neighbors[arcs[i][1]] = [arcs[i][0]]
        else:
            neighbors[arcs[i][1]].append(arcs[i][0])

    return neighbors

def insideInterval(
    val:            "The value to be compared with the interval" = None,
    interval:   "List, in the format of [start, end], None if no constraint on that side" = [None, None]    
    ) -> "Given a value `val`, returns true if `val` is inside the interval (or on the edge), false else wise.":

    # Initialize ==============================================================
    insideFlag = True
    [s, e] = interval

    # Check left side =========================================================
    if (s != None):
        if (val < s):
            insideFlag = False
            return insideFlag

    # Check right side ========================================================
    if (e != None):
        if (val > e):
            insideFlag = False
            return insideFlag

    return insideFlag
