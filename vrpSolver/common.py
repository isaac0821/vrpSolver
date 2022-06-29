import math
import random
import pickle5 as pickle

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

def rndPick(
    coefficients: "A list of float numbers as the weight"
    ) -> "Given a list of coefficients, randomly return an index according to that coefficient.":
    totalSum = sum(coefficients)
    tmpSum = 0
    rnd = random.uniform(0, totalSum)
    idx = 0
    for i in range(len(coefficients)):
        tmpSum += coefficients[i]
        if rnd <= tmpSum:
            idx = i
            break
    return idx

def iterSeq(seqL, i, direction):
    if (direction == 'next'):
        return i + 1 if i < seqL - 1 else 0
    elif (direction == 'prev'):
        return i - i if i > 0 else seqL - 1
    else:
        return None

def insideInterval(
    val:        "The value to be compared with the interval" = None,
    interval:   "List, in the format of [start, end], None if no constraint on that side" = [None, None]    
    ) -> "Given a value `val`, returns true if `val` is inside the interval (or on the edge), false else wise.":
    [s, e] = interval
    if (s != None and val < s):
        return False
    if (e != None and val > e):
        return False
    return True

def listSetMinus(a, b):
    return [v for v in a if v not in b]

def listSetIntersect(a, b):
    return [v for v in a if v in b]

def listSetUnion(a, b):
    l = [v for v in a]
    l.extend([v for v in b if v not in a])
    return l

def list2String(l):
    listString = "["
    listString += ', '.join([list2String(elem) if type(elem) == list else str(elem) for elem in l.copy()])
    listString += "]"
    return listString

def list2Tuple(l):
    sortedList = [i for i in l]
    sortedList.sort()
    tp = tuple(sortedList)
    return tp

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
            return None

    if (closeFlag):
        if ((seq[-1], seq[0]) in tau):
            cost += tau[seq[-1], seq[0]]
        else:
            return None
        
    return cost

