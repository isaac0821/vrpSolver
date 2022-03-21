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

def insideInterval(
    val:            "The value to be compared with the interval" = None,
    interval:   "List, in the format of [start, end], None if no constraint on that side" = [None, None]    
    ) -> "Given a value `val`, returns true if `val` is inside the interval (or on the edge), false else wise.":

    # Initialize ==============================================================
    [s, e] = interval

    # Check left side =========================================================
    if (s != None and val < s):
        return False

    # Check right side ========================================================
    if (e != None) and val > e:
        return False

    return True

def listSetMinus(a, b):
    return list(set(a) - set(b))

def listSetIntersect(a, b):
    return [value for value in a if value in b]

def listSetUnion(a, b):
    return list(set(a).union(b))

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
