import random

try:
    import pickle5 as pickle
except(ImportError):
    import pickle

from .const import *

# Type alias
pt = list[float] | tuple[float, float]
poly = list[pt]
polys = list[poly]
line = list[pt]

def saveDictionary(obj, name: str) -> None:
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def loadDictionary(name: str) -> None:
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def rndPick(coefficients: list[int | float]) -> int:
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

def binary2StartEndPair(l, trueValue = True):
    p = []
    start = None
    if (l[0] == trueValue):
        start = 0
    for i in range(len(l)):
        if (l[i] == trueValue):
            if (start == None):
                start = i
        else:
            if (start != None):
                p.append([start, i - 1])
                start = None
    if (l[-1] == trueValue):
        if (start != None):
            p.append([start, len(l) - 1])
    return p

def iterSeq(seqL, i, direction):
    if (direction == 'next'):
        return i + 1 if i < seqL - 1 else 0
    elif (direction == 'prev'):
        return i - i if i > 0 else seqL - 1
    else:
        return None

def insideInterval(val: float, interval: list[float | None]) -> bool:
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
    
def hyphenStr(s, length=75, sym='-'):
    lenMidS = len(s)
    if (s == ""):
        return length * sym
    elif (lenMidS + 2 < length):
        lenLeftS = (int)((length - lenMidS - 2) / 2)
        lenRightS = length - lenMidS - lenLeftS - 2
        return (lenLeftS * sym) + " " + s + " " + (lenRightS * sym)
    else:
        return s
