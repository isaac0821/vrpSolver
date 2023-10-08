import random
import math

try:
    import pickle5 as pickle
except(ImportError):
    import pickle

from .const import *

# Type alias
pt = list[float] | tuple[float, float]
poly = list[list[float]] | list[tuple[float, float]]
polys = list[list[list[float]]] | list[list[tuple[float, float]]]
circle = tuple[pt, float]
arcSeg = list[pt] | tuple[pt, pt, pt|None]
arcPoly = list[arcSeg]
arcPolys = list[arcPoly]
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

def splitList(inputList, binNum):
    listLength = len(inputList)
    perFloor = math.floor(listLength / binNum)
    sizePerBin = [perFloor for i in range(binNum)]
    residual = listLength - sum(sizePerBin)
    for i in range(residual):
        sizePerBin[i] += 1
    bins = []
    acc = 0
    for i in range(len(sizePerBin)):
        bins.append([])
        for k in range(acc, acc + sizePerBin[i]):
            bins[i].append(inputList[k])
        acc += sizePerBin[i]
    return bins