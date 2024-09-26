import random
import math
import datetime
import functools

try:
    import pickle5 as pickle
except(ImportError):
    import pickle

CONST_EPSILON = 0.00001
CONST_EARTH_RADIUS_MILES = 3958.8
CONST_EARTH_RADIUS_METERS = 6378137.0

# Type alias
pt = list[float] | tuple[float, float]
pt3D = list[float] | tuple[float, float, float]
poly = list[list[float]] | list[tuple[float, float]]
polys = list[list[list[float]]] | list[list[tuple[float, float]]]
circle = tuple[pt, float]
line = list[pt]

DEBUG_WRITE_LOG = True
DEBUG_PRINT_LOG = False
DEBUG_LOG_PATH = "log.log"

class UnsupportedInputError(Exception):
    pass

class ZeroVectorError(Exception):
    pass

class InvalidPolygonError(Exception):
    pass

class EmptyError(Exception):
    pass

class MissingParameterError(Exception):
    pass

class KeyExistError(Exception):
    pass

class KeyNotExistError(Exception):
    pass

class OutOfRangeError(Exception):
    pass

class VrpSolverNotAvailableError(Exception):
    pass

def saveDictionary(obj, name: str) -> None:
    """
    Save the dictionary to local file as `.pkl`

    Parameters
    ----------

    obj: dict, required
        The dictionary to be saved
    name: str, required
        The name of local file without `.pkl`
    """
    saveName = name + '.pkl'
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def loadDictionary(name: str) -> None:
    """
    Load a dictionary from a local `.pkl` file

    Parameters
    ----------

    name: str, required
        The name of local file with `.pkl`

    Returns
    -------

    dict
        The dictionary loaded from local

    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
 
def rndPick(coefficients: list[int|float]) -> int:
    """
    Given a list of coefficients, randomly returns an index of the list by coefficient

    Parameters
    ----------

    coefficient: list, required
        A list of probabilities.

    Returns
    -------

    int
        An index randomly selected

    """

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

def list2String(l, noCommaFlag=False):
    listString = "["
    if (noCommaFlag==False):
        listString += ', '.join([list2String(elem) if type(elem) == list else str(elem) for elem in l.copy()])
    else:
        listString += ''.join([list2String(elem) if type(elem) == list else str(elem) for elem in l.copy()])
    listString += "]"
    return listString

def list2Tuple(l):
    sortedList = [i for i in l]
    sortedList.sort()
    tp = tuple(sortedList)
    return tp

def hyphenStr(s="", length=75, sym='-'):
    if (s == ""):
        return length * sym
    lenMidS = len(s)
    if (lenMidS + 2 < length):
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

def writeLog(string, logPath = None):
    if (DEBUG_WRITE_LOG):
        if (logPath == None):
            logPath = DEBUG_LOG_PATH
        f = open(logPath, "a")
        f.write(string + "\n")
        f.close()
    if (DEBUG_PRINT_LOG):
        print(string)
    return

def splitIntoSubSeq(inputList, selectFlag):
    if (len(inputList) != len(selectFlag)):
        raise UnsupportedInputError("ERROR: The length of `inputList` should be the same as the length of `selectFlag`")
    splitSub = []
    sub = []
    for i in range(len(selectFlag)):
        if (selectFlag[i] == True):
            sub.append(inputList[i])
        elif (selectFlag[i] == False):
            if (len(sub) > 0):
                splitSub.append([k for k in sub])
                sub = []
    if (len(sub) > 0):
        splitSub.append([k for k in sub])
    return splitSub
 
def findBoundingBox(
    boundingBox = (None, None, None, None),
    pts: list[pt] = None,
    nodes: dict = None,
    locFieldName = 'loc',
    arcs: dict = None,
    arcFieldName = 'arc',
    arcStartLocFieldName = 'startLoc',
    arcEndLocFieldName = 'endLoc',
    poly: poly = None,
    polys: polys = None, 
    polygons: dict = None,
    anchorFieldName: str = 'anchor',
    polyFieldName: str = 'poly',
    xyReverseFlag: bool = False,
    edgeWidth: float = 0.1):

    """
    Given a list of objects, returns a bounding box of all given objects.

    Parameters
    ----------
    boundingBox: list|tuple, optional, default as None
        An existing bounding box
    pts: list of pts, optional, default as None
        A list of pts
    nodes: dict, optional, default as None
        A `nodes` dictionary
    locFieldName: str, optional, default as 'loc'
        The field in `nodes` indicates locations of nodes
    arcs: dict, optional, default as None
        An `arcs` dictionary
    arcFieldName: str, optional, default as 'arc'
        The field in `arcs` indicates locations of arcs
    poly: poly, optional, default as None
        A poly
    polys: polys, optional, default as None
        A list of polys
    polygons: dict, optional, default as None
        A `polygons` dictionary
    anchorFieldName: str, optional, default as `anchor`
        The field in `polygons` indicates anchor of each polygon
    polyFieldName: str, optional, default as `poly`
        The field in `polygons` indicates polygons
    xyReverseFlag: bool, optional, default as True
        True if x, y is reversed.
    edgeWidth: float, optional, default as 0.1
        The extra space around bounding box

    Returns
    -------
    (float, float, float, float)
        A bounding box

    """

    (xMin, xMax, yMin, yMax) = boundingBox
    allX = []
    allY = []
    if (xMin != None):
        allX.append(xMin)
    if (xMax != None):
        allX.append(xMax)
    if (yMin != None):
        allY.append(yMin)
    if (yMax != None):
        allY.append(yMax)

    if (pts != None):
        for pt in pts:
            allX.append(pt[0])
            allY.append(pt[1])
    if (nodes != None):
        for i in nodes:
            allX.append(nodes[i][locFieldName][0])
            allY.append(nodes[i][locFieldName][1])
    if (arcs != None):
        for i in arcs:
            if (arcFieldName in arcs[i]):
                allX.append(arcs[i][arcFieldName][0][0])
                allX.append(arcs[i][arcFieldName][1][0])
                allY.append(arcs[i][arcFieldName][0][1])
                allY.append(arcs[i][arcFieldName][1][1])
            elif (arcStartLocFieldName in arcs[i] and arcEndLocFieldName in arcs[i]):
                allX.append(arcs[i][arcStartLocFieldName][0])
                allY.append(arcs[i][arcStartLocFieldName][1])
                allX.append(arcs[i][arcEndLocFieldName][0])
                allY.append(arcs[i][arcEndLocFieldName][1])     
    if (poly != None):
        for pt in poly:
            allX.append(pt[0])
            allY.append(pt[1])
    if (polys != None):
        for poly in polys:
            for pt in poly:
                allX.append(pt[0])
                allY.append(pt[1])
    if (polygons != None):
        for p in polygons:
            for pt in polygons[p][polyFieldName]:
                allX.append(pt[0])
                allY.append(pt[1])

    xMin = min(allX) - edgeWidth * abs(max(allX) - min(allX))
    xMax = max(allX) + edgeWidth * abs(max(allX) - min(allX))
    yMin = min(allY) - edgeWidth * abs(max(allY) - min(allY))
    yMax = max(allY) + edgeWidth * abs(max(allY) - min(allY))

    if (xyReverseFlag):
        xMin, xMax, yMin, yMax = yMin, yMax, xMin, xMax

    return (xMin, xMax, yMin, yMax)
 
def findBoundingBox3D(
    boundingBox3D = (None, None, None, None, None, None),
    pts3D: list[pt] = None,
    cone: dict = None,
    arcs: dict = None,
    arcFieldName = 'arc',
    arcStartLocFieldName = 'startLoc',
    arcEndLocFieldName = 'endLoc',
    xyReverseFlag: bool = False,
    edgeWidth: float = 0.1):

    (xMin, xMax, yMin, yMax, zMin, zMax) = boundingBox3D
    allX = []
    allY = []
    allZ = []
    if (xMin != None):
        allX.append(xMin)
    if (xMax != None):
        allX.append(xMax)
    if (yMin != None):
        allY.append(yMin)
    if (yMax != None):
        allY.append(yMax)
    if (zMin != None):
        allZ.append(zMin)
    if (zMax != None):
        allZ.append(zMax)

    if (pts3D != None):
        for pt in pts3D:
            allX.append(pt[0])
            allY.append(pt[1])
            allZ.append(pt[2])
    if (arcs != None):
        for i in arcs:
            if (arcFieldName in arcs[i]):
                allX.append(arcs[i][arcFieldName][0][0])
                allX.append(arcs[i][arcFieldName][1][0])
                allY.append(arcs[i][arcFieldName][0][1])
                allY.append(arcs[i][arcFieldName][1][1])
                allZ.append(arcs[i][arcFieldName][0][2])
                allZ.append(arcs[i][arcFieldName][1][2])
            elif (arcStartLocFieldName in arcs[i] and arcEndLocFieldName in arcs[i]):
                allX.append(arcs[i][arcStartLocFieldName][0])
                allY.append(arcs[i][arcStartLocFieldName][1])
                allZ.append(arcs[i][arcStartLocFieldName][2])
                allX.append(arcs[i][arcEndLocFieldName][0])
                allY.append(arcs[i][arcEndLocFieldName][1])
                allZ.append(arcs[i][arcEndLocFieldName][2])
    if (cone != None):
        # 顶点
        allX.append(cone['center'][0])
        allY.append(cone['center'][1])
        allZ.append(cone['center'][2] if len(cone['center']) >= 3 else 0)
        # 上边缘四个点
        rMax = cone['maxHeight'] * cone['tanAlpha']
        allX.append(cone['center'][0] + rMax)
        allX.append(cone['center'][0] - rMax)
        allY.append(cone['center'][1] + rMax)
        allY.append(cone['center'][1] - rMax)
        allZ.append(cone['center'][2] if len(cone['center']) >= 3 else 0 + cone['maxHeight'])

    xMin = min(allX) - edgeWidth * abs(max(allX) - min(allX))
    xMax = max(allX) + edgeWidth * abs(max(allX) - min(allX))
    yMin = min(allY) - edgeWidth * abs(max(allY) - min(allY))
    yMax = max(allY) + edgeWidth * abs(max(allY) - min(allY))
    zMin = min(allZ) - edgeWidth * abs(max(allZ) - min(allZ))
    zMax = max(allZ) + edgeWidth * abs(max(allZ) - min(allZ))

    if (xyReverseFlag):
        xMin, xMax, yMin, yMax = yMin, yMax, xMin, xMax

    return (xMin, xMax, yMin, yMax, zMin, zMax)


def findFigSize(boundingBox, width = None, height = None):
    """
    Given a bounding box, a width(or height), returns the height(or width) of the figure

    Parameters
    ----------

    boundingBox: 4-tuple, required
        The bounding box of the figure
    width: float|None, optional, default as None
        The desired width of the figure
    height: float|None, optional, default as None
        The desired height of the figure

    Returns
    -------
    float, float
        The (width, height) proportional to bounding box

    """
    (xMin, xMax, yMin, yMax) = boundingBox
    w = None
    h = None
    if (width == None and height == None):
        if (xMax - xMin > yMax - yMin):
            w = 5
            h = 5 * ((yMax - yMin) / (xMax - xMin))
        else:
            w = 5 * ((xMax - xMin) / (yMax - yMin))
            h = 5
    elif (width != None and height == None):
        w = width
        h = width * ((yMax - yMin) / (xMax - xMin))
    elif (width == None and height != None):
        w = height * ((xMax - xMin) / (yMax - yMin))
        h = height
    else:
        w = width
        h = height
    return w, h

globalRuntimeAnalysis = {}
# Support runtime tracking and store in dictionary of at most three level
def runtime(key1, key2=None, key3=None):
    def deco(func):
        @functools.wraps(func)
        def fun(*args, **kwargs):
            t = datetime.datetime.now()
            result = func(*args, **kwargs)
            dt = (datetime.datetime.now() - t).total_seconds()
            global globalRuntimeAnalysis
            if (key1 != None and key2 != None and key3 != None):
                # Store in globalRuntimeAnalysis[key1][key2][key3]
                if (key1 not in globalRuntimeAnalysis):
                    globalRuntimeAnalysis[key1] = {}
                if (key2 not in globalRuntimeAnalysis[key1]):
                    globalRuntimeAnalysis[key1][key2] = {}
                if (key3 not in globalRuntimeAnalysis[key1][key2]):
                    globalRuntimeAnalysis[key1][key2][key3] = [dt, 1]
                else:
                    globalRuntimeAnalysis[key1][key2][key3][0] += dt
                    globalRuntimeAnalysis[key1][key2][key3][1] += 1
            elif (key1 != None and key2 != None and key3 == None):
                # Store in globalRuntimeAnalysis[key1][key2]
                if (key1 not in globalRuntimeAnalysis):
                    globalRuntimeAnalysis[key1] = {}
                if (key2 not in globalRuntimeAnalysis[key1]):
                    globalRuntimeAnalysis[key1][key2] = [dt, 1]
                else:
                    globalRuntimeAnalysis[key1][key2][0] += dt
                    globalRuntimeAnalysis[key1][key2][1] += 1
            elif (key1 != None and key2 == None and key3 == None):
                # Store in globalRuntimeAnalysis[key1]
                if (key1 not in globalRuntimeAnalysis):
                    globalRuntimeAnalysis[key1] = [dt, 1]
                else:
                    globalRuntimeAnalysis[key1][0] += dt
                    globalRuntimeAnalysis[key1][1] += 1
            return result
        return fun
    return deco