import math
import numpy as np
import random
import pickle

from .const import *

def saveDictionary(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def loadDictionary(name):
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

def getTauEuclidean(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    speed:      "Ratio: Dist / Time" = 1
    ) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":
    tau = {}
    lstNodeID = list(nodes.keys())
    for i in lstNodeID:
        for j in lstNodeID:
            if (i != j):
                t = distEuclidean2D(nodes[i]['loc'], nodes[j]['loc']) / speed
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = CONST_EPSILON
        tau[None, i] = CONST_EPSILON
        tau[i, None] = CONST_EPSILON

    return tau

def getTauSphereEuclidean(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    speed:      "Ratio: Dist / Time" = 1
    ) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":
    tau = {}
    lstNodeID = list(nodes.keys())
    for i in lstNodeID:
        for j in lstNodeID:
            if (i != j):
                t = distSphereEuclidean2D(nodes[i]['loc'], nodes[j]['loc']) / speed
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = CONST_EPSILON
        tau[None, i] = CONST_EPSILON
        tau[i, None] = CONST_EPSILON

    return tau

def calTriangleAreaByEdges(a, b, c):
    # Using Heron's Formula
    s = (a + b + c) / 2
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return area

def distEuclidean2D(
    coord1:     "First coordinate, in (x, y)", 
    coord2:     "Second coordinate, in (x, y)"
    ) -> "Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number":
    if (coord1 != None and coord2 != None):
        return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2)
    else:
        return 0

def distSphereEuclidean2D(
    coord1:     "First coordinate, in (lat, lon)", 
    coord2:     "Second coordinate, in (lat, lon)",
    distUnit:   "Unit of distance\
                 1) String (default) 'mile'\
                 2) String 'meter'\
                 3) String 'kilometer'" = 'mile'
    ) -> "Gives a Euclidean distance based on two lat/lon coords, if two coordinates are the same, return a small number":
    R = None
    if (distUnit == 'mile'):
        R = 3958.8 # Earth radius in miles
    elif (distUnit == 'meter'):
        R = 6371000
    elif (distUnit == 'kilometer'):
        R = 6371
    else:
        print("ERROR: Unrecognized distance unit, options are 'mile', 'meter', 'kilometer'")
        return

    if (coord1 != None and coord2 != None):
        lat1, lon1 = coord1
        lat2, lon2 = coord2
        phi1, phi2 = math.radians(lat1), math.radians(lat2) 
        dphi       = math.radians(lat2 - lat1)
        dlambda    = math.radians(lon2 - lon1)
        a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
        return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    else:
        return CONST_EPSILON

def findComponentsUndirected(
    arcs:       "A list of 2-tuples", 
    ) -> "A list of components":
    # Create adj list, each vertex start with an empty list
    adjList = arcs2AdjList(arcs)

    # Initialize
    found = {}
    for node in adjList:
        found[node] = 0
    components = []

    # Main algorithm, mark neighbors
    for i in adjList:
        comp = []
        q = []
        if (found[i] == 0):
            found[i] = 1
            comp.append(i)
            q.append(i)
            while (q):
                v = q.pop(0)
                for u in adjList[v]:
                    if (found[u] == 0):
                        found[u] = 1
                        comp.append(u)
                        q.append(u)
            components.append(comp)
    return components

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

def valDeg2Vec(
    vVal:       "Norm of the vector",
    vDeg:       "Degree of the vector, 0 as North, in [0, 360)"
    ) -> "Given vector's norm and its degree to North, convert it into a 2-tuple vector":

    vX = None
    vY = None

    while(vDeg < 0):
        vDeg = vDeg + 360

    while(vDeg >= 360):
        vDeg = vDeg - 360

    vX = vVal * math.sin(math.radians(vDeg))
    vY = vVal * math.cos(math.radians(vDeg))

    return (vX, vY)

def vec2ValDeg(
    vec:        "2-tuple, as the vector"
    ) -> "Given a 2-tuple, convert it into a norm and a direction in degree":
    
    (vX, vY) = vec
    
    vDeg = None
    vVal = None
    if (abs(vX) <= 0.0001):
        if (vY >= 0):
            vDeg = 0
            vVal = vY
        elif (vY < 0):
            vDeg = 180
            vVal = -vY
    elif (abs(vY) <= 0.0001):
        if (vX >= 0):
            vVal = vX
            vDeg = 90
        elif (vX < 0):
            vVal = -vX
            vDeg = 270
    else:
        vVal = math.sqrt(vX**2 + vY**2)
        # 1st quad
        if (vX > 0 and vY >= 0):
            vDeg = math.degrees(math.atan(vX / vY))
        # 2nd quad
        elif (vX > 0 and vY < 0):
            vDeg = 180 + math.degrees(math.atan(vX / vY))
        # 3rd quad
        elif (vX < 0 and vY < 0):
            vDeg = 180 + math.degrees(math.atan(vX / vY))
        # 4th quad
        elif (vX < 0 and vY >= 0):
            vDeg = 360 + math.degrees(math.atan(vX / vY))

    return vVal, vDeg

def calVecAddition(
    v1Val:      "Norm of vector 1", 
    v1Deg:      "Degree of the vector 1, 0 as North, in [0, 360)", 
    v2Val:      "Norm of vector 2", 
    v2Deg:      "Degree of the vector 2, 0 as North, in [0, 360)",
    ) -> "Given two vectors' norm and their meteorological degrees, get v3 that v3 = v1 + v2":
    # Change to 2-tuple vector ================================================
    (v1X, v1Y) = valDeg2Vec(v1Val, v1Deg)
    (v2X, v2Y) = valDeg2Vec(v2Val, v2Deg)

    # Get v3 ==================================================================
    (v3X, v3Y) = (v1X + v2X, v1Y + v2Y)

    # Get vector norm and direction ===========================================
    v3Val, v3Deg = vec2ValDeg((v3X, v3Y))

    return v3Val, v3Deg

def calVecSubtract(
    v1Val:      "Norm of vector 1", 
    v1Deg:      "Degree of the vector 1, 0 as North, in [0, 360)", 
    v2Val:      "Norm of vector 2", 
    v2Deg:      "Degree of the vector 2, 0 as North, in [0, 360)", 
    ) -> "Given two vectors' norm and their meteorological degrees, get v3 that v1 = v2 + v3":
    # Change to 2-tuple vector ================================================
    (v1X, v1Y) = valDeg2Vec(v1Val, v1Deg)
    (v2X, v2Y) = valDeg2Vec(v2Val, v2Deg)

    # Get v3 ==================================================================
    (v3X, v3Y) = (v1X - v2X, v1Y - v2Y)

    # Get vector norm and direction ===========================================
    v3Val, v3Deg = vec2ValDeg((v3X, v3Y))

    return v3Val, v3Deg

def getHeading(
    currentLoc: "Current location", 
    goalLoc:    "Targeted location"
    ) -> "Given current location and a goal location, calculate the heading. North is 0-degrees, east is 90-degrees.":
    
    radPreLat = np.radians(currentLoc[0])
    radNextLat = np.radians(goalLoc[0])
    deltaLon = np.radians(goalLoc[1] - currentLoc[1])
    x = np.sin(deltaLon) * np.cos(radNextLat)
    y = (np.cos(radPreLat) * np.sin(radNextLat) - (np.sin(radPreLat) * np.cos(radNextLat) * np.cos(deltaLon)))
    bearingInDegree = np.degrees(np.arctan2(x, y))
    if bearingInDegree < 0:
        bearingInDegree += 360

    return bearingInDegree
