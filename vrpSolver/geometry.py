import math
import heapq
import geopy.distance

from .const import *
from .common import *
from .matrices import *

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
    distUnit:   "Unit of distance\
                 1) String (default) 'mile'\
                 2) String 'meter'\
                 3) String 'kilometer'" = 'mile',
    speed:      "Ratio: Dist / Time" = 1
    ) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":
    tau = {}
    lstNodeID = list(nodes.keys())
    for i in lstNodeID:
        for j in lstNodeID:
            if (i != j):
                t = distSphereEuclidean2D(nodes[i]['loc'], nodes[j]['loc'], distUnit) / speed
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = CONST_EPSILON
        tau[None, i] = CONST_EPSILON
        tau[i, None] = CONST_EPSILON

    return tau

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
        dphi = math.radians(lat2 - lat1)
        dlambda = math.radians(lon2 - lon1)
        a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
        return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    else:
        return CONST_EPSILON

def calTriangleAreaByEdges(a, b, c):
    # Using Heron's Formula
    s = (a + b + c) / 2
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return area

def calTriangleAreaByCoords(loc1, loc2, loc3):
    x1 = loc1[0]
    x2 = loc2[0]
    x3 = loc3[0]
    y1 = loc1[1]
    y2 = loc2[1]
    y3 = loc3[1]
    val = (x2 * y3 + x3 * y1 + x1 * y2) - (x2 * y1 + x3 * y2 + x1 * y3)
    area = abs(val)

    return area

def getHeadingXY(
    coord1:     "Current location", 
    coord2:     "Targeted location"
    ) -> "Given current location and a goal location, calculate the heading. Up is 0-degrees, clock-wise":
    
    vec = (coord2[0] - coord1[0], coord2[1] - coord1[1])
    vVal, vDeg = vec2ValDeg(vec)

    return vDeg

def getHeadingLatLon(
    coord1:     "Current location as in [lat, lon]", 
    coord2:     "Target location as in [lat, lon]"
    ) -> "Given current location and a goal location, calculate the heading. Up is 0-degrees, clock-wise":
    # Ref: https://github.com/manuelbieh/Geolib/issues/28
    
    lat1, lon1 = coord1
    lat2, lon2 = coord2
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)
    
    dLon = lon2 - lon1
    if abs(dLon) > math.pi:
        if dLon > 0.0:
            dLon = -(2.0 * math.pi - dLon)
        else:
            dLon = (2.0 * math.pi + dLon)

    tan1 = math.tan(lat1 / 2.0 + math.pi / 4.0)
    tan2 = math.tan(lat2 / 2.0 + math.pi / 4.0)
    phi = math.log(tan2 / tan1)
    deg = (math.degrees(math.atan2(dLon, phi)) + 360.0) % 360.0
    
    return deg

def pointInDistXY(
    loc:        "Starting location" = None, 
    direction:  "Direction towards the destination, Up is 0-degrees, clock-wise" = None, 
    dist:       "Distance from origin location to destination location" = None
    ) -> "A location in distance with given direction, in [lat, lon] form.":

    x = loc[0] + dist * math.sin(math.radians(direction))
    y = loc[1] + dist * math.cos(math.radians(direction))

    return (x, y)
def pointInDistLatLon(
    loc:        "Starting location" = None, 
    direction:  "Direction towards the destination, North is 0-degree, East is 90-degrees" = None, 
    distMeters: "Distance from origin location to destination location" = None
    ) -> "A location in distance with given direction, in [lat, lon] form.":
    newLoc = list(geopy.distance.distance(meters=distMeters).destination(point=loc, bearing=direction))
    return newLoc

def getSweepSeq(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    excludeIDs: "List of nodeIDs that will be excluded in sweeping, usually includes depot" = [],
    centerLoc:  "List, [x, y], the center point" = None,
    isClockwise: "True if the sweeping direction is clock-wise, False otherwise" = True,
    initDeg:    "Starting direction of the sweeping, 0 as North" = 0
    ) -> "Given a set of locations, and a center point, gets the sequence from sweeping":

    # Initialize heap =========================================================
    degHeap = []
    centerLocNodes = []
    
    # Build heap ==============================================================
    for nodeID in nodes:
        if (nodeID not in excludeIDs):
            dist = distEuclidean2D(nodes[nodeID]['loc'], centerLoc)
            # If the nodes are too close, separate it/them
            if (dist <= CONST_EPSILON):
                centerLocNodes.append(nodeID)
            else:
                dx = nodes[nodeID]['loc'][0] - centerLoc[0]
                dy = nodes[nodeID]['loc'][1] - centerLoc[1]
                (_, deg) = vec2ValDeg([dx, dy])
                # Calculate angles
                evalDeg = None
                if (isClockwise):
                    evalDeg = deg - initDeg
                else:
                    evalDeg = initDeg - deg
                while(evalDeg > 360):
                    evalDeg -= 360
                while(evalDeg < 0):
                    evalDeg += 360
                heapq.heappush(degHeap, (evalDeg, nodeID))

    # Sweep ===================================================================
    sweepSeq = []
    while (len(degHeap)):
        sweepSeq.append(heapq.heappop(degHeap)[1])
    sweepSeq.extend(centerLocNodes)

    return sweepSeq

def getNodesConvexHull(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    algo:       "1) String, 'Jarvis', O(nH) (current implementation is O(nHlog n)) or\
                 2) String, (not available) 'DNC', O(nlog n) or\
                 3) String, (not available) 'Graham', O(nlog n) or\
                 4) String, (not available) 'Melkman'" = "Jarvis"
    ) -> "Given a set of node locations, return a list of nodeID which construct the convex hull":

    # Initialize ==============================================================
    chSeq = None

    # Some extreme cases ======================================================
    if (len(nodes) == 0):
        return None
    elif (len(nodes) <= 3):
        chSeq = []
        for n in nodes:
            chSeq.append(n)
        return chSeq

    # Call subroutines ========================================================
    if (algo == "Jarvis"):
        chSeq = _getNodesConvexHullJavis(nodes)
    else:
        return None
    
    return chSeq

def _getNodesConvexHullJavis(nodes):
    # References ==============================================================
    # 1. https://blog.csdn.net/Bone_ACE/article/details/46239187
    # 2. chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/viewer.html?pdfurl=http%3A%2F%2Fwww.ams.sunysb.edu%2F~jsbm%2Fcourses%2F345%2F13%2Fmelkman.pdf&clen=46562&chunk=true
    # 3. https://en.wikipedia.org/wiki/Convex_hull_algorithms

    # Initialize ==============================================================
    chSeq = []

    # Get the location of the left-most nodeID ================================
    # Find an initial point which guaranteed to be in convex hull
    leftMostID = None
    leftMostX = None
    for n in nodes:
        if (leftMostID == None or nodes[n]['loc'][0] < leftMostX):
            leftMostID = n
            leftMostX = nodes[n]['loc'][0]

    # Jarvis march ============================================================
    curNodeID = leftMostID
    curDirection = 0
    marchFlag = True
    while (marchFlag):
        sweepSeq = getSweepSeq(
            nodes = nodes,
            excludeIDs = chSeq,
            centerLoc = nodes[curNodeID]['loc'],
            initDeg = curDirection)
        if (sweepSeq[0] == leftMostID):
            marchFlag = False
        chSeq.append(sweepSeq[0])
        curDirection = getHeadingXY(nodes[curNodeID]['loc'], nodes[sweepSeq[0]]['loc'])    
        curNodeID = sweepSeq[0]

    return chSeq