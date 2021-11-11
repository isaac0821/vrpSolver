import math
import numpy as np
import heapq
import geopy.distance

from .const import *
from .msg import *
from .common import *
from .vector import *
from .relation import *

def ptXY2LatLon(
    ptXY:       "Point in (x, y) coordinates",
    ratio:      "Distance ratio between euclidean traveling matrix (`edgeInEuc`) and geological distance (as in meters)" = 1,
    refLatLon:  "(Lat, Lon) of (0, 0)" = None
    ) -> "Given a point in (x, y), try to map it with a (lat, lon) coordinate with necessary inputs":

    

    return ptLatLon

def ptLatLon2XY(

    ):

    return ptXY

def twMovingPtInsidePolyXY(
    pt:         "Point that are moving" = None,
    poly:       "Polygon that are moving" = None,
    vecXYPt:    "Moving speed vector of the point in XY" = None,
    vecXYPoly:  "Moving speed vector of the polygon in XY" = None,
    ) -> "Given a moving point, a moving polygon, returns \
                1) the time window of when the point will be inside the polygon, or \
                2) None if not going to be inside the polygon":

    # Initialize ==============================================================
    tw = []
    vecXY = calXYVecSubtract(vecXYPt, vecXYPoly)

    # Ray of point ============================================================
    ray = (pt, calXYVecAddition(pt, vecXY))

    # Check if the point is going to go into the polygon ======================
    intersectFlag = isRayCrossPoly(ray, poly)
    if (not intersectFlag):
        return []

    # Find the intersecting pts of the ray and all edges ======================
    # Polygon might not be convex, there could be multiple time windows
    intPts = []
    for i in range(len(poly) - 1):
        ptInt = intSeg2Ray([poly[i], poly[i + 1]], ray)
        if (ptInt != None):
            intPts.append(ptInt)
    ptInt = intSeg2Ray([poly[0], poly[-1]], ray)
    if (ptInt != None):
        intPts.append(ptInt)

    # Sort intersection points by dist ========================================
    # Convert intPts to nodes
    nodes = {}
    nodes[0] = {'loc': (pt[0], pt[1])}
    for i in range(len(intPts)):
        nodes[i + 1] = {'loc': (intPts[i][0], intPts[i][1])}
    sortedSeq = getSortedNodesByDist(
        nodes = nodes,
        edges = 'Euclidean',
        refNodeID = 0)

    # Calculate relative distances to pt ======================================
    dist = []
    for i in sortedSeq:
        dist.append(distEuclidean2D(pt, nodes[sortedSeq[i]]['loc']))
    timeStamp = []
    absVXY = math.sqrt(vecXY[0]**2 + vecXY[1]**2)
    for i in dist:
        timeStamp.append(i / absVXY)

    # Confirm time windows ====================================================
    for i in range(len(dist) - 1):
        mid = ((nodes[sortedSeq[i]]['loc'][0] + nodes[sortedSeq[i + 1]]['loc'][0]) / 2, 
            (nodes[sortedSeq[i]]['loc'][1] + nodes[sortedSeq[i + 1]]['loc'][1]) / 2)
        if (isPtInsidePoly(mid, poly)):
            tw.append([timeStamp[i], timeStamp[i + 1]])

    return tw

def twMovingPtInsidePolyLatLon(
    pt:         "Point that are moving" = None,
    poly:       "Polygon that are moving" = None,
    vecPolyPt:    "Moving speed vector of the point in XY" = None,
    vecPolyPoly:  "Moving speed vector of the polygon in XY" = None,
    ) -> "Given a moving point, a moving polygon, returns \
                1) the time window of when the point will be inside the polygon, or \
                2) None if not going to be inside the polygon":

    # Initialize ==============================================================
    tw = []
    vecPoly = calPolyVecSubtract(vecPolyPoly, vecPolyPoly)

    # Need to convert the speed vector into a proper coordinate system ========
    

    return tw

def getTauEuclidean(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    speed:      "Ratio: Dist / Time" = 1
    ) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Get tau =================================================================
    tau = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                t = distEuclidean2D(nodes[i]['loc'], nodes[j]['loc']) / speed
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = CONST_EPSILON
    return tau

def getTauLatLon(
    nodes:      "Dictionary, returns the coordinate of given nodeIDs, \
                    {\
                        nodeIDs1: {'loc': (x, y)}, \
                        nodeIDs2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',                    
    distUnit:   "Unit of distance\
                 1) String (default) 'mile'\
                 2) String 'meter'\
                 3) String 'kilometer'" = 'mile',
    speed:      "Ratio: Dist / Time" = 1
    ) -> "Dictionary, {(nodeID1, nodeID2): dist, ...}":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Get tau =================================================================
    tau = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                t = distLatLon(nodes[i]['loc'], nodes[j]['loc'], distUnit) / speed
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = CONST_EPSILON
    return tau

def getPerpendicularLine(
    pt:         "Point which the line will go through",
    vec:        "The vector that perpendicular to the line"
    ) -> "Given a point, a vector, returns a line that is going through this point and perpendicular to the vector":

    # Get the direction =======================================================
    heading = getHeadingXY(pt, (pt[0] + vec[0], pt[1] + vec[1]))
    newHeading1 = heading + 90
    newHeading2 = heading - 90
    while (newHeading1 > 360):
        newHeading1 -= 360
    while (newHeading2 > 360):
        newHeading2 -= 360
    while (newHeading1 < 0):
        newHeading1 += 360
    while (newHeading2 < 0):
        newHeading2 += 360

    # Get line ================================================================
    line = [pointInDistXY(pt, newHeading1, 10), pointInDistXY(pt, newHeading2, 10)]

    return line

def distEuclidean2D(
    coord1:     "First coordinate, in (x, y)", 
    coord2:     "Second coordinate, in (x, y)"
    ) -> "Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number":
    if (coord1 != None and coord2 != None):
        return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2)
    else:
        return 0

def distLatLon(
    coord1:     "First coordinate, in (lat, lon)", 
    coord2:     "Second coordinate, in (lat, lon)",
    distUnit:   "Unit of distance\
                 1) String (default) 'mile'\
                 2) String 'meter'\
                 3) String 'kilometer'" = 'mile'
    ) -> "Gives a Euclidean distance based on two lat/lon coords, if two coordinates are the same, return a small number":
    R = None
    if (distUnit == 'mile'):
        R = CONST_EARTH_RADIUS_MILES
    elif (distUnit == 'meter'):
        R = CONST_EARTH_RADIUS_METERS
    elif (distUnit == 'kilometer'):
        R = CONST_EARTH_RADIUS_METERS / 1000
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

def distPt2Line(
    pt:         "2-tuple of coordinate (x, y)",
    line:       "List of 2 pts which defines a line"
    ) -> "Distance from a given point to a given line":

    # Validation ==============================================================
    if (is2PtsSame(line[0], line[1])):
        print(ERROR_ZERO_VECTOR)
        return None
    
    # Initialize ==============================================================
    ptA = line[0]
    ptB = line[1]
    ptS = pt

    # Check if the pt is on line, if so the distance is 0 =====================
    if (isPtOnLine(pt, line)):
        return 0.0

    # Calculate dist to line ==================================================
    areaSAB = calTriangleAreaByCoords(ptS, ptA, ptB)
    bottom = distEuclidean2D(ptA, ptB)
    height = 2 * areaSAB / bottom
    dist = height

    return dist

def distPt2Seg(
    pt:         "2-tuple of coordinate (x, y)",
    seg:        "List of 2 pts which defines a line segment"
    ) -> "Distance from a given point to a given line segment":
    
    # Validation ==============================================================
    if (is2PtsSame(line[0], line[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Initialize ==============================================================
    ptA = line[0]
    ptB = line[1]
    ptS = pt

    # Check if the pt is on line, if so the distance is 0 =====================
    if (isPtOnLine(pt, line)):
        return 0.0

    # Rays ====================================================================
    rayAS = [ptS[0] - ptA[0], ptS[1] - ptA[1]]
    rayAB = [ptB[0] - ptA[0], ptB[1] - ptA[1]]
    rayBS = [ptS[0] - ptB[0], ptS[1] - ptB[1]]
    rayBA = [ptA[0] - ptB[0], ptA[1] - ptB[1]]

    # cos value for A angle and B angle =======================================
    cosSAB = cosRay2Ray(rayAS, rayAB)
    cosSBA = cosRay2Ray(rayBS, rayBA)

    # Calculate dist to line ==================================================
    # if both angles are sharp, the closest point will be in the line, otherwise the closest point is at the edge
    if (cosSAB >= 0 and cosSBA >= 0):
        areaSAB = calTriangleAreaByCoords(ptS, ptA, ptB)
        bottom = distEuclidean2D(ptA, ptB)
        height = 2 * areaSAB / bottom
        dist = height
    else:
        distAS = distEuclidean2D(ptS, ptA)
        distBS = distEuclidean2D(ptS, ptB)
        dist = min(distAS, distBS)

    return dist

def calTriangleAreaByEdges(
    a:          "Length of edge", 
    b:          "Length of edge", 
    c:          "Length of edge"
    ) -> "Given the length of three edges, calculates the area":
    # Using Heron's Formula
    s = (a + b + c) / 2
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return area

def calTriangleAreaByCoords(
    pt1:       "Coordinate of point", 
    pt2:       "Coordinate of point", 
    pt3:       "Coordinate of point"
    ) -> "Given the coordinates of three points, calculates the area":
    [x1, y1] = pt1
    [x2, y2] = pt2
    [x3, y3] = pt3
    val = (x2 * y3 + x3 * y1 + x1 * y2) - (x2 * y1 + x3 * y2 + x1 * y3)
    area = abs(val)
    return area

def getHeadingXY(
    coord1:     "Current location", 
    coord2:     "Targeted location"
    ) -> "Given current location and a goal location, calculate the heading. Up is 0-degrees, clock-wise":
    
    vec = (coord2[0] - coord1[0], coord2[1] - coord1[1])
    (_, vDeg) = vecXY2Polar(vec)

    return vDeg

def getHeadingLatLon(
    coord1:     "Current location as in [lat, lon]", 
    coord2:     "Target location as in [lat, lon]"
    ) -> "Given current location and a goal location, calculate the heading. Up is 0-degrees, clock-wise":
    # Ref: https://github.com/manuelbieh/Geolib/issues/28
    
    [lat1, lon1] = coord1
    [lat2, lon2] = coord2
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
    pt:         "Starting location" = None, 
    direction:  "Direction towards the destination, Up is 0-degrees, clock-wise" = None, 
    dist:       "Distance from origin location to destination location" = None
    ) -> "A location in distance with given direction, in [lat, lon] form.":
    x = pt[0] + dist * math.sin(math.radians(direction))
    y = pt[1] + dist * math.cos(math.radians(direction))
    return (x, y)

def pointInDistLatLon(
    pt:         "Starting location" = None, 
    direction:  "Direction towards the destination, North is 0-degree, East is 90-degrees" = None, 
    distMeters: "Distance from origin location to destination location" = None
    ) -> "A location in distance with given direction, in [lat, lon] form.":
    newLoc = list(geopy.distance.distance(meters=distMeters).destination(point=pt, bearing=direction))
    return newLoc

def getSortedNodesByDist(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    refNodeID:  "List, [x, y], the reference location to calculate distance" = None    
    ) -> "Given a set of locations, and a referencing node, sort the nodes by distance to this referencing node":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'LatLon'):
            edges = getTauLatLon(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Sort distance ===========================================================
    sortedSeq = []
    sortedSeqHeap = []
    for n in nodeIDs:
        dist = edges[refNodeID, n]
        heapq.heappush(sortedSeqHeap, (dist, n))
    while (len(sortedSeqHeap) > 0):
        sortedSeq.append(heapq.heappop(sortedSeqHeap)[1])  

    return sortedSeq

def getSweepSeq(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    centerLoc:  "List, [x, y], the center point" = None,
    isClockwise: "True if the sweeping direction is clock-wise, False otherwise" = True,
    initDeg:    "Starting direction of the sweeping, 0 as North" = 0
    ) -> "Given a set of locations, and a center point, gets the sequence from sweeping":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Initialize heap =========================================================
    degHeap = []
    centerLocNodes = []
    
    # Build heap ==============================================================
    for n in nodeIDs:
        dist = distEuclidean2D(nodes[n]['loc'], centerLoc)
        # If the nodes are too close, separate it/them
        if (dist <= CONST_EPSILON):
            centerLocNodes.append(n)
        else:
            dx = nodes[n]['loc'][0] - centerLoc[0]
            dy = nodes[n]['loc'][1] - centerLoc[1]
            (_, deg) = vecXY2Polar((dx, dy))
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
            heapq.heappush(degHeap, (evalDeg, n))

    # Sweep ===================================================================
    sweepSeq = []
    while (len(degHeap)):
        sweepSeq.append(heapq.heappop(degHeap)[1])
    sweepSeq.extend(centerLocNodes)

    return sweepSeq

# [Testing]
def getScan(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    direction:  "Direction of scanning" = 0
    ) -> "Scan nodes from one direction":

    # Initialize ==============================================================
    baseline = []
    maxDist = 0
    firstNodeID = list(nodes.keys())[0]
    for n in nodes:
        d = distEuclidean2D(nodes[n]['loc'], nodes[firstNodeID]['loc'])
        if (maxDist == None or d > maxDist):
            maxDist = d
    basePt = pointInDistXY(nodes[firstNodeID]['loc'], direction, d)
    baseline = getPerpendicularLine(basePt, vecPolar2XY([10, direction]))

    # Distance to the baseline ================================================
    distHeap = []
    for n in nodes:
        dist2Baseline = distPt2Line(nodes[n]['loc'], baseline)
        heapq.heappush(distHeap, (dist2Baseline, n))

    # Scan ====================================================================
    scanSeq = []
    while (len(distHeap)):
        scanSeq.append(heapq.heappop(distHeap)[1])

    return scanSeq

def getConvexHull(
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

    # Subroutines for finding convex hull =====================================
    def _getConvexHullJavis():
        # References ----------------------------------------------------------
        # 1. https://blog.csdn.net/Bone_ACE/article/details/46239187
        # 2. chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/viewer.html?pdfurl=http%3A%2F%2Fwww.ams.sunysb.edu%2F~jsbm%2Fcourses%2F345%2F13%2Fmelkman.pdf&clen=46562&chunk=true
        # 3. https://en.wikipedia.org/wiki/Convex_hull_algorithms

        # Initialize ----------------------------------------------------------
        chSeq = []

        # Get the location of the left-most nodeID ----------------------------
        # Find an initial point which guaranteed to be in convex hull
        leftMostID = None
        leftMostX = None
        for n in nodes:
            if (leftMostID == None or nodes[n]['loc'][0] < leftMostX):
                leftMostID = n
                leftMostX = nodes[n]['loc'][0]

        # Jarvis march --------------------------------------------------------
        curNodeID = leftMostID
        curDirection = 0
        marchFlag = True
        while (marchFlag):
            sweepSeq = getSweepSeq(
                nodes = nodes,
                nodeIDs = [i for i in nodes if i not in chSeq],
                centerLoc = nodes[curNodeID]['loc'],
                initDeg = curDirection)
            if (sweepSeq[0] == leftMostID):
                marchFlag = False
            chSeq.append(sweepSeq[0])
            curDirection = getHeadingXY(nodes[curNodeID]['loc'], nodes[sweepSeq[0]]['loc'])    
            curNodeID = sweepSeq[0]
        return chSeq

    # Call subroutines ========================================================
    if (algo == "Jarvis"):
        chSeq = _getConvexHullJavis()
    else:
        return None
    
    return chSeq
