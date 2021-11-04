import math
import numpy as np
import heapq
import geopy.distance

from .const import *
from .msg import *
from .common import *
from .vector import *

def isPtSame(
    pt1:    "2-tuple of coordinate (x, y)",
    pt2:    "2-tuple of coordinate (x, y)",
    ) -> "Check if two points are the same":

    d = distEuclidean2D(pt1, pt2)
    if (d <= CONST_EPSILON):
        return True
    else:
        return False

def isPtOnSeg(
    pt:     "2-tuple of coordinate (x, y)", 
    seg:    "List of 2 pts as two ends"
    ) -> "Return true if the point is on the line segment, including initial pts, with error less than CONST_EPSILON":

    # Distance of the vector ==================================================
    if (isPtSame(seg[0], seg[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Get pts =================================================================
    [x1, y1] = [seg[0][0], seg[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [seg[1][0], seg[1][1]]

    # Calculate distance, needs to be in between ==============================
    S = calTriangleAreaByCoords(pt, seg[0], seg[1])
    d = distEuclidean2D(seg[0], seg[1])
    onSeg = (abs(2 * S / d) <= CONST_EPSILON 
        and x2 >= min(x1, x3) 
        and x2 <= max(x1, x3) 
        and y2 >= min(y1, y3) 
        and y2 <= max(y1, y3))
    
    return onSeg

def isPtInsideSeg(
    pt:     "2-tuple of coordinate (x, y)", 
    seg:    "List of 2 pts as two ends"
    ) -> "Return true if the point is inside the line segment, not including initial pts, with error less than CONST_EPSILON":

    # Distance of the vector ==================================================
    if (isPtSame(seg[0], seg[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Get pts =================================================================
    [x1, y1] = [seg[0][0], seg[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [seg[1][0], seg[1][1]]

    # Calculate distance, needs to be in between ==============================
    S = calTriangleAreaByCoords(pt, seg[0], seg[1])
    d = distEuclidean2D(seg[0], seg[1])
    insideSeg = (abs(2 * S / d) <= CONST_EPSILON 
        and x2 > min(x1, x3) 
        and x2 < max(x1, x3) 
        and y2 > min(y1, y3) 
        and y2 < max(y1, y3))
    
    return insideSeg

def isPtOnRay(
    pt:     "2-tuple of coordinate (x, y)",
    ray:    "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return true if the point is on the ray, with error less than CONST_EPSILON":

    # Distance of the vector ==================================================
    if (isPtSame(ray[0], ray[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Get pts =================================================================
    [x1, y1] = [ray[0][0], ray[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [ray[1][0], ray[1][1]]

    # Calculate distance, needs to be in between ==============================
    S = (x2 * y3 + x3 * y1 + x1 * y2) - (x2 * y1 + x3 * y2 + x1 * y3)
    d = distEuclidean2D(ray[0], ray[1])
    onRay = (abs(2 * S / d) < CONST_EPSILON
        and (
                (
                    x2 >= min(x1, x3) 
                    and x2 <= max(x1, x3) 
                    and y2 >= min(y1, y3) 
                    and y2 <= max(y1, y3)
                )
                or
                (
                    x3 >= min(x1, x2) 
                    and x3 <= max(x1, x2) 
                    and y3 >= min(y1, y2) 
                    and y3 <= max(y1, y2)
                )
            ))

    return onRay

def isPtOnPolyEdge(
    pt:     "2-tuple of coordinate (x, y)",
    poly:   "List of pts, form a close area"
    ) -> "Return true if the point is on the edge of polygon, with error less than CONST_EPSILON":

    # Check if the pt is on any of the edge segment ===========================
    onEdge = False
    for i in range(len(poly) - 1):
        currSeg = [poly[i], poly[i + 1]]
        onEdge = isPtOnSeg(pt, currSeg)
        if (onEdge):
            break
    if (not onEdge):
        onEdge = isPtOnSeg(pt, [poly[0], poly[-1]])

    return onEdge

def isPtInPoly(
    pt:     "2-tuple of coordinate (x, y)",
    poly:   "List of pts, form a close area"
    ) -> "Return true if the pt is in the interior of polygon":

    if (pt in poly):
        return False

    x = pt[1]
    y = pt[0]
    inside = False
    j = len(poly) - 1
    for i in range(len(poly)):
        xi = poly[i][1]
        yi = poly[i][0]
        xj = poly[j][1]
        yj = poly[j][0]
        intersect = (yi > y) != (yj > y)
        if (intersect):
            intersect = (x < (xj - xi) * (y - yi) / float(yj - yi) + xi)
        if (intersect):
            inside = not inside
        j = i

    return inside

def is3PtsClockWise(
    pt1: "2-tuple of coordinate (x, y)",
    pt2: "2-tuple of coordinate (x, y)",
    pt3: "2-tuple of coordinate (x, y)"
    ) -> "True if three given points are clockWise, False otherwise (could be collinear)":

    # Use Determinant to determine ============================================
    [x1, y1] = [pt1[0], pt1[1]]
    [x2, y2] = [pt2[0], pt2[1]]
    [x3, y3] = [pt3[0], pt3[1]]
    clockWise = ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) < 0)

    return clockWise

def calAngleBetweenRays(
    vec1:   "Ending pt, the origin pt is (0, 0)",
    vec2:   "Ending pt, the origin pt is (0, 0)"
    ) -> "Return angle between two rays, in radius":

    v1_u = vec1 / np.linalg.norm(vec1)
    v2_u = vec2 / np.linalg.norm(vec2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def isSegCrossSeg(
    seg1:   "List of 2 pts as two bounds ",
    seg2:   "List of 2 pts as two bounds "
    ) -> "Return true if segs are crossing in the interior of both segs, false otherwise":
    
    intPt = intSeg2Seg(seg1, seg2)
    if (intPt == None):
        return False

    if (isPtSame(intPt, seg1[0]) or isPtSame(intPt, seg1[1]) or isPtSame(intPt, seg2[0]) or isPtSame(intPt, seg2[1])):
        return False

    return True

def isSegIntSeg(
    seg1:   "List of 2 pts as two bounds ",
    seg2:   "List of 2 pts as two bounds "
    ) -> "Return true if segs are intersecting, false otherwise":

    [p, q] = seg1
    [u, w] = seg2

    loopPQU = is3PtsClockWise(p, q, u)
    loopPQW = is3PtsClockWise(p, q, w)
    loopUWP = is3PtsClockWise(u, w, p)
    loopUWQ = is3PtsClockWise(u, w, q)

    intersectFlag = False
    if (loopPQU != loopPQW and loopUWP != loopUWQ):
        intersectFlag = True

    if (not intersectFlag):
        if (loopPQU != loopPQW and loopUWP == loopUWQ):
            if (isPtOnPolyEdge(p, seg2) 
                or isPtOnPolyEdge(q, seg2)):
                intersectFlag = True

    if (not intersectFlag):
        if (loopPQU == loopPQW and loopUWP != loopUWQ):
            if (isPtOnPolyEdge(u, seg1) 
                or isPtOnPolyEdge(w, seg1)):
                intersectFlag = True

    if (not intersectFlag):
        if (loopPQU == loopPQW and loopUWP == loopUWQ):
            if (isPtOnPolyEdge(u, seg1) 
                or isPtOnPolyEdge(w, seg1) 
                or isPtOnPolyEdge(p, seg2) 
                or isPtOnPolyEdge(q, seg2)):
                intersectFlag = True
    return intersectFlag

def isSegCrossRay(
    seg:    "List of 2 pts as two bounds ",
    ray:    "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return true if ray is crossing interior of seg, false otherwise":

    crossFlag = (intSeg2Ray(seg, ray) != None)

    return crossFlag

def isRayCrossRay(
    ray1:   "List of 2 pts, the first pt defines the initial point, the second pt is on the ray",
    ray2:   "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return true if ray is crossing interior of ray, false otherwise":

    crossFlag = (intRay2Ray(ray1, ray2) != None)

    return crossFlag

def isSegIntRay(
    seg:    "List of 2 pts as two bounds ",
    ray:    "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return true if seg is intersecting with ray, false otherwise":

    [p, q] = seg
    [u, w] = ray

    loopPQU = is3PtsClockWise(p, q, u)
    loopPQW = is3PtsClockWise(p, q, w)
    loopUWP = is3PtsClockWise(u, w, p)
    loopUWQ = is3PtsClockWise(u, w, q)

    intersectFlag = False
    if (isSegIntSeg(seg, ray)):
        intersectFlag = True

    if (not intersectFlag):
        if (loopPQU == loopPQW
            and (isPtInPoly(w, [p, q, u])
                or isPtOnPolyEdge(w, [p, q, u]))):
            intersectFlag = True

    return intersectFlag

def isSegCrossPoly(
    seg:    "List of 2 pts as two bounds ",
    poly:   "List of pts, form a close area"
    ) -> "Return true if the line segment is crossing the interior of the poly, false otherwise":

    crossFlag = False
    if (isPtInPoly(seg[0], poly) or isPtInPoly(seg[1], poly)):
        crossFlag = True
    if (not crossFlag):
        for i in range(len(poly) - 1):
            edge = [poly[i], poly[i + 1]]
            if (isSegCrossSeg(edge, seg)):
                crossFlag = True
                break
    if (not crossFlag):
        edge = [poly[0], poly[-1]]
        if (isSegCrossSeg(edge, seg)):
            crossFlag = True

    return crossFlag

def isRayCrossPoly(
    ray:    "List of 2 pts, the first pt defines the initial point, the second pt is on the ray",
    poly:   "List of pts, form a close area"
    ) -> "Return true if the ray is crossing the interior of the poly, false otherwise":

    crossFlag = False
    if (isPtInPoly(ray[0], poly) or isPtInPoly(ray[1], poly)):
        crossFlag = True
    if (not crossFlag):
        for i in range(len(poly) - 1):
            edge = [poly[i], poly[i + 1]]
            if (isSegCrossRay(edge, ray)):
                crossFlag = True
                break
    if (not crossFlag):
        edge = [poly[0], poly[-1]]
        if (isSegCrossRay(edge, ray)):
            crossFlag = True

    return crossFlag

def intLine2Line(
    line1:  "List of 2 pts that defines a line",
    line2:  "List of 2 pts that defines a line"
    ) -> "Return 1) None if parallel or overlapped, or \
                 2) The intersect point if two lines intersect":

    intersectLoc = None

    # Get Ax + By + C = 0 =====================================================
    def abc(pt1, pt2):
        x1, y1 = pt1
        x2, y2 = pt2
        a = y1 - y2
        b = x2 - x1
        c = x1 * y2 - x2 * y1
        return a, b, c

    # Calculate intersection ==================================================
    a1, b1, c1 = abc(line1[0], line1[1])
    a2, b2, c2 = abc(line2[0], line2[1])
    D = a1 * b2 - a2 * b1

    # Check if parallel =======================================================
    if (D == 0):
        return None

    # Intersection ============================================================
    x = (b1 * c2 - b2 * c1) / D
    y = (a2 * c1 - a1 * c2) / D
    ptInt = (x, y)

    return ptInt

def intSeg2Seg(
    seg1:   "List of 2 pts as two bounds ",
    seg2:   "List of 2 pts as two bounds "
    ) -> "Return 1) None if no intersection, or \
                 2) The intersect point if two segments intersect":

    # Calculate intersection
    ptInt = intLine2Line(seg1, seg2)

    # Check if it is on segs
    if (ptInt != None and isPtOnSeg(ptInt, seg1) and isPtOnSeg(ptInt, seg2)):
        return ptInt
    else:
        return None

def intSeg2Ray(
    seg:    "List of 2 pts as two bounds ",
    ray:    "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return 1) None if no intersection, or \
                 2) The intersect point if two segments intersect":

    # Calculate intersection
    ptInt = intLine2Line(seg, ray)

    # Check if it is on segs
    if (ptInt != None and isPtOnSeg(ptInt, seg) and isPtOnRay(ptInt, ray)):
        return ptInt
    else:
        return None

def intRay2Ray(
    ray1:   "List of 2 pts, the first pt defines the initial point, the second pt is on the ray",
    ray2:   "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return 1) None if no intersection, or \
                 2) The intersect point if two segments intersect":

    # Calculate intersection
    ptInt = intLine2Line(ray1, ray2)

    # Check if it is on segs
    if (ptInt != None and isPtOnRay(ptInt, ray1) and isPtOnRay(ptInt, ray2)):
        return ptInt
    else:
        return None

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
    sortedNode = sortNodesByDist(
        nodes = nodes,
        edges = 'Euclidean',
        refNodeID = 0)

    # Calculate relative distances to pt ======================================
    dist = []
    for i in sortedNode:
        dist.append(distEuclidean2D(pt, nodes[sortedNode[i]]['loc']))
    timeStamp = []
    absVXY = math.sqrt(vecXY[0]**2 + vecXY[1]**2)
    for i in dist:
        timeStamp.append(i / absVXY)

    # Confirm time windows ====================================================
    for i in range(len(dist) - 1):
        mid = ((nodes[sortedNode[i]]['loc'][0] + nodes[sortedNode[i + 1]]['loc'][0]) / 2, 
            (nodes[sortedNode[i]]['loc'][1] + nodes[sortedNode[i + 1]]['loc'][1]) / 2)
        if (isPtInPoly(mid, poly)):
            tw.append([timeStamp[i], timeStamp[i + 1]])

    return tw

def twMovingPtInsidePolyLatLon(
    pt:         "Point that are moving" = None,
    poly:       "Polygon that are moving" = None,
    vecXYPt:    "Moving speed vector of the point in XY" = None,
    vecXYPoly:  "Moving speed vector of the polygon in XY" = None,
    ) -> "Given a moving point, a moving polygon, returns \
                1) the time window of when the point will be inside the polygon, or \
                2) None if not going to be inside the polygon":

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
    loc1:       "Coordinate of point", 
    loc2:       "Coordinate of point", 
    loc3:       "Coordinate of point"
    ) -> "Given the coordinates of three points, calculates the area":
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

def sortNodesByDist(
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
    sortedNode = []
    sortedNodeHeap = []
    for n in nodeIDs:
        dist = edges[refNodeID, n]
        heapq.heappush(sortedNodeHeap, (dist, n))
    while (len(sortedNodeHeap) > 0):
        sortedNode.append(heapq.heappop(sortedNodeHeap)[1])  

    return sortedNode

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
            heapq.heappush(degHeap, (evalDeg, n))

    # Sweep ===================================================================
    sweepSeq = []
    while (len(degHeap)):
        sweepSeq.append(heapq.heappop(degHeap)[1])
    sweepSeq.extend(centerLocNodes)

    return sweepSeq

def getScan(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    direction:  "Direction of scanning" = 0
    ) -> "Scan nodes from one direction":

    

    return seq

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

    # Subroutines for finding convex hull =====================================
    def getNodesConvexHullJavis():
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
        chSeq = getNodesConvexHullJavis()
    else:
        return None
    
    return chSeq
