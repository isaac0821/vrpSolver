import geopy.distance
import heapq
import math
import numpy as np
from pyproj import Geod
from shapely.geometry import Polygon
import tripy

from .common import *
from .const import *
from .graph import *
from .msg import *
from .relation import *
from .vector import *

def ptXY2LatLonMercator(
    ptXY:       "Point in (x, y) coordinates"
    ) -> "Given a point in (x, y), try to map it with a (lat, lon) coordinate with necessary inputs":
    # Ref: https://wiki.openstreetmap.org/wiki/Mercator#Python
    (x, y) = ptXY
    lon = math.degrees(y / CONST_EARTH_RADIUS_METERS)
    lat = math.degrees(2 * math.atan(math.exp(x / CONST_EARTH_RADIUS_METERS)) - math.pi / 2.0)
    ptLatLon = (lat, lon)
    return ptLatLon

def ptLatLon2XYMercator(
    ptLatLon:   "Point in (lat, lon) coordinates" = None
    ) -> "Given a point int (lat, lon), try to map it with a (x, y) coordinates with necessary inputs":
    # Ref: https://wiki.openstreetmap.org/wiki/Mercator#Python
    (lat, lon) = ptLatLon
    y = math.radians(lon) * CONST_EARTH_RADIUS_METERS
    x = math.log(math.tan(math.pi / 4 + math.radians(lat) / 2)) * CONST_EARTH_RADIUS_METERS
    ptXY = (x, y)
    return ptXY

def mapLatLon2LatLon(
    locsAreaA:    "A list of (lat, lon) in area A" = None,
    locAnchorA:   "The anchor point of area A" = None,
    locAnchorB:   "The anchor point of area B" = None
    ) -> "Given a set of lat/lon coords in area A, a relative location A of area A, and a relative \
          location B of area B. Map the lat/lon coords to area B":

    # Initialize ==============================================================
    locsAreaB = []

    # Mapping, using directions and distance ==================================
    for locA in locsAreaA:
        deg = getHeadingLatLon(locAnchorA, locA)
        dist = distLatLon(locAnchorA, locA)
        locB = pointInDistLatLon(locAnchorB, deg, dist)
        locsAreaB.append(locB)

    return locsAreaB

def htMovingPtTowardsLineSegXY(
    ptXY:       "Point that is moving" = None,
    segXY:      "Line segment that is moving - no rotation" = None,
    vecXYPt:    "Moving speed vector of the point in XY" = None,
    vecXYSeg:   "Moving speed vector of the line segment in XY" = None
    ) -> "Given a moving point and a moving line segment, returns\
        1) when the point is going to hit the line segment, or\
        2) None if the point is not going to hit the line segment":

    # Initialize ==============================================================
    # Convert the speed of both objects into the speed of pt
    vecXY = calXYVecSubtract(vecXYPt, vecXYSeg)

    # Move axis, so that ptXY = (0, 0) ========================================
    segOffSet = [[segXY[0][0] - ptXY[0], segXY[0][1] - ptXY[1]], [segXY[1][0] - ptXY[0], segXY[1][1] - ptXY[1]]]
    rayOffset = [(0, 0), vecXY]

    # Find intersecting point of the ray and the line seg =====================
    ptInt = intSeg2Ray(segOffSet, rayOffset)
    if (ptInt == None):
        return None

    # Calculate the distance from ptXY to ptInt and calculate the time needed =
    dist = distEuclidean2D((0, 0), ptInt)
    absVXY = math.sqrt(vecXY[0]**2 + vecXY[1]**2)
    hitTime = dist / absVXY

    return hitTime

def htMovingPtTowardsLineSegLatLon(
    ptLatLon:   "Point that is moving" = None,
    segLatLon:  "Line segment that is moving - no rotation" = None,
    vecPolarPt: "Moving speed vector of the point in XY" = None,
    vecPolarSeg: "Moving speed vector of the line segment in XY" = None
    ) -> "Lat/lon version of `hitTimeMovingPtTowardsLineSegXY()`":

    # Initialize ==============================================================
    ptXY = ptLatLon2XYMercator(ptLatLon)
    segXY = [ptLatLon2XYMercator(segLatLon[0]), ptLatLon2XYMercator(segLatLon[1])]
    segOffset = [[segXY[0][0] - ptXY[0], segXY[0][1] - ptXY[1]], [segXY[1][0] - ptXY[0], segXY[1][1] - ptXY[1]]]
    vecXY = calXYVecSubtract(
        vecPolar2XY([vecPolarPt[0], 90 - vecPolarPt[1]]), 
        vecPolar2XY([vecPolarSeg[0], 90 - vecPolarSeg[1]]))
    rayOffset = [(0, 0), vecXY]

    # Find intersection in XY space ===========================================
    ptIntXY = intSeg2Ray(segOffset, rayOffset)
    if (ptIntXY == None):
        return None
    ptIntLatLon = ptXY2LatLonMercator([ptIntXY[0] + ptXY[0], ptIntXY[1] + ptXY[1]])

    # Calculate dist and time =================================================
    dist = distLatLon(ptLatLon, ptIntLatLon)
    absVXY = vecXY2Polar(vecXY)[0]
    hitTime = dist / absVXY

    return hitTime

def htMovingPtTowardsLineLatLon(
    ptLatLon:   "Point that is moving" = None,
    lineLatLon:  "Line segment that is moving - no rotation" = None,
    vecPolarPt: "Moving speed vector of the point in XY" = None,
    vecPolarSeg: "Moving speed vector of the line segment in XY" = None
    ) -> "Lat/lon version of `hitTimeMovingPtTowardsLineSegXY()`":

    # Initialize ==============================================================
    ptXY = ptLatLon2XYMercator(ptLatLon)
    segXY = [ptLatLon2XYMercator(lineLatLon[0]), ptLatLon2XYMercator(lineLatLon[1])]
    segOffset = [[segXY[0][0] - ptXY[0], segXY[0][1] - ptXY[1]], [segXY[1][0] - ptXY[0], segXY[1][1] - ptXY[1]]]
    vecXY = calXYVecSubtract(
        vecPolar2XY([vecPolarPt[0], 90 - vecPolarPt[1]]), 
        vecPolar2XY([vecPolarSeg[0], 90 - vecPolarSeg[1]]))
    rayOffset = [(0, 0), vecXY]

    # Find intersection in XY space ===========================================
    ptIntXY = intRay2Line(segOffset, rayOffset)
    if (ptIntXY == None):
        return None
    ptIntLatLon = ptXY2LatLonMercator([ptIntXY[0] + ptXY[0], ptIntXY[1] + ptXY[1]])

    # Calculate dist and time =================================================
    dist = distLatLon(ptLatLon, ptIntLatLon)
    absVXY = vecXY2Polar(vecXY)[0]
    hitTime = dist / absVXY

    return hitTime

def twMovingPtInsidePolyXY(
    ptXY:       "Point that is moving" = None,
    polyXY:     "Polygon that are moving" = None,
    vecXYPt:    "Moving speed vector of the point in XY" = None,
    vecXYPoly:  "Moving speed vector of the polygon in XY" = None,
    ) -> "Given a moving point in (x, y), a moving polygon in Euclidean2D, returns \
        1) the time window of when the point is inside the polygon, or \
        2) empty list if the point is not going to be inside the polygon":

    # Initialize ==============================================================
    tw = []
    # Convert the speed of both objects into the speed of pt
    vecXY = calXYVecSubtract(vecXYPt, vecXYPoly)

    # Move axis, so that ptXY = (0, 0) ========================================
    polyOffSet = []
    for p in polyXY:
        polyOffSet.append((p[0] - ptXY[0], p[1] - ptXY[1]))
    rayOffset = [(0, 0), vecXY]

    # Check if the point is going to go into the polygon ======================
    intersectFlag = isRayCrossPoly(rayOffset, polyOffSet)
    if (not intersectFlag):
        return []

    # Find the intersecting pts of the ray and all edges ======================
    # Polygon might not be convex, there could be multiple time windows
    intPts = []
    for i in range(-1, len(polyOffSet) - 1):
        ptInt = intSeg2Ray([polyOffSet[i], polyOffSet[i + 1]], rayOffset)
        if (ptInt != None):
            intPts.append(ptInt)

    # Sort intersection points by dist ========================================
    # Convert intPts to nodes
    nodes = {}
    nodes[0] = {'loc': (0, 0), 'color': 'blue'}
    for i in range(len(intPts)):
        nodes[i + 1] = {'loc': (intPts[i][0], intPts[i][1]), 'color': 'red'}
    sortedSeq = getSortedNodesByDist(
        nodes = nodes,
        edges = 'Euclidean',
        refNodeID = 0)

    # Calculate relative distances to ptXY ====================================
    dist = []
    for i in range(len(sortedSeq)):
        dist.append(distEuclidean2D((0, 0), nodes[sortedSeq[i]]['loc']))
    timeStamp = []
    absVXY = math.sqrt(vecXY[0]**2 + vecXY[1]**2)
    for i in dist:
        timeStamp.append(i / absVXY)

    # Confirm time windows ====================================================
    for i in range(len(dist) - 1):
        mid = ((nodes[sortedSeq[i]]['loc'][0] / 2 + nodes[sortedSeq[i + 1]]['loc'][0] / 2), 
            (nodes[sortedSeq[i]]['loc'][1] / 2 + nodes[sortedSeq[i + 1]]['loc'][1] / 2))
        if (isPtInsidePoly(mid, polyOffSet)):
            tw.append([timeStamp[i], timeStamp[i + 1]])

    return tw

def twMovingPtInsidePolyLatLon(
    ptLatLon:      "Point that are moving" = None,
    polyLatLon:    "Polygon that are moving" = None,
    vecPolarPt:    "Moving speed vector of the point in XY, in the format of (val, deg), 0 deg as north" = None,
    vecPolarPoly:  "Moving speed vector of the polygon in XY, in the format of (val, deg), 0 deg as north" = None,
    ) -> "Given a moving point in (lat, lon), a moving polygon, returns \
        1) the time window of when the point will be inside the polygon, or \
        2) None if not going to be inside the polygon":

    # Steps ===================================================================
    # 1. Project the Lat/Lon into XY space, using Mercator Projection
    # 2. Calculate the `twMovingPtInsidePolyXY()`, find the intersections
    # 3. Convert the intersections back to (lat, lon)
    # 4. Calculated the time windows based on speed vector in LatLon

    # Initialize ==============================================================
    tw = []
    ptXY = ptLatLon2XYMercator(ptLatLon)
    polyXYOffset = []
    for p in polyLatLon:
        pXY = ptLatLon2XYMercator(p)
        polyXYOffset.append((pXY[0] - ptXY[0], pXY[1] - ptXY[1]))
    vecXY = calXYVecSubtract(
        vecPolar2XY([vecPolarPt[0], 90 - vecPolarPt[1]]), 
        vecPolar2XY([vecPolarPoly[0], 90 - vecPolarPoly[1]]))
    rayXYOffset = [(0, 0), vecXY]

    # Find intersection in XY space ===========================================
    intersectFlag = isRayCrossPoly(rayXYOffset, polyXYOffset)
    # print(intersectFlag)
    if (not intersectFlag):
        return []
    intPtsXY = []
    for i in range(-1, len(polyXYOffset) - 1):
        ptIntXY = intSeg2Ray([polyXYOffset[i], polyXYOffset[i + 1]], rayXYOffset)
        if (ptIntXY != None):
            intPtsXY.append(ptIntXY)

    # Convert those points back to (lat, lon) =================================
    intPtsLatLon = []
    for p in intPtsXY:
        pXY = (p[0] + ptXY[0], p[1] + ptXY[1])
        intPtsLatLon.append(ptXY2LatLonMercator(pXY))

    # Sort intersection points by dist ========================================
    # Convert intPts to nodes
    intPtsLatLonDict = {}
    intPtsLatLonDict[0] = {'loc': (ptLatLon[0], ptLatLon[1])}
    for i in range(len(intPtsLatLon)):
        intPtsLatLonDict[i + 1] = {'loc': (intPtsLatLon[i][0], intPtsLatLon[i][1])}
    sortedSeq = getSortedNodesByDist(
        nodes = intPtsLatLonDict,
        edges = 'LatLon',
        refNodeID = 0)

    # Calculate relative distances to ptLatLon ======================================
    dist = []
    for i in range(len(sortedSeq)):
        dist.append(distLatLon(ptLatLon, intPtsLatLonDict[sortedSeq[i]]['loc'], distUnit='meter'))
    timeStamp = []
    vecPolar = vecXY2Polar(vecXY)
    absVXY = vecPolar[0]
    for i in dist:
        timeStamp.append(i / absVXY)

    # Confirm time windows ====================================================
    for i in range(len(dist) - 1):
        mid = ((intPtsLatLonDict[sortedSeq[i]]['loc'][0] / 2 + intPtsLatLonDict[sortedSeq[i + 1]]['loc'][0] / 2), 
            (intPtsLatLonDict[sortedSeq[i]]['loc'][1] / 2 + intPtsLatLonDict[sortedSeq[i + 1]]['loc'][1] / 2))
        if (isPtInsidePoly(mid, polyLatLon)):
            tw.append([timeStamp[i], timeStamp[i + 1]])

    return tw

def projSeg2LineXY(
    seg:        "Line segment to be projected to line" = None,
    line:       "Line to project to" = None,
    vec:        "Vector of projecting" = None
    ) -> "Given a line segment, project it to a line using a given vector":

    line4SegEnd1 = [seg[0], [seg[0][0] + vec[0], seg[0][1] + vec[1]]]
    line4SegEnd2 = [seg[1], [seg[1][0] + vec[0], seg[1][1] + vec[1]]]

    projPt1 = intLine2Line(line4SegEnd1, line)
    projPt2 = intLine2Line(line4SegEnd2, line)

    return {
        'shadowSegOnLine': [projPt1, projPt2]
    }

def distEuclidean2D(
    pt1:        "First coordinate, in (x, y)", 
    pt2:        "Second coordinate, in (x, y)"
    ) -> "Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number":

    # 2-Norm ==================================================================
    if (pt1 != None and pt2 != None):
        return math.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)
    else:
        return 0

def distLatLon(
    pt1:        "First coordinate, in (lat, lon)", 
    pt2:        "Second coordinate, in (lat, lon)",
    distUnit:   "Unit of distance\
                 1) String 'mile'\
                 2) String (default) 'meter'\
                 3) String 'kilometer'" = 'meter'
    ) -> "Gives a Euclidean distance based on two lat/lon coords, if two coordinates are the same, return a small number":
    
    # Get radius as in distUnit ===============================================
    R = None
    if (distUnit in ['mile', 'mi']):
        R = CONST_EARTH_RADIUS_MILES
    elif (distUnit in ['meter', 'm']):
        R = CONST_EARTH_RADIUS_METERS
    elif (distUnit in ['kilometer', 'km']):
        R = CONST_EARTH_RADIUS_METERS / 1000
    else:
        print("ERROR: Unrecognized distance unit, options are 'mile', 'meter', 'kilometer'")
        return

    # Calculate distance ======================================================
    if (pt1 != None and pt2 != None):
        (lat1, lon1) = pt1
        (lat2, lon2) = pt2
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
    # Using Heron's Formula ===================================================
    s = (a / 2 + b / 2 + c / 2)
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return area

def calTriangleAreaByCoords(
    pt1:       "Coordinate of point", 
    pt2:       "Coordinate of point", 
    pt3:       "Coordinate of point"
    ) -> "Given the coordinates of three points, calculates the area":
    # Using determinant =======================================================
    [x1, y1] = pt1
    [x2, y2] = pt2
    [x3, y3] = pt3
    val = (x2 * y3 + x3 * y1 + x1 * y2) - (x2 * y1 + x3 * y2 + x1 * y3)
    area = abs(val)
    return area

def calPolygonAreaLatLon(
    polyLatLon: "List of [lat, lon] coordinates"
    ) -> "Returns the area surrounded by polyLatLon on the Earth":

    # Create polygon ==========================================================
    # NOTE: shapely is in [lon, lat] format
    rev = []
    for p in polyLatLon:
        rev.append((p[1], p[0]))
    polygon = Polygon(rev)

    # Using pyproj to calculate ===============================================
    # Ref: https://hypc.github.io/2020/03/16/python-geo-area/
    geod = Geod(ellps = "WGS84")
    area = abs(geod.geometry_area_perimeter(polygon)[0])

    return area

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
                tau[j, i] = CONST_EPSILON
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
                tau[j, i] = CONST_EPSILON
    return tau

def getTauGrid(
    nodes:      "Dictionary, returns the coordinate of given nodeIDs, \
                    {\
                        nodeIDs1: {'loc': (x, y)}, \
                        nodeIDs2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All', 
    colRow:     "Number of columns, and number of rows" = (None, None),
    barriers:   "List of blocking grids" = []
    ) -> "Given a grid that has barrier, returns tau matrix":

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
                t = gridPathFinding(
                    gridColRow = colRow,
                    barriers = barriers,
                    startGridCoord = nodes[i]['loc'],
                    endGridCoord = nodes[j]['loc'])['dist']
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = 0
                tau[j, i] = 0
    return tau

def getTauNetwork(
    nodes:      "Dictionary, returns the coordinate of given nodeIDs, \
                    {\
                        nodeIDs1: {'loc': (x, y)}, \
                        nodeIDs2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All', 
    network:    "List of dictionary of networks in the format of \
                    {\
                        'start': (lat, lon),\
                        'end': (lat, lon),\
                    }" = None,
    ) -> "Given a network, e.g., road network, and a set of locations on network (could be on arc), returns distance matrix":

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

def getHeadingXY(
    pt1:        "Current location", 
    pt2:        "Targeted location"
    ) -> "Given current location and a goal location, calculate the heading. North is 0-degrees, clock-wise":
    
    vec = (pt2[0] - pt1[0], pt2[1] - pt1[1])
    (_, vDeg) = vecXY2Polar(vec)

    return vDeg

def getHeadingLatLon(
    pt1:        "Current location as in [lat, lon]", 
    pt2:        "Target location as in [lat, lon]"
    ) -> "Given current location and a goal location, calculate the heading. North is 0-degrees, clock-wise":
    # Ref: https://github.com/manuelbieh/Geolib/issues/28
    (lat1, lon1) = pt1
    (lat2, lon2) = pt2
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

def getMileageInPathLatLon(
    path:       "A list of coordinates in [lat, lon]" = None,
    distMeters: "Distance from starting point of the path, in [m]" = None
    ) -> "Given a list of lat/lon coordinates, and a traveling mileage, returns the coordinate":

    # Initialize ==============================================================
    inPathFlag = False
    accDist = 0
    preLoc = []
    nextLoc = []

    # Find segment ============================================================
    for i in range(0, len(path) - 1):
        accDist += distLatLon(path[i], path[i + 1])
        if (accDist > distMeters):
            preLoc = path[i]
            nextLoc = path[i + 1]
            inPathFlag = True
            break

    if (inPathFlag == False):
        return None

    # Find location on the segment ============================================
    remainDist = accDist - distMeters
    segDist = distLatLon(preLoc, nextLoc)
    lat = nextLoc[0] + (remainDist / segDist) * (preLoc[0] - nextLoc[0])
    lon = nextLoc[1] + (remainDist / segDist) * (preLoc[1] - nextLoc[1])
    locInMileage = [lat, lon]

    return locInMileage

def getRndPtRoadNetworkPoly(
    N:          "Number of nodes" = 1,
    roadNetwork: "Dictionary of road network in the format of \
                {\
                    roadID: {\
                        'line': [[lat, lon], [lat, lon], ...]\
                    }\
                }" = None,
    poly:       "Nodes should also within this polygon" = None,
    ) -> "Given a road network, generate customers that locates on the road network":
    
    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for road in roadNetwork:
        roadLength = 0
        includedFlag = False
        if (poly == None):
            includedFlag = True
        else:
            for i in range(len(roadNetwork[road]['line'])):
                if (isPtOnPoly(roadNetwork[road]['line'][i], poly)):
                    includedFlag = True
                    break

        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(roadNetwork[road]['line']) - 1):
                roadLength += distLatLon(roadNetwork[road]['line'][i], roadNetwork[road]['line'][i + 1])
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(road)

    # Check if there are roads included =======================================
    if (sum(lengths) == 0):
        return None

    # Use accept-denial to test if the node is within poly ====================
    # FIXME: Inefficient approach, will need to be rewritten
    nodeLocs = []
    for i in range(N):
        lat = None
        lon = None
        if (poly == None):
            idx = rndPick(lengths)
            edgeLength = lengths[idx]
            edgeDist = random.uniform(0, 1) * edgeLength
            (lat, lon) = getMileageInPathLatLon(roadNetwork[roadIDs[idx]]['line'], edgeDist)
        else:
            insideFlag = False
            while (not insideFlag):
                idx = rndPick(lengths)
                edgeLength = lengths[idx]
                edgeDist = random.uniform(0, 1) * edgeLength
                (lat, lon) = getMileageInPathLatLon(roadNetwork[roadIDs[idx]]['line'], edgeDist)
                if (isPtOnPoly([lat, lon], poly)):
                    insideFlag = True
        nodeLocs.append((lat, lon))
    return nodeLocs

def getRndPtRoadNetworkCircle(
    N:          "Number of nodes" = 1,
    roadNetwork: "Dictionary of road network in the format of \
                {\
                    roadID: {\
                        'line': [[lat, lon], [lat, lon], ...]\
                    }\
                }" = None,
    circle:     "Dictionary, {'centerLoc': centerLoc; 'radius': radius in [m]}" = None,
    ) -> "Given a road network, generate customers that locates on the road network":
    
    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for road in roadNetwork:
        roadLength = 0
        includedFlag = False
        if (circle == None):
            includedFlag = True
        else:
            for i in range(len(roadNetwork[road]['line'])):
                if (distLatLon(roadNetwork[road]['line'][i], circle['centerLoc']) <= circle['radius']):
                    includedFlag = True
                    break

        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(roadNetwork[road]['line']) - 1):
                roadLength += distLatLon(roadNetwork[road]['line'][i], roadNetwork[road]['line'][i + 1])
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(road)

    # Check if there are roads included =======================================
    if (sum(lengths) == 0):
        return None

    # Use accept-denial to test if the node is within poly ====================
    # FIXME: Inefficient approach, will need to be rewritten
    nodeLocs = []
    for i in range(N):
        lat = None
        lon = None
        insideFlag = False
        while (not insideFlag):
            idx = rndPick(lengths)
            edgeLength = lengths[idx]
            edgeDist = random.uniform(0, 1) * edgeLength
            (lat, lon) = getMileageInPathLatLon(roadNetwork[roadIDs[idx]]['line'], edgeDist)
            if (distLatLon([lat, lon], circle['centerLoc']) <= circle['radius']):
                insideFlag = True
        nodeLocs.append((lat, lon))
    return nodeLocs

def getRndPtUniformSquare(
    xRange:    "The range of x coordinates" = (0, 100),
    yRange:    "The range of y coordinates" = (0, 100)
    ) -> "Given the range of x, y, returns a random point in the square defined by the ranges":
    x = random.randrange(xRange[0], xRange[1])
    y = random.randrange(yRange[0], yRange[1])
    return (x, y)

def getRndPtUniformTriangle(
    triangle:   "The triangle for generating random points" = None
    ) -> "Given a triangle, generate a random point in the triangle uniformly":
    
    # Get three extreme points ================================================
    [x1, y1] = triangle[0]
    [x2, y2] = triangle[1]
    [x3, y3] = triangle[2]

    # Generate random points ==================================================
    rndR1 = random.uniform(0, 1)
    rndR2 = random.uniform(0, 1)
    x = (1 - math.sqrt(rndR1)) * x1 + math.sqrt(rndR1) * (1 - rndR2) * x2 + math.sqrt(rndR1) * rndR2 * x3
    y = (1 - math.sqrt(rndR1)) * y1 + math.sqrt(rndR1) * (1 - rndR2) * y2 + math.sqrt(rndR1) * rndR2 * y3

    return (x, y)

def getRndPtUniformPoly(
    poly:       "The polygon for generating random points" = None
    ) -> "Given a polygon, generate a random point in the polygons uniformly":

    # Get list of triangles ===================================================
    lstTriangle = tripy.earclip(poly)

    # Weight them and make draws ==============================================
    lstWeight = []
    for i in range(len(lstTriangle)):
        lstWeight.append(calTriangleAreaByCoords(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2]))

    # Select a triangle and randomize a point in the triangle =================
    idx = rndPick(lstWeight)
    (x, y) = getRndPtUniformTriangle(lstTriangle[idx])

    return (x, y)

def getRndPtUniformPolys(
    polys:       "A list of polygons for generating random points" = None
    ) -> "Given a list of polygons, generate a random point in the polygons uniformly":

    # Get all triangulated triangles ==========================================
    lstTriangle = []
    for p in polys:
        lstTriangle.extend(tripy.earclip(p))

    # Weight them and make draws ==============================================
    lstWeight = []
    for i in range(len(lstTriangle)):
        lstWeight.append(calTriangleAreaByCoords(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2]))

    # Select a triangle and randomize a point in the triangle =================
    idx = rndPick(lstWeight)
    (x, y) = getRndPtUniformTriangle(lstTriangle[idx])

    return (x, y)

def getRndPtClusterXY(
    centroidLocs: "A list of center locs of clusters" = None,
    clusterDiameter: "Diameter of cluster" = None
    ):
    idx = random.randint(0, len(centroidLocs) - 1)
    ctrLoc = centroidLocs[idx]
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, clusterDiameter))
    x = ctrLoc[0] + r * math.cos(theta)
    y = ctrLoc[1] + r * math.sin(theta)
    return (x, y)

def getRndPtClusterLatLon(
    centroidLocs: "A list of center locs of clusters" = None,
    clusterDiameterInMeters: "Diameter of cluster" = None
    ):
    idx = random.randint(0, len(centroidLocs) - 1)
    ctrLoc = centroidLocs[idx]
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, clusterDiameterInMeters))
    (lat, lon) = pointInDistLatLon(ctrLoc, theta, r)
    return (lat, lon)

def getRndPtUniformCircleXY(
    radius:     "Radius of the circle" = None,
    centerLoc:  "Center location of the circle" = None
    ):
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, radius ** 2))
    x = centerLoc[0] + r * math.cos(theta)
    y = centerLoc[1] + r * math.sin(theta)
    return (x, y)

def getRndPtUniformCircleLatLon(
    radius:     "Radius of the circle" = None,
    centerLoc:  "Center location of the circle" = None
    ):
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, radius ** 2))
    (lat, lon) = pointInDistLatLon(centerLoc, theta, r)
    return (lat, lon)

def rectInWidthLengthOrientationXY(
    centroidXY: "Centroid of rectangular in (x, y)" = None,
    width:      "Width of the rectangular" = None,
    length:     "Length of the rectangular" = None,
    oriDeg:     "Orientation of the rectangular" = None
    ) -> "Given args for the width, length, and orientation of rectangular, returns the coordinates in (x, y)":

    # Create four corner points ===============================================
    ptTemp1 = pointInDistXY(centroidXY, oriDeg, length / 2)
    pt1 = pointInDistXY(ptTemp1, oriDeg + 90, width / 2)
    pt2 = pointInDistXY(ptTemp1, oriDeg - 90, width / 2)
    ptTemp2 = pointInDistXY(centroidXY, oriDeg + 180, length / 2)
    pt3 = pointInDistXY(ptTemp2, oriDeg + 90, width / 2)
    pt4 = pointInDistXY(ptTemp2, oriDeg - 90, width / 2)

    # Get the rectangular =====================================================
    rect = [pt1, pt2, pt4, pt3]

    return rect

def rectInWidthLengthOrientationLatLon(
    centroidLatLon: "Centroid of rectangular in (lat, lon)" = None,
    widthInMeter: "Width of the rectangular in meters" = None,
    lengthInMeter: "Length of the rectangular in meters" = None,
    oriDeg:     "Orientation of the rectangular" = None
    ) -> "Given args for the width, length, and orientation of rectangular, returns the coordinates in (lat, lon)":

    # Create four corner points ===============================================
    ptTemp1 = pointInDistLatLon(centroidLatLon, oriDeg, lengthInMeter / 2)
    pt1 = pointInDistLatLon(ptTemp1, oriDeg + 90, widthInMeter / 2)
    pt2 = pointInDistLatLon(ptTemp1, oriDeg - 90, widthInMeter / 2)
    ptTemp2 = pointInDistLatLon(centroidLatLon, oriDeg + 180, lengthInMeter / 2)
    pt3 = pointInDistLatLon(ptTemp2, oriDeg + 90, widthInMeter / 2)
    pt4 = pointInDistLatLon(ptTemp2, oriDeg - 90, widthInMeter / 2)

    # Get the rectangular =====================================================
    rect = [pt1, pt2, pt4, pt3]

    return rect

def getRndPolyXY(

    ) -> "Give a centroid of poly and some arguments, return a randomized polygon in (x, y)":

    return rndPolyXY

def getRndPolyLatLon(

    ) -> "Give a centroid of poly and some arguments, return a randomized polygon in (lat, lon)":

    return rndPolyLatLon

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
    newLoc = list(geopy.distance.distance(meters=distMeters).destination(point=pt, bearing=direction))[:2]
    return newLoc

def getNeighborCluster(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon'" = "Euclidean",
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    diameter:   "The diameter of the neighborhood" = None,
    maxSize:    "Maximum number of nodes in each cluster" = 5
    ) -> "Given a dictionary of locations, an a radius, return the sets that any two locations are within the distance":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define edges ============================================================
    if (edges == 'Euclidean'):
        edges = getTauEuclidean(nodes)
    elif (edges == 'LatLon'):
        edges = getTauLatLon(nodes)
    else:
        print("Error: Incorrect type `edges`")
        return None

    # Initialize ==============================================================
    seedClique = []
    # For each node get the neighbors within diameter
    neighbor = {}
    for n1 in nodeIDs:
        for n2 in nodeIDs:
            if (n1 != n2 and edges[n1, n2] <= diameter):
                if (n1 not in neighbor):
                    neighbor[n1] = []
                neighbor[n1].append(n2)

    # Find seed cliques in the neighbor graph =================================
    # FIXME: Now using stupid method, will be replaced by clique searching algorithm, or should I?
    # 3-clique
    for n1 in neighbor:
        for n2 in neighbor[n1]:
            for n3 in neighbor[n2]:
                if (n1 < n2 < n3 and n1 in neighbor[n3]):
                    seedClique.append([n1, n2, n3])

    # Try to union existing seed-clique to find larger clique =================
    canUnionFlag = True
    while (canUnionFlag):
        canUnionFlag = False
        for c1 in range(len(seedClique) - 1):
            for c2 in range(c1 + 1, len(seedClique)):
                # Try to merge two cliques
                clique1 = [i for i in seedClique[c1]]
                clique2 = [i for i in seedClique[c2]]
                intersect = listSetIntersect(clique1, clique2)                
                # Two cliques can be merged iff they have intersection
                if (len(intersect) > 0):
                    diff1 = listSetMinus(clique1, intersect)
                    diff2 = listSetMinus(clique2, intersect)                    
                    # Try to see if the nodes that are not in the intersection are close
                    mergeFlag = True
                    for n1 in diff1:
                        for n2 in diff2:
                            if (n1 not in neighbor[n2]):
                                mergeFlag = False
                                break
                        if (not mergeFlag):
                            break
                    if (mergeFlag):
                        newClique = listSetUnion(clique1, clique2)
                        if (maxSize == None or len(newClique) <= maxSize):
                            newClique.sort()
                            seedClique.remove(clique1)
                            seedClique.remove(clique2)
                            seedClique.append(newClique)
                            canUnionFlag = True
                            break
            if (canUnionFlag):
                break
    return seedClique

def getSortedNodesByDist(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon'" = "Euclidean",
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
    centerLoc:  "1) (Default) String, 'Centroid' the centroid of nodes, or\
                 2) List, [x, y], the center point" = None,
    isClockwise: "True if the sweeping direction is clock-wise, False otherwise" = True,
    initDeg:    "Starting direction of the sweeping, 0 as North" = 0
    ) -> "Given a set of locations, and a center point, gets the sequence from sweeping":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Initialize centroid =====================================================
    if (centerLoc == 'Centroid'):
        centerLoc = getCentroid(nodes = nodes, nodeIDs = nodeIDs)

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

def getCentroid(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    algo:       "1) String, 'Weiszfeld'" = 'Weiszfeld'
    ) -> "Get centroid location for given nodes":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Subroutine for finding centroid =========================================
    def _getCentroidWeiszfeld(nodes, nodeIDs):
        # Initialize ----------------------------------------------------------
        q = [1 for i in range(len(nodeIDs))]
        a = [nodes[nodeIDs[i]]['loc'][0] for i in range(len(nodeIDs))]
        b = [nodes[nodeIDs[i]]['loc'][1] for i in range(len(nodeIDs))]
        x = sum(a) / len(a)
        y = sum(b) / len(b)
        f = 0
        for i in range(len(nodeIDs)):
            f += math.sqrt((x - a[i])**2 + (y - b[i])**2)

        # Iterations ----------------------------------------------------------
        canGoFlag = True
        while (canGoFlag):
            # update q
            q = []
            for i in range(len(nodeIDs)):
                q.append(1 / math.sqrt((x - a[i])**2 + (y - b[i])**2))

            # update x, y
            x = 0
            y = 0
            for i in range(len(nodeIDs)):
                x += q[i] * a[i]
                y += q[i] * b[i]
            x /= sum(q)
            y /= sum(q)

            # update f
            newF = 0
            for i in range(len(nodeIDs)):
                newF += math.sqrt((x - a[i])**2 + (y - b[i])**2)
            if (abs(newF - f) < CONST_EPSILON):
                canGoFlag = False
            f = newF

        # Output --------------------------------------------------------------
        centroid = (x, y)

        return centroid

    # Call subroutines ========================================================
    centroid = None
    if (algo == "Weiszfeld"):
        centroid = _getCentroidWeiszfeld(nodes, nodeIDs)
    else:
        return None

    return centroid

def getScan(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    direction:  "Direction of scanning, 0 as North" = 0
    ) -> "Scan nodes from one direction":

    # Initialize ==============================================================
    baseline = []
    maxDist = 0
    centroid = getCentroid(nodes)
    for n in nodes:
        d = distEuclidean2D(nodes[n]['loc'], centroid)
        if (maxDist == None or d > maxDist):
            maxDist = 1.2 * d
    basePt = pointInDistXY(centroid, direction, maxDist)
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
