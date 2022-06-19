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
    locsInAreaA:  "A list of (lat, lon) in area A" = None,
    locAnchorA:   "The anchor point of area A" = None,
    locAnchorB:   "The anchor point of area B" = None
    ) -> "Given a set of lat/lon coords in area A, a relative location A of area A, and a relative \
          location B of area B. Map the lat/lon coords to area B":

    # Initialize ==============================================================
    locsInAreaB = []

    # Mapping, using directions and distance ==================================
    for locA in locsInAreaA:
        deg = headingLatLon(locAnchorA, locA)
        dist = distLatLon(locAnchorA, locA)
        locB = ptInDistLatLon(locAnchorB, deg, dist)
        locsInAreaB.append(locB)

    return locsInAreaB

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
        msgError(ERROR_ZERO_VECTOR)
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
        msgError(ERROR_ZERO_VECTOR)
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

def linePerpendicularLine(
    pt:         "Point which the line will go through",
    vec:        "The vector that perpendicular to the line"
    ) -> "Given a point, a vector, returns a line that is going through this point and perpendicular to the vector":

    # Get the direction =======================================================
    heading = headingXY(pt, (pt[0] + vec[0], pt[1] + vec[1]))
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
    line = [ptInDistXY(pt, newHeading1, 10), ptInDistXY(pt, newHeading2, 10)]

    return line

def headingXY(
    pt1:        "Current location", 
    pt2:        "Targeted location"
    ) -> "Given current location and a goal location, calculate the heading. North is 0-degrees, clock-wise":
    
    vec = (pt2[0] - pt1[0], pt2[1] - pt1[1])
    (_, vDeg) = vecXY2Polar(vec)

    return vDeg

def headingLatLon(
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

def mileageInPathLatLon(
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

def rectInWidthLengthOrientationXY(
    centroidXY: "Centroid of rectangular in (x, y)" = None,
    width:      "Width of the rectangular" = None,
    length:     "Length of the rectangular" = None,
    oriDeg:     "Orientation of the rectangular" = None
    ) -> "Given args for the width, length, and orientation of rectangular, returns the coordinates in (x, y)":

    # Create four corner points ===============================================
    ptTemp1 = ptInDistXY(centroidXY, oriDeg, length / 2)
    pt1 = ptInDistXY(ptTemp1, oriDeg + 90, width / 2)
    pt2 = ptInDistXY(ptTemp1, oriDeg - 90, width / 2)
    ptTemp2 = ptInDistXY(centroidXY, oriDeg + 180, length / 2)
    pt3 = ptInDistXY(ptTemp2, oriDeg + 90, width / 2)
    pt4 = ptInDistXY(ptTemp2, oriDeg - 90, width / 2)

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
    ptTemp1 = ptInDistLatLon(centroidLatLon, oriDeg, lengthInMeter / 2)
    pt1 = ptInDistLatLon(ptTemp1, oriDeg + 90, widthInMeter / 2)
    pt2 = ptInDistLatLon(ptTemp1, oriDeg - 90, widthInMeter / 2)
    ptTemp2 = ptInDistLatLon(centroidLatLon, oriDeg + 180, lengthInMeter / 2)
    pt3 = ptInDistLatLon(ptTemp2, oriDeg + 90, widthInMeter / 2)
    pt4 = ptInDistLatLon(ptTemp2, oriDeg - 90, widthInMeter / 2)

    # Get the rectangular =====================================================
    rect = [pt1, pt2, pt4, pt3]

    return rect

def ptInDistXY(
    pt:         "Starting location" = None, 
    direction:  "Direction towards the destination, Up is 0-degrees, clock-wise" = None, 
    dist:       "Distance from origin location to destination location" = None
    ) -> "A location in distance with given direction, in [lat, lon] form.":
    x = pt[0] + dist * math.sin(math.radians(direction))
    y = pt[1] + dist * math.cos(math.radians(direction))
    return (x, y)

def ptInDistLatLon(
    pt:         "Starting location" = None, 
    direction:  "Direction towards the destination, North is 0-degree, East is 90-degrees" = None, 
    distMeters: "Distance from origin location to destination location" = None
    ) -> "A location in distance with given direction, in [lat, lon] form.":
    newLoc = list(geopy.distance.distance(meters=distMeters).destination(point=pt, bearing=direction))[:2]
    return newLoc
