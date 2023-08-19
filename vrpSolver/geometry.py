import geopy.distance
import heapq
import math

# Shapely can do a lot of things, but it seems that shapely does not 
#     support (infinite) lines and rays.
import shapely

from .error import *
from .common import *
from .const import *
from .graph import *
from .msg import *


import math

from .const import *
from .msg import *

# Point versus Objects ========================================================
def is2PtsSame(pt1: pt, pt2: pt) -> bool:
    """Are two points at the 'same' location"""
    if (abs(pt1[0] - pt2[0]) >= CONST_EPSILON):
        return False
    if (abs(pt1[1] - pt2[1]) >= CONST_EPSILON):
        return False
    return True

def is3PtsClockWise(pt1: pt, pt2: pt, pt3: pt) -> bool | None:
    """Are three given pts in a clock-wise order, None as they are colliner"""
    if (is2PtsSame(pt1, pt2) or is2PtsSame(pt2, pt3) or is2PtsSame(pt1, pt3)):
        # If points are overlapped, return None as collinear
        return None
    [x1, y1] = [pt1[0], pt1[1]]
    [x2, y2] = [pt2[0], pt2[1]]
    [x3, y3] = [pt3[0], pt3[1]]
    ori = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
    if (abs(ori) <= CONST_EPSILON):
        # collinear
        return None
    elif (ori < 0):
        # clockwise
        return True
    else:
        # counter-clockwise
        return False

def isPtOnLine(pt: pt, line: line) -> bool:
    """Is a pt on the line."""
    if (is2PtsSame(line[0], line[1])):
        raise ZeroVectorError()
    if (is3PtsClockWise(pt, line[0], line[1]) == None):
        return True
    else:
        return False

def isPtOnSeg(pt: pt, seg: line, interiorOnly: bool=False) -> bool:
    """Is a pt on the segment"""
    onLine = isPtOnLine(pt, seg)
    if (onLine == False):
        return False
    # Get pts =================================================================
    [x1, y1] = [seg[0][0], seg[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [seg[1][0], seg[1][1]]
    onSeg = (
        abs(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) 
        + math.sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2) 
        - math.sqrt((x1 - x3) ** 2 + (y1 - y3) ** 2)) <= CONST_EPSILON)
    # Check if the intertion is in the interior ===============================
    if (interiorOnly):
        return onSeg and not is2PtsSame(pt, seg[0]) and not is2PtsSame(pt, seg[1])
    else:
        return onSeg

def isPtOnRay(pt: pt, ray: line, interiorOnly: bool=False) -> bool:
    """Is a pt on the ray, could be at the end."""
    onLine = isPtOnLine(pt, ray)
    if (onLine == False):
        return False
    # Get pts =================================================================
    [x1, y1] = [ray[0][0], ray[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [ray[1][0], ray[1][1]]
    onRay = (
        (abs(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) 
            + math.sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2) 
            - math.sqrt((x1 - x3) ** 2 + (y1 - y3) ** 2)) <= CONST_EPSILON)
        or (math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) >= math.sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2)))
    # Check if the intertion is in the interior ===============================
    if (interiorOnly):
        return onRay and not is2PtsSame(pt, ray[0])
    else:
        return onRay

def isPtOnPolyEdge(pt: pt, poly: poly) -> bool:
    """Is a point on any edge of a polygon"""
    # Check if the pt is on any of the edge segment ===========================
    for i in range(-1, len(poly) - 1):
        if (isPtOnSeg(pt, [poly[i], poly[i + 1]])):
            return True
    return False

def isPtInPoly(pt: pt, poly: poly, interiorOnly: bool=False) -> bool:
    """Is a pt in the polygon, could be on the edge"""
    x = pt[1]
    y = pt[0]
    inPoly = False
    # NOTE: 忘记这段魔法是从哪里来的了...
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
            inPoly = not inPoly
        j = i
    # Check if the intertion is in the interior ===============================
    if (interiorOnly):
        return inPoly and not isPtOnPolyEdge(pt, poly)
    else:
        return inPoly

# Line-shape versus Line-shape ===================================================
def isLineIntLine(line1: line, line2: line) -> bool:
    """Is two line intersect with each other"""
    intPt = intLine2Line(line1, line2)
    if (intPt == None):
        return False
    elif (intPt == "Collinear"):
        return True
    else:
        return True

def isLineIntSeg(line: line, seg: line, interiorOnly:bool=False) -> bool:
    """Is a line intersect with a segment"""
    intPt = intLine2Line(line, seg)
    if (intPt == None):
        return False
    elif (intPt == "Collinear"):
        return True
    else:
        if (isPtOnSeg(intPt, seg, interiorOnly)):
            return True
        else:
            return False

def isLineIntRay(line: line, ray: line, interiorOnly: bool=False) -> bool:
    """Is a line intersect with a ray"""
    intPt = intLine2Line(line, ray)
    if (intPt == None):
        return False
    elif (intPt == "Collinear"):
        return True
    else:
        if (isPtOnRay(intPt, ray, interiorOnly)):
            return True
        else:
            return False

def isSegIntLine(seg: line, line: line, interiorOnly: bool=False) -> bool:
    """Is a segment intersect with a line"""
    return isLineIntSeg(line, seg, interiorOnly)

def isSegIntSeg(seg1: line, seg2: line, interiorOnly: bool=False) -> bool:
    """Are two segment intersect with each other, could be at the end of a segment"""
    intPt = intLine2Line(seg1, seg2)
    if (intPt == None):
        return False
    elif (intPt == "Collinear"):
        if (isPtOnSeg(seg1[0], seg2, interiorOnly=False) 
                or isPtOnSeg(seg1[1], seg2, interiorOnly=False)
                or isPtOnSeg(seg2[0], seg1, interiorOnly=False)
                or isPtOnSeg(seg2[1], seg1, interiorOnly=False)):
            if (interiorOnly
                    and ((is2PtsSame(seg1[0], seg2[0]) and not is2PtsSame(seg1[1], seg2[1]))
                        or (is2PtsSame(seg1[0], seg2[1]) and not is2PtsSame(seg1[1], seg2[0]))
                        or (is2PtsSame(seg1[1], seg2[1]) and not is2PtsSame(seg1[0], seg2[0]))
                        or (is2PtsSame(seg1[1], seg2[0]) and not is2PtsSame(seg1[0], seg2[1])))):
                return False
            return True
        else:
            return False
    else:
        if (isPtOnSeg(intPt, seg1, interiorOnly) and isPtOnSeg(intPt, seg2, interiorOnly)):
            return True
        else:
            return False

def isSegIntRay(seg: line, ray: line, interiorOnly: bool=False) -> bool:
    """Is a segment intersect with a ray, could be at the end of the segment/ray"""
    intPt = intLine2Line(seg, ray)
    if (intPt == None):
        return False
    elif (intPt == "Collinear"):
        if (isPtOnSeg(ray[0], seg, interiorOnly=False)):
            if (interiorOnly
                    and ((is2PtsSame(seg[0], ray[0]) and not isPtOnRay(seg[1], ray))
                        or (is2PtsSame(seg[1], ray[0]) and not isPtOnRay(seg[0], ray)))):
                return False
            return True
        else:
            return False
    else:
        if (isPtOnSeg(intPt, seg, interiorOnly) and isPtOnRay(intPt, ray, interiorOnly)):
            return True
        else:
            return False

def isRayIntLine(ray: line, line: line, interiorOnly: bool=False) -> bool:
    """Is a ray intersect with a line"""
    return isLineIntRay(line, ray, interiorOnly)

def isRayIntSeg(ray: line, seg: line, interiorOnly: bool=False) -> bool:
    """Is a ray intersect with a segment"""
    return isSegIntRay(seg, ray, interiorOnly)

def isRayIntRay(ray1: line, ray2: line, interiorOnly: bool=False) -> bool:
    """Is a segment intersect with a ray, could be at the end of the segment/ray"""
    intPt = intLine2Line(ray1, ray2)
    if (intPt == None):
        return False
    elif (intPt == "Collinear"):
        if (isPtOnRay(ray1[0], ray2, interiorOnly=False) or isPtOnRay(ray2[0], ray1, interiorOnly=False)):
            if (interiorOnly
                    and (is2PtsSame(ray1[0], ray2[0]) and not isPtOnRay(ray1[1], ray2))):
                return False
            return True
        else:
            return False
    else:
        if (isPtOnRay(intPt, ray1, interiorOnly) and isPtOnRay(intPt, ray2, interiorOnly)):
            return True
        else:
            return False

# Line-shape versus polygon ===================================================
def isLineIntPoly(line: line, poly: poly, interiorOnly: bool=False) -> bool:
    """Is a line intersect with a polygon"""
    for i in range(-1, len(poly) - 1):
        edge = [poly[i], poly[i + 1]]
        if (isLineIntSeg(line, edge, interiorOnly)):
            return True

def isSegIntPoly(seg: line, poly: poly, interiorOnly: bool=False) -> bool:
    """Is a segment intersect with a polygon"""
    # 若线段中的一个端点妥妥的在poly的interior里，interiorOnly True
    if (isPtInPoly(seg[0], poly, interiorOnly) or isPtInPoly(seg[1], poly, interiorOnly)):
        return True

    # 若线段和多边形的至少一条边妥妥的相交，interiorOnly True
    for i in range(-1, len(poly) - 1):
        edge = [poly[i], poly[i + 1]]
        if (not is2PtsSame(poly[i], poly[i + 1])):
            if (isSegIntSeg(edge, seg, interiorOnly)):
                return True

    # 若线段只和多边形的boundary相交： 线段的一个顶点在poly edge，或一个顶点在线段上
    if (interiorOnly == False):
        if (isPtOnPolyEdge(seg[0], poly) or isPtOnPolyEdge(seq[1], poly)):
            return True
        for pt in poly:
            if (isPtOnSeg(pt, seg, interiorOnly=False)):
                return True

    # 最后最麻烦的情况：如果线段的两端都在poly edge上，是否相交呢？
    # NOTE: 思路是找线段interior上的任何一点，若在poly内，就相交，否则一定不相交（interior-wise）
    midPt = [seg[0][0] + (seg[1][0] - seg[0][0])/2, seg[0][1] + (seg[1][1] - seg[0][1])/2]
    if (isPtInPoly(midPt, poly, interiorOnly=True)):
        return True
    else:
        if (interiorOnly == True):
            return False
        else:
            return True
    return False

def isRayIntPoly(ray: line, poly: poly, interiorOnly: bool=False) -> bool:
    """Is a ray intersect with a polygon"""
    # 若射线的端点妥妥的在poly的interior里，interiorOnly True
    if (isPtInPoly(ray[0], poly, interiorOnly)):
        return True

    # Check number of intersection between ray and each line segment
    intList = []
    for i in range(-1, len(poly) - 1):
        edge = [poly[i], poly[i + 1]]
        sxr = intSeg2Ray(edge, ray)
        if (sxr != None):
            intList.append(sxr)

    if (len(intList) >= 3):
        return True
    elif (len(intList) == 2):
        if (interiorOnly == True):
            return False
        else:
            return True
    elif (len(intList) <= 1):
        return False

    return False


# Distance from Point to Object ===============================================
def distPt2Line(pt: pt, line: line) -> float:
    return

def distPt2Seg(pt: pt, seg: line) -> float:
    return

def distPt2Seq(pt: pt, seq: list[pt]) -> float:
    return


# Intersection calculation ====================================================
def intLine2Line(line1: line, line2: line) -> pt | None | str:
    # Validation ==============================================================
    if (is2PtsSame(line1[0], line1[1])):
        raise ZeroVectorError(line1)
    if (is2PtsSame(line2[0], line2[1])):
        raise ZeroVectorError(line2)

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
        if (is3PtsClockWise(line1[0], line1[1], line2[0]) == None):
            return {
                'status': 'Collinear',
                'intersect': line1,
                'intersectType': 'Line'
            }
        else:
            return {
                'status': 'NoCross',
                'intersect': None,
                'intersectType': None
            }

    # Intersection ============================================================
    x = (b1 * c2 - b2 * c1) / D
    y = (a2 * c1 - a1 * c2) / D
    return {
        'status': 'Cross',
        'intersect': (x, y),
        'intersectType': 'Point'
    }

def intLine2Seg(line: line, seg: line) -> pt | None | line:
    intPt = intLine2Line(line, seg)

    if (intPt['status'] == 'NoCross'):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None
        }
    elif (intPt['status'] == 'Collinear'):
        return {
            'status': 'Collinear',
            'intersect': seg,
            'intersectType': 'Segment'
        }
    else:
        if (isPtOnSeg(intPt['intersect'], seg, interiorOnly=False)):
            return intPt
        else:
            return {
                'status': 'NoCross',
                'intersect': None,
                'intersectType': None
            }

def ptLine2Ray(line: line, ray: line) -> pt | None | str:
    intPt = intLine2Line(line, ray)

    if (intPt['status'] == 'NoCross'):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None
        }
    elif (intPt == "Collinear"):
        return {
            'status': 'Collinear',
            'intersect': ray,
            'intersectType': 'Ray'
        }
    else:
        if (isPtOnRay(intPt['intersect'], ray, interiorOnly=False)):
            return intPt
        else:
            return {
                'status': 'NoCross',
                'intersect': None,
                'intersectType': None
            }

def ptSeg2Line(seg: line, line: line) -> pt | None | line:
    return intLine2Seg(line, seg)

def ptSeg2Seg(seg1: line, seg2: line) -> pt | None | line:
    intPt = intLine2Line(seg1, seg2)
    if (intPt == None):
        return None
    elif (intPt == "Collinear"):
        if (isPtOnSeg(seg1[0], seg2, interiorOnly=False) 
                or isPtOnSeg(seg1[1], seg2, interiorOnly=False)
                or isPtOnSeg(seg2[0], seg1, interiorOnly=False)
                or isPtOnSeg(seg2[1], seg1, interiorOnly=False)):
            if (interiorOnly
                    and ((is2PtsSame(seg1[0], seg2[0]) and not is2PtsSame(seg1[1], seg2[1]))
                        or (is2PtsSame(seg1[0], seg2[1]) and not is2PtsSame(seg1[1], seg2[0]))
                        or (is2PtsSame(seg1[1], seg2[1]) and not is2PtsSame(seg1[0], seg2[0]))
                        or (is2PtsSame(seg1[1], seg2[0]) and not is2PtsSame(seg1[0], seg2[1])))):
                return False
            return True
        else:
            return False

        # 
        if (isPtOnSeg(seg1[0], seg2, interiorOnly=False)
                and (seg1[1], seg2)):

    else:
        if (isPtOnSeg(intPt, seg1, interiorOnly=False)
                and isPtOnSeg(intPt, seg2, interiorOnly=False)):
            return intPt
        else:
            return None

def ptSeg2Ray(seg: line, ray: line) -> pt | None | line:
    return

def ptRay2Line(ray: line, line: line) -> pt | None | line:
    return

def ptRay2Seg(ray: line, seg: line) -> pt | None | line:
    return

def ptRay2Ray(ray1: line, ray2: line) -> pt | None | line:
    return

# 
def vecPolar2XY(vecPolar: pt) -> pt:

    """Given vector's norm and its degree to North, convert it into a 2-tuple vector

    Parameters
    ----------
    vecPolar: tuple[float|int, float|int], required
        2-tuple (vVal, vDeg), `vVal` is the norm and `vDeg` is the direction, 0 as North, clockwise, in [0, 360)

    Returns
    -------

    tuple[float|int, float|int]

    """

    # Initialize ==============================================================
    (vVal, vDeg) = vecPolar

    vX = 0
    vY = 0

    while(vDeg < 0):
        vDeg = vDeg + 360

    while(vDeg >= 360):
        vDeg = vDeg - 360

    vX = vVal * math.sin(math.radians(vDeg))
    vY = vVal * math.cos(math.radians(vDeg))

    return (vX, vY)

def vecXY2Polar(vecXY: pt) -> pt:

    """Given a 2-tuple, convert it into a norm and a direction in degree

    Parameters
    ----------

    vecXY: tuple[float|int, float|int], required
        2-tuple (vX, vY), the coordinate of vector


    Returns
    -------

    tuple[float|int, float|int]

    """    
    (vX, vY) = vecXY    
    vDeg = 0
    vVal = 0
    if (abs(vX) <= CONST_EPSILON):
        if (vY >= 0):
            vDeg = 0
            vVal = vY
        elif (vY < 0):
            vDeg = 180
            vVal = -vY
    elif (abs(vY) <= CONST_EPSILON):
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

    return (vVal, vDeg)

def ptXY2LatLonMercator(ptXY: pt) -> pt:
    """Given a point in (x, y), try to map it with a (lat, lon) coordinate with necessary inputs"""

    # Ref: https://wiki.openstreetmap.org/wiki/Mercator#Python
    (x, y) = ptXY
    lon = math.degrees(y / CONST_EARTH_RADIUS_METERS)
    lat = math.degrees(2 * math.atan(math.exp(x / CONST_EARTH_RADIUS_METERS)) - math.pi / 2.0)
    ptLatLon = (lat, lon)
    return ptLatLon

def ptLatLon2XYMercator(ptLatLon: pt) -> pt:
    """Given a point int (lat, lon), try to map it with a (x, y) coordinates with necessary inputs"""

    # Ref: https://wiki.openstreetmap.org/wiki/Mercator#Python
    (lat, lon) = ptLatLon
    y = math.radians(lon) * CONST_EARTH_RADIUS_METERS
    x = math.log(math.tan(math.pi / 4 + math.radians(lat) / 2)) * CONST_EARTH_RADIUS_METERS
    ptXY = (x, y)
    return ptXY

def distEuclideanXY(pt1: pt, pt2: pt) -> float:
    """Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number"""
    return math.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)

def distManhattenXY(pt1: pt, pt2: pt) -> float:
    """Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number"""
    return abs(pt1[0] - pt2[0]) + abs(pt1[1] - pt2[1])

def distLatLon(pt1: pt, pt2: pt, distUnit: str = 'meter') -> float:
    """Gives a Euclidean distance based on two lat/lon coords, if two coordinates are the same, return a small number"""
    
    # Get radius as in distUnit ===============================================
    R = None
    if (distUnit in ['mile', 'mi']):
        R = CONST_EARTH_RADIUS_MILES
    elif (distUnit in ['meter', 'm']):
        R = CONST_EARTH_RADIUS_METERS
    elif (distUnit in ['kilometer', 'km']):
        R = CONST_EARTH_RADIUS_METERS / 1000
    else:
        raise UnsupportedInputError("ERROR: Unrecognized distance unit, options are 'mile', 'meter', 'kilometer'")

    # Calculate distance ======================================================
    (lat1, lon1) = pt1
    (lat2, lon2) = pt2
    phi1, phi2 = math.radians(lat1), math.radians(lat2) 
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))

# Time seq related ============================================================
def locInTimedSeq(seq: list[pt], timeStamp: list[float], t: float) -> pt:
    if (len(seq) != len(timeStamp)):
        raise UnsupportedInputError("ERROR: `timeStamp` does not match with `seq`.")
    for i in range(len(timeStamp) - 1):
        if (timeStamp[i] > timeStamp[i + 1]):
            raise UnsupportedInputError("ERROR: `timeStamp` should be a non-descending sequence.")

    if (t <= timeStamp[0]):
        return seq[0]
    if (t >= timeStamp[-1]):
        return seq[-1]

    for i in range(len(timeStamp) - 1):
        if (timeStamp[i] <= t < timeStamp[i + 1]):
            if (timeStamp[i] == timeStamp[i + 1]):
                raise UnsupportedInputError("ERROR: an object cannot be two places at the same time.")
            curLocX = seq[i][0] + (seq[i + 1][0] - seq[i][0]) * (t - timeStamp[i]) / (timeStamp[i + 1] - timeStamp[i])
            curLocY = seq[i][1] + (seq[i + 1][1] - seq[i][1]) * (t - timeStamp[i]) / (timeStamp[i + 1] - timeStamp[i])
            return [curLocX, curLocY]
    raise UnsupportedInputError("ERROR: cannot find time stamp")

def speedInTimedSeq(seq: list[pt], timeStamp: list[float], t: float) -> float:
    spd = 0

    if (len(seq) != len(timeStamp)):
        raise UnsupportedInputError("ERROR: `timeStamp` does not match with `seq`.")
    for i in range(len(timeStamp) - 1):
        if (timeStamp[i] > timeStamp[i + 1]):
            raise UnsupportedInputError("ERROR: `timeStamp` should be a non-descending sequence.")

    if (t <= timeStamp[0]):
        return 0
    if (t >= timeStamp[-1]):
        return 0

    for i in range(len(timeStamp) - 1):
        if (timeStamp[i] <= t < timeStamp[i + 1]):
            if (timeStamp[i] == timeStamp[i + 1]):
                raise UnsupportedInputError("ERROR: an object cannot be two places at the same time.")
            dist = distEuclideanXY(seq[i], seq[i + 1])
            spd = dist / (timeStamp[i + 1] - timeStamp[i])

    return spd

def traceInTimedSeq(seq: list[pt], timeStamp: list[float], ts: float, te: float) -> list[pt]:
    trace = []

    if (len(seq) != len(timeStamp)):
        raise UnsupportedInputError("ERROR: `timeStamp` does not match with `seq`.")
    for i in range(len(timeStamp) - 1):
        if (timeStamp[i] > timeStamp[i + 1]):
            raise UnsupportedInputError("ERROR: `timeStamp` should be a non-descending sequence.")
    if (ts >= te):
        raise UnsupportedInputError("ERROR: `ts` should be earlier than `te`")

    if (ts <= timeStamp[0] and te >= timeStamp[-1]):
        return [pt for pt in seq]
    if (ts >= timeStamp[-1]):
        return []
    if (te <= timeStamp[0]):
        return []

    tsIndex = -1
    teIndex = -1
    tsLoc = []
    teLoc = []

    if (ts <= timeStamp[0]):
        tsIndex = 0
        ts = timeStamp[0]
        tsLoc = seq[0]
    if (te >= timeStamp[-1]):
        teIndex = len(timeStamp) - 1
        te = timeStamp[-1]
        teLoc = seq[-1]

    for i in range(len(timeStamp) - 1):
        if (timeStamp[i] <= ts < timeStamp[i + 1]):
            if (timeStamp[i] == timeStamp[i + 1]):
                raise UnsupportedInputError("ERROR: an object cannot be two places at the same time.")
            tsIndex = i
            tsX = seq[i][0] + (seq[i + 1][0] - seq[i][0]) * (ts - timeStamp[i]) / (timeStamp[i + 1] - timeStamp[i])
            tsY = seq[i][1] + (seq[i + 1][1] - seq[i][1]) * (ts - timeStamp[i]) / (timeStamp[i + 1] - timeStamp[i])
            tsLoc = [tsX, tsY]
            for j in range(tsIndex, len(timeStamp) - 1):
                if (timeStamp[j] <= te < timeStamp[j + 1]):
                    if (timeStamp[j] == timeStamp[j + 1]):
                        raise UnsupportedInputError("ERROR: an object cannot be two places at the same time.")
                    teIndex = j
                    teX = seq[j][0] + (seq[j + 1][0] - seq[j][0]) * (te - timeStamp[j]) / (timeStamp[j + 1] - timeStamp[j])
                    teY = seq[j][1] + (seq[j + 1][1] - seq[j][1]) * (te - timeStamp[j]) / (timeStamp[j + 1] - timeStamp[j])
                    teLoc = [teX, teY]

    if (tsIndex == teIndex):
        trace = [tsLoc, teLoc]
    elif (tsIndex + 1 == teIndex):
        trace.append(tsLoc)
        trace.append(seq[tsIndex + 1])
        trace.append(teLoc)
    else:
        trace.append(tsLoc)
        for i in range(tsIndex + 1, teIndex + 1):
            trace.append(seq[i])
        trace.append(teLoc)

    return trace


def locOnPolyExt(poly: poly, startPt: pt, dist: int|float, reverseFlag: bool=False, dimension: str = 'XY') -> pt:
    perimeter = calPolygonPerimeter(poly)
    while(dist > perimeter):
        dist -= perimeter

    seq = [startPt]
    startID = len(poly) + 1

    if (dimension == 'XY'):
        for i in range(len(poly) - 1):
            if (abs(distEuclideanXY(poly[i], poly[i + 1]) - distEuclideanXY(poly[i], startPt) - distEuclideanXY(startPt, poly[i + 1])) <= CONST_EPSILON):
                startID = i + 1
                break
        if (startID == len(poly) + 1):
            if (abs(distEuclideanXY(poly[0], poly[-1]) - distEuclideanXY(poly[0], startPt) - distEuclideanXY(startPt, poly[-1])) <= CONST_EPSILON):
                startID = 0
            else:
                raise UnsupportedInputError("ERROR: `startPt` should be on the boundary of `poly`")
    elif (dimension == 'LatLon'):
        for i in range(len(poly) - 1):
            if (abs(distLatLon(poly[i], poly[i + 1]) - distLatLon(poly[i], startPt) - distLatLon(startPt, poly[i + 1])) <= CONST_EPSILON):
                startID = i + 1
                break
        if (startID == len(poly) + 1):
            if (abs(distLatLon(poly[0], poly[-1]) - distLatLon(poly[0], startPt) - distLatLon(startPt, poly[-1])) <= CONST_EPSILON):
                startID = 0
            else:
                raise UnsupportedInputError("ERROR: `startPt` should be on the boundary of `poly`")
    else:
        raise UnsupportedInputError("ERROR: options for parameter `dimension` includes ['XY', 'LatLon']")

    if (not reverseFlag):
        for i in range(startID, len(poly)):
            seq.append(poly[i])
        for i in range(startID):
            seq.append(poly[i])
        seq.append(startPt)
    else:
        for i in range(startID, len(poly)):
            seq.insert(0, poly[i])
        for i in range(startID):
            seq.insert(0, poly[i])
        seq.insert(0, startPt) 

    loc = locInSeq(seq, dist, dimension)

    return loc

def traceOnPolyExt(poly: poly, startPt: pt, endPt:pt, dist: int|float, reverseFlag: bool=False) -> list[pt]:
    return trace

def locInSeq(seq: list[pt], dist: int|float, dimension: str = 'XY') -> pt:
    """Given a list of lat/lon coordinates, and a traveling mileage, returns the coordinate"""

    # Initialize ==============================================================
    inPathFlag = False
    accDist = 0
    preLoc = []
    nextLoc = []

    # Find segment ============================================================
    for i in range(0, len(seq) - 1):
        if (dimension == 'LatLon'):
            accDist += distLatLon(seq[i], seq[i + 1])
        elif (dimension == 'XY'):
            accDist += distEuclideanXY(seq[i], seq[i + 1])
        if (accDist > dist):
            preLoc = seq[i]
            nextLoc = seq[i + 1]
            inPathFlag = True
            break

    if (inPathFlag == False):
        raise UnsupportedInputError("ERROR: `dist` is longer than the length of `seq`")

    # Find location on the segment ============================================
    remainDist = accDist - dist
    segDist = 0
    if (dimension == 'LatLon'):
        segDist = distLatLon(preLoc, nextLoc)
    elif (dimension == 'XY'):
        segDist = distEuclideanXY(preLoc, nextLoc)
    if (segDist <= CONST_EPSILON):
        raise ZeroDivisionError
    lat = nextLoc[0] + (remainDist / segDist) * (preLoc[0] - nextLoc[0])
    lon = nextLoc[1] + (remainDist / segDist) * (preLoc[1] - nextLoc[1])
    return (lat, lon)


def calPolygonPerimeter(poly: poly) -> float:
    p = 0
    for i in range(-1, len(poly) - 1):
        p += distEuclideanXY(poly[i], poly[i + 1])
    return p

def calTriangleAreaEdge(a: float, b: float, c: float) -> float:
    # Using Heron's Formula ===================================================
    s = (a / 2 + b / 2 + c / 2)
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return area

def calTriangleAreaXY(pt1: pt, pt2: pt, pt3: pt) -> float:
    # Using determinant =======================================================
    [x1, y1] = pt1
    [x2, y2] = pt2
    [x3, y3] = pt3
    val = (x2 * y3 + x3 * y1 + x1 * y2) - (x2 * y1 + x3 * y2 + x1 * y3)
    area = abs(val)
    return area

def calPolygonAreaLatLon(polyLatLon: poly) -> float:
    """Returns the area surrounded by polyLatLon on the Earth"""

    # Additional packages =====================================================
    from pyproj import Geod

    # Create polygon ==========================================================
    # NOTE: shapely is in [lon, lat] format
    rev = []
    for p in polyLatLon:
        rev.append((p[1], p[0]))
    polygon = shapely.Polygon(rev)

    # Using pyproj to calculate ===============================================
    # Ref: https://hypc.github.io/2020/03/16/python-geo-area/
    geod = Geod(ellps = "WGS84")
    area = abs(geod.geometry_area_perimeter(polygon)[0])

    return area


def headingXY(pt1: pt, pt2: pt) -> float:
    
    vec = (pt2[0] - pt1[0], pt2[1] - pt1[1])
    (_, vDeg) = vecXY2Polar(vec)

    return vDeg

def headingLatLon(pt1: pt, pt2: pt) -> float:

    """Given current location and a goal location, calculate the heading. North is 0-degrees, clock-wise"""

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

def ptInDistXY(pt: pt, direction: int|float, dist: int|float):
    """A location in distance with given direction, in [lat, lon] form."""
    x = pt[0] + dist * math.sin(math.radians(direction))
    y = pt[1] + dist * math.cos(math.radians(direction))
    return (x, y)

def ptInDistLatLon(pt: pt, direction: int|float, distMeters: int|float):
    """A location in distance with given direction, in [lat, lon] form."""
    # Bearing in degrees: 0 – North, 90 – East, 180 – South, 270 or -90 – West.
    newLoc = list(geopy.distance.distance(meters=distMeters).destination(point=pt, bearing=direction))[:2]
    return newLoc


def getTau(
    nodes: dict, 
    edges: dict,
    depotID: int|str = 0,
    nodeIDs: list|str = 'All',
    serviceTime: float = 0
    ) -> dict:

    # Define tau ==============================================================
    tau = {}
    if (type(edges) != dict or 'method' not in edges):
        raise MissingParameterError(ERROR_MISSING_EDGES)

    if (edges['method'] == 'Euclidean'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        tau = _getTauEuclidean(nodes, nodeIDs, ratio)
    elif (edges['method'] == 'LatLon'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        tau = _getTauLatLon(nodes, nodeIDs, speed=ratio)
    elif (edges['method'] == 'Manhatten'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        tau = _getTauManhatten(nodes, nodeIDs, ratio)
    elif (edges['method'] == 'Dictionary'):
        if ('dictionary' not in edges or edges['dictionary'] == None):
            raise MissingParameterError("'dictionary' is not specified")
        for p in edges['dictionary']:
            ratio = 1 if 'ratio' not in edges else edges['ratio']
            tau[p] = edges['dictionary'][p] * ratio
    elif (edges['method'] == 'Grid'):
        if ('grid' not in edges or edges['grid'] == None):
            raise MissingParameterError("'grid' is not specified")
        if ('column' not in edges['grid'] or 'row' not in edges['grid']):
            raise MissingParameterError("'column' and 'row' need to be specified in 'grid'")
        tau = _getTauGrid(nodes, nodeIDs, edges['grid'])
    else:
        raise UnsupportedInputError(ERROR_MISSING_EDGES)        

    # Service time ============================================================
    if (depotID != None and serviceTime != None and serviceTime > 0):
        for (i, j) in tau:
            if (i != depotID and j != depotID and i != j):
                tau[i, j] += serviceTime
            elif (i == depotID or j == depotID and i != j):
                tau[i, j] += serviceTime / 2 

    return tau

def _getTauEuclidean(nodes: dict, nodeIDs: list|str, speed = 1):
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
                t = distEuclideanXY(nodes[i]['loc'], nodes[j]['loc']) / speed
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = CONST_EPSILON
                tau[j, i] = CONST_EPSILON
    return tau

def _getTauManhatten(nodes: dict, nodeIDs: list|str, speed = 1):
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
                t = distManhattenXY(nodes[i]['loc'], nodes[j]['loc']) / speed
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = CONST_EPSILON
                tau[j, i] = CONST_EPSILON
    return tau

def _getTauLatLon(nodes: dict, nodeIDs: list|str, distUnit = 'meter', speed = 1):
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

def _getTauGrid(nodes: dict, nodeIDs: list|str, grid: dict):
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
                t = gridPathFinding(grid = grid, startCoord = nodes[i]['loc'], endCoord = nodes[j]['loc'])['dist']
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = 0
                tau[j, i] = 0
    return tau

def getSortedNodesByDist(nodes: dict, refLoc: pt, nodeIDs: list|str = 'All') -> list:
    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Sort distance ===========================================================
    sortedSeq = []
    sortedSeqHeap = []
    for n in nodeIDs:
        dist = distEuclideanXY(refLoc, nodes[n]['loc'])
        heapq.heappush(sortedSeqHeap, (dist, n))
    while (len(sortedSeqHeap) > 0):
        sortedSeq.append(heapq.heappop(sortedSeqHeap)[1])  

    return sortedSeq

def getSweepSeq(nodes: dict, nodeIDs: list|str = 'All', centerLoc: None|pt = None, isClockwise: bool = True, initDeg: float = 0) -> list:
    """Given a set of locations, and a center point, gets the sequence from sweeping"""
    
    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Initialize centroid =====================================================
    if (centerLoc == None):
        lstNodeLoc = []
        for n in nodeIDs:
            lstNodeLoc.append(shapely.Point(nodes[n]['loc'][0], nodes[n]['loc'][1]))
        centerLoc = list(shapely.centroid(shapely.MultiPoint(points = lstNodeLoc)))

    # Initialize heap =========================================================
    degHeap = []
    centerLocNodes = []
    
    # Build heap ==============================================================
    for n in nodeIDs:
        dist = distEuclideanXY(nodes[n]['loc'], centerLoc)
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
            heapq.heappush(degHeap, (evalDeg, dist, n))

    # Sweep ===================================================================
    sweepSeq = []
    while (len(degHeap)):
        sweepSeq.append(heapq.heappop(degHeap)[2])
    sweepSeq.extend(centerLocNodes)

    return sweepSeq
