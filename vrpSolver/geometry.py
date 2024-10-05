import geopy.distance
import heapq
import math
import tripy
import warnings
import random
import gurobipy as grb

import shapely
from shapely.geometry import mapping
from shapely.ops import nearest_points
import networkx as nx

from .common import *
from .msg import *
from .ds import *

# Point versus Objects ========================================================
def is2PtsSame(pt1: pt, pt2: pt, error: float = CONST_EPSILON) -> bool:
    """
    Are two points at the 'same' location?

    Parameters
    ----------
    pt1: pt, required
        Coordinate of the first point
    pt2: pt, required
        Coordinate of the second point
    error: float, optional, default as CONST_EPSILON
        Error tolerance

    Return
    ------
    bool
        True if two points are at the same location, False else-wise

    """
    if (abs(pt1[0] - pt2[0]) >= error):
        return False
    if (abs(pt1[1] - pt2[1]) >= error):
        return False
    return True

def is3PtsClockWise(pt1: pt, pt2: pt, pt3: pt, error: float = CONST_EPSILON) -> bool | None:
    """
    Are three given pts in a clock-wise order, None as they are collinear

    Parameters
    ----------
    pt1: pt, required
        Coordinate of the first point
    pt2: pt, required
        Coordinate of the second point
    pt3: pt, required
        Coordinate of the third point
    error: float, optional, default as CONST_EPSILON
        Error tolerance

    Return
    ------
    bool | None
        True if three given points are in a clock-wise order, False if they are in counter-clock-wise order, None if they are collinear
    """
    if (is2PtsSame(pt1, pt2) or is2PtsSame(pt2, pt3) or is2PtsSame(pt1, pt3)):
        # If points are overlapped, return None as collinear
        return None
    [x1, y1] = [pt1[0], pt1[1]]
    [x2, y2] = [pt2[0], pt2[1]]
    [x3, y3] = [pt3[0], pt3[1]]
    ori = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
    # collinear
    if (abs(ori) <= error):
        return None
    # clockwise 
    elif (ori < 0):        
        return True
    # counter-clockwise
    else:        
        return False

def isPtOnLine(pt: pt, line: line, error: float = CONST_EPSILON) -> bool:
    """
    Is a pt on the line?

    Parameters
    ----------
    pt: pt, required
        Coordinate of the point
    line: line, required
        Two coordinates to form a line
    error: float, optional, default as CONST_EPSILON
        Error tolerance

    Return
    ------
    bool
        True if the point is on the line, False else-wise

    """
    if (is2PtsSame(line[0], line[1], error = error)):
        raise ZeroVectorError()
    if (is3PtsClockWise(pt, line[0], line[1], error = error) == None):
        return True
    else:
        return False

def isPtOnSeg(pt: pt, seg: line, interiorOnly: bool=False, error:float = CONST_EPSILON) -> bool:
    """
    Is a pt on the lines segment?

    Parameters
    ----------
    pt: pt, required
        Coordinate of the point
    seg: line, required
        Two coordinates to form a line segment
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
    error: float, optional, default as CONST_EPSILON
        Error tolerance

    Return
    ------
    bool
        True if the point is on the line segment, False else-wise

    """
    onLine = isPtOnLine(pt, seg, error = error)
    if (onLine == False):
        return False
    # Get pts =================================================================
    [x1, y1] = [seg[0][0], seg[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [seg[1][0], seg[1][1]]
    # NOTE: 判断在线段上的标准为到两端点的距离的和等于端点间距离
    onSeg = (
        abs(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) 
        + math.sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2) 
        - math.sqrt((x1 - x3) ** 2 + (y1 - y3) ** 2)) <= error)
    # Check if the intersection is in the interior ============================
    if (interiorOnly):
        return onSeg and not is2PtsSame(pt, seg[0], error = error) and not is2PtsSame(pt, seg[1], error = error)
    else:
        return onSeg

def isPtOnRay(pt: pt, ray: line, interiorOnly: bool=False, error = CONST_EPSILON) -> bool:
    """
    Is a pt on the ray?

    Parameters
    ----------
    pt: pt, required
        Coordinate of the point
    ray: line, required
        Two coordinates to form a line segment
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
    error: float, optional, default as CONST_EPSILON
        Error tolerance

    Return
    ------
    bool
        True if the point is on the line segment, False else-wise

    """
    onLine = isPtOnLine(pt, ray, error = error)
    if (onLine == False):
        return False
    # Get pts =================================================================
    [x1, y1] = [ray[0][0], ray[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [ray[1][0], ray[1][1]]
    onRay = (
        (abs(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) 
            + math.sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2) 
            - math.sqrt((x1 - x3) ** 2 + (y1 - y3) ** 2)) <= error)
        or (math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) >= math.sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2)))
    # Check if the intertion is in the interior ===============================
    if (interiorOnly):
        return onRay and not is2PtsSame(pt, ray[0], error = error)
    else:
        return onRay

def isPtOnPolyEdge(pt: pt, poly: poly, error: float = CONST_EPSILON) -> bool:
    """
    Is a pt on the edge of the polygon?

    Parameters
    ----------
    pt: pt, required
        Coordinate of the point
    poly: poly, required
        The polygon
    error: float, optional, default as CONST_EPSILON
        Error tolerance

    Return
    ------
    bool
        True if the point is on the line segment, False else-wise

    """

    # Check if the pt is on any of the edge segment ===========================
    for i in range(-1, len(poly) - 1):
        if (isPtOnSeg(pt, [poly[i], poly[i + 1]], error)):
            return True
    return False

def isPtInPoly(pt: pt, poly: poly, interiorOnly: bool=False, error: float = CONST_EPSILON) -> bool:
    """
    Is a pt in the polygon?

    Parameters
    ----------
    pt: pt, required
        Coordinate of the point
    poly: poly, required
        The polygon
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
    error: float, optional, default as CONST_EPSILON
        Error tolerance

    Return
    ------
    bool
        True if the point is on the line segment, False else-wise

    """

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
        return inPoly and not isPtOnPolyEdge(pt, poly, error)
    else:
        return inPoly

# Object to line-shape ========================================================
def rayPerp2Line(pt: pt, line: line) -> line:
    """
    Given a point and a line, return a ray from that point and perpendicular to the given line

    Parameters
    ----------
    pt: pt, required
        The point where the ray will go through
    line: line Required
        The line which the ray is perpendicular to

    Return
    ------
    line
        A ray that is perpendicular to the line

    """

    if (isPtOnLine(pt, line)):
        raise ZeroVectorError()

    # Heading of line
    heading = headingXY(line[0], line[1])
    heading += 90
    while (heading > 360):
        heading -= 360
    while (heading < 0):
        heading += 360

    # Perp line    
    perp = [pt, ptInDistXY(pt, heading, 10)]
    intPt = intLine2Line(line, perp)
    ray = [pt, intPt['intersect']]

    return ray

# Line-shape intersection =====================================================
def intLine2Line(line1: line, line2: line) -> dict:
    """
    The intersection of a line to another line

    Parameters
    ----------
    line1: line, required
        The first line
    line2: line, required
        The second line

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """
    
    if (is2PtsSame(line1[0], line1[1])):
        raise ZeroVectorError(line1)
    if (is2PtsSame(line2[0], line2[1])):
        raise ZeroVectorError(line2)

    # Get Ax + By + C = 0
    def abc(pt1, pt2):
        x1, y1 = pt1
        x2, y2 = pt2
        a = y1 - y2
        b = x2 - x1
        c = x1 * y2 - x2 * y1
        return a, b, c

    # Calculate intersection
    a1, b1, c1 = abc(line1[0], line1[1])
    a2, b2, c2 = abc(line2[0], line2[1])
    D = a1 * b2 - a2 * b1

    # Check if parallel
    # 共线情形
    if (D == 0 and is3PtsClockWise(line1[0], line1[1], line2[0]) == None):
        return {
            'status': 'Collinear',
            'intersect': line1,
            'intersectType': 'Line',
            'interiorFlag': True
        }
    # 平行情形
    elif (D == 0):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 相交情形
    else:
        x = (b1 * c2 - b2 * c1) / D
        y = (a2 * c1 - a1 * c2) / D
        return {
            'status': 'Cross',
            'intersect': (x, y),
            'intersectType': 'Point',
            'interiorFlag': True
        }

def intLine2Seg(line: line, seg: line) -> dict:
    """
    The intersection of a line to another line segment

    Parameters
    ----------
    line: line, required
        The first line
    seg: line, required
        The second line segment

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    intPt = intLine2Line(line, seg)

    # 如果直线不相交，返回不相交
    if (intPt['status'] == 'NoCross'):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': True
        }
    # 如果直线和线段共线，返回线段
    elif (intPt['status'] == 'Collinear'):
        return {
            'status': 'Collinear',
            'intersect': seg,
            'intersectType': 'Segment',
            'interiorFlag': True
        }
    # 若相交点不在线段上，返回不相交
    elif (not isPtOnSeg(intPt['intersect'], seg, interiorOnly=False)):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若相交点在线段端点上，返回相交且不在interior
    elif (is2PtsSame(intPt['intersect'], seg[0]) or is2PtsSame(intPt['intersect'], seg[1])):
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': False
        }
    # 若相交点在线段内，返回相交且在interior内
    else:
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': True
        }

def intLine2Ray(line: line, ray: line) -> dict:
    """
    The intersection of a line to a ray

    Parameters
    ----------
    line: line, required
        The first line
    ray: line, required
        The second ray

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    intPt = intLine2Line(line, ray)

    # 若直线不相交，返回不相交
    if (intPt['status'] == 'NoCross'):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若直线与射线共线，返回射线
    elif (intPt == "Collinear"):
        return {
            'status': 'Collinear',
            'intersect': ray,
            'intersectType': 'Ray',
            'interiorFlag': True
        }
    # 若相交点不在射线内，返回不相交
    elif (not isPtOnRay(intPt['intersect'], ray, interiorOnly=False)):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若相交点在射线端点上，返回相交且不在interior
    elif (is2PtsSame(intPt['intersect'], ray[0])):
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': False
        }
    # 若相交点在射线内，返回相交且在interior内
    else:
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': True
        }        

def intSeg2Line(seg: line, line: line) -> dict:
    """
    The intersection of a line segment to another line

    Parameters
    ----------
    seg: line, required
        The first line segment
    line: line, required
        The second line

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """
    return intLine2Seg(line, seg)

def intSeg2Seg(seg1: line, seg2: line) -> dict:
    """
    The intersection of a line segment to another line segment

    Parameters
    ----------
    seg1: line, required
        The first line segment
    seg2: line, required
        The second line segment

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    intPt = intLine2Line(seg1, seg2)

    # 若直线不相交，返回不相交
    if (intPt['status'] == 'NoCross'):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若直线共线，进一步分情况讨论
    elif (intPt['status'] == "Collinear"):
        # 线段端点是否在另一条线段上
        seg1aInSeg2Flag = isPtOnSeg(seg1[0], seg2, interiorOnly=False)
        seg1bInSeg2Flag = isPtOnSeg(seg1[1], seg2, interiorOnly=False)
        seg2aInSeg1Flag = isPtOnSeg(seg2[0], seg1, interiorOnly=False)
        seg2bInSeg1Flag = isPtOnSeg(seg2[1], seg1, interiorOnly=False)

        # Case 1: seg1完全在seg2内
        if (seg1aInSeg2Flag and seg1bInSeg2Flag):
            return {
                'status': 'Collinear',
                'intersect': seg1,
                'intersectType': 'Segment',
                'interiorFlag': True
            }
        # Case 2: seg2完全在seg1内
        elif (seg2aInSeg1Flag and seg2bInSeg1Flag):
            return {
                'status': 'Collinear',
                'intersect': seg2,
                'intersectType': 'Segment',
                'interiorFlag': True
            }
        # Case 3: seg1b在seg2内，seg2a在seg1内
        elif (seg1bInSeg2Flag and seg2aInSeg1Flag):
            # Case 3.1: 相交在端点上
            if (is2PtsSame(seg1[1], seg2[0])):
                return {
                    'status': 'Collinear',
                    'intersect': seg1[1],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            # Case 3.2: 有重合线段
            else:
                return {
                    'status': 'Collinear',
                    'intersect': [seg1[1], seg2[0]],
                    'intersectType': 'Segment',
                    'interiorFlag': True
                }
        # Case 4: seg1a在seg2内，seg2a在seg1内
        elif (seg1aInSeg2Flag and seg2aInSeg1Flag):
            # Case 4.1: 相交在端点上
            if (is2PtsSame(seg1[0], seg2[0])):
                return {
                    'status': 'Collinear',
                    'intersect': seg1[0],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            # Case 4.2: 有重合线段
            else:
                return {
                    'status': 'Collinear',
                    'intersect': [seg1[0], seg2[0]],
                    'intersectType': 'Segment',
                    'interiorFlag': True
                }
        # Case 5: seg1b在seg2内，seg2b在seg1内
        elif (seg1bInSeg2Flag and seg2bInSeg1Flag):
            # Case 5.1: 相交在端点上
            if (is2PtsSame(seg1[1], seg2[1])):
                return {
                    'status': 'Collinear',
                    'intersect': seg1[1],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            # Case 5.2: 有重合线段
            else:
                return {
                    'status': 'Collinear',
                    'intersect': [seg1[1], seg2[1]],
                    'intersectType': 'Segment',
                    'interiorFlag': True
                }
        # Case 6: seg1a在seg2内，seg2b在seg1内
        elif (seg1aInSeg2Flag and seg2bInSeg1Flag):
            # Case 6.1: 相交在端点上
            if (is2PtsSame(seg1[0], seg2[1])):
                return {
                    'status': 'Collinear',
                    'intersect': seg1[0],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            # Case 6.2: 有重合线段
            else:
                return {
                    'status': 'Collinear',
                    'intersect': [seg1[0], seg2[1]],
                    'intersectType': 'Segment',
                    'interiorFlag': True
                }
        # Case 7: 无相交
        else:
            return {
                'status': 'Collinear',
                'intersect': None,
                'intersectType': None,
                'interiorFlag': None
            }
    # 若交点不在任何一个线段上，不相交
    elif (not isPtOnSeg(intPt['intersect'], seg1, interiorOnly=False)
            or not isPtOnSeg(intPt['intersect'], seg2, interiorOnly=False)):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若交点在任意一个端点上，交于boundary
    if (is2PtsSame(intPt['intersect'], seg1[0])
            or is2PtsSame(intPt['intersect'], seg1[1])
            or is2PtsSame(intPt['intersect'], seg2[0])
            or is2PtsSame(intPt['intersect'], seg2[1])):
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': False
        }
    # 若交点在两个线段内，交于interior
    else:
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': True
        }

def intSeg2Ray(seg: line, ray: line) -> dict:
    """
    The intersection of a line segment to a ray

    Parameters
    ----------
    line: line, required
        The first line
    ray: line, required
        The second line

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    intPt = intLine2Line(seg, ray)

    # 若直线不相交，返回不相交
    if (intPt['status'] == 'NoCross'):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若直线共线，进一步分情况讨论
    elif (intPt['status'] == "Collinear"):
        # 线段端点是否在射线内
        segAInRayFlag = isPtOnRay(seg[0], ray, interiorOnly=False)
        segBInRayFlag = isPtOnRay(seg[1], ray, interiorOnly=False)
        # Case 1: seg整个在ray内
        if (segAInRayFlag and segBInRayFlag):
            return {
                'status': 'Collinear',
                'intersect': seg,
                'intersectType': 'Segment',
                'interiorFlag': True
            }
        # Case 2: segA在ray内，segB不在ray内
        elif (segAInRayFlag and not segBInRayFlag):
            # Case 2.1: 相交于射线的端点
            if (is2PtsSame(seg[1], ray[0])):
                return {
                    'status': 'Collinear',
                    'intersect': seg[1],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            # Case 2.2: 相较于射线内部
            else:
                return {
                    'status': 'Collinear',
                    'intersect': [ray[0], seg[1]],
                    'intersectType': 'Segment',
                    'interiorFlag': True
                }
        # Case 3: segB在ray内，segA不在ray内
        elif (segBInRayFlag and not segAInRayFlag):
            # Case 3.1: 相交于射线的端点
            if (is2PtsSame(seg[0], ray[0])):
                return {
                    'status': 'Collinear',
                    'intersect': seg[0],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            # Case 3.2: 相较于射线内部
            else:
                return {
                    'status': 'Collinear',
                    'intersect': [ray[0], seg[0]],
                    'intersectType': 'Segment',
                    'interiorFlag': True
                }
        # Case 4: 无相交
        else:
            return {
                'status': 'Collinear',
                'intersect': None,
                'intersectType': None,
                'interiorFlag': None
            }
    # 若交点不在线段内或不在射线内，不相交
    elif (not isPtOnSeg(intPt['intersect'], seg, interiorOnly=False)
            or not isPtOnRay(intPt['intersect'], ray, interiorOnly=False)):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若交点在线段端点上或者在射线端点上，交于boundary
    elif (is2PtsSame(intPt['intersect'], seg[0])
            or is2PtsSame(intPt['intersect'], seg[1])
            or is2PtsSame(intPt['intersect'], ray[0])):
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': False
        }
    # 若交点在线段内且在射线内，交于interior
    else:
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': True
        }

def intRay2Line(ray: line, line: line) -> dict:
    """
    The intersection of a ray to another line

    Parameters
    ----------
    ray: line, required
        The first line
    line: line, required
        The second line

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """
    return intLine2Ray(line, ray)

def intRay2Seg(ray: line, seg: line) -> dict:
    """
    The intersection of a ray to another line segment

    Parameters
    ----------
    ray: line, required
        The first line
    seg: line, required
        The second line

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """
    return intSeg2Ray(seg, ray)

def intRay2Ray(ray1: line, ray2: line) -> dict:
    """
    The intersection of a ray to another ray

    Parameters
    ----------
    ray1: line, required
        The first line
    ray2: line, required
        The second line

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    intPt = intLine2Line(ray1, ray2)

    # 若直线不相交，返回不相交
    if (intPt['status'] == 'NoCross'):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若直线共线，进一步情况讨论
    elif (intPt['status'] == "Collinear"):
        # 射线端点是否在其他射线内
        ray1AInRay2Flag = isPtOnRay(ray1[0], ray2, interiorOnly=False)
        ray2AInRay1Flag = isPtOnRay(ray2[0], ray1, interiorOnly=False)
        # Case 1: ray2a在ray1内，且ray1a在ray2内
        if (ray1AInRay2Flag and ray2AInRay1Flag):
            # Case 1.1: 射线对头相交于端点
            if (is2PtsSame(ray1[0], ray2[0])):
                return {
                    'status': 'Collinear',
                    'intersect': ray1[0],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            # Case 1.2: 射线有重合线段
            else:
                return {
                    'status': 'Collinear',
                    'intersect': [ray1[0], ray2[0]],
                    'intersectType': 'Segment',
                    'interiorFlag': True
                }
        # Case 2: ray2a在ray1内，且ray1a不在ray2内
        elif (ray2AInRay1Flag and not ray1AInRay2Flag):
            return {
                'status': 'Collinear',
                'intersect': ray2,
                'intersectType': 'Ray',
                'interiorFlag': True
            }
        # Case 3: ray1a在ray2内，且ray2a不在ray1内
        elif (ray1AInRay2Flag and not ray2AInRay1Flag):
            return {
                'status': 'Collinear',
                'intersect': ray1,
                'intersectType': 'Ray',
                'interiorFlag': True
            }
        # Case 4: 无相交
        else:
            return {
                'status': 'Collinear',
                'intersect': None,
                'intersectType': None,
                'interiorFlag': None
            }
    # 若直线相交点不在任何一个射线上，不相交
    elif (not isPtOnRay(intPt['intersect'], ray1, interiorOnly=False)
            or not isPtOnRay(intPt['intersect'], ray2, interiorOnly=False)):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若相交点在射线的端点上，交于boundary
    elif (is2PtsSame(intPt['intersect'], ray1[0])
            or is2PtsSame(intPt['intersect'], ray2[0])):
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': False
        }
    # 若相交点在射线内，交于interior
    else:
        return {
            'status': 'Cross',
            'intersect': intPt['intersect'],
            'intersectType': 'Point',
            'interiorFlag': True
        }

# Line-shape versus Line-shape ===================================================
def isLineIntLine(line1: line, line2: line) -> bool:
    """
    Is two line intersect with each other?

    Parameters
    ----------
    line1: line, required
        The first line
    line2: line, required
        The second line

    Return
    ------
    bool
        True if intersects
    """
    intPt = intLine2Line(line1, line2)
    if (intPt['intersect'] == None):
        return False
    else:
        return True

def isLineIntSeg(line: line, seg: line, interiorOnly:bool=False) -> bool:
    """
    Is a line intersect with a line segment?

    Parameters
    ----------
    line: line, required
        The first line
    seg: line, required
        The second line segment
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior

    Return
    ------
    bool
        True if intersects
    """
    intPt = intLine2Seg(line, seg)
    # 无相交点
    if (intPt['intersect'] == None):
        return False
    # 需要交于interior但未能交于interiror
    elif (interiorOnly and not intPt['interiorFlag']):
        return False
    # 不需要交于interior或需要且交于interior
    else:
        return True

def isLineIntRay(line: line, ray: line, interiorOnly: bool=False) -> bool:
    """
    Is a line intersect with a ray?

    Parameters
    ----------
    line: line, required
        The first line
    ray: line, required
        The second ray
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior

    Return
    ------
    bool
        True if intersects
    """    
    intPt = intLine2Ray(line, ray)
    # 无相交点
    if (intPt['intersect'] == None):
        return False
    # 需要交于interior但未能交于interiror
    elif (interiorOnly and not intPt['interiorFlag']):
        return False
    # 不需要交于interior或需要且交于interior
    else:
        return True

def isSegIntLine(seg: line, line: line, interiorOnly: bool=False) -> bool:
    """
    Is a line segment intersect with a line?

    Parameters
    ----------
    seg: line, required
        The first line segment
    line: line, required
        The second line
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
        
    Return
    ------
    bool
        True if intersects
    """   
    return isLineIntSeg(line, seg, interiorOnly)

def isSegIntSeg(seg1: line, seg2: line, interiorOnly: bool=False) -> bool:
    """
    Is a line segment intersect with another line segment?

    Parameters
    ----------
    seg1: line, required
        The first line segment
    seg2: line, required
        The second line
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
        
    Return
    ------
    bool
        True if intersects
    """ 
    intPt = intSeg2Seg(seg1, seg2)
    # 无相交点
    if (intPt['intersect'] == None):
        return False
    # 需要交于interior但未能交于interiror
    elif (interiorOnly and not intPt['interiorFlag']):
        return False
    # 不需要交于interior或需要且交于interior
    else:
        return True

def isSegIntRay(seg: line, ray: line, interiorOnly: bool=False) -> bool:
    """
    Is a line segment intersect with a ray?

    Parameters
    ----------
    seg: line, required
        The first line segment
    ray: line, required
        The second line
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
        
    Return
    ------
    bool
        True if intersects
    """ 
    intPt = intSeg2Ray(seg, ray)
    # 无相交点
    if (intPt['intersect'] == None):
        return False
    # 需要交于interior但未能交于interiror
    elif (interiorOnly and not intPt['interiorFlag']):
        return False
    # 不需要交于interior或需要且交于interior
    else:
        return True

def isSegIntBoundingbox(seg: line, boundingBox: list) -> bool:
    """
    Is a line segment intersect with a given bounding box?

    Parameters
    ----------
    seg: line, required
        The first line segment
    boundingBox: list, required
        A list, in the following format: [minX, maxX, minY, maxY]
        
    Return
    ------
    bool
        True if intersects
    """ 
    if (boundingBox[0] <= seg[0][0] <= boundingBox[2] and boundingBox[1] <= seg[0][1] <= boundingBox[3]):
        return True
    if (boundingBox[0] <= seg[1][0] <= boundingBox[2] and boundingBox[1] <= seg[1][1] <= boundingBox[3]):
        return True

    # If both end is in the same side of bounding box
    if (seg[0][0] <= boundingBox[0] and seg[1][0] <= boundingBox[0]):
        return False
    if (seg[0][0] >= boundingBox[2] and seg[1][0] >= boundingBox[2]):
        return False
    if (seg[0][1] <= boundingBox[1] and seg[1][1] <= boundingBox[1]):
        return False
    if (seg[0][1] >= boundingBox[3] and seg[1][1] >= boundingBox[3]):
        return False

    # Clockwise check
    c1 = is3PtsClockWise(seg[0], seg[1], [boundingBox[0], boundingBox[1]])
    c2 = is3PtsClockWise(seg[0], seg[1], [boundingBox[2], boundingBox[1]])
    if (c1 != c2):
        return True
    c3 = is3PtsClockWise(seg[0], seg[1], [boundingBox[0], boundingBox[3]])
    if (c1 != c3 or c2 != c3):
        return True
    c4 = is3PtsClockWise(seg[0], seg[1], [boundingBox[2], boundingBox[3]])
    if (c1 != c4 or c2 != c4 or c3 != c4):
        return True

    return False

def isRayIntLine(ray: line, line: line, interiorOnly: bool=False) -> bool:
    """
    Is a ray intersect with a line?

    Parameters
    ----------
    ray: line, required
        The first line segment
    line: line, required
        The second line
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
        
    Return
    ------
    bool
        True if intersects
    """ 
    return isLineIntRay(line, ray, interiorOnly)

def isRayIntSeg(ray: line, seg: line, interiorOnly: bool=False) -> bool:
    """
    Is a ray intersect with a line segment?

    Parameters
    ----------
    ray: line, required
        The first line segment
    seg: line, required
        The second line
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
        
    Return
    ------
    bool
        True if intersects
    """ 
    return isSegIntRay(seg, ray, interiorOnly)

def isRayIntRay(ray1: line, ray2: line, interiorOnly: bool=False) -> bool:
    """
    Is a ray intersect with a ray?

    Parameters
    ----------
    ray1: line, required
        The first line segment
    ray2: line, required
        The second line
    interiorOnly: bool, optional, default as False
        True if only consider intersecting in the interior
        
    Return
    ------
    bool
        True if intersects
    """ 
    intPt = intRay2Ray(ray1, ray2)
    # 无相交点
    if (intPt['intersect'] == None):
        return False
    # 需要交于interior但未能交于interiror
    elif (interiorOnly and not intPt['interiorFlag']):
        return False
    # 不需要交于interior或需要且交于interior
    else:
        return True

# Line-shape intersect with polygon ===========================================
def intLine2Poly(line: line, poly: poly=None, polyShapely: shapely.Polygon=None, returnShaplelyObj: bool=False) -> dict | list[dict] | shapely.Point | shapely.Polygon | shapely.GeometryCollection:
    """
    The intersection of a line to a polygon

    Parameters
    ----------
    line: line, required
        The first line
    poly: poly, optional, default as None
        The second polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object        

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    # Sanity check
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: `poly` and `polyShapely` cannot be None at the same time.")

    # Projct points to the line
    projPts = []    
    for pt in poly:
        projPts.append(ptFoot2Line(pt, line))

    # Find two pts on the line that are the farthest away from each other =====
    projPts.sort()
    seg = [projPts[0], projPts[-1]]

    return intSeg2Poly(seg, poly, polyShapely, returnShaplelyObj)

def intSeg2Poly(seg: line, poly: poly=None, polyShapely: shapely.Polygon=None, returnShaplelyObj: bool=False) -> dict | list[dict] | shapely.Point | shapely.Polygon | shapely.GeometryCollection:
    """
    The intersection of a line segment to a polygon

    Parameters
    ----------
    seg: line, required
        The first line segment
    poly: poly, optional, default as None
        The second polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object        

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    # WARNING: Results may not be reliable on both ends of the seg.
    # NOTE: 20231103 暂时用了一个stupid way来处理了

    # Sanity check ============================================================
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: `poly` and `polyShapely` cannot be None at the same time.")

    # get shapely objects =====================================================
    segShapely = shapely.LineString([seg[0], seg[1]])
    if (polyShapely == None):
        polyShapely = shapely.Polygon(poly)
    intShape = shapely.intersection(segShapely, polyShapely)

    # If return shapely objects no processing needed ==========================
    if (returnShaplelyObj):
        return intShape

    # 若不相交，返回不相交
    # FIXME: 现在的精度可能有问题，需要计算两者间距离
    if (intShape.is_empty):
        dist = shapely.distance(segShapely, polyShapely)
        if (dist <= CONST_EPSILON):
            # Case 1: 若两个端点足够近
            end1Dist = distPt2Poly(seg[0], polyShapely = polyShapely)
            if (end1Dist <= CONST_EPSILON):
                return {
                    'status': 'Cross',
                    'intersect': seg[0],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            end2Dist = distPt2Poly(seg[1], polyShapely = polyShapely)
            if (end2Dist <= CONST_EPSILON):
                return {
                    'status': 'Cross',
                    'intersect': seg[1],
                    'intersectType': 'Point',
                    'interiorFlag': False
                }
            # Case 2: 相切的情形
            for pt in poly:
                ptDist = distPt2Seg(pt, seg)
                if (ptDist <= CONST_EPSILON):
                    return {
                        'status': 'Cross',
                        'intersect': pt,
                        'intersectType': 'Point',
                        'interiorFlag': False
                    }
        else: 
            return {
                'status': 'NoCross',
                'intersect': None,
                'intersectType': None,
                'interiorFlag': None
            }
    elif (isinstance(intShape, shapely.Point)):
        return {
            'status': 'Cross',
            'intersect': (intShape.x, intShape.y),
            'intersectType': 'Point',
            'interiorFlag': False
        }
    elif (isinstance(intShape, shapely.LineString)):
        seg = [tuple(intShape.coords[0]), tuple(intShape.coords[1])]
        midPt = (seg[0][0] + (seg[1][0] - seg[0][0]) / 2,
                 seg[0][1] + (seg[1][1] - seg[0][1]) / 2)
        interiorFlag = shapely.contains(polyShapely, shapely.Point(midPt))
        # NOTE: 由于精度的问题，实际上是Point的情况可能会返回Segment
        if (distEuclideanXY(seg[0], seg[1])['dist'] <= CONST_EPSILON):
            return {
                'status': 'Cross',
                'intersect': seg[0],
                'intersectType': 'Point',
                'interiorFlag': interiorFlag
            }
        else:
            return {
                'status': 'Cross',
                'intersect': seg,
                'intersectType': 'Segment',
                'interiorFlag': interiorFlag
            }
    else:
        intSp = []
        for obj in intShape.geoms:
            if (isinstance(obj, shapely.Point)):
                intSp.append({
                    'status': 'Cross',
                    'intersect': (obj.x, obj.y),
                    'intersectType': 'Point',
                    'interiorFlag': False
                })
            elif (isinstance(obj, shapely.LineString)):
                seg = [tuple(obj.coords[0]), tuple(obj.coords[1])]
                midPt = (seg[0][0] + (seg[1][0] - seg[0][0]) / 2,
                         seg[0][1] + (seg[1][1] - seg[0][1]) / 2)
                interiorFlag = shapely.contains(polyShapely, shapely.Point(midPt))
                # NOTE: 由于精度的问题，实际上是Point的情况可能会返回Segment
                if (distEuclideanXY(seg[0], seg[1])['dist'] <= CONST_EPSILON):
                    intSp.append({
                        'status': 'Cross',
                        'intersect': seg[0],
                        'intersectType': 'Point',
                        'interiorFlag': interiorFlag
                    })
                else:
                    intSp.append({
                        'status': 'Cross',
                        'intersect': seg,
                        'intersectType': 'Segment',
                        'interiorFlag': interiorFlag
                    })
        return intSp

def intRay2Poly(ray: line, poly: poly=None, polyShapely: shapely.Polygon=None, returnShaplelyObj: bool=False) -> dict | list[dict] | shapely.Point | shapely.Polygon | shapely.GeometryCollection:
    """
    The intersection of a ray to a polygon

    Parameters
    ----------
    ray: line, required
        The first ray
    poly: poly, optional, default as None
        The second polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object        

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    # Sanity check
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: `poly` and `polyShapely` cannot be None at the same time.")

    # Project points to the line
    projPts = []    
    for pt in poly:
        projPts.append(ptFoot2Line(pt, ray))

    projPts.sort()
    minPt = projPts[0]
    maxPt = projPts[-1]
    
    isMinPtOnRay = isPtOnRay(minPt, ray)
    isMaxPtOnRay = isPtOnRay(maxPt, ray)
    # 若minPt和maxPt均不能投影到射线上，肯定不相交
    if (not isMinPtOnRay and not isMaxPtOnRay):
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    # 若两个都在射线上，射线可能的与之相交部分在两点间
    elif (isMinPtOnRay and isMaxPtOnRay):
        seg = [minPt, maxPt]
        return intSeg2Poly(seg, poly, polyShapely, returnShaplelyObj)
    # 若minPt在射线上，maxPt不在
    elif (isMinPtOnRay):
        seg = [ray[0], minPt]
        return intSeg2Poly(seg, poly, polyShapely, returnShaplelyObj)
    # 若maxPt在射线上，minPt不在
    else:
        seg = [ray[0], maxPt]
        return intSeg2Poly(seg, poly, polyShapely, returnShaplelyObj)

def intSeq2Poly(seq: list[pt], poly: poly, seqShapely: shapely.LineString=None, polyShapely: shapely.Polygon=None, returnShaplelyObj: bool=False):
    """
    The intersection of a sequence to a polygon

    Parameters
    ----------
    seq: list of pt, required
        The first sequence of points
    poly: poly, optional, default as None
        The second polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object        

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    # Sanity check
    if (seq == None and seqShapely == None):
        raise MissingParameterError("ERROR: `seq` and `seqShapely` cannot be None at the same time.")
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: `poly` and `polyShapely` cannot be None at the same time.")

    # get shapely objects
    if (seqShapely == None):
        seqShapely = shapely.LineString(seq)
    if (polyShapely == None):
        polyShapely = shapely.Polygon(poly)
    intShape = shapely.intersection(seqShapely, polyShapely)

    # If return shapely objects no processing needed
    if (returnShaplelyObj):
        return intShape 
    # 若不相交，返回不相交
    if (intShape.is_empty):
        # FIXME: 误差的部分之后加，现在只考虑用在clipRoadNetworkByPoly()里
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }
    elif (isinstance(intShape, shapely.Point)):
        return {
            'status': 'Cross',
            'intersect': (intShape.x, intShape.y),
            'intersectType': 'Point',
            'interiorFlag': False
        }
    elif (isinstance(intShape, shapely.LineString)):
        intSeq = [list(pt) for pt in intShape.coords]
        return {
            'status': 'Cross',
            'intersect': intSeq,
            'intersectType': 'Segment',
            'interiorFlag': True
        }
    else:
        intSp = []
        for obj in intShape.geoms:
            if (isinstance(obj, shapely.Point)):
                intSp.append({
                    'status': 'Cross',
                    'intersect': (obj.x, obj.y),
                    'intersectType': 'Point',
                    'interiorFlag': False
                })
            elif (isinstance(obj, shapely.LineString)):
                intSeq = [list(pt) for pt in obj.coords]
                intSp.append({
                    'status': 'Cross',
                    'intersect': intSeq,
                    'intersectType': 'Segment',
                    'interiorFlag': True
                })
        return intSp

# Line-shape versus polygon ===================================================
def isLineIntPoly(line: line, poly: poly=None, polyShapely: shapely.Polygon=None, interiorOnly: bool=False) -> bool:
    """
    Is a line intersects to a polygon?

    Parameters
    ----------
    line: line, required
        The first line
    poly: poly, optional, default as None
        The second polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object        

    Return
    ------
    bool
        True if intersects
    """
    intSp = intLine2Poly(line, poly, polyShapely)
    # 若只输出了一个字典，按字典判断
    if (isinstance(intSp, dict)):
        return (intSp['status'] == 'Cross' 
            and not (interiorOnly and not intSp['interiorFlag']))
    elif (isinstance(intSp, list)):
        for intPt in intSp:
            trueWhen = (intPt['status'] == 'Cross' 
                and not (interiorOnly and not intPt['interiorFlag']))
            if (trueWhen):
                return True    
        return False

def isSegIntPoly(seg: line, poly: poly=None, polyShapely: shapely.Polygon=None, interiorOnly: bool=False) -> bool:
    """
    Is a line segment intersects to a polygon?

    Parameters
    ----------
    seg: line, required
        The first line
    poly: poly, optional, default as None
        The second polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object        

    Return
    ------
    bool
        True if intersects
    """

    # WARNING: results may not be reliable if `interiorFlag` == False.
    
    """Is a segment intersect with a polygon"""
    intSp = intSeg2Poly(seg, poly, polyShapely)
    # 若只输出了一个字典，按字典判断
    if (isinstance(intSp, dict)):
        return (intSp['status'] == 'Cross' 
            and not (interiorOnly and not intSp['interiorFlag']))
    elif (isinstance(intSp, list)):
        for intPt in intSp:
            trueWhen = (intPt['status'] == 'Cross' 
                and not (interiorOnly and not intPt['interiorFlag']))
            if (trueWhen):
                return True    
        return False

def isRayIntPoly(ray: line, poly: poly=None, polyShapely: shapely.Polygon=None, interiorOnly: bool=False) -> bool:
    """
    Is a ray intersects to a polygon?

    Parameters
    ----------
    ray: line, required
        The first line
    poly: poly, optional, default as None
        The second polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object        

    Return
    ------
    bool
        True if intersects
    """

    intSp = intRay2Poly(ray, poly, polyShapely)
    # 若只输出了一个字典，按字典判断
    if (isinstance(intSp, dict)):
        return (intSp['status'] == 'Cross' 
            and not (interiorOnly and not intSp['interiorFlag']))
    elif (isinstance(intSp, list)):
        for intPt in intSp:
            trueWhen = (intPt['status'] == 'Cross' 
                and not (interiorOnly and not intPt['interiorFlag']))
            if (trueWhen):
                return True    
        return False

# Poly vs poly ================================================================
def intPoly2Poly(poly1: poly=None, poly2: poly=None, poly1Shapely: shapely.Polygon=None, poly2Shapely: shapely.Polygon=None):
    """
    Intersection between two polygons

    Parameters
    ----------
    poly1: poly, optional, default as None
        The first polygon
    poly2: poly, optional, default as None
        The second polygon
    poly1Shapely: shapely.Polygon, optional, default as None
        The correspond shapely object for the first polygon. Need to provide one of the following fields: [`poly1`, `poly1Shapely`]
    poly2Shapely: shapely.Polygon, optional, default as None
        The correspond shapely object for the second polygon. Need to provide one of the following fields: [`poly2`, `poly2Shapely`]

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    if (poly1Shapely == None):
        poly1Shapely = shapely.Polygon([p for p in poly1])
    if (poly2Shapely == None):
        poly2Shapely = shapely.Polygon([p for p in poly2])

    intShape = shapely.intersection(poly1Shapely, poly2Shapely)
    intType = shapely.get_type_id(intShape)

    # Point
    if (intType == 0):
        return {
            'status': 'Cross',
            'intersect': (intShape.x, intShape.y),
            'intersectType': 'Point',
            'interiorFlag': False
        }
    # LineString
    elif (intType == 1):
        return {
            'status': 'Cross',
            'intersect': [intShape.coords[0], intShape.coords[1]],
            'intersectType': 'Segment',
            'interiorFlag': False
        }
    # Polygon
    elif (intType == 3 and not intShape.is_empty):
        return {
            'status': 'Cross',
            'intersect': [i for i in mapping(intShape)['coordinates'][0]],
            'intersectType': 'Polygon',
            'interiorFlag': True
        }
    # MultiPoint/MultiLineString/MultiPolygon/GeometryCollection
    elif (intType in [4, 5, 6, 7]):
        intSp = []
        for g in intShape.geoms:
            intG = shapely.get_type_id(g)
            # Point
            if (intG == 0):
                intSp.append({
                    'status': 'Cross',
                    'intersect': (g.x, g.y),
                    'intersectType': 'Point',
                    'interiorFlag': False
                })
            # LineString
            elif (intG == 1):
                intSp.append({
                    'status': 'Cross',
                    'intersect': [g.coords[0], g.coords[1]],
                    'intersectType': 'Segment',
                    'interiorFlag': False
                })
            # Polygon
            elif (intG == 3):
                intSp.append({
                    'status': 'Cross',
                    'intersect': [i for i in mapping(g)['coordinates'][0]],
                    'intersectType': 'Polygon',
                    'interiorFlag': True
                })
        return intSp
    else:
        return {
            'status': 'NoCross',
            'intersect': None,
            'intersectType': None,
            'interiorFlag': None
        }

def isPolyIntPoly(poly1: poly=None, poly2: poly=None, poly1Shapely: shapely.Polygon=None, poly2Shapely: shapely.Polygon=None, interiorOnly: bool=False) -> bool:
    """
    Is a polygon intersect to another polygon

    Parameters
    ----------
    poly1: poly, optional, default as None
        The first polygon
    poly2: poly, optional, default as None
        The second polygon
    poly1Shapely: shapely.Polygon, optional, default as None
        The correspond shapely object for the first polygon. Need to provide one of the following fields: [`poly1`, `poly1Shapely`]
    poly2Shapely: shapely.Polygon, optional, default as None
        The correspond shapely object for the second polygon. Need to provide one of the following fields: [`poly2`, `poly2Shapely`]

    Return
    ------
    dict
        >>> {
        ...     'status': 'Cross', # Relations between two objects,
        ...     'intersect': pt, # The intersection,
        ...     'intersectType': 'Point', # Type of intersection
        ...     'interiorFlag': False, # True if the intersection is at the boundary
        ... }
    """

    intSp = intPoly2Poly(poly1, poly2, poly1Shapely, poly2Shapely)
    # 若只输出了一个字典，按字典判断
    if (isinstance(intSp, dict)):
        return (intSp['status'] == 'Cross' 
            and not (interiorOnly and not intSp['interiorFlag']))
    elif (isinstance(intSp, list)):
        for intPt in intSp:
            trueWhen = (intPt['status'] == 'Cross' 
                and not (interiorOnly and not intPt['interiorFlag']))
            if (trueWhen):
                return True    
        return False

# Distance from Point to Object ===============================================
def distPt2Line(pt: pt, line: line) -> float:
    """
    The distance between a point and a line

    Parameters
    ----------
    pt: pt, required
        The point
    line: line, required
        The line

    Return
    ------
    float
        The distance between two objects

    """

    area = calTriangleAreaXY(pt, line[0], line[1])
    a = distEuclideanXY(line[0], line[1])['dist']
    h = 2 * area / a
    return h

def distPt2Seg(pt: pt, seg: line) -> float:
    """
    The distance between a point and a line segment

    Parameters
    ----------
    pt: pt, required
        The point
    seg: line, required
        The line segment

    Return
    ------
    float
        The distance between two objects

    """

    foot = ptFoot2Line(pt, seg)
    if (isPtOnSeg(foot, seg)):
        return distEuclideanXY(pt, foot)['dist']
    else:
        return min(distEuclideanXY(pt, seg[0])['dist'], distEuclideanXY(pt, seg[1])['dist'])

def distPt2Ray(pt: pt, ray: line) -> float:
    """
    The distance between a point and a ray

    Parameters
    ----------
    pt: pt, required
        The point
    ray: line, required
        The ray

    Return
    ------
    float
        The distance between two objects

    """

    foot = ptFoot2Line(pt, ray)
    if (isPtOnRay(foot)):
        return distEuclideanXY(pt, foot)
    else:
        return distEuclideanXY(pt, ray[0])

def distPt2Seq(pt: pt, seq: list[pt], closedFlag = False) -> float:
    """
    The distance between a point and a sequence of points

    Parameters
    ----------
    pt: pt, required
        The point
    seq: list of pt, required
        A sequence of points
    closedFlag: bool, optional, default False
        True if the sequence is closed

    Return
    ------
    float
        The distance between two objects

    """

    # FIXME: stupid way, needs improvement
    if (len(seq) == 2):
        return distPt2Seg(pt, seq)

    dist2Seg = []
    for p in seq:
        dist2Seg.append(distEuclideanXY(pt, p)['dist'])
    minIndex = dist2Seg.index(min(dist2Seg))
    if (minIndex == 0):
        if (closedFlag == False):
            return distPt2Seg(pt, [seq[0], seq[1]])
        else:
            return min(distPt2Seg(pt, [seq[0], seq[1]]),
                       distPt2Seg(pt, [seq[0], seq[-1]]))
    elif (minIndex == len(dist2Seg) - 1):
        if (closedFlag == False):
            return distPt2Seg(pt, [seq[-2], seq[-1]])
        else:
            return min(distPt2Seg(pt, [seq[-2], seq[-1]]),
                       distPt2Seg(pt, [seq[0], seq[-1]]))
    else:
        return min(distPt2Seg(pt, [seq[minIndex], seq[minIndex + 1]]),
                   distPt2Seg(pt, [seq[minIndex], seq[minIndex - 1]]))

def distPt2Poly(pt: pt, poly: poly=None, polyShapely: shapely.Polygon=None) -> float:
    """
    The distance between a point and a polygon

    Parameters
    ----------
    pt: pt, required
        The point
    poly: poly, optional, default as None
        The polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]

    Return
    ------
    float
        The distance between two objects

    """

    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: Missing required field `poly` or `polyShapely`")

    if (polyShapely == None):
        polyShapely = shapely.Polygon([p for p in poly])
    ptShapely = shapely.Point(pt)
    return shapely.distance(ptShapely, polyShapely)

def distPoly2Poly(poly1: poly=None, poly2: poly=None, poly1Shapely: shapely.Polygon=None, poly2Shapely: shapely.Polygon=None) -> float:
    """
    The distance between two polys

    Parameters
    ----------
    poly1: poly, optional, default as None
        The first polygon
    poly1Shapely: shapely.Polygon, optional, default as None
        The correspond shapely object for the first polygon. Need to provide one of the following fields: [`poly1`, `poly1Shapely`]
    poly2: poly, optional, default as None
        The second polygon
    poly2Shapely: shapely.Polygon, optional, default as None
        The correspond shapely object for the second polygon. Need to provide one of the following fields: [`poly2`, `poly2Shapely`]

    Return
    ------
    float
        The distance between two objects

    """

    if (poly1 == None and poly1Shapely == None):
        raise MissingParameterError("ERROR: Missing required field `poly1` or `poly1Shapely`")
    if (poly2 == None and poly2Shapely == None):
        raise MissingParameterError("ERROR: Missing required field `poly2` or `poly2Shapely`")

    if (poly1Shapely == None):
        poly1Shapely = shapely.Polygon([p for p in poly1])
    if (poly2Shapely == None):
        poly2Shapely = shapely.Polygon([p for p in poly2])
    return shapely.distance(poly1Shapely, poly2Shapely)

def nearestPtLine2Poly(line: line, poly: poly=None, polyShapely: shapely.Polygon=None) -> dict:
    """
    Find the nearest point between a line and a polygon

    Parameters
    ----------
    line: line, required
        The line
    poly: poly, optional, default as None
        The polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]

    Return
    ------
    dict
        >>> {
        ...     'ptOnLine': pt, # The point from the line,
        ...     'ptOnPoly': pt, # The point from the polygon,
        ... }

    """

    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: Missing required field `poly` or `polyShapely`")

    if (polyShapely == None):
        polyShapely = shapely.Polygon([p for p in poly])
    lineShapely = shapely.LineString(line)
    nearestPts = nearest_points(lineShapely, polyShapely)
    
    return {
        'ptOnLine': (nearestPts[0].x, nearestPts[0].y),
        'ptOnPoly': (nearestPts[1].x, nearestPts[1].y),
    }

# Distance between timed-objects ==============================================
# [Constructing]
def distTimedSeg2TimedSeg(timedSeg1: list[tuple[pt, float]], timedSeg2: list[tuple[pt, float]], startTime: float = 0):
    # 先区分出timedSeg哪段在前，哪段在后
    if (timedSeg1[0][1] < timedSeg2[0][1]):
        (sx1, sy1) = timedSeg1[0][0]
        (ex1, ey1) = timedSeg1[1][0]
        st1 = timedSeg1[0][1]
        et1 = timedSeg1[1][1]
        (sx2, sy2) = timedSeg2[0][0]
        (ex2, ey2) = timedSeg2[1][0]
        st2 = timedSeg2[0][1]
        et2 = timedSeg2[1][1]
    else:
        (sx1, sy1) = timedSeg2[0][0]
        (ex1, ey1) = timedSeg2[1][0]
        st1 = timedSeg2[0][1]
        et1 = timedSeg2[1][1]
        (sx2, sy2) = timedSeg1[0][0]
        (ex2, ey2) = timedSeg1[1][0]
        st2 = timedSeg1[0][1]
        et2 = timedSeg1[1][1]

    # 分三类情况讨论求解
    minDist = None
    if (st1 <= et1 <= st2 <= et2):
        # t \in [st1, et1] 
        # 1 move, 2 stay
        pt2 = (sx2, sy2)
        seg1 = [(sx1, sy1), (ex1, ey1)]
        minDist1 = distPt2Seg(pt2, seg1)

        # t \in [et1, st2]
        # 1 stay, 2 stay
        minDist2 = distEuclideanXY((ex1, ey1), (ex2, ey2))['dist']

        # t \in [st2, et2]
        # 1 stay, 2 move
        pt1 = (ex1, ey1)
        seg2 = [(sx2, sy2), (ex2, ey2)]
        minDist3 = distPt2Seg(pt1, seg2)

        return min([minDist1, minDist2, minDist3])

    elif (st1 <= st2 <= et1 <= et2):
        # t \in [st1, st2]
        # 1 move, 2 stay
        seg1AtS2 = snapInTimedSeq([((sx1, sy1), st1), ((ex1, ey1), et1)], st2)['loc']
        seg1 = [(sx1, ex1), seg1AtS2]
        pt2 = (sx2, sy2)
        minDist1 = distPt2Seg(pt2, seg1)

        # t \in [st2, et1]
        seg2AtE1 = snapInTimedSeq([((sx2, sy2), st2), ((ex2, ey2), et2)], et1)['loc']
        seg1 = [seg1AtS2, (ex1, ey1)]
        seg2 = None
        # 1 move, 2 move

        # t \in [et1, et2]
        # 1 stay, 2 move

        return min([minDist1, minDist2, minDist3])

    elif (st1 <= st2 <= et2 <= et1):
        # t \in [st1, st2]
        # 1 move, 2 stay

        # t \in [st2, et2]
        # 1 move, 2 move

        # t \in [et2, et1]
        # 1 move, 2 stay

        return min([minDist1, minDist2, minDist3])

def timedSeg2Vec(timedSeg):
    dt = timedSeg[1][1] - timedSeg[0][1]
    dx = timedSeg[1][0][0] - timedSeg[0][0][0]
    dy = timedSeg[1][0][1] - timedSeg[0][0][1]
    l = math.sqrt(dx ** 2 + dy ** 2)
    lx = dx * dt / l 
    ly = dy * dt / l
    return timedSeg[0][0], (lx, ly)

def distVec2Vec(pt1, vec1, pt2, vec2, earliest:None|float = None, latest:None|float = None):
    # NOTE: 俩点同时都在动
    # NOTE: 这段代码目前先用gurobi来做，之后要换成O(1)代入公式
    x1, y1 = pt1
    x2, y2 = pt2
    vx1, vy1 = vec1
    vx2, vy2 = vec2

    model = grb.Model("SOCP")
    model.setParam("OutputFlag", 0)

    # Decision variables ======================================================
    d = model.addVar(vtype=grb.GRB.CONTINUOUS, obj=1)
    t = model.addVar(vtype=grb.GRB.CONTINUOUS)
    dx = model.addVar(vtype=grb.GRB.CONTINUOUS, lb = -float('inf'))
    dy = model.addVar(vtype=grb.GRB.CONTINUOUS, lb = -float('inf'))

    # Constraints =============================================================
    model.addConstr(dx == (x1 - x2) + t * (vx1 - vx2))
    model.addConstr(dy == (y1 - y2) + t * (vy1 - vy2))
    model.addConstr(d ** 2 >= dx ** 2 + dy ** 2)
    if (earliest != None):
        model.addConstr(t >= earliest)
    if (latest != None):
        model.addConstr(t <= latest)

    # Optimize ================================================================
    model.modelSense = grb.GRB.MINIMIZE
    model.optimize()
    minDist = model.getObjective().getValue()

    return minDist

def travelVec2Vec(pt1, vec1, pt2, vec2, speed, earliest:None|float = None, latest:None|float = None):
    # NOTE: 俩点同时都在动
    # NOTE: 这段代码目前先用gurobi来做，之后要换成O(1)代入公式
    x1, y1 = pt1
    x2, y2 = pt2
    vx1, vy1 = vec1
    vx2, vy2 = vec2

    model = grb.Model("SOCP")
    model.setParam("OutputFlag", 0)

    # Decision variables ======================================================
    d = model.addVar(vtype=grb.GRB.CONTINUOUS, obj=1)

    t1 = model.addVar(vtype=grb.GRB.CONTINUOUS)
    t2 = model.addVar(vtype=grb.GRB.CONTINUOUS)

    dx = model.addVar(vtype=grb.GRB.CONTINUOUS, lb = -float('inf'))
    dy = model.addVar(vtype=grb.GRB.CONTINUOUS, lb = -float('inf'))

    # Constraints =============================================================
    model.addConstr(dx == (x1 - x2) + t1 * (vx1 - vx2))
    model.addConstr(dy == (y1 - y2) + t2 * (vy1 - vy2))
    model.addConstr(d ** 2 >= dx ** 2 + dy ** 2)
    model.addConstr(t2 == t1 + d * (1 / speed))

    if (earliest != None):
        model.addConstr(t1 >= earliest)
    if (latest != None):
        model.addConstr(t2 <= latest)

    # Optimize ================================================================
    model.modelSense = grb.GRB.MINIMIZE
    model.optimize()

    if (model.status == grb.GRB.status.OPTIMAL):
        minDist = model.getObjective().getValue()
        timedSeg = [
            ((x1 + t1.x * vx1, y1 + t1.x * vy1), t1.x), 
            ((x2 + t2.x * vx2, y2 + t2.x * vy2), t2.x)
        ]
        return {
            'minDist': minDist,
            'minTime': (t2.x - t1.x),
            'timedSeg': timedSeg
        }
    else:
        return None

def earliest2Vec(oriPt, speed, movPt, movVec):
    # NOTE: 这俩同时都在动
    oriPtX, oriPtY = oriPt
    movPtX, movPtY = movPt
    movVX, movVY = movVec

    A = (movVX ** 2 + movVY ** 2 - speed ** 2)
    B = 2 * (movVX * (movPtX - oriPtX) + movVY * (movPtY - oriPtY))
    C = (movPtX - oriPtX) ** 2 + (movPtY - oriPtY) ** 2

    # 求根，如有两个解，均输出
    delta = B ** 2 - 4 * A * C

    if (delta < 0):
        return None
    else:
        t1 = (- B + math.sqrt(delta)) / (2 * A)
        t2 = (- B - math.sqrt(delta)) / (2 * A)

        if (max([t1, t2]) < 0):
            return None
        elif (min([t1, t2]) < 0 and max([t1, t2]) >= 0):
            meetTime = max([t1, t2])
            meetPtX = movPtX + meetTime * movVX
            meetPtY = movPtY + meetTime * movVY
            meetPt = (meetPtX, meetPtY)
            return (meetPt, meetTime)
        else: 
            meetTime = min([t1, t2])
            meetPtX = movPtX + meetTime * movVX
            meetPtY = movPtY + meetTime * movVY
            meetPt = (meetPtX, meetPtY)
            return (meetPt, meetTime)

# Dimension mapping ===========================================================
def vecPolar2XY(vecPolar: pt) -> pt:

    """Given vector's norm and its degree to North, convert it into a 2-tuple vector

    Parameters
    ----------
    vecPolar: tuple[float|int, float|int], required
        A 2-tuple (vVal, vDeg), `vVal` is the norm and `vDeg` is the direction, 0 as North, clockwise, in [0, 360)

    Returns
    -------
    tuple[float|int, float|int]
        The XY vector
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
        A 2-tuple (vX, vY), the coordinate of vector

    Returns
    -------
    tuple[float|int, float|int]
        The polar vector
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
    """
    Given a point in (x, y), map it to a (lat, lon) coordinate using Mercator projection

    Parameters
    ----------
    ptXY: pt, required
        The coordinate in (x, y)

    Return
    ------
    tuple
        The coordinate in (lat, lon)

    """

    # Ref: https://wiki.openstreetmap.org/wiki/Mercator#Python
    (x, y) = ptXY
    lon = math.degrees(y / CONST_EARTH_RADIUS_METERS)
    lat = math.degrees(2 * math.atan(math.exp(x / CONST_EARTH_RADIUS_METERS)) - math.pi / 2.0)
    ptLatLon = (lat, lon)
    return ptLatLon

def ptLatLon2XYMercator(ptLatLon: pt) -> pt:
    """
    Given a point in (lat, lon), map it to a (x, y) coordinate using Mercator projection

    Parameters
    ----------
    ptLatLon: pt, required
        The coordinate in (lat, lon)

    Return
    ------
    tuple
        The coordinate in (x, y)

    """

    # Ref: https://wiki.openstreetmap.org/wiki/Mercator#Python
    (lat, lon) = ptLatLon
    y = math.radians(lon) * CONST_EARTH_RADIUS_METERS
    x = math.log(math.tan(math.pi / 4 + math.radians(lat) / 2)) * CONST_EARTH_RADIUS_METERS
    ptXY = (x, y)
    return ptXY

# Get pt ======================================================================
def ptFoot2Line(pt: pt, line: line) -> pt:
    """
    Given a point and a line, return the foot of that point on the line

    Parameters
    ----------
    pt: pt, required
        The point where the foot will go through
    line: line Required
        The line which the foot is perpendicular to

    Return
    ------
    pt
        The foot of that point on the line

    """
    if (isPtOnLine(pt, line)):
        return tuple(pt)
    else:
        ray = rayPerp2Line(pt, line)
        return ray[1]

def ptInSeqMileage(seq: list[pt], dist: int|float, dimension: str = 'XY') -> pt:
    """Given a list of lat/lon coordinates, and a traveling mileage, returns the coordinate"""
    # Initialize ==============================================================
    inPathFlag = False
    accDist = 0
    preLoc = []
    nextLoc = []
    # Find segment ============================================================
    for i in range(0, len(seq) - 1):
        if (dimension == 'LatLon'):
            accDist += distLatLon(seq[i], seq[i + 1])['dist']
        elif (dimension == 'XY'):
            accDist += distEuclideanXY(seq[i], seq[i + 1])['dist']
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
        segDist = distLatLon(preLoc, nextLoc)['dist']
    elif (dimension == 'XY'):
        segDist = distEuclideanXY(preLoc, nextLoc)['dist']
    if (segDist <= CONST_EPSILON):
        raise ZeroDivisionError
    x = nextLoc[0] + (remainDist / segDist) * (preLoc[0] - nextLoc[0])
    y = nextLoc[1] + (remainDist / segDist) * (preLoc[1] - nextLoc[1])
    return (x, y)

def ptPolyCenter(poly: poly=None, polyShapely: shapely.Polygon=None) -> pt:
    """
    Given a poly, returns the centroid of the poly

    Parameters
    ----------
    poly: poly, optional, default as None
        The polygon
    polyShapely: shapely.Polygon, optional, default as None
        The correspond shapely object for polygon. Need to provide one of the following fields: [`poly`, `polyShapely`]

    Returns
    -------
    pt
        The centroid

    """

    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: Missing required field 'poly' or 'polyShapely'.")
    if (polyShapely == None):
        polyShapely = shapely.Polygon(poly)

    ptShapely = shapely.centroid(polyShapely)
    center = (ptShapely.x, ptShapely.y)
    return center

# Vectors =====================================================================
def rndVec(norm: float = 1):
    deg = random.random() * 360
    vec = ptInDistXY((0, 0), deg, norm)
    return vec

# Polys =======================================================================
def polysUnion(polys:polys=None, polysShapely:list[shapely.Polygon]=None, returnShaplelyObj:bool=False) -> list:
    """
    Given a list of polygons which could be intersecting to each other, return unioned polygons that are not intersecting

    Parameters
    ----------
    polys: poly, optional, default as None
        A list of polygons
    polysShapely: shapely.Polygon, optional, default as None
        The correspond shapely objects for polygons. Need to provide one of the following fields: [`polys`, `polysShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object

    Return
    ------
    list of polys
        A list of polygons
    """
    if (polys == None and polysShapely == None):
        raise MissingParameterError("ERROR: Missing required field 'polys' or 'polysShapely'.")
    if (polysShapely == None):
        polysShapely = []
        for p in polys:
            polysShapely.append(shapely.Polygon(p))
    unionAll = shapely.union_all(polysShapely)
    if (returnShaplelyObj):
        return unionAll

    unionPolys = []
    if (isinstance(unionAll, shapely.geometry.polygon.Polygon)):
        unionPolys = [[[i[0], i[1]] for i in list(unionAll.exterior.coords)]]
    elif (isinstance(unionAll, shapely.geometry.multipolygon.MultiPolygon)):
        for p in unionAll.geoms:
            unionPolys.append([[i[0], i[1]] for i in list(p.exterior.coords)])
    for k in range(len(unionPolys)):
        unionPolys[k] = [unionPolys[k][i] for i in range(len(unionPolys[k])) if distEuclideanXY(unionPolys[k][i], unionPolys[k][i - 1])['dist'] > CONST_EPSILON]
        if (distEuclideanXY(unionPolys[k][0], unionPolys[k][-1])['dist'] <= CONST_EPSILON):
            unionPolys[k] = unionPolys[k][:-1]
    return unionPolys

def polysSubtract(polys:polys=None, polysShapely:list[shapely.Polygon]=None, subPolys:polys=None, subPolysShapely: list[shapely.Polygon]=None, returnShaplelyObj:bool=False) -> list:
    """
    Given a list of polygons, subtract a list of polygon from the first list of polygons

    Parameters
    ----------
    polys: poly, optional, default as None
        A list of polygons
    polysShapely: shapely.Polygon, optional, default as None
        The correspond shapely objects for polygons. Need to provide one of the following fields: [`polys`, `polysShapely`]
    subPolys: poly, optional, default as None
        A list of polygons to substract
    subPolysShapely: shapely.Polygon, optional, default as None
        The correspond shapely objects for polygons. Need to provide one of the following fields: [`subPolys`, `subPolysShapely`]    

    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object

    Return
    ------
    list of polys
        A list of polygons
    """

    if (polys == None and polysShapely == None):
        raise MissingParameterError("ERROR: Missing required field 'polys' or 'polysShapely'.")
    if (subPolys == None and subPolysShapely == None):
        raise MissingParameterError("ERROR: Missing required field 'subPolys' or 'subPolysShapely'.")
    if (polysShapely == None):
        polysShapely = []
        for p in polys:
            polysShapely.append(shapely.Polygon(p))
    unionAll = shapely.union_all(polysShapely)

    if (subPolysShapely == None):
        subPolysShapely = []
        for p in subPolys:
            subPolysShapely.append(shapely.Polygon(p))
    unionSub = shapely.union_all(subPolysShapely)

    diffShapely = shapely.difference(unionAll, unionSub)
    if (returnShaplelyObj):
        return diffShapely

    diffPolys = []
    if (isinstance(diffShapely, shapely.geometry.polygon.Polygon)):
        diffPolys = [[[i[0], i[1]] for i in list(diffShapely.exterior.coords)]]
    elif (isinstance(diffShapely, shapely.geometry.multipolygon.MultiPolygon)):
        for p in diffShapely.geoms:
            diffPolys.append([[i[0], i[1]] for i in list(p.exterior.coords)])
    for k in range(len(diffPolys)):
        diffPolys[k] = [diffPolys[k][i] for i in range(len(diffPolys[k])) if distEuclideanXY(diffPolys[k][i], diffPolys[k][i - 1])['dist'] > CONST_EPSILON]
        if (distEuclideanXY(diffPolys[k][0], diffPolys[k][-1])['dist'] <= CONST_EPSILON):
            diffPolys[k] = diffPolys[k][:-1]
    return diffPolys

def polysIntersect(polys: polys=None, polysShapely:list[shapely.Polygon]=None, returnShaplelyObj:bool=False) -> list:
    """
    Given a list of polygons which could be intersecting to each other, return the intersecting polygons

    Parameters
    ----------
    polys: poly, optional, default as None
        A list of polygons
    polysShapely: shapely.Polygon, optional, default as None
        The correspond shapely objects for polygons. Need to provide one of the following fields: [`polys`, `polysShapely`]
    returnShaplelyObj: bool, optional, default as False
        True if alter the result to be a shapely object

    Return
    ------
    list of polys
        A list of polygons
    """

    if (polys == None and polysShapely == None):
        raise MissingParameterError("ERROR: Missing required field 'polys' or 'polysShapely'.")
    if (polysShapely == None):
        polysShapely = []
        for p in polys:
            polysShapely.append(shapely.Polygon(p))
    intersectionAll = shapely.intersection_all(polysShapely)

    if (returnShaplelyObj):
        return intersectionAll

    intersectionPoly = []
    if (isinstance(intersectionAll, shapely.geometry.polygon.Polygon)):
        intersectionPoly = [[[i[0], i[1]] for i in list(intersectionAll.exterior.coords)]]
    elif (isinstance(intersectionAll, shapely.geometry.multipolygon.MultiPolygon)):
        for p in intersectionAll.geoms:
            intersectionPoly.append([[i[0], i[1]] for i in list(p.exterior.coords)])
    for k in range(len(intersectionPoly)):
        intersectionPoly[k] = [intersectionPoly[k][i] for i in range(len(intersectionPoly[k])) if distEuclideanXY(intersectionPoly[k][i], intersectionPoly[k][i - 1])['dist'] > CONST_EPSILON]
        if (distEuclideanXY(intersectionPoly[k][0], intersectionPoly[k][-1])['dist'] <= CONST_EPSILON):
            intersectionPoly[k] = intersectionPoly[k][:-1]
    return intersectionPoly

def polysSteinerZone(polys: dict, order: int|None = None) -> list[dict]:
    """Given a node dictionary, returns a list of Steiner zones

    Warning
    -------
    This function needs to be rewritten. It's not traceable now.

    Parameters
    ----------

    polys: dictionary, required
        The polys dictionary with neighborhood.
    order: int, optional, default None
        Maximum order of Steiner zone

    Returns
    -------

    list[dict]
        A list of Steiner zone dictionaries, each in the following format::
            >>> SteinerZone = {
            ...     'poly': poly,
            ...     'repPt': centroid,
            ...     'nodeID': []
            ... }

    """

    # List of Steiner Zones
    lstSteinerZone = []
    lstSteinerZoneShape = []

    # Check overlapping
    overlapMatrix = {}
    registeredSZ = []

    # First check by any two pairs of neighbor
    for i in polys:
        for j in polys:
            if (i < j):
                neiI = None
                if ('poly' in polys[i]):
                    neiI = shapely.Polygon([[p[0], p[1]] for p in polys[i]['poly']])
                else:
                    neiI = shapely.Point([polys[i]['loc'][0], polys[i]['loc'][1]])
                neiJ = None
                if ('poly' in polys[j]):
                    neiJ = shapely.Polygon([[p[0], p[1]] for p in polys[j]['poly']])
                else:
                    neiJ = shapely.Point([polys[j]['loc'][0], polys[j]['loc'][1]])

                intersectIJ = shapely.intersection(neiI, neiJ)
                if (not intersectIJ.is_empty):
                    overlapMatrix[i, j] = 1
                    overlapMatrix[j, i] = 1                    
                    lstSteinerZoneShape.append({
                            'polyShape': intersectIJ,
                            'repPtShape': intersectIJ.centroid,
                            'nodeIDs': [i, j]
                        })
                else:
                    overlapMatrix[i, j] = 0
                    overlapMatrix[j, i] = 0

    # If no two neighbor are overlapped, every neighborhood is a Steiner Zone of order 1
    if (sum(overlapMatrix.values()) == 0 or order == 1):
        return [
            {
                'poly': polys[n]['poly'] if 'poly' in polys[n] else [polys[n]['loc']],
                'repPt': list(shapely.Polygon([[p[0], p[1]] for p in polys[n]['poly']]).centroid.coords[0]) if 'poly' in polys[n] else polys[n]['loc'],
                'nodeIDs': [n]
            } for n in polys]

    pointer = 0
    checkOrder = 2
    while (checkOrder <= (order if order != None else len(polys))):
        # Q: How many SZ needs to be checked? A: From `pointer` to `endPointer`
        endPointer = len(lstSteinerZoneShape)

        # Check each SZ in this order, to see if it can be increase by order 1
        for p in range(pointer, endPointer):
            # Get a SZ, see if there is a node n intersect with this SZ
            SZShape = lstSteinerZoneShape[p]
            for n in polys:
                # First, n should not be a SZ member
                if (n not in SZShape['nodeIDs']):
                    # Assume all neighbor in SZShape is intersected with n
                    overlapAllFlag = True
                    # Check one by one, if one of neighbor is not intersected with n, skip
                    for i in SZShape['nodeIDs']:
                        if (overlapMatrix[i, n] == 0):
                            overlapAllFlag = False
                            break
                    if (overlapAllFlag):
                        newNodeIDs = [k for k in SZShape['nodeIDs']]
                        newNodeIDs.append(n)
                        if (list2Tuple(newNodeIDs) not in registeredSZ):
                            # The neighbor of node n could be either a Polygon or a Point
                            neiN = shapely.Polygon([[p[0], p[1]] for p in polys[n]['poly']]) if 'poly' in polys[n] else shapely.Point(polys[n]['loc'])
                            newIntersect = shapely.intersection(SZShape['polyShape'], neiN)
                            if (not newIntersect.is_empty):
                                registeredSZ.append(list2Tuple(newNodeIDs))
                                lstSteinerZoneShape.append({
                                        'polyShape': newIntersect,
                                        'repPtShape': newIntersect.centroid,
                                        'nodeIDs': newNodeIDs
                                    })

        # Set `pointer` to be `endPointer`
        pointer = endPointer
        checkOrder += 1

    lstSteinerZone = [{
        'poly': polys[n]['poly'] if 'poly' in polys[n] else [polys[n]['loc']],
        'repPt': list(shapely.Polygon([[p[0], p[1]] for p in polys[n]['poly']]).centroid.coords[0]) if 'poly' in polys[n] else polys[n]['loc'],
        'nodeIDs': [n]
    } for n in polys]

    for n in lstSteinerZoneShape:
        if (type(n['polyShape']) == shapely.Point):
            lstSteinerZone.append({
                'poly': [[n['polyShape'].x, n['polyShape'].y]],
                'repPt': list(n['repPtShape'].coords[0]),
                'nodeIDs': [i for i in n['nodeIDs']]
            })
        else:
            lstSteinerZone.append({
                'poly': [i for i in mapping(n['polyShape'])['coordinates'][0]],
                'repPt': list(n['repPtShape'].coords[0]),
                'nodeIDs': [i for i in n['nodeIDs']]
            })    

    return lstSteinerZone

def polygonsAlongLocSeq(seq, polygons:dict, polyFieldName = 'poly'):
     # First, for each leg in the seq, find the individual polygons intersect with the leg
    actions = {}
    actionIDOnLeg = 0
    for i in range(len(seq) - 1):
        seg = [seq[i], seq[i + 1]]

        # NOTE: 准备用segment tree
        # NOTE: For now use the naive way - checking the bounding box
        for pID in polygons:
            # 根据Seg和poly的相交情况做判断
            segIntPoly = intSeg2Poly(seg = seg, poly = polygons[pID][polyFieldName])
            # if (type(segIntPoly) == list or segIntPoly['status'] != 'NoCross'):
                # print(segIntPoly, pID)
            # 如果相交得到多个部分，则分别进行处理
            if (type(segIntPoly) == list):
                for intPart in segIntPoly:
                    if (intPart['status'] == 'Cross' and intPart['intersectType'] == 'Point'):
                        intPt = segIntPoly['intersect']
                        actions[actionIDOnLeg] = {
                            'loc': intPt,
                            'action': 'touch',
                            'polyID': pID,                       
                        }
                        actionIDOnLeg += 1
                    elif (intPart['status'] == 'Cross' and intPart['intersectType'] == 'Segment'):
                        intPt1 = segIntPoly['intersect'][0]
                        intPt2 = segIntPoly['intersect'][1]
                        intPt1InnerFlag = False
                        intPt2InnerFlag = False
                        if (isPtInPoly(intPt1, polygons[pID][polyFieldName], interiorOnly=True)):
                            intPt1InnerFlag = True
                        if (isPtInPoly(intPt2, polygons[pID][polyFieldName], interiorOnly=True)):
                            intPt2InnerFlag = True
                        if (distEuclideanXY(seq[i], intPt1)['dist'] <
                            distEuclideanXY(seq[i], intPt2)['dist']):
                            if (not intPt1InnerFlag):
                                actions[actionIDOnLeg] = {
                                    'loc': intPt1,
                                    'action': 'enter',
                                    'polyID': pID,
                                }
                                actionIDOnLeg += 1
                            if (not intPt2InnerFlag):
                                actions[actionIDOnLeg] = {
                                    'loc': intPt2,
                                    'action': 'leave',
                                    'polyID': pID,
                                }
                                actionIDOnLeg += 1
                        else:
                            if (not intPt1InnerFlag):
                                actions[actionIDOnLeg] = {
                                    'loc': intPt1,
                                    'action': 'leave',
                                    'polyID': pID,
                                }
                                actionIDOnLeg += 1
                            if (not intPt2InnerFlag):
                                actions[actionIDOnLeg] = {
                                    'loc': intPt2,
                                    'action': 'enter',
                                    'polyID': pID,
                                }
                                actionIDOnLeg += 1
            else:
                if (segIntPoly['status'] == 'NoCross'):
                    # No intersection pass
                    pass
                elif (segIntPoly['status'] == 'Cross' and segIntPoly['intersectType'] == 'Point'):
                    intPt = segIntPoly['intersect']
                    actions[actionIDOnLeg] = {
                        'loc': intPt,
                        'action': 'touch',
                        'polyID': pID,
                    }
                    actionIDOnLeg += 1
                elif (segIntPoly['status'] == 'Cross' and segIntPoly['intersectType'] == 'Segment'):
                    intPt1 = segIntPoly['intersect'][0]
                    intPt2 = segIntPoly['intersect'][1]
                    intPt1InnerFlag = False
                    intPt2InnerFlag = False
                    if (isPtInPoly(intPt1, polygons[pID][polyFieldName], interiorOnly=True)):
                        intPt1InnerFlag = True
                    if (isPtInPoly(intPt2, polygons[pID][polyFieldName], interiorOnly=True)):
                        intPt2InnerFlag = True
                    if (distEuclideanXY(seq[i], intPt1)['dist'] <
                        distEuclideanXY(seq[i], intPt2)['dist']):
                        if (not intPt1InnerFlag):
                            actions[actionIDOnLeg] = {
                                'loc': intPt1,
                                'action': 'enter',
                                'polyID': pID,
                            }
                            actionIDOnLeg += 1
                        if (not intPt2InnerFlag):
                            actions[actionIDOnLeg] = {
                                'loc': intPt2,
                                'action': 'leave',
                                'polyID': pID,
                            }
                            actionIDOnLeg += 1
                    else:
                        if (not intPt1InnerFlag):
                            actions[actionIDOnLeg] = {
                                'loc': intPt1,
                                'action': 'leave',
                                'polyID': pID,
                            }
                            actionIDOnLeg += 1
                        if (not intPt2InnerFlag):
                            actions[actionIDOnLeg] = {
                                'loc': intPt2,
                                'action': 'enter',
                                'polyID': pID,
                            }
                            actionIDOnLeg += 1
    sortedIDOnLeg = locSeqSortNodes(seq, actions, allowNotIncludeFlag = False)['sortedNodeIDs']

    sortedActions = [actions[i] for i in sortedIDOnLeg]
    sortedActions = [sortedActions[i] for i in range(len(sortedActions)) if (
        i == 0
        or not is2PtsSame(sortedActions[i - 1]['loc'], sortedActions[i]['loc']))]

    return sortedActions

def polyClockWise(poly) -> bool:
    """
    Given a poly, return True if the poly is clockwise, False otherwise

    Parameters
    ----------
    poly: poly, required
        The poly

    Returns
    -------
    bool
        True if clockwise, false otherwise

    """
    numCW = 0
    numCCW = 0

    for i in range(len(poly)):
        pt1 = poly[i]
        pt2 = None
        pt3 = None
        if (i == len(poly) - 1):
            pt2 = poly[0]
            pt3 = poly[1]
        elif (i == len(poly) - 2):
            pt2 = poly[i + 1]
            pt3 = poly[0]
        else:
            pt2 = poly[i + 1]
            pt3 = poly[i + 2]
        cwFlag = is3PtsClockWise(pt1, pt2, pt3)
        if (cwFlag == True):
            numCW += 1
        elif (cwFlag == False):
            numCCW += 1
    if (numCW > numCCW):
        return True
    else:
        return False

def poly2CW(poly) -> poly:
    if (polyClockWise(poly)):
        return poly
    else:
        return [poly[len(poly) - 1 - i] for i in range(len(poly))]

def poly2CCW(poly) -> poly:
    if (polyClockWise(poly)):
        return [poly[len(poly) - 1 - i] for i in range(len(poly))]
    else:
        return poly

# Visibility check ============================================================
def polysVisibleGraph(polys:polys) -> dict:
    """
    Create a visual graph for given polys.

    Parameters
    ----------
    polys: polys, required
        The polys to create visual graph

    Return
    ------
    dict
        Each polygon has a index p, each point in the polygon has a index e, therefore, a (p, e) pair defines a location. The visual graph returns use (p, e) as keys, collects the location of (p, e) in 'loc', and finds the set of visible (p, e) in 'visible'

    """

    vg = {}
    for p in range(len(polys)):
        for e in range(len(polys[p])):
            vg[(p, e)] = {'loc': polys[p][e], 'visible': []}
            W = _visPtAmongPolys((p, e), polys, knownVG=vg)
            for w in W:
                vg[(p, e)]['visible'].append(w)
    return vg

def _visPtAmongPolys(v:int|str|tuple, polys:polys, standalonePts:dict|None=None, knownVG:dict={}) -> list:
    # NOTE: 该函数不需要使用shapely
    vertices = {}
    polyVertices = []
    for p in range(len(polys)):
        for e in range(len(polys[p])):
            vertices[(p, e)] = {
                'loc': polys[p][e],
                'visible': []
            }
            polyVertices.append((p, e))
    if (standalonePts != None):
        if (v not in standalonePts):
            raise MissingParameterError("ERROR: Cannot find `v` in `polys` or `standalonePts`")
        else:
            vertices[v] = {
                'loc': standalonePts[v]['loc'],
                'visible': []
            }

    # 把所有的poly vertices按到v的距离排序，从而给每个点得到一个可排序的唯一编码
    verticeDistIndex = {}
    sortedSeq = nodeSeqByDist(
        nodes = vertices,
        refLoc = vertices[v]['loc'],
        nodeIDs = polyVertices)
    for i in range(len(sortedSeq)):
        verticeDistIndex[sortedSeq[i]] = i
    # 把所有的poly vertices按到v的角度排序，从x轴正方向开始，从而得到后续可视性检查的顺序
    sweepSeq = nodeSeqBySweeping(
        nodes = vertices,
        nodeIDs = polyVertices,
        refLoc = vertices[v]['loc'],
        initDeg = 90)

    # 用一个红黑树来维护射线通过的边，边的键值用(closer-index, further-index)，这样排序的时候无论如何都能保持距离关系
    def rayIntersectPolyEdges(ray, polys) -> RedBlackTree:
        T = RedBlackTree()
        # FIXME: 这个显然需要用线段树来优化，在这里就先实现再说吧
        for p in range(len(polys)):
            for i in range(-1, len(polys[p]) - 1):
                edge = [polys[p][i], polys[p][i + 1]]
                # NOTE: 注意，这里射线可以与线段交于端点上
                if (isSegIntRay(edge, ray, interiorOnly=True)):
                    # NOTE: 防止数组溢出
                    k = i
                    if (i == -1):
                        k = len(polys[p]) - 1
                    vIdx1 = verticeDistIndex[(p, k)]
                    vIdx2 = verticeDistIndex[(p, i + 1)]
                    T.insert(RedBlackTreeNode(
                        key = (min(vIdx1, vIdx2), max(vIdx1, vIdx2)), 
                        value = [(p, k), (p, i + 1)]))
        return T

    # 初始射线方向为x-轴正方向
    xAxisRay = [vertices[v]['loc'], (vertices[v]['loc'][0] + 1, vertices[v]['loc'][1])]
    T = rayIntersectPolyEdges(xAxisRay, polys)

    # 给定一个障碍物的顶点w_i，确定v是否可以见到w_i
    # NOTE: 这里的wi, wim（也就是w_{i-1})得是polys的顶点
    def visible(wi, wim, polys):
        # 如果wiv已经确定可视，直接返回True
        if (wi in knownVG and v in knownVG[wi]['visible']):
            # print("Time saved")
            return True

        # 所在的polygon的编号
        polyV = polys[v[0]] if len(v) == 2 else None
        polyW = polys[wi[0]]
        vwi = [vertices[v]['loc'], vertices[wi]['loc']]

        vNext = None
        vPrev = None
        if (len(v) == 2):
            vNext = (v[0], v[1] + 1 if v[1] < len(polys[v[0]]) - 1 else 0)
            vPrev = (v[0], v[1] - 1 if v[1] > 0 else len(polys[v[0]]) - 1)
        wiNext = (wi[0], wi[1] + 1 if wi[1] < len(polys[wi[0]]) - 1 else 0)
        wiPrev = (wi[0], wi[1] - 1 if wi[1] > 0 else len(polys[wi[0]]) - 1)        

        # 判断是否是相邻节点，相邻节点直接返回可见
        if (polyV == polyW and (wi == vNext or wi == vPrev)):
            return True

        # 需要w_{i-1}不存在，或者w_{i-1}不在线段vwi上
        notOnlineFlag = False
        if (wim == None or not isPtOnSeg(wim, vwi)):
            notOnlineFlag = True

        # 若T非空，查有没有阻拦线段
        noSegBlockFlag = True
        if (notOnlineFlag and not T.isEmpty and not T.min(T.root).isNil):
            for edge in T.traverse():
                if (isSegIntSeg(
                        seg1 = [vertices[edge.value[0]]['loc'], vertices[edge.value[1]]['loc']], 
                        seg2 = vwi, 
                        interiorOnly = True)):
                    noSegBlockFlag = False
                    break

        # 如果visibleFlag，查交点上是否相切
        int2PolyVTangenFlag = False
        int2PolyWTangenFlag = False
        if (noSegBlockFlag):
            # Check point V
            if (polyV == None):
                int2PolyVTangenFlag = True
            else:
                segV = [vertices[vNext]['loc'], vertices[vPrev]['loc']]
                int2PolyVTangenFlag = not isLineIntSeg(vwi, segV, interiorOnly=True)
        if (int2PolyVTangenFlag):
            # Check point W
            segW = [vertices[wiNext]['loc'], vertices[wiPrev]['loc']]
            int2PolyWTangenFlag = not isLineIntSeg(vwi, segW, interiorOnly=True)
        bothEndTangenFlag = int2PolyVTangenFlag and int2PolyWTangenFlag

        # 如果没有阻挡线段，查是不是在多边形内部
        # NOTE: 似乎可以通过查中点确定是不是在多边形内部
        visibleFlag = False
        if (noSegBlockFlag and bothEndTangenFlag and not isSegIntPoly(
                seg = vwi,
                poly = polyW,
                interiorOnly = True)):
            visibleFlag = True        


        if (visibleFlag):
            return True

        # 前面任何一关不通过，则不可视
        return False

    W = []
    for i in range(len(sweepSeq)):
        wi = sweepSeq[i]
        wim = sweepSeq[i - 1] if i > 1 else None
        # print(T)
        if (not is2PtsSame(vertices[v]['loc'], vertices[wi]['loc']) and visible(wi, wim, polys)):
            W.append(wi)

        # 将wi的边加入/移出平衡树T
        # wi所在的poly编号
        polyIdx = wi[0]
        # wi在poly内的编号
        idInPoly = wi[1]

        # Edge 1: wi -> wi.next
        wiNext = (wi[0], wi[1] + 1 if wi[1] < len(polys[wi[0]]) - 1 else 0)
        # Edge 2: wi.prev -> wi
        wiPrev = (wi[0], wi[1] - 1 if wi[1] > 0 else len(polys[wi[0]]) - 1)

        # 判断edgeNext
        vIdxWi = verticeDistIndex[wi]
        vIdxWiNext = verticeDistIndex[wiNext]
        vIdxWiPrev = verticeDistIndex[wiPrev]
        if (is3PtsClockWise(vertices[v]['loc'], vertices[wi]['loc'], vertices[wiNext]['loc'])):
            if (T.query((min(vIdxWi, vIdxWiNext), max(vIdxWi, vIdxWiNext))).isNil):
                T.insert(RedBlackTreeNode(
                    key = (min(vIdxWi, vIdxWiNext), max(vIdxWi, vIdxWiNext)), 
                    value = [wi, wiNext]))
        else:
            if (not T.query((min(vIdxWi, vIdxWiNext), max(vIdxWi, vIdxWiNext))).isNil):
                T.delete((min(vIdxWi, vIdxWiNext), max(vIdxWi, vIdxWiNext)))        
        if (is3PtsClockWise(vertices[v]['loc'], vertices[wi]['loc'], vertices[wiPrev]['loc'])):
            if (T.query((min(vIdxWiPrev, vIdxWi), max(vIdxWiPrev, vIdxWi))).isNil):
                T.insert(RedBlackTreeNode(
                    key = (min(vIdxWiPrev, vIdxWi), max(vIdxWiPrev, vIdxWi)), 
                    value = [wiPrev, wi]))
        else:
            if (not T.query((min(vIdxWiPrev, vIdxWi), max(vIdxWiPrev, vIdxWi))).isNil):
                T.delete((min(vIdxWiPrev, vIdxWi), max(vIdxWiPrev, vIdxWi)))
    return W

# Time seq related ============================================================
def snapInTimedSeq(timedSeq: list[tuple[pt, float]], t: float) -> dict:
    """
    Given a timedSeq, return the location, speed, and trajectory at time t

    Parameters
    ----------

    timedSeq: timedSeq, required
        The timed sequence
    t: float, required
        The snapshot timestamp

    Return
    ------
    dict
        >>> {
        ...     'loc': location,
        ...     'speed': speed,
        ...     'trajectory': trajectory
        ... }

    """

    for i in range(len(timedSeq) - 1):
        if (timedSeq[i][1] > timedSeq[i + 1][1]):
            raise UnsupportedInputError("ERROR: `timedSeq` should be a non-descending sequence.")
        if (timedSeq[i][1] == timedSeq[i + 1][1] 
            and (abs(timedSeq[i][0][0] - timedSeq[i + 1][0][0]) >= CONST_EPSILON
                or abs(timedSeq[i][0][1] - timedSeq[i + 1][0][1]) >= CONST_EPSILON)):
            raise UnsupportedInputError("ERROR: an object cannot be two places at the same time.")

    curLocX = None
    curLocY = None
    curSpeed = None
    curTrajectory = None

    if (t <= timedSeq[0][1]):
        return {
            'loc': timedSeq[0][0],
            'speed': 0,
            'trajectory': None
        }
    if (t >= timedSeq[-1][1]):
        return {
            'loc': timedSeq[-1][0],
            'speed': 0,
            'trajectory': None
        }

    for i in range(len(timedSeq) - 1):
        if (timedSeq[i][1] <= t < timedSeq[i + 1][1]):
            dist = distEuclideanXY(timedSeq[i][0], timedSeq[i + 1][0])['dist']
            if (dist > 0):
                dt = (t - timedSeq[i][1]) / (timedSeq[i + 1][1] - timedSeq[i][1])                
                curLocX = timedSeq[i][0][0] + (timedSeq[i + 1][0][0] - timedSeq[i][0][0]) * dt
                curLocY = timedSeq[i][0][1] + (timedSeq[i + 1][0][1] - timedSeq[i][0][1]) * dt
                curSpeed = dist / (timedSeq[i + 1][1] - timedSeq[i][1])
                trajX = (timedSeq[i + 1][0][0] - timedSeq[i][0][0]) / dist
                trajY = (timedSeq[i + 1][0][1] - timedSeq[i][0][1]) / dist
                curTrajectory = [(curLocX, curLocY), (curLocX + trajX, curLocY + trajY)]
            else:
                curLocX = timedSeq[i][0][0]
                curLocY = timedSeq[i][0][1]
                curSpeed = 0
                curTrajectory = None
            break
    return {
        'loc': [curLocX, curLocY],
        'speed': curSpeed,
        'trajectory': curTrajectory
    }

def traceInTimedSeq(timedSeq: list[tuple[pt, float]], ts: float, te: float) -> list[pt]:
    """
    Given a timedSeq, a start time and an end time, returns the trace between start time and end time.

    Parameters
    ----------

    timedSeq: timedSeq, required
        The timed sequence
    ts: float, required
        The start time
    te: float, required
        The end time

    Return
    ------
    list
        The line segment sequence between start time and end time in the timedSeq

    """

    for i in range(len(timedSeq) - 1):
        if (timedSeq[i][1] > timedSeq[i + 1][1]):
            raise UnsupportedInputError("ERROR: `timedSeq` should be a non-descending sequence.")
        if (timedSeq[i][1] == timedSeq[i + 1][1] 
            and (abs(timedSeq[i][0][0] - timedSeq[i + 1][0][0]) >= CONST_EPSILON
                or abs(timedSeq[i][0][1] - timedSeq[i + 1][0][1]) >= CONST_EPSILON)):
            raise UnsupportedInputError("ERROR: an object cannot be two places at the same time.")
    
    trace = []

    if (ts >= te):
        raise UnsupportedInputError("ERROR: `ts` should be earlier than `te`")
    if (ts <= timedSeq[0][1] and te >= timedSeq[-1][1]):
        return [timedSeq[i][0] for i in range(len(timedSeq))]
    if (ts >= timedSeq[-1][1]):
        return []
    if (te <= timedSeq[0][1]):
        return []

    tsIndex = -1
    teIndex = -1
    tsLoc = []
    teLoc = []

    if (ts <= timedSeq[0][1]):
        tsIndex = 0
        ts = timedSeq[0][1]
        tsLoc = timedSeq[0][0]
    if (te >= timedSeq[-1][1]):
        teIndex = len(timedSeq) - 1
        te = timedSeq[-1][1]
        teLoc = timedSeq[-1][0]

    for i in range(len(timedSeq) - 1):
        if (timedSeq[i][1] <= ts < timedSeq[i + 1][1]):
            tsIndex = i
            dTs = (ts - timedSeq[i][1]) / (timedSeq[i + 1][1] - timedSeq[i][1])
            tsX = timedSeq[i][0][0] + (timedSeq[i + 1][0][0] - timedSeq[i][0][0]) * dTs
            tsY = timedSeq[i][0][1] + (timedSeq[i + 1][0][1] - timedSeq[i][0][1]) * dTs
            tsLoc = [tsX, tsY]
            for j in range(tsIndex, len(timedSeq) - 1):
                if (timedSeq[j][1] <= te < timedSeq[j + 1][1]):
                    teIndex = j
                    dTe = (te - timedSeq[j][1]) / (timedSeq[j + 1][1] - timedSeq[j][1])
                    teX = timedSeq[j][0][0] + (timedSeq[j + 1][0][0] - timedSeq[j][0][0]) * dTe
                    teY = timedSeq[j][0][1] + (timedSeq[j + 1][0][1] - timedSeq[j][0][1]) * dTe
                    teLoc = [teX, teY]

    if (tsIndex == teIndex):
        trace = [tsLoc, teLoc]
    elif (tsIndex + 1 == teIndex):
        trace.append(tsLoc)
        trace.append(timedSeq[tsIndex + 1][0])
        trace.append(teLoc)
    else:
        trace.append(tsLoc)
        for i in range(tsIndex + 1, teIndex + 1):
            trace.append(timedSeq[i][0])
        trace.append(teLoc)

    return trace

def locSeq2TimedSeq(seq: list[pt], vehSpeed: float, timeEachStop: float = 0, startTime: float = 0):
    timedSeq = []
    accTime = startTime
    for i in range(len(seq) - 1):
        timedSeq.append((seq[i], accTime))
        if (timeEachStop > 0):
            accTime += timeEachStop
            timedSeq.append((seq[i], accTime))
        accTime += (distEuclideanXY(seq[i], seq[i + 1])['dist']) / vehSpeed
    timedSeq.append((seq[-1], accTime))
    return timedSeq

# Loc seq related =============================================================
def locSeqSortPts(seq: list[pt], pts: list[pt], allowNotIncludeFlag = True) -> list:
    # First, calculate accumulated dist since start for each turning point of locSeq
    acc = [0]
    sofar = 0
    for i in range(len(seq) - 1):
        sofar += distEuclideanXY(seq[i], seq[i + 1])['dist']
        acc.append(sofar)
    ptHeap = []
    notIncludedPts = []

    # 看看每个节点落在哪个区间
    for pt in pts:
        isOnSegFlag = False
        for i in range(len(seq) - 1):
            # 若和多个线段相交，只考虑第一次相交的情况，否则turnpoint上的始终会被记录两次
            if (isPtOnSeg(pt, [seq[i], seq[i + 1]])):
                isOnSegFlag = True
                addDist = distEuclideanXY(seq[i], pt)['dist']
                heapq.heappush(ptHeap, (acc[i] + addDist, pt))
                break
        if (not isOnSegFlag):
            if (allowNotIncludeFlag):
                notIncludedPts.append(pt)
            else:
                raise UnsupportedInputError("ERROR: %s does not belong to given loc sequence" % str(pt))

    sortedPts = []
    while (len(ptHeap) > 0):
        sortedPts.append(heapq.heappop(ptHeap)[1])

    return {
        'sortedPts': sortedPts,
        'notIncludedPts': notIncludedPts
    }

def locSeqSortNodes(seq: list[pt], nodes: dict, locFieldName: str = 'loc', allowNotIncludeFlag: bool = True) -> list: 
    acc = [0]
    sofar = 0
    for i in range(len(seq) - 1):
        sofar += distEuclideanXY(seq[i], seq[i + 1])['dist']
        acc.append(sofar)
    ptHeap = []
    notIncludedNodeIDs = []

    for nodeID in nodes:
        isOnSegFlag = False
        for i in range(len(seq) - 1):
            if (isPtOnSeg(nodes[nodeID][locFieldName], [seq[i], seq[i + 1]])):
                isOnSegFlag = True
                addDist = distEuclideanXY(seq[i], nodes[nodeID][locFieldName])['dist']
                print(acc[i] + addDist, nodeID)
                heapq.heappush(ptHeap, (acc[i] + addDist, nodeID))
                break
        if (not isOnSegFlag):
            if (allowNotIncludeFlag):
                notIncludedNodeIDs.append(nodeID)
            else:
                raise UnsupportedInputError("ERROR: %s does not belong to given loc sequence" % str(nodeID))
    
    sortedNodeIDs = []
    while (len(ptHeap) > 0):
        sortedNodeIDs.append(heapq.heappop(ptHeap)[1])

    return {
        'sortedNodeIDs': sortedNodeIDs,
        'notIncludedNodeIDs': notIncludedNodeIDs
    }

# Area calculation ============================================================
def calTriangleAreaEdge(a: float, b: float, c: float) -> float:
    
    s = (a / 2 + b / 2 + c / 2)
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return area

def calTriangleAreaXY(pt1: pt, pt2: pt, pt3: pt) -> float:
    
    [x1, y1] = pt1
    [x2, y2] = pt2
    [x3, y3] = pt3
    val = (x2 * y3 + x3 * y1 + x1 * y2) - (x2 * y1 + x3 * y2 + x1 * y3)
    area = abs(val)
    return area

def calPolyAreaXY(poly: poly) -> float:
    lstTriangle = tripy.earclip(poly)
    area = 0
    for i in range(len(lstTriangle)):
        area += calTriangleAreaXY(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2])
    return area

def calPolyAreaLatLon(polyLatLon: poly) -> float:
    """Returns the area surrounded by polyLatLon on the Earth"""
    try:
        from pyproj import Geod
    except:
        raise ImportError

    # NOTE: shapely is in [lon, lat] format
    rev = []
    for p in polyLatLon:
        rev.append((p[1], p[0]))
    polygon = shapely.Polygon(rev)

    # Using pyproj to calculate
    # Ref: https://hypc.github.io/2020/03/16/python-geo-area/
    geod = Geod(ellps = "WGS84")
    area = abs(geod.geometry_area_perimeter(polygon)[0])

    return area

def calPolyPerimeterXY(poly: poly) -> float:
    p = 0
    for i in range(-1, len(poly) - 1):
        p += distEuclideanXY(poly[i], poly[i + 1])['dist']
    return p

# Location of points ==========================================================
def headingXY(pt1: pt, pt2: pt) -> float:
    """
    Given two points, returns the direction from the first point to the second point

    Parameters
    ----------
    pt1: pt, required
        The first point
    pt2: pt, required
        The second point

    Return
    ------
    float
        The direction from pt1 to pt2, north as 0, clock-wise

    """
    vec = (pt2[0] - pt1[0], pt2[1] - pt1[1])
    (_, vDeg) = vecXY2Polar(vec)
    return vDeg

def headingLatLon(pt1: pt, pt2: pt) -> float:
    """
    Given current location and a goal location, calculate the heading. North is 0-degrees, clock-wise

    Parameters
    ----------
    pt1: pt, required
        The first point
    pt2: pt, required
        The second point

    Return
    ------
    float
        The direction from pt1 to pt2, north as 0, clock-wise

    """

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

def ptInDistXY(pt: pt, direction: int|float, dist: int|float) -> pt:
    """
    A location in distance with given direction.

    Parameters
    ----------
    pt: pt, required
        The origin point
    direction: float, required
        The direction from origin point
    dist: float, required

    Returns
    -------
    pt
        The new location

    """
    x = pt[0] + dist * math.sin(math.radians(direction))
    y = pt[1] + dist * math.cos(math.radians(direction))
    return (x, y)

def ptInDistLatLon(pt: pt, direction: int|float, distMeters: int|float) -> pt:
    """
    A location in distance with given direction in lat/lon.

    Parameters
    ----------
    pt: pt, required
        The origin point
    direction: float, required
        The direction from origin point
    dist: float, required

    Returns
    -------
    pt
        The new location

    """
    # Bearing in degrees: 0 – North, 90 – East, 180 – South, 270 or -90 – West.
    newLoc = list(geopy.distance.distance(meters=distMeters).destination(point=pt, bearing=direction))[:2]
    return newLoc

def circleByCenterLatLon(center: pt, radius: int|float, lod: int = 30) -> poly:
    """
    Create a circle by the center and a radius. The circle is approximated by a x-gon, e.g., 30-gon polygon.

    Parameters
    ----------
    center: pt, required
        The center of the circle
    radius: float, required
        The radius of the circle
    lod: int, optional, default as 30
        Level of details. The circle is approximated as a x-gon. E.g. 30-gon

    Return
    ------
    poly
        A polygon that approximates the circle.

    """
    circle = []
    for i in range(lod):
        deg = float(360 * i) / float(lod)
        pt = ptInDistLatLon(pt = center, direction = deg, distMeters = radius)
        circle.append(pt)
    return circle

def circleByCenterXY(center: pt, radius: int|float, lod: int = 30) -> poly:
    """
    Create a circle by the center and a radius. The circle is approximated by a x-gon, e.g., 30-gon polygon.

    Parameters
    ----------
    center: pt, required
        The center of the circle
    radius: float, required
        The radius of the circle
    lod: int, optional, default as 30
        Level of details. The circle is approximated as a x-gon. E.g. 30-gon

    Return
    ------
    poly
        A polygon that approximates the circle.

    """
    circle = []
    for i in range(lod):
        deg = float(360 * i) / float(lod)
        pt = ptInDistXY(pt = center, direction = deg, dist = radius)
        circle.append(pt)
    return circle

# Sort nodes ==================================================================
def nodeSeqByDist(nodes: dict, nodeIDs: list|str = 'All', locFieldName = 'loc', refLoc: pt|None = None, refNodeID: int|str|None = None) -> list:
    # Define nodeIDs
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define refLoc
    if (refLoc == None):
        if (refNodeID == None):
            raise MissingParameterError("ERROR: Missing reference location")
        elif (refNodeID not in nodeIDs):
            raise OutOfRangeError("ERROR: `refNodeID` cannot be found.")
        else:
            refLoc = nodes[refNodeID][locFieldName]

    # Sort distance
    sortedSeq = []
    sortedSeqHeap = []
    for n in nodeIDs:
        dist = scaleDist(loc1 = refLoc, loc2 = nodes[n][locFieldName], edges = 'Euclidean')['dist']
        heapq.heappush(sortedSeqHeap, (dist, n))
    while (len(sortedSeqHeap) > 0):
        sortedSeq.append(heapq.heappop(sortedSeqHeap)[1])  

    return sortedSeq

def nodeSeqBySweeping(nodes: dict, nodeIDs: list|str = 'All', locFieldName = 'loc', refLoc: None|pt = None, refNodeID: int|str|None = None, isClockwise: bool = True, initDeg: float = 0) -> list:
    """Given a set of locations, and a center point, gets the sequence from sweeping"""
    # Define nodeIDs
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Initialize centroid
    if (refLoc == None):
        lstNodeLoc = []
        for n in nodeIDs:
            lstNodeLoc.append(shapely.Point(nodes[n][locFieldName][0], nodes[n][locFieldName][1]))
        refLoc = list(shapely.centroid(shapely.MultiPoint(points = lstNodeLoc)))

    # Initialize heap
    degHeap = []
    refLocNodes = []
    
    # Build heap
    for n in nodeIDs:
        dist = distEuclideanXY(nodes[n][locFieldName], refLoc)['dist']
        # If the nodes are too close, separate it/them
        if (dist <= CONST_EPSILON):
            refLocNodes.append(n)
        else:
            dx = nodes[n][locFieldName][0] - refLoc[0]
            dy = nodes[n][locFieldName][1] - refLoc[1]
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

    # Sweep
    sweepSeq = []
    while (len(degHeap)):
        sweepSeq.append(heapq.heappop(degHeap)[2])
    sweepSeq.extend(refLocNodes)

    return sweepSeq

def nodesInIsochrone(nodes: dict, nodeIDs: list|str = 'All', locFieldName = 'loc', refLoc: pt|None = None, refNodeID: int|str|None = None, isoRange: float = None, sortFlag: bool = False) -> list: 
    # FIXME: Need an algorithm to filter out locations that are clearly too far from refLoc
    # Define nodeIDs
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define refLoc
    if (refLoc == None):
        if (refNodeID == None):
            raise MissingParameterError("ERROR: Missing reference location")
        elif (refNodeID not in nodeIDs):
            raise OutOfRangeError("ERROR: `refNodeID` cannot be found.")
        else:
            refLoc = nodes[refNodeID][locFieldName]

    # Sort distance
    if (isoRange == None):
        raise MissingParameterError("ERROR: Missing required field `isoRange`.")
    nearSet = []
    nearest = None
    nearestDist = float('inf')
    if (sortFlag):        
        nearSetHeap = []
        for n in nodeIDs:
            dist = scaleDist(loc1 = refLoc, loc2 = nodes[n][locFieldName], edges = 'Euclidean')['dist']
            if (dist <= isoRange):
                heapq.heappush(nearSetHeap, (dist, n))
        while (len(nearSetHeap) > 0):
            nearSet.append(heapq.heappop(nearSetHeap)[1])
        nearest = nearSet[0]
    else:
        for n in nodeIDs:
            dist = scaleDist(loc1 = refLoc, loc2 = nodes[n][locFieldName], edges = 'Euclidean')['dist']
            if (dist <= isoRange):
                nearSet.append(n)
            if (dist <= nearestDist):
                nearest = n
                nearestDist = dist
    return {
        'nearSet': nearSet,
        'nearest': nearest
    }

# Distance calculation ========================================================
def distEuclideanXY(pt1: pt, pt2: pt) -> dict:
    """
    Gives a Euclidean distance based on two coords.

    Parameters
    ----------
    pt1: pt, required
        The first location
    pt2: pt, required
        The second location

    Returns
    -------
    dict
        A dictionary, with the distance in 'dist', and the path in 'path'
    """
    return {
        'dist': math.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2),
        'path': [pt1, pt2]
    }

def distManhattenXY(pt1: pt, pt2: pt) -> dict:
    """
    Gives a Manhatten distance based on two coords.

    Parameters
    ----------
    pt1: pt, required
        The first location
    pt2: pt, required
        The second location

    Returns
    -------
    dict
        A dictionary, with the distance in 'dist', and the path in 'path'
    """
    return {
        'dist': abs(pt1[0] - pt2[0]) + abs(pt1[1] - pt2[1]),
        'path': [pt1, (pt1[0], pt2[1]), pt2]
    }

def distBtwPolysXY(pt1:pt, pt2:pt, polys:polys, polyVG: dict = None) -> dict:
    """
    Gives a Manhatten distance based on two coords.

    Parameters
    ----------
    pt1: pt, required
        The first location
    pt2: pt, required
        The second location
    polys: polys, required
        The polygons as barriers.
    polyVG: dict, optional, default as None
        The pre-calculated visual-graph using :func:`~vrpSolver.polysVisibleGraph()`. To avoid repeated calculation

    Returns
    -------
    dict
        A dictionary, with the distance in 'dist', and the path in 'path'
    """

    # Reference: Computational Geometry: Algorithms and Applications Third Edition
    # By Mark de Berg et al. Page 326 - 330
    # With some modifications

    # First check if start pt or end pt is in one of the polygons =============
    for poly in polys:
        if (isPtInPoly(pt1, poly, interiorOnly = True)):
            raise OutOfRangeError("Point (%s, %s) is inside `polys` when it is not suppose to." % (pt1[0], pt1[1]))
        if (isPtInPoly(pt2, poly, interiorOnly = True)):
            raise OutOfRangeError("Point (%s, %s) is inside `polys` when it is not suppose to." % (pt2[0], pt2[1]))

    # Quick checkout ==========================================================
    visibleDirectly = True
    for poly in polys:
        if (isSegIntPoly([pt1, pt2], poly)):
            visibleDirectly = False
            break
    if (visibleDirectly):
        return {
            'dist': distEuclideanXY(pt1, pt2)['dist'],
            'path': [pt1, pt2]
        }

    # Create visible graph for polys ==========================================
    if (polyVG == None):      
        for p in range(len(polys)):
            polys[p] = [polys[p][i] for i in range(len(polys[p])) if distEuclideanXY(polys[p][i], polys[p][i - 1])['dist'] > CONST_EPSILON]
        polyVG = polysVisibleGraph(polys)

    # Create a visible graph ==================================================
    # NOTE: startPt可视的vertices将不需要测试是不是相互之间可视，同样地，可视endPt的vertices之间也不需要可视
    vertices = {}
    for p in polyVG:
        vertices[p] = {
            'loc': polyVG[p]['loc'],
            'visible': [i for i in polyVG[p]['visible']]
        }
    vertices['s'] = {'loc': pt1, 'visible': []}
    Ws = _visPtAmongPolys('s', polys, {'s': {'loc': pt1, 'visible': []}})
    vertices['s']['visible'] = Ws
    vertices['e'] = {'loc': pt2, 'visible': []}
    We = _visPtAmongPolys('e', polys, {'e': {'loc': pt2, 'visible': []}})
    vertices['e']['visible'] = We

    # Find shortest path ======================================================
    vg = nx.Graph()
    for v in vertices:
        vg.add_node(v)
    for v in vertices:
        for e in vertices[v]['visible']:
            vg.add_edge(v, e, weight=distEuclideanXY(vertices[v]['loc'], vertices[e]['loc'])['dist'])
    sp = nx.dijkstra_path(vg, 's', 'e')

    dist = 0
    for i in range(len(sp) - 1):
        dist += distEuclideanXY(vertices[sp[i]]['loc'], vertices[sp[i + 1]]['loc'])['dist']

    return {
        'dist': dist,
        'path': [vertices[wp]['loc'] for wp in sp]
    }

def distLatLon(pt1: pt, pt2: pt, distUnit: str = 'meter') -> dict:
    """
    Gives a distance based on two lat/lon coords.

    Parameters
    ----------
    pt1: pt, required
        The first location
    pt2: pt, required
        The second location

    Returns
    -------
    dict
        A dictionary, with the distance in 'dist', and the path in 'path'
    """
    
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
    return {
        'dist': 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a)),
        'path': [pt1, pt2]
    }

def distOnGrid(pt1: pt, pt2: pt, column, row, barriers = [], algo: str = 'A*', **kwargs) -> dict:
    """
    Given two coordinates on the grid, finds the 'shortest' path to travel

    Parameters
    ----------

    pt1: pt, required
        Starting location on the grid
    pt2: pt, required
        Ending location on the grid 
    column: int, required
        Number of columns
    row: int, required
        Number of rows
    barriers: list[pt], optional, default as []
        A list of coordinates as barriers on the grid.
    algo: dict, required, default as 'A*'
        The algorithm configuration. For example

        1) A*, use the A star algorithm, additional information needed is as follows
            - measure: str, optional, default as 'Manhatten'            
    **kwargs: optional
        Provide additional inputs for different `algo` options

    Returns
    -------
    dict
        A dictionary, with the distance in 'dist', and the path in 'path'
    """

    res = None
    if (algo == 'A*'):
        measure = None
        if ('measure' not in kwargs or kwargs['measure'] not in ['Manhatten', 'Euclidean']):
            warnings.warn("WARNING: Set distance measurement to be default as 'Manhatten")
            measure = 'Manhatten'
        else:
            measure = kwargs['measure']
        res = _distOnGridAStar(column, row, barriers, pt1, pt2, measure)
    else:
        raise UnsupportedInputError("Error: Incorrect or not available grid path finding option!")
    return res

def _distOnGridAStar(column, row, barriers, pt1, pt2, distMeasure):
    # Heuristic measure ==================================================-
    def _calManhattenDist(coord1, coord2):
        return abs(coord1[0] - coord2[0]) + abs(coord1[1] - coord2[1])
    def _calEuclideanDist(coord1, coord2):
        return math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2)

    # Initialize grid ====================================================-
    # Evaluate value f(n) = g(n) + h(n)
    gridStatus = {}
    for col in range(column):
        for ro in range(row):
            if ((col, ro) not in barriers):
                # Content in the dictionary (g(n), h(n), fromCoord)
                # At this stage, no need to calculate h(n) 
                gridStatus[(col, ro)] = (None, None, None)
            else:
                gridStatus[(col, ro)] = 'block'
    if (distMeasure == 'Manhatten'):
        gridStatus[pt1] = (0, _calManhattenDist(pt1, pt2), None)
    elif (distMeasure == 'Euclidean'):
        gridStatus[pt1] = (0, _calEuclideanDist(pt1, pt2), None)
    gridStatus[pt2] = (None, 0, None)

    # Open/close set ======================================================
    openList = [pt1]
    closeList = [i for i in barriers]

    # Find smallest Fn ====================================================
    def _findSmallestFnGrid():
        bestFn = None
        bestCoord = None
        for coord in openList:
            if (gridStatus[coord] != None
                and gridStatus[coord] != 'block' 
                and (bestFn == None or gridStatus[coord][0] + gridStatus[coord][1] < bestFn)):
                bestFn = gridStatus[coord][0] + gridStatus[coord][1]
                bestCoord = coord
        if (bestCoord != None):
            return bestCoord
        else:
            raise

    # For each grid in open set, update g(n) ==============================
    while (len(openList) > 0):
        tmpOpenList = []
        coord = _findSmallestFnGrid()
        # Up
        upCoord = (coord[0], coord[1] + 1)
        if (coord[1] + 1 < row and gridStatus[upCoord] != None and gridStatus[upCoord] != 'block' and upCoord not in closeList):
            if (gridStatus[upCoord][0] == None or gridStatus[upCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[upCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(upCoord, pt2), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[upCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(upCoord, pt2), coord)
                if (upCoord == pt2):
                    break
                else:
                    tmpOpenList.append(upCoord)
        # Down
        downCoord = (coord[0], coord[1] - 1)
        if (coord[1] - 1 >= 0 and gridStatus[downCoord] != None and gridStatus[downCoord] != 'block' and downCoord not in closeList):
            if (gridStatus[downCoord][0] == None or gridStatus[downCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[downCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(downCoord, pt2), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[downCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(downCoord, pt2), coord)
                if (downCoord == pt2):
                    break
                else:
                    tmpOpenList.append(downCoord)
        # Left
        leftCoord = (coord[0] - 1, coord[1])
        if (coord[0] - 1 >= 0 and gridStatus[leftCoord] != None and gridStatus[leftCoord] != 'block' and leftCoord not in closeList):
            if (gridStatus[leftCoord][0] == None or gridStatus[leftCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[leftCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(leftCoord, pt2), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[leftCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(leftCoord, pt2), coord)
                if (leftCoord == pt2):
                    break
                else:
                    tmpOpenList.append(leftCoord)
        # Right
        rightCoord = (coord[0] + 1, coord[1])
        if (coord[0] + 1 < column and gridStatus[rightCoord] != None and gridStatus[rightCoord] != 'block' and rightCoord not in closeList):
            if (gridStatus[rightCoord][0] == None or gridStatus[rightCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[rightCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(rightCoord, pt2), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[rightCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(rightCoord, pt2), coord)
                if (rightCoord == pt2):
                    break
                else:
                    tmpOpenList.append(rightCoord)
        openList.remove(coord)
        openList.extend(tmpOpenList)
        closeList.append(coord)

    # Recover path ========================================================
    path = []
    curCoord = pt2
    finishReconstructFlag = True
    while (finishReconstructFlag):
        finishReconstructFlag = False
        path.insert(0, curCoord)
        curCoord = gridStatus[curCoord][2]
        if (curCoord != None):
            finishReconstructFlag = True
    return {
        'dist': len(path) - 1,
        'path': path
    }

# Path touring through polygons ===============================================
def locSeqRemoveDegen(seq: list[pt], error:float=CONST_EPSILON):
    """
    Given a sequence of points, returns a subset of points that only includes turning points of the sequence. If there are multiple points overlapped at the same location, keeps one of those points.

    Parameters
    ----------
    seq: list[pt], required
        The coordinates of a sequence of points.
    error: float, optional, default as CONST_EPSILON
        Error tolerance.

    Returns
    -------
    dict
        A new dictionary, in the format of 
        
        >>> {
        ...     'newSeq': newSeq, 
        ...     'aggNodeList': aggNodeList, 
        ...     'removedFlag': removedFlag, 
        ...     'locatedSeg': locatedSeg
        ... }

        - The 'newSeq' returns a new sequence which only has turn points of the origin sequence
        - The 'aggNodeList' is a list of lists, if a point overlaps with its previous/next point, the index of both points will be aggregated into the same list.
        - The 'removedFlag' indicates whether a point is removed as a non-turning point, true if the point is not included in the 'newSeq'
        - The 'locatedSeq' returns a list of line segments, for each point removed, it returns the line segment it belongs to 

    Examples
    --------
    For the following inputs
        >>> seq = [[-1, 0], [0, 0], [0, 0], [1, 0], [2, 0], [2, 1], [1, 1], [1, 0], [1, -1]]
        >>> res = locSeqRemoveDegen(seq)
    The result is as follows
        >>> res = {
        ...     'newSeq': [[-1, 0], [2, 0], [2, 1], [1, 1], [1, -1]],
        ...     'aggNodeList': [[0], [1, 2], [3], [4], [5], [6], [7], [8]],
        ...     'removedFlag': [False, True, True, False, False, False, True, False],
        ...     'locatedSeg': [None,
        ...         [[-1, 0], [2, 0]],
        ...         [[-1, 0], [2, 0]],
        ...         None,
        ...         None,
        ...         None,
        ...         [[1, 1], [1, -1]],
        ...         None]
        ... }
    The result shows that, the new sequence is [[-1, 0], [2, 0], [2, 1], [1, 1], [1, -1]], in the new sequence, seq[1], seq[2], seq[6] are not included since they are not turn points.
    seq[1] and seq[2] are aggregated due to overlaps. Although seq[3] and seq[6] are overlapped, they are not aggregated because they are not neighboring. 
    For the removed points, 'locatedSeq' finds the segment they located.

    """

    # Step 1: 先按是否重合对点进行聚合  
    curLocList = [seq[0]]
    curAgg = [0]
    aggNodeList = []
    # 挤香肠算法
    for i in range(1, len(seq)):
        # 如果当前点和任意一个挤出来的点足够近，则计入
        samePtFlag = False

        for pt in curLocList:
            if (is2PtsSame(pt, seq[i], error)):
                curAgg.append(i)
                curLocList.append(seq[i])
                samePtFlag = True
                break

        # 若不重复，了结
        if (not samePtFlag):
            aggNodeList.append([k for k in curAgg])
            curAgg = [i]
            curLocList = [seq[i]]

    aggNodeList.append([k for k in curAgg])

    # Step 2: 对聚合后的点，判断是否为转折点
    # NOTE: removeFlag的长度和aggNodeList一致
    removedFlag = [False]
    for i in range(1, len(aggNodeList) - 1):
        preLoc = seq[aggNodeList[i - 1] if type(aggNodeList[i - 1]) != list else aggNodeList[i - 1][0]]
        curLoc = seq[aggNodeList[i] if type(aggNodeList[i]) != list else aggNodeList[i][0]]
        sucLoc = seq[aggNodeList[i + 1] if type(aggNodeList[i + 1]) != list else aggNodeList[i + 1][0]]

        dev = distPt2Seg(curLoc, [preLoc, sucLoc])
        if (dev <= error):
            removedFlag.append(True)
        else:
            removedFlag.append(False)
    removedFlag.append(False)

    # 得到去掉共线和重合点后的折线
    newSeq = []
    for i in range(len(aggNodeList)):
        if (removedFlag[i] == False):
            if (type(aggNodeList[i]) == list):
                newSeq.append(seq[aggNodeList[i][0]])
            else:
                newSeq.append(seq[aggNodeList[i]])

    # 对于被移除的共线点，找到其所在的线段
    locatedSeg = []
    # 把没有移除的点的序号记一下
    seqPre = []
    seqSuc = []
    # 查找移除点之前和之后一个removeFlag为False的对应aggNode，得到对应线段
    for i in range(len(removedFlag)):
        # Log the prev that is not removed
        if (removedFlag[i] == False):
            seqPre.append(i)
        else:
            seqPre.append(seqPre[-1])
        # Log the next that is not removed
        if (removedFlag[len(removedFlag) - 1 - i] == False):
            seqSuc.insert(0, len(removedFlag) - 1 - i)
        else:
            seqSuc.insert(0, seqSuc[0])

    for i in range(len(removedFlag)):
        if (removedFlag[i] == False):
            locatedSeg.append(None)
            para = []
        else:
            startLoc = None
            if (type(aggNodeList[seqPre[i]]) == list):
                startLoc = seq[aggNodeList[seqPre[i]][0]]
            else:
                startLoc = seq[aggNodeList[seqPre[i]]]
            endLoc = None
            if (type(aggNodeList[seqSuc[i]]) == list):
                endLoc = seq[aggNodeList[seqSuc[i]][0]]
            else:
                endLoc = seq[aggNodeList[seqSuc[i]]]
            locatedSeg.append([startLoc, endLoc])

    return {
        'newSeq': newSeq,
        'aggNodeList': aggNodeList,
        'removedFlag': removedFlag,
        'locatedSeg': locatedSeg
    }

