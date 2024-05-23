import geopy.distance
import heapq
import math
import tripy

import shapely
from shapely.geometry import mapping
from shapely.ops import nearest_points
import networkx as nx

from .common import *
from .msg import *
from .ds import *

# History =====================================================================
# 20231219 - A stable version from the past
# 20240117 - Visible area visFromPt()
# =============================================================================

# Point versus Objects ========================================================
def is2PtsSame(pt1: pt, pt2: pt, error: float = CONST_EPSILON) -> bool:
    """Are two points at the 'same' location"""
    if (abs(pt1[0] - pt2[0]) >= error):
        return False
    if (abs(pt1[1] - pt2[1]) >= error):
        return False
    return True

def is3PtsClockWise(pt1: pt, pt2: pt, pt3: pt, error: float = CONST_EPSILON) -> bool | None:
    """Are three given pts in a clock-wise order, None as they are collinear"""
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
    """Is a pt on the line."""
    if (is2PtsSame(line[0], line[1], error = error)):
        raise ZeroVectorError()
    if (is3PtsClockWise(pt, line[0], line[1], error = error) == None):
        return True
    else:
        return False

def isPtOnSeg(pt: pt, seg: line, interiorOnly: bool=False, error:float = CONST_EPSILON) -> bool:
    """Is a pt on the segment"""
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
    """Is a pt on the ray, could be at the end."""
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

# Point to line-shape =========================================================
def rayPerp2Line(pt: pt, line: line) -> line:
    """Given a point and a line, return a ray from that point and perpendicular to the givne line"""
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

def ptFoot2Line(pt: pt, line: line) -> pt:
    """Given a point and a line, return the foot of that point on the line"""
    if (isPtOnLine(pt, line)):
        return tuple(pt)
    else:
        ray = rayPerp2Line(pt, line)
        return ray[1]

# Line-shape intersection =====================================================
def intLine2Line(line1: line, line2: line) -> dict:
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
    return intLine2Seg(line, seg)

def intSeg2Seg(seg1: line, seg2: line) -> dict:
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
    return intLine2Ray(line, ray)

def intRay2Seg(ray: line, seg: line) -> dict:
    return intSeg2Ray(seg, ray)

def intRay2Ray(ray1: line, ray2: line) -> dict:
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

def crossSeg2Seg(seg1: line, seg2: line) -> list[list[pt]] | None:
    """Given two line segments return the line segments that are crossed by each other"""
    intSeg = intSeg2Seg(seg1, seg2)
    if (intSeg['status'] == 'Cross'):
        segs = []
        if (intSeg['interiorFlag'] == True):
            segs.append([seg1[0], intSeg['intersect']])
            segs.append([intSeg['intersect'], seg1[1]])
            segs.append([seg2[0], intSeg['intersect']])
            segs.append([intSeg['intersect'], seg2[1]])
        else:
            if (is2PtsSame(seg1[0], intSeg['intersect'])):
                segs.append(seg1)
            elif (is2PtsSame(seg1[1], intSeg['intersect'])):
                segs.append(seg1)
            else:
                segs.append([seg1[0], intSeg['intersect']])
                segs.append([intSeg['intersect'], seg1[1]])
            if (is2PtsSame(seg2[0], intSeg['intersect'])):
                segs.append(seg2)
            elif (is2PtsSame(seg2[1], intSeg['intersect'])):
                segs.append(seg2)
            else:
                segs.append([seg2[0], intSeg['intersect']])
                segs.append([intSeg['intersect'], seg2[1]])
        return segs
    else:
        return [seg1, seg2]

# Line-shape versus Line-shape ===================================================
def isLineIntLine(line1: line, line2: line) -> bool:
    """Is two line intersect with each other"""
    intPt = intLine2Line(line1, line2)
    if (intPt['intersect'] == None):
        return False
    else:
        return True

def isLineIntSeg(line: line, seg: line, interiorOnly:bool=False) -> bool:
    """Is a line intersect with a segment"""
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
    """Is a line intersect with a ray"""
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
    """Is a segment intersect with a line"""
    return isLineIntSeg(line, seg, interiorOnly)

def isSegIntSeg(seg1: line, seg2: line, interiorOnly: bool=False) -> bool:
    """Are two segment intersect with each other, could be at the end of a segment"""
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
    """Is a segment intersect with a ray, could be at the end of the segment/ray"""
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
    # If any end is inside bounding box
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
    """Is a ray intersect with a line"""
    return isLineIntRay(line, ray, interiorOnly)

def isRayIntSeg(ray: line, seg: line, interiorOnly: bool=False) -> bool:
    """Is a ray intersect with a segment"""
    return isSegIntRay(seg, ray, interiorOnly)

def isRayIntRay(ray1: line, ray2: line, interiorOnly: bool=False) -> bool:
    """Is a segment intersect with a ray, could be at the end of the segment/ray"""
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
    # Sanity check ============================================================
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: `poly` and `polyShapely` cannot be None at the same time.")

    # Projct points to the line ===============================================
    projPts = []    
    for pt in poly:
        projPts.append(ptFoot2Line(pt, line))

    # Find two pts on the line that are the farthest away from each other =====
    projPts.sort()
    seg = [projPts[0], projPts[-1]]

    return intSeg2Poly(seg, poly, polyShapely, returnShaplelyObj)

def intSeg2Poly(seg: line, poly: poly=None, polyShapely: shapely.Polygon=None, returnShaplelyObj: bool=False) -> dict | list[dict] | shapely.Point | shapely.Polygon | shapely.GeometryCollection:
    # WARNING: Results are not reliable on both ends of the seg.
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
            for pt in polys:
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
    # Convert ray to a proper segment
    # Sanity check ============================================================
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: `poly` and `polyShapely` cannot be None at the same time.")

    # Projct points to the line ===============================================
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
    # Sanity check ============================================================
    if (seq == None and seqShapely == None):
        raise MissingParameterError("ERROR: `seq` and `seqShapely` cannot be None at the same time.")
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: `poly` and `polyShapely` cannot be None at the same time.")

    # get shapely objects =====================================================
    if (seqShapely == None):
        seqShapely = shapely.LineString(seq)
    if (polyShapely == None):
        polyShapely = shapely.Polygon(poly)
    intShape = shapely.intersection(seqShapely, polyShapely)

    # If return shapely objects no processing needed ==========================
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
    """Is a line intersect with a polygon"""
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
    # WARNING: results are not reliable if `interiorFlag` == False.
    
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
    """Is a ray intersect with a polygon"""
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
    # MultiPoint/MultiLineString/MultiPolygon/Geometrycollection
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
    area = calTriangleAreaXY(pt, line[0], line[1])
    a = distEuclideanXY(line[0], line[1])['dist']
    h = 2 * area / a
    return h

def distPt2Seg(pt: pt, seg: line) -> float:
    foot = ptFoot2Line(pt, seg)
    if (isPtOnSeg(foot, seg)):
        return distEuclideanXY(pt, foot)['dist']
    else:
        return min(distEuclideanXY(pt, seg[0])['dist'], distEuclideanXY(pt, seg[1])['dist'])

def distPt2Ray(pt: pt, ray: line) -> float:
    foot = ptFoot2Line(pt, ray)
    if (isPtOnRay(foot)):
        return distEuclideanXY(pt, foot)
    else:
        return distEuclideanXY(pt, ray[0])

def distPt2Seq(pt: pt, seq: list[pt]) -> float:
    # FIXME: stupid way, needs improvement
    if (len(seq) == 2):
        return distPt2Seg(pt, seq)

    dist2Seg = []
    for p in seq:
        dist2Seg.append(distEuclideanXY(pt, p))
    minIndex = dist2Seg.index(min(dist2Seg))
    if (minIndex == 0):
        return distPt2Seg(pt, [seq[0], seq[1]])
    elif (minIndex == len(dist2Seg) - 1):
        return distPt2Seg(pt, [seq[-2], seq[-1]])
    else:
        return min(distPt2Seg(pt, [seq[minIndex], seq[minIndex + 1]]),
                   dist2Seg(pt, [seq[minIndex], seq[minIndex - 1]]))

def distPt2Poly(pt: pt, poly: poly=None, polyShapely: shapely.Polygon=None) -> float:
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: Missing required field `poly` or `polyShapely`")

    if (polyShapely == None):
        polyShapely = shapely.Polygon([p for p in poly])
    ptShapely = shapely.Point(pt)
    return shapely.distance(ptShapely, polyShapely)

def nearestPtLine2Poly(line: line, poly: poly=None, polyShapely: shapely.Polygon=None) -> dict:
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

# Dimension mapping ===========================================================
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

# Polys =======================================================================
def polysUnion(polys:polys=None, polysShapely:list[shapely.Polygon]=None, returnShaplelyObj:bool=False) -> list:
    """Given a list of polygons which could be intersecting to each other, return unioned polygons that are not intersecting"""
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

    Warning!!!
    ----------
    This function needs to be rewritten. It's not trackable now

    Parameters
    ----------

    polys: dictionary, required
        The polys dictionary with neighborhood.
    order: int, optional, defalut None
        Maximum order of Steiner zone

    Returns
    -------

    list[dict]
        A list of Steiner zone dictionaris, each in the following format::
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
 
def polysBoundingBox(polys:polys=None, polysShapely:list[shapely.Polygon]=None):
    if (polys == None and polysShapely == None):
        raise MissingParameterError("ERROR: Missing required field 'polys' or 'polysShapely'.")
    if (polys == None):
        polys = []
        for p in polyShapely:
            polys.append([i for i in mapping(p)['coordinates'][0]])
    polysBox = []
    for p in polys:
        x = []
        y = []
        for pt in p:
            x.append(pt[0])
            y.append(pt[1])
        polysBox.append([min(x), max(x), min(y), max(y)])

    return polysBox

def polygonsAlongLocSeq(seq, polygons:dict, polyFieldName = 'poly'):
    if (polygons == None):
        raise MissingParameterError("ERROR: Missing required field 'polygons'.")

    # polysBox = polysBoundingBox(polys, polysShapely)

    # First, for each leg in the seq, find the individual polygons intersect with the leg
    actions = {}
    actionIDOnLeg = 0
    for i in range(len(seq) - 1):
        seg = [seq[i], seq[i + 1]]

        # NOTE: 准备用segment tree
        # NOTE: For now use the naive way - checking the boundingbox
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

def ptPolyCenter(poly: poly=None, polyShapely: shapely.Polygon=None) -> pt:
    if (poly == None and polyShapely == None):
        raise MissingParameterError("ERROR: Missing required field 'poly' or 'polyShapely'.")
    if (polyShapely == None):
        polyShapely = shapely.Polygon(poly)

    ptShapely = shapely.centroid(polyShapely)
    center = (ptShapely.x, ptShapely.y)
    return center

# Visibility check ============================================================
def polysVisibleGraph(polys:polys) -> dict:
    vg = {}
    for p in range(len(polys)):
        for e in range(len(polys[p])):
            vg[(p, e)] = {'loc': polys[p][e], 'visible': []}
            W = visPtAmongPolys((p, e), polys, knownVG=vg)
            for w in W:
                vg[(p, e)]['visible'].append(w)
    return vg

def visPtAmongPolys(v:int|str|tuple, polys:polys, standalonePts:dict|None=None, knownVG:dict={}) -> list:
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

def visPolyFromPt(v:int|str|tuple, range:float, polys:polys, standalonePts:dict|None=None, knownVG:dict={}, lod:int = 36) -> list:
    
    # Step 1: fi


    return poly

def circleShadowByPolys(circle: dict = None, polys: polys = None, lod: int = 30):
    standalonePts = {'c': {'loc': circle['center']}}
    cirR = circle['radius']

    remain = [[
        circle['center'][0] + circle['radius'] * math.sin(2 * d * math.pi / lod),
        circle['center'][1] + circle['radius'] * math.cos(2 * d * math.pi / lod),
    ] for d in range(lod + 1)]

    # 分别计算每个poly的阴影遮挡
    for poly in polys:
        # 首先计算poly从

        visExpt = visPtAmongPolys('c', [poly], standalonePts)



    visExpt = visPtAmongPolys('c', polys, standalonePts)
    visLocs = []
    for v in visExpt:
        visLocs.append(polys[v[0]][v[1]])
    return visCircle

def polyShadowByPolys(poly: poly = None, anchor: pt = None, polys: polys = None):
    standalonePts = {'c': {'loc': circle['center']}}

    # FIXME: 不得已的办法，把polygons拿去triangulate，然后逐个计算遮挡阴影

    # 分别计算每个poly的遮挡
    for p in polys:
        # 首先计算有几个遮挡的转折点，若只有两个，则可以简单的划分为阴面和阳面
        pass

    return visPoly

# Time seq related ============================================================
def snapInTimedSeq(timedSeq: list[tuple[pt, float]], t: float) -> pt:
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

def locSeqMileage(seq: list[pt], dist: int|float, dimension: str = 'XY') -> pt:
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

# Area calculation ============================================================
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

def calPolyAreaXY(poly: poly) -> float:
    lstTriangle = tripy.earclip(poly)

    # Weight them and make draws ==============================================
    area = 0
    for i in range(len(lstTriangle)):
        area += calTriangleAreaXY(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2])
    return area

def calPolyAreaLatLon(polyLatLon: poly) -> float:
    """Returns the area surrounded by polyLatLon on the Earth"""

    # Additional packages =====================================================
    try:
        from pyproj import Geod
    except:
        raise ImportError

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

def calPolyPerimeterXY(poly: poly) -> float:
    p = 0
    for i in range(-1, len(poly) - 1):
        p += distEuclideanXY(poly[i], poly[i + 1])['dist']
    return p

# Location of points ==========================================================
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

def circleByCenterLatLon(center: pt, radius: int|float, lod: int = 30):
    circle = []
    for i in range(lod):
        deg = float(360 * i) / float(lod)
        pt = ptInDistLatLon(pt = center, direction = deg, distMeters = radius)
        circle.append(pt)
    return circle

def circleByCenterXY(center: pt, radius: int|float, lod: int = 30):
    circle = []
    for i in range(lod):
        deg = float(360 * i) / float(lod)
        pt = ptInDistXY(pt = center, direction = deg, distMeters = radius)
        circle.append(pt)
    return circle

# Sort nodes ==================================================================
def nodeSeqByDist(nodes: dict, nodeIDs: list|str = 'All', locFieldName = 'loc', edges: dict = {'method': 'Euclidean'}, refLoc: pt|None = None, refNodeID: int|str|None = None) -> list:
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
        dist = scaleDist(loc1 = refLoc, loc2 = nodes[n][locFieldName], edges = edges)['dist']
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

def nodesInIsochrone(nodes: dict, nodeIDs: list|str = 'All', locFieldName = 'loc', edges: dict = {'method': 'Euclidean'}, refLoc: pt|None = None, refNodeID: int|str|None = None, isoRange: float = None, sortFlag: bool = False) -> list: 
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
            dist = scaleDist(loc1 = refLoc, loc2 = nodes[n][locFieldName], edges = edges)['dist']
            if (dist <= isoRange):
                heapq.heappush(nearSetHeap, (dist, n))
        while (len(nearSetHeap) > 0):
            nearSet.append(heapq.heappop(nearSetHeap)[1])
        nearest = nearSet[0]
    else:
        for n in nodeIDs:
            dist = scaleDist(loc1 = refLoc, loc2 = nodes[n][locFieldName], edges = edges)['dist']
            if (dist <= isoRange):
                nearSet.append(n)
            if (dist <= nearestDist):
                nearest = n
                nearestDist = dist
    return {
        'nearSet': nearSet,
        'nearest': nearest
    }

# Create distance matrix ======================================================
def matrixDist(nodes: dict, edges: dict = {'method': 'Euclidean'}, nodeIDs: list|str = 'All', locFieldName: str = 'loc') -> dict:
    # Define tau
    tau = {}
    pathLoc = {}

    if (type(edges) != dict or 'method' not in edges):
        raise MissingParameterError(ERROR_MISSING_EDGES)
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    if (edges['method'] == 'Euclidean'):
        speed = 1 if 'speed' not in edges else edges['speed']
        tau, pathLoc = _matrixDistEuclideanXY(
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            speed = speed, 
            locFieldName = locFieldName)
    elif (edges['method'] == 'EuclideanBarrier'):
        if ('polys' not in edges or edges['polys'] == None):
            warnings.warning("WARNING: No barrier provided.")
            tau, pathLoc = _matrixDistEuclideanXY(
                nodes = nodes, 
                nodeIDs = nodeIDs, 
                locFieldName = locFieldName)
        else:
            tau, pathLoc = _matrixDistBtwPolysXY(
                nodes = nodes, 
                nodeIDs = nodeIDs, 
                polys = edges['polys'], 
                locFieldName = locFieldName)
    elif (edges['method'] == 'LatLon'):
        speed = 1 if 'speed' not in edges else edges['speed']
        tau, pathLoc = _matrixDistLatLon(
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            speed=speed, 
            locFieldName = locFieldName)
    elif (edges['method'] == 'Manhatten'):
        speed = 1 if 'speed' not in edges else edges['speed']
        tau, pathLoc = _matrixDistManhattenXY(
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            speed = speed, 
            locFieldName = locFieldName)
    elif (edges['method'] == 'Dictionary'):
        if ('tau' not in edges or edges['tau'] == None):
            raise MissingParameterError("'tau' is not specified")
        for p in edges['tau']:
            speed = 1 if 'speed' not in edges else edges['speed']
            tau[p] = edges['tau'][p] / speed
            pathLoc[p] = edges['path'][p] if 'path' in edges else [nodes[p[0]][locFieldName], nodes[p[1]][locFieldName]]
    elif (edges['method'] == 'Grid'):
        if ('grid' not in edges or edges['grid'] == None):
            raise MissingParameterError("'grid' is not specified")
        if ('column' not in edges['grid'] or 'row' not in edges['grid']):
            raise MissingParameterError("'column' and 'row' need to be specified in 'grid'")
        tau, pathLoc = _matrixDistGrid(
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            grids = edges['grid'], 
            locFieldName = locFieldName)
    else:
        raise UnsupportedInputError(ERROR_MISSING_EDGES)        

    return tau, pathLoc

def _matrixDistEuclideanXY(nodes: dict, nodeIDs: list, speed = 1, locFieldName = 'loc'):
    tau = {}
    pathLoc = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                d = distEuclideanXY(nodes[i][locFieldName], nodes[j][locFieldName])
                tau[i, j] = d['dist'] / speed
                tau[j, i] = d['dist'] / speed
                pathLoc[i, j] = [nodes[i][locFieldName], nodes[j][locFieldName]]
                pathLoc[j, i] = [nodes[j][locFieldName], nodes[i][locFieldName]]
            else:
                tau[i, j] = CONST_EPSILON
                tau[j, i] = CONST_EPSILON
                pathLoc[i, j] = []
                pathLoc[j, i] = []
    return tau, pathLoc

def _matrixDistManhattenXY(nodes: dict, nodeIDs: list, speed = 1, locFieldName = 'loc'):
    tau = {}
    pathLoc = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                d = distManhattenXY(nodes[i][locFieldName], nodes[j][locFieldName])
                tau[i, j] = d['dist'] / speed
                tau[j, i] = d['dist'] / speed
                pathLoc[i, j] = d['path']
                pathLoc[j, i] = [d['path'][len(d['path']) - 1 - i] for i in range(len(d['path']))]
            else:
                tau[i, j] = CONST_EPSILON
                tau[j, i] = CONST_EPSILON
                pathLoc[i, j] = []
                pathLoc[j, i] = []
    return tau, pathLoc

def _matrixDistLatLon(nodes: dict, nodeIDs: list, distUnit = 'meter', speed = 1, locFieldName = 'loc'):
    tau = {}
    pathLoc = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                d = distLatLon(nodes[i][locFieldName], nodes[j][locFieldName], distUnit)
                tau[i, j] = d['dist'] / speed
                tau[j, i] = d['dist'] / speed
                pathLoc[i, j] = [nodes[i][locFieldName], nodes[j][locFieldName]]
                pathLoc[j, i] = [nodes[j][locFieldName], nodes[i][locFieldName]]
            else:
                tau[i, j] = CONST_EPSILON
                tau[j, i] = CONST_EPSILON
                pathLoc[i, j] = []
                pathLoc[j, i] = []
    return tau, pathLoc

def _matrixDistGrid(nodes: dict, nodeIDs: list, grid: dict, locFieldName = 'loc'):
    tau = {}
    pathLoc = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                d = distOnGrid(pt1 = nodes[i][locFieldName], pt2 = nodes[j][locFieldName], grid = grid)
                tau[i, j] = d['dist']
                tau[j, i] = d['dist']
                pathLoc[i, j] = d['path']
                pathLoc[j, i] = [d['path'][len(d['path']) - 1 - i] for i in range(len(d['path']))]
            else:
                tau[i, j] = CONST_EPSILON
                tau[j, i] = CONST_EPSILON
                pathLoc[i, j] = []
                pathLoc[j, i] = []
    return tau, pathLoc

def _matrixDistBtwPolysXY(nodes: dict, nodeIDs: list, polys: polys, polyVG = None, locFieldName = 'loc'):
    tau = {}
    pathLoc = {}
    
    if (polyVG == None):
        polyVG = polysVisibleGraph(polys)

    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                d = distBtwPolysXY(pt1 = nodes[i][locFieldName], pt2 = nodes[j][locFieldName], polys = polys, polyVG = polyVG)
                tau[i, j] = d['dist']
                tau[j, i] = d['dist']
                pathLoc[i, j] = d['path']
                pathLoc[j, i] = [d['path'][len(d['path']) - 1 - i] for i in range(len(d['path']))]
            else:
                tau[i, j] = CONST_EPSILON
                tau[j, i] = CONST_EPSILON
                pathLoc[i, j] = []
                pathLoc[j, i] = []
    return tau, pathLoc

def vectorDist(loc: pt, nodes: dict, edges: dict = {'method': 'Euclidean'}, nodeIDs: list|str = 'All', locFieldName: str = 'loc') -> dict:
    # Define tau
    tau = {}
    revTau = {}
    pathLoc = {}
    revPathLoc = {}

    if (type(edges) != dict or 'method' not in edges):
        raise MissingParameterError(ERROR_MISSING_EDGES)
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    if (edges['method'] == 'Euclidean'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        tau, revTau, pathLoc, revPath = _vectorDistEuclideanXY(
            loc = loc,
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            ratio = ratio, 
            locFieldName = locFieldName)
    elif (edges['method'] == 'EuclideanBarrier'):
        if ('polys' not in edges or edges['polys'] == None):
            warnings.warning("WARNING: No barrier provided.")
            tau, revTau, pathLoc, revPath = _vectorDistEuclideanXY(
                loc = loc,
                nodes = nodes, 
                nodeIDs = nodeIDs, 
                locFieldName = locFieldName)
        else:
            tau, revTau, pathLoc, revPath = _vectorDistBtwPolysXY(
                loc = loc,
                nodes = nodes, 
                nodeIDs = nodeIDs, 
                polys = edges['polys'], 
                locFieldName = locFieldName)
    elif (edges['method'] == 'LatLon'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        tau, revTau, pathLoc, revPath = _vectorDistLatLon(
            loc = loc,
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            speed=ratio, 
            locFieldName = locFieldName)
    elif (edges['method'] == 'Manhatten'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        tau, revTau, pathLoc, revPath = _vectorDistManhattenXY(
            loc = loc,
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            ratio = ratio, 
            locFieldName = locFieldName)
    elif (edges['method'] == 'Grid'):
        if ('grid' not in edges or edges['grid'] == None):
            raise MissingParameterError("'grid' is not specified")
        if ('column' not in edges['grid'] or 'row' not in edges['grid']):
            raise MissingParameterError("'column' and 'row' need to be specified in 'grid'")
        tau, revTau, pathLoc, revPath = _vectorDistGrid(
            loc = loc,
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            grids = edges['grid'], 
            locFieldName = locFieldName)
    else:
        raise UnsupportedInputError(ERROR_MISSING_EDGES)        

    return tau, revTau, pathLoc, revPathLoc

def _vectorDistEuclideanXY(loc: pt, nodes: dict, nodeIDs: list, speed = 1, locFieldName = 'loc'):
    tau = {}
    revTau = {}
    pathLoc = {}
    revPathLoc = {}
    for i in nodeIDs:
        d = distEuclideanXY(loc, nodes[i][locFieldName])
        tau[i] = d['dist'] / speed
        revTau[i] = d['dist'] / speed
        pathLoc[i] = [loc, nodes[i][locFieldName]]
        revPathLoc[i] = [nodes[i][locFieldName], loc]
    return tau, revTau, pathLoc, revPathLoc

def _vectorDistManhattenXY(loc: pt, nodes: dict, nodeIDs: list, speed = 1, locFieldName = 'loc'):
    tau = {}
    revTau = {}
    pathLoc = {}
    revPathLoc = {}
    for i in nodeIDs:
        d = distManhattenXY(loc, nodes[i][locFieldName])
        tau[i] = d['dist'] / speed
        revTau[i] = d['dist'] / speed
        pathLoc[i] = d['path']
        revPathLoc[i] = [d['path'][len(d['path']) - 1 - i] for i in range(len(d['path']))]
    return tau, revTau, pathLoc, revPathLoc

def _vectorDistLatLon(loc: pt, nodes: dict, nodeIDs: list, distUnit = 'meter', speed = 1, locFieldName = 'loc'):
    tau = {}
    revTau = {}
    pathLoc = {}
    revPathLoc = {}
    for i in nodeIDs:
        d = distLatLon(loc, nodes[i][locFieldName], distUnit)
        tau[i] = d['dist'] / speed
        revTau[i] = d['dist'] / speed
        pathLoc[i] = [loc, nodes[i][locFieldName]]
        revPathLoc[i] = [nodes[i][locFieldName], loc]
    return tau, revTau, pathLoc, revPathLoc

def _vectorDistGrid(loc: pt, nodes: dict, nodeIDs: list, grid: dict, locFieldName = 'loc'):
    tau = {}
    revTau = {}
    pathLoc = {}
    revPathLoc = {}
    for i in nodeIDs:
        d = distOnGrid(pt1 = loc, pt2 = nodes[i][locFieldName], grid = grid)
        tau[i] = d['dist'] / speed
        revTau[i] = d['dist'] / speed
        pathLoc[i] = d['path']
        revPathLoc[i] = [d['path'][len(d['path']) - 1 - i] for i in range(len(d['path']))]
    return tau, revTau, pathLoc, revPathLoc

def _vectorDistBtwPolysXY(loc: pt, nodes: dict, nodeIDs: list, polys: polys, polyVG = None, locFieldName = 'loc'):
    tau = {}
    revTau = {}
    pathLoc = {}
    revPathLoc = {}

    if (polyVG == None):
        polyVG = polysVisibleGraph(polys)

    for i in nodeIDs:
        d = distBtwPolysXY(pt1 = loc, pt2 = nodes[i][locFieldName], polys = polys, polyVG = polyVG)
        tau[i] = d['dist'] / speed
        revTau[i] = d['dist'] / speed
        pathLoc[i] = d['path']
        revPathLoc[i] = [d['path'][len(d['path']) - 1 - i] for i in range(len(d['path']))]
    return tau, revTau, pathLoc, revPathLoc

def scaleDist(loc1: pt, loc2: pt, edges: dict = {'method': 'Euclidean'}) -> dict:
    # Define tau
    dist = None
    revDist = None
    pathLoc = []
    revPathLoc = []

    if (type(edges) != dict or 'method' not in edges):
        raise MissingParameterError(ERROR_MISSING_EDGES)

    if (edges['method'] == 'Euclidean'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        dist = distEuclideanXY(loc1, loc2)['dist']
        revDist = dist
        pathLoc = [loc1, loc2]
        revPathLoc = [loc2, loc1]
    elif (edges['method'] == 'EuclideanBarrier'):
        if ('polys' not in edges or edges['polys'] == None):
            warnings.warning("WARNING: No barrier provided.")
            dist = distEuclideanXY(loc1, loc2)['dist']
            revDist = dist
            pathLoc = [loc1, loc2]
            revPathLoc = [loc2, loc1]
        else:
            res = distBtwPolysXY(pt1, pt2, edges['polys'])
            dist = res['dist']
            revDist = dist
            pathLoc = [i for i in res['path']]
            revPathLoc = [pathLoc[len(pathLoc) - i - 1] for i in range(len(pathLoc))]
    elif (edges['method'] == 'LatLon'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        dist = distLatLon(loc1, loc2)['dist']
        revDist = dist
        pathLoc = [loc1, loc2]
        revPathLoc = [loc2, loc1]
    elif (edges['method'] == 'Manhatten'):
        ratio = 1 if 'ratio' not in edges else edges['ratio']
        dist = distManhattenXY(loc1, loc2)['dist']
        revDist = dist
        pathLoc = [loc1, (pt1[0], pt2[1]), loc2]
        revPathLoc = [loc2, (pt1[0], pt2[1]), loc1]
    elif (edges['method'] == 'Grid'):
        if ('grid' not in edges or edges['grid'] == None):
            raise MissingParameterError("'grid' is not specified")
        if ('column' not in edges['grid'] or 'row' not in edges['grid']):
            raise MissingParameterError("'column' and 'row' need to be specified in 'grid'")
        res = distOnGrid(loc1, loc2, edges['grid'])
        dist = res['dist']
        revDist = dist
        pathLoc = [i for i in res['path']]
        revPathLoc = [pathLoc[len(pathLoc) - i - 1] for i in range(len(pathLoc))]
    else:
        raise UnsupportedInputError(ERROR_MISSING_EDGES)        

    return {
        'dist': dist,
        'revDist': revDist,
        'pathLoc': pathLoc,
        'revPathLoc': revPathLoc
    }

# Distance calculation ========================================================
def distEuclideanXY(pt1: pt, pt2: pt) -> dict:
    """Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number"""
    return {
        'dist': math.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2),
        'path': [pt1, pt2]
    }

def distManhattenXY(pt1: pt, pt2: pt) -> dict:
    """Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number"""
    return {
        'dist': abs(pt1[0] - pt2[0]) + abs(pt1[1] - pt2[1]),
        'path': [pt1, (pt1[0], pt2[1]), pt2]
    }

def distBtwPolysXY(pt1:pt, pt2:pt, polys:polys, polyVG:dict=None) -> dict:
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
    Ws = visPtAmongPolys('s', polys, {'s': {'loc': pt1, 'visible': []}})
    vertices['s']['visible'] = Ws
    vertices['e'] = {'loc': pt2, 'visible': []}
    We = visPtAmongPolys('e', polys, {'e': {'loc': pt2, 'visible': []}})
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
    return {
        'dist': 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a)),
        'path': [pt1, pt2]
    }

def distOnGrid(pt1: pt, pt2: pt, grid: dict, algo: dict = {'method': 'A*', 'measure': 'Manhatten'}) -> dict:
    """Given two coordinates on the grid, finds the 'shortest' path to travel

    Parameters
    ----------

    grid: dictionary, required, default as None
        The environment of a grid area, in the following format:
            >>> grid = {
            ...     'column': col, # Number of columns,
            ...     'row': row, # Number of rows,
            ...     'barriers': barriers, # A list of coordinates,
            ... }
    pt1: 2-tuple|2-list, required
        Starting location on the grid
    pt2: 2-tuple|2-list, required
        Ending location on the grid 
    algo: dictionary, required, default as {'method': 'A*', 'distMeasure': 'Manhatten'}
        The algorithm configuration. For example
        1) A*
            >>> algo = {
            ...     'method': A*,
            ...     'measure': 'Manhatten', # Options: 'Manhatten', 'Euclidean'
            ... }

    Returns
    -------
    dictionary
        A path on the given grid, in the following formatt::
            >>> res = {
            ...     'dist': dist,
            ...     'path': path,
            ... }

    """

    # Decode ==================================================================
    column = grid['column']
    row = grid['row']
    barriers = grid['barriers']
    res = None

    # Call path finding =======================================================
    if (algo['method'] == 'A*'):
        if ('measure' not in algo or algo['measure'] not in ['Manhatten', 'Euclidean']):
            warnings.warn("WARNING: Set distance measurement to be default as 'Manhatten")
        res = _distOnGridAStar(column, row, barriers, pt1, pt2, algo['measure'])
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
