import math

from .const import *
from .msg import *

# Description =================================================================
# This script is comparing the relative relation between basic geometry objects,
#     including: pt (point), line, seg (line segment), ray, and poly (polygon).
# Notice:
# 1. All the geometry objects are in Euclidean space
# 2. The data structure of the geometry objects are all lists (or tuples) and list of lists
# 3. "On" means the object (point) is in the closure of another object
#    "Inside" means the object (point) is in the interior of another object
#    "Cross" means the intersection of two objects needs to be in the interior of both objects
#    "Int" (Intersect) means the intersection, which could be in the interior of either objects.

def is2PtsSame(
    pt1:        "2-tuple of coordinate (x, y)",
    pt2:        "2-tuple of coordinate (x, y)"
    ) -> "Check if two points are the same":
    # Check if two points are very close to each other ========================
    if (abs(pt1[0] - pt2[0]) >= CONST_EPSILON):
        return False
    if (abs(pt1[1] - pt2[1]) >= CONST_EPSILON):
        return False
    return True

def isPtOnLine(
    pt:         "2-tuple of coordinate (x, y)",
    line:       "List of 2 pts which defines a line"
    ) -> "Return True if the point is on the line, with error less than CONST_EPSILON":
    # Validation ==============================================================
    if (is2PtsSame(line[0], line[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Calculate the distance between pt and the line ==========================
    if (is3PtsClockWise(pt, line[0], line[1]) == None):
        return True
    else:
        return False

def isPtOnSeg(
    pt:         "2-tuple of coordinate (x, y)", 
    seg:        "List of 2 pts as two ends"
    ) -> "Return true if the point is on the line segment, including initial pts, with error less than CONST_EPSILON":
    # Check if is on the line seg =============================================
    onLine = isPtOnLine(pt, seg)
    if (onLine != True):
        return onLine

    # Get pts =================================================================
    [x1, y1] = [seg[0][0], seg[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [seg[1][0], seg[1][1]]

    # Relative location =======================================================
    return (x2 >= min(x1, x3) 
        and x2 <= max(x1, x3) 
        and y2 >= min(y1, y3) 
        and y2 <= max(y1, y3))

def isPtInsideSeg(
    pt:         "2-tuple of coordinate (x, y)", 
    seg:        "List of 2 pts as two ends"
    ) -> "Return true if the point is inside the line segment, not including initial pts, with error less than CONST_EPSILON":
    # Check if is on the line seg =============================================
    onSeg = isPtOnSeg(pt, seg)
    if (onSeg != True):
        return onSeg

    # Check if is at any one the ends =========================================
    if (is2PtsSame(pt, seg[0])):
        return False
    if (is2PtsSame(pt, seg[1])):
        return False

    return True

def isPtOnRay(
    pt:         "2-tuple of coordinate (x, y)",
    ray:        "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return true if the point is on the ray, with error less than CONST_EPSILON":
    # Check if is on the line seg =============================================
    onLine = isPtOnLine(pt, ray)
    if (onLine != True):
        return onLine

    # Get pts =================================================================
    [x1, y1] = [ray[0][0], ray[0][1]]
    [x2, y2] = [pt[0], pt[1]]
    [x3, y3] = [ray[1][0], ray[1][1]]

    # Relative location =======================================================
    return (
        x2 >= min(x1, x3) 
        and x2 <= max(x1, x3) 
        and y2 >= min(y1, y3) 
        and y2 <= max(y1, y3)
    ) or (
        x3 >= min(x1, x2) 
        and x3 <= max(x1, x2) 
        and y3 >= min(y1, y2) 
        and y3 <= max(y1, y2)
    )

def isPtOnPolyEdge(
    pt:         "2-tuple of coordinate (x, y)",
    poly:       "List of pts, form a close area"
    ) -> "Return true if the point is on the edge of polygon, with error less than CONST_EPSILON":
    # Check if the pt is on any of the edge segment ===========================
    for i in range(-1, len(poly) - 1):
        if (isPtOnSeg(pt, [poly[i], poly[i + 1]])):
            return True
    return False

def isPtOnPoly(
    pt:         "2-tuple of coordinate (x, y)",
    poly:       "List of pts, form a close area"
    ) -> "Return true if the pt is on the polygon, including the edge":

    x = pt[1]
    y = pt[0]
    onPoly = False

    # [Need Rewriting]
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
            onPoly = not onPoly
        j = i
    return onPoly

def isPtInsidePoly(
    pt:         "2-tuple of coordinate (x, y)",
    poly:       "List of pts, form a close area"
    ) -> "Return true if the pt is in the interior of polygon":
    # Check if the pt is on the edge of the polygon ===========================
    if (isPtOnPoly(pt, poly) and not isPtOnPolyEdge(pt, poly)):
        return True
    return False

def is3PtsClockWise(
    pt1:        "2-tuple of coordinate (x, y)",
    pt2:        "2-tuple of coordinate (x, y)",
    pt3:        "2-tuple of coordinate (x, y)"
    ) -> "True if three given points are clockWise, False otherwise (could be collinear)":
    # Use Determinant to determine ============================================
    [x1, y1] = [pt1[0], pt1[1]]
    [x2, y2] = [pt2[0], pt2[1]]
    [x3, y3] = [pt3[0], pt3[1]]
    ori = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)

    if (abs(ori) <= CONST_EPSILON):
        return None
    elif (ori < 0):
        return True
    else:
        return False

def isSegCrossSeg(
    seg1:       "List of 2 pts as two bounds ",
    seg2:       "List of 2 pts as two bounds "
    ) -> "Return true if segs are crossing in the interior of both segs, false otherwise":    
    # Check isSegIntSeg =======================================================
    segIntSeg = isSegIntSeg(seg1, seg2)
    if (segIntSeg != True):
        return segIntSeg

    # Check for collinear =====================================================
    if (isPtOnSeg(seg1[0], seg2) or isPtOnSeg(seg1[1], seg2) or isPtOnSeg(seg2[0], seg1) or isPtOnSeg(seg2[1], seg2)):
        return False

    return True

def isSegIntSeg(
    seg1:       "List of 2 pts as two bounds ",
    seg2:       "List of 2 pts as two bounds "
    ) -> "Return true if segs are intersecting, including the case where the intersection is at one end of line segment, false otherwise":
    # Validation ==============================================================
    if (is2PtsSame(seg1[0], seg1[1])):
        print(ERROR_ZERO_VECTOR)
        return None
    if (is2PtsSame(seg2[0], seg2[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Initialize ==============================================================
    [p, q] = seg1
    [u, w] = seg2

    # Clockwise checking ======================================================
    loopPQU = is3PtsClockWise(p, q, u)
    loopPQW = is3PtsClockWise(p, q, w)
    loopUWP = is3PtsClockWise(u, w, p)
    loopUWQ = is3PtsClockWise(u, w, q)

    # Check for collinear =====================================================
    if (loopPQU == None or loopPQW == None or loopUWP == None or loopUWQ == None):
        return True

    # Else check intersection =================================================
    if ((loopPQU != loopPQW and loopUWP != loopUWQ)
        or (loopPQU != loopPQW and loopUWP == loopUWQ)
        or (loopPQU == loopPQW and loopUWP != loopUWQ)
        or (loopPQU == loopPQW and loopUWP == loopUWQ)):
        return True
    return False

def isSegCrossRay(
    seg:        "List of 2 pts as two bounds ",
    ray:        "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return true if ray is crossing interior of seg, false otherwise":
    crossFlag = (intSeg2Ray(seg, ray) != None)
    return crossFlag

def isRayCrossRay(
    ray1:       "List of 2 pts, the first pt defines the initial point, the second pt is on the ray",
    ray2:       "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return true if ray is crossing interior of ray, false otherwise":
    crossFlag = (intRay2Ray(ray1, ray2) != None)
    return crossFlag

def isSegIntRay(
    seg:        "List of 2 pts as two bounds ",
    ray:        "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
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
            and (isPtInsidePoly(w, [p, q, u])
                or isPtOnPolyEdge(w, [p, q, u]))):
            intersectFlag = True

    return intersectFlag

def isSegCrossPoly(
    seg:        "List of 2 pts as two bounds ",
    poly:       "List of pts, form a close area"
    ) -> "Return true if the line segment is crossing the interior of the poly, false otherwise":

    if (isPtInsidePoly(seg[0], poly) or isPtInsidePoly(seg[1], poly)):
        return True
    for i in range(-1, len(poly) - 1):
        edge = [poly[i], poly[i + 1]]
        if (isSegCrossSeg(edge, seg)):
            return True
    return False

def isRayCrossPoly(
    ray:        "List of 2 pts, the first pt defines the initial point, the second pt is on the ray",
    poly:       "List of pts, form a close area"
    ) -> "Return true if the ray is crossing the interior of the poly, false otherwise":
    if (isPtInsidePoly(ray[0], poly)):
        return True
    for i in range(-1, len(poly) - 1):
        edge = [poly[i], poly[i + 1]]
        if (isSegCrossRay(edge, ray)):
            return True
    return False

def intLine2Line(
    line1:      "List of 2 pts that defines a line",
    line2:      "List of 2 pts that defines a line"
    ) -> "Return 1) None if parallel or overlapped, or \
                 2) The intersect point if two lines intersect":
    # Validation ==============================================================
    if (is2PtsSame(line1[0], line1[1])):
        print(ERROR_ZERO_VECTOR)
        return None
    if (is2PtsSame(line2[0], line2[1])):
        print(ERROR_ZERO_VECTOR)
        return None

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
    seg1:       "List of 2 pts as two bounds ",
    seg2:       "List of 2 pts as two bounds "
    ) -> "Return 1) None if no intersection, or \
                 2) The intersect point if two segments intersect":
    # Calculate intersection for two lines ====================================
    ptInt = intLine2Line(seg1, seg2)

    # Check if it is on segs ==================================================
    if (ptInt != None and isPtOnSeg(ptInt, seg1) and isPtOnSeg(ptInt, seg2)):
        return ptInt
    else:
        return None

def intSeg2Ray(
    seg:        "List of 2 pts as two bounds ",
    ray:        "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return 1) None if no intersection, or \
                 2) The intersect point if two segments intersect":
    # Calculate intersection for two lines ====================================
    ptInt = intLine2Line(seg, ray)
    # print("Intersect: ", ptInt, isPtOnSeg(ptInt, seg), isPtOnRay(ptInt, ray))

    # Check if it is on segs
    if (ptInt != None and isPtOnSeg(ptInt, seg) and isPtOnRay(ptInt, ray)):
        return ptInt
    else:
        return None

def intRay2Line(
    ray:        "List of 2 pts, the first pt defines the initial point, the second pt is on the ray",
    line:       "List of 2 pts that defines a line"
    ) -> "Return 1) None if no intersection, or \
                 2) The intersect point if ray and line intersect":
    # Validation ==============================================================
    if (is2PtsSame(ray[0], ray[1])):
        print(ERROR_ZERO_VECTOR)
        return None
    if (is2PtsSame(line[0], line[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Calculate intersection ==================================================
    ptInt = intLine2Line(ray, line)

    # Check if it is on segs ==================================================
    if (ptInt != None and isPtOnRay(ptInt, ray)):
        return ptInt
    else:
        return None

def intRay2Ray(
    ray1:       "List of 2 pts, the first pt defines the initial point, the second pt is on the ray",
    ray2:       "List of 2 pts, the first pt defines the initial point, the second pt is on the ray"
    ) -> "Return 1) None if no intersection, or \
                 2) The intersect point if two rays intersect":
    # Validation ==============================================================
    if (is2PtsSame(ray1[0], ray1[1])):
        print(ERROR_ZERO_VECTOR)
        return None
    if (is2PtsSame(ray2[0], ray2[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Calculate intersection ==================================================
    ptInt = intLine2Line(ray1, ray2)

    # Check if it is on segs ==================================================
    if (ptInt != None and isPtOnRay(ptInt, ray1) and isPtOnRay(ptInt, ray2)):
        return ptInt
    else:
        return None

def cosRay2Ray(
    ray1:       "Ending pt, the origin pt is (0, 0)",
    ray2:       "Ending pt, the origin pt is (0, 0)"
    ) -> "Calculate cosine between two vectors":
    # Validation ==============================================================
    if (is2PtsSame(ray1[0], ray1[1])):
        print(ERROR_ZERO_VECTOR)
        return None
    if (is2PtsSame(ray2[0], ray2[1])):
        print(ERROR_ZERO_VECTOR)
        return None

    # Calculate cosine ========================================================
    cosAngle = (ray1[0] * ray2[0] + ray1[1] * ray2[1]) / (math.sqrt(ray1[0] * ray1[0] + ray1[1] * ray1[1]) * math.sqrt(ray2[0] * ray2[0] + ray2[1] * ray2[1]))
    return cosAngle

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
