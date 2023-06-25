import geopy.distance
import heapq
import math
import shapely

from .error import *
from .common import *
from .const import *
from .graph import *
from .msg import *
# from .relation import *

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

def distEuclidean2D(pt1: pt, pt2: pt) -> float:
    """Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number"""

    if (pt1 != None and pt2 != None):
        return math.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2)
    else:
        return 0

def distManhattenXY(pt1: pt, pt2: pt) -> float:
    """Gives a Euclidean distance based on two coords, if two coordinates are the same, return a small number"""

    if (pt1 != None and pt2 != None):
        return abs(pt1[0] - pt2[0]) + abs(pt1[1] - pt2[1])
    else:
        return 0

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

def spdInTimedSeq(seq: list[pt], timeStamp: list[float], t: float) -> float:
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
            dist = distEuclidean2D(seq[i], seq[i + 1])
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

def locInPoly(poly: poly, startPt: pt, dist: int|float, reverseFlag: bool=False, dimension: str = 'XY') -> pt:

    perimeter = calPolygonPerimeter(poly)
    while(dist > perimeter):
        dist -= perimeter

    seq = [startPt]
    startID = len(poly) + 1

    if (dimension == 'XY'):
        for i in range(len(poly) - 1):
            if (abs(distEuclidean2D(poly[i], poly[i + 1]) - distEuclidean2D(poly[i], startPt) - distEuclidean2D(startPt, poly[i + 1])) <= CONST_EPSILON):
                startID = i + 1
                break
        if (startID == len(poly) + 1):
            if (abs(distEuclidean2D(poly[0], poly[-1]) - distEuclidean2D(poly[0], startPt) - distEuclidean2D(startPt, poly[-1])) <= CONST_EPSILON):
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

def traceInPoly(poly: poly, startPt: pt, endPt:pt, dist: int|float, reverseFlag: bool=False) -> list[pt]:

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
            accDist += distEuclidean2D(seq[i], seq[i + 1])
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
        segDist = distEuclidean2D(preLoc, nextLoc)
    if (segDist <= CONST_EPSILON):
        raise ZeroDivisionError
    lat = nextLoc[0] + (remainDist / segDist) * (preLoc[0] - nextLoc[0])
    lon = nextLoc[1] + (remainDist / segDist) * (preLoc[1] - nextLoc[1])
    return (lat, lon)

def calPolygonPerimeter(poly: poly) -> float:
    p = 0
    for i in range(-1, len(poly) - 1):
        p += distEuclidean2D(poly[i], poly[i + 1])
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
    nodeIDs: list[int|str]|str = 'All',
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

def _getTauEuclidean(nodes: dict, nodeIDs: list[int|str]|str, speed = 1):
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

def _getTauManhatten(nodes: dict, nodeIDs: list[int|str]|str, speed = 1):
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

def _getTauLatLon(nodes: dict, nodeIDs: list[int|str]|str, distUnit = 'meter', speed = 1):
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

def _getTauGrid(nodes: dict, nodeIDs: list[int|str]|str, grid: dict):
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

def getSortedNodesByDist(nodes: dict, edges: dict, refNodeID: int|str, nodeIDs: list[int|str]|str = 'All'):

    """Given a set of locations, and a referencing node, sort the nodes by distance to this referencing node"""

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define tau ==============================================================
    tau = getTau(nodes, edges, refNodeID, nodeIDs, 0)

    # Sort distance ===========================================================
    sortedSeq = []
    sortedSeqHeap = []
    for n in nodeIDs:
        dist = tau[refNodeID, n]
        heapq.heappush(sortedSeqHeap, (dist, n))
    while (len(sortedSeqHeap) > 0):
        sortedSeq.append(heapq.heappop(sortedSeqHeap)[1])  

    return sortedSeq

def getSweepSeq( nodes: dict, nodeIDs: list[int|str]|str = 'All', centerLoc: None|pt = None, isClockwise: bool = True, initDeg: float = 0) -> list:
    
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

def gridPathFinding(grid: dict, startCoord: pt, endCoord: pt, algo: dict = {'method': 'A*', 'measure': 'Manhatten'}) -> dict:

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
    startCoord: 2-tuple|2-list, required
        Starting location on the grid
    endCoord: 2-tuple|2-list, required
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
    res = {}

    # Call path finding =======================================================
    if (algo['method'] == 'A*'):
        if ('measure' not in algo or algo['measure'] not in ['Manhatten', 'Euclidean']):
            warnings.warn("WARNING: Set distance measurement to be default as 'Manhatten")
        res = _gridPathFindingAStar(column, row, barriers, startCoord, endCoord, algo['measure'])
    else:
        print("Error: Incorrect or not available grid path finding option!")
    return res

def _gridPathFindingAStar(column, row, barriers, startCoord, endCoord, distMeasure):
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
        gridStatus[startCoord] = (0, _calManhattenDist(startCoord, endCoord), None)
    elif (distMeasure == 'Euclidean'):
        gridStatus[startCoord] = (0, _calEuclideanDist(startCoord, endCoord), None)
    gridStatus[endCoord] = (None, 0, None)

    # Open/close set ======================================================
    openList = [startCoord]
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
                    gridStatus[upCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(upCoord, endCoord), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[upCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(upCoord, endCoord), coord)
                if (upCoord == endCoord):
                    break
                else:
                    tmpOpenList.append(upCoord)
        # Down
        downCoord = (coord[0], coord[1] - 1)
        if (coord[1] - 1 >= 0 and gridStatus[downCoord] != None and gridStatus[downCoord] != 'block' and downCoord not in closeList):
            if (gridStatus[downCoord][0] == None or gridStatus[downCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[downCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(downCoord, endCoord), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[downCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(downCoord, endCoord), coord)
                if (downCoord == endCoord):
                    break
                else:
                    tmpOpenList.append(downCoord)
        # Left
        leftCoord = (coord[0] - 1, coord[1])
        if (coord[0] - 1 >= 0 and gridStatus[leftCoord] != None and gridStatus[leftCoord] != 'block' and leftCoord not in closeList):
            if (gridStatus[leftCoord][0] == None or gridStatus[leftCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[leftCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(leftCoord, endCoord), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[leftCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(leftCoord, endCoord), coord)
                if (leftCoord == endCoord):
                    break
                else:
                    tmpOpenList.append(leftCoord)
        # Right
        rightCoord = (coord[0] + 1, coord[1])
        if (coord[0] + 1 < column and gridStatus[rightCoord] != None and gridStatus[rightCoord] != 'block' and rightCoord not in closeList):
            if (gridStatus[rightCoord][0] == None or gridStatus[rightCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[rightCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(rightCoord, endCoord), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[rightCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(rightCoord, endCoord), coord)
                if (rightCoord == endCoord):
                    break
                else:
                    tmpOpenList.append(rightCoord)
        openList.remove(coord)
        openList.extend(tmpOpenList)
        closeList.append(coord)

    # Recover path ========================================================
    path = []
    curCoord = endCoord
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
