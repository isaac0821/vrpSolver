import math
import random
import tripy
import warnings

from .common import *
from .const import *
from .geometry import *
from .msg import *
from .error import *

# History =====================================================================
# 20230510 - Cleaned up for v0.0.55, including the following available nodes 
#            distributions: 'UniformSquareXY', 'UniformPolyXY', 'UniformCircleXY', 
#            'UniformCircleLatLon', 'RoadNetworkPolyLatLon', and 'RoadNetworkCircleLatLon'
# 20230515 - Revise the parameter annotation according to PEP 3107
# 20230624 - Add `rndPlainArcs()`
# =============================================================================

def rndPlainNodes(
    N: int|None = None, 
    nodeIDs: list[int|str] = [], 
    distr: dict = {
            'method': 'UniformSquareXY', 
            'xRange': (0, 100), 
            'yRange': (0, 100)
        }
    ) -> dict:

    """Randomly create a set of locations

    Parameters
    ----------

    N: integer, optional, default None
        Number of locations/vertices/customers to be randomly created
    nodeIDs: list, optional, default None
        Alternative input parameter of `N`. A list of node IDs, `N` will be overwritten if `nodeIDs` is given
    distr: dictionary, optional, default {'method': 'UniformSquareXY', 'xRange': (0, 100), 'yRange': (0, 100)}
        Spatial distribution of nodes, options are as following:
            1) (default) Uniformly sample from a square on the Euclidean space
            >>> distr = {
            ...     'method': 'UniformSquareXY', 
            ...     'xRange': (0, 100), # A 2-tuple with minimum/maximum range of x, default as (0, 100), 
            ...     'yRange': (0, 100), # A 2-tuple with minimum/maximum range of y, default as (0, 100), 
            ... }
            2) Uniformly sample from a given polygon on the Euclidean space
            >>> distr = {
            ...     'method': 'UniformPolyXY', 
            ...     'polyXY': poly, # polygon of the area, (no holes)
            ...     'polyXYs': polys, # alternative option for 'polyXY', as a list of polygons 
            ... }
            3) Uniformly sample from a circle on the Euclidean space
            >>> distr = {
            ...     'method': 'UniformCircleXY',
            ...     'centerXY': (0, 0), # centering location, default as (0, 0), 
            ...     'radius': 100, # radius of the circle , default as 100
            ... }
            4) Uniformly sample from a given polygon by lat/lon
            >>> distr = {
            ...     'method': 'UniformPolyLatLon', 
            ...     'polyLatLon': polygon of the area, (no holes)
            ...     'polyLatLons': alternative option for 'polyLatLon', as a list of polygons 
            ... }
            5) Uniformly sample from a given circle by lat/lon,
            >>> distr = {
            ...     'method': 'UniformCircleLatLon', 
            ...     'centerLatLon': required, centering location in lat/lon, 
            ...     'radiusInMeters': radius of the circle in meters 
            ... }
            6) Uniformly generate from a given polygon on a road network
            >>> distr = {
            ...     'method': 'RoadNetworkPolyLatLon'
            ...     'RoadNetwork': list of arcs that can be sampled 
            ...     'polyLatLon': nodes should generated within the polygon, if not provided, will consider the entire network, 
            ...     'roadClass': list of classes that can be sampled 
            ... }
            7) Uniformly generate from a given circle on a road network
            >>> distr = {
            ...     'method': 'RoadNetworkCircleLatLon', 
            ...     'RoadNetwork': list of arcs that can be sampled 
            ...     'centerLatLon': [lat, lon], 
            ...     'radiusInMeters': radius in [m] 
            ...     'roadClass': list of classes that can be sampled
            ... }

    Returns
    -------
    dictionray
        A set of randomly created locations, as in the following format::
        >>> nodes[nodeID] = {
        ...     'loc': (lat, lon)
        ... }

    Raises
    ------
    MissingParameterError
        Missing `distr` field or values in `distr`.
    UnsupportedInputError
        Option is not supported for `distr['method']`
    NotAvailableError
        Functions/options that are not ready yet.
    EmptyError
        The sample area is empty.
    """

    # Sanity checks ===========================================================
    if (distr == None or 'method' not in distr):
        raise MissingParameterError(ERROR_MISSING_NODES_DISTR)

    nodes = {}
    if (nodeIDs == [] and N == None):
        raise MissingParameterError(ERROR_MISSING_N)
    elif (nodeIDs == [] and N != None):
        nodeIDs = [i for i in range(N)]

    # Generate instance =======================================================
    # Uniformly sample from a square on the Euclidean space
    if (distr['method'] == 'UniformSquareXY'):
        xRange = None
        yRange = None
        if ('xRange' not in distr or 'yRange' not in distr):
            xRange = [0, 100]
            yRange = [0, 100]
            warnings.warn("WARNING: Set sampled area to be default as a (0, 100) x (0, 100) square")
        else:
            xRange = [float(distr['xRange'][0]), float(distr['xRange'][1])]
            yRange = [float(distr['yRange'][0]), float(distr['yRange'][1])]
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformSquareXY(xRange, yRange)
            }

    # Uniformly sample from a polygon/a list of polygons on the Euclidean space
    elif (distr['method'] == 'UniformPolyXY'):
        if ('polyXY' not in distr and 'polyXYs' not in distr):
            raise MissingParameterError("ERROR: Missing required key 'polyXY' or 'polyXYs' in field `distr`, which indicates a polygon / a list of polygons in the Euclidean space")
        if ('polyXY' in distr):
            for n in nodeIDs:
                nodes[n] = {
                    'loc': _rndPtUniformPolyXY(distr['polyXY'])
                }
        elif ('polyXYs' in distr):
            for n in nodeIDs:
                nodes[n] = {
                    'loc': _rndPtUniformPolyXYs(distr['polyXYs'])
                }

    # Uniformly sample from a circle on the Euclidean space
    elif (distr['method'] == 'UniformCircleXY'):
        centerXY = None
        radius = None
        if ('centerXY' not in distr or 'radius' not in distr):
            centerXY = (0, 0)
            radius = 100
            warnings.warn("WARNING: Set sample area to be default as a circle with radius of 100 centering at (0, 0)")
        else:
            centerXY = distr['centerXY']
            radius = distr['radius']
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformCircleXY(radius, centerXY)
            }

    # Uniformly sample from a polygon by lat/lon
    elif (distr['method'] == 'UniformPolyLatLon'):
        if ('polyLatLon' not in distr and 'polyLatLons' not in distr):
            raise MissingParameterError("ERROR: Missing required key 'polyXY' or 'polyXYs' in field `distr`, which indicates a polygon / a list of polygons in the Euclidean space")
        # TODO: Mercator projection
        raise VrpSolverNotAvailableError("ERROR: 'UniformPolyLatLon' is not available yet, please stay tune.")

    # Uniformly sample from a circle by lat/lon
    elif (distr['method'] == 'UniformCircleLatLon'):
        if ('centerLatLon' not in distr or 'radiusInMeters' not in distr):
            raise MissingParameterError("ERROR: Missing required key 'centerLatLon' or 'radiusInMeters' in field `distr`.")
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformCircleLatLon(distr['radiusInMeters'], distr['centerLatLon'])
            }

    # Uniformly sample from the roads/streets within a polygon/a list of polygons from given road networks
    elif (distr['method'] == 'RoadNetworkPolyLatLon'):
        if ('polyLatLon' not in distr):
            raise MissingParameterError("ERROR: Missing required key 'polyXY' or 'polyXYs' in field `distr`, which indicates a polygon / a list of polygons in the Euclidean space")
        elif ('roadNetwork' not in distr):
            raise MissingParameterError("ERROR: Missing required key 'RoadNetwork' in field `distr`. Need to provide the road network where the nodes are generated.")
        elif ('roadClass' not in distr):
            warnings.warn("WARNING: Set 'roadClass' to be default as ['residential']")
        nodeLocs = _rndPtRoadNetworkPolyLatLon(
            N if N != None else len(nodeIDs),
            distr['roadNetwork'], 
            distr['polyLatLon'],
            distr['roadClass'] if 'roadClass' in distr else ['residential'])
        for n in range(len(nodeIDs)):
            nodes[nodeIDs[n]] = {
                'loc': nodeLocs[n]
            }

    # Uniformly sample from the roads/streets within a circle from given road network
    elif (distr['method'] == 'RoadNetworkCircleLatLon'):
        if ('centerLatLon' not in distr or 'radiusInMeters' not in distr):
            raise MissingParameterError("ERROR: Missing required key 'centerLatLon' or 'radiusInMeters' in field `distr`.")
        elif ('roadNetwork' not in distr):
            raise MissingParameterError("ERROR: Missing required key 'RoadNetwork' in field `distr`. Need to provide the road network where the nodes are generated.")
        elif ('roadClass' not in distr):
            warnings.warn("WARNING: Set 'roadClass' to be default as ['residential']")
        nodeLocs = _rndPtRoadNetworkCircleLatLon(
            N if N != None else len(nodeIDs),
            distr['roadNetwork'], 
            distr['radiusInMeters'],
            distr['centerLatLon'],
            distr['roadClass'] if 'roadClass' in distr else ['residential'])
        for n in range(len(nodeIDs)):
            nodes[nodeIDs[n]] = {
                'loc': nodeLocs[n]
            }
    else:
        raise UnsupportedInputError(ERROR_MISSING_NODES_DISTR)

    return nodes

def _rndPtUniformSquareXY(xRange: list[int]|list[float], yRange: list[int]|list[float]) -> pt:
    x = random.uniform(xRange[0], xRange[1])
    y = random.uniform(yRange[0], yRange[1])
    return (x, y)

def _rndPtUniformTriangleXY(triangle: poly) -> pt:
    
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

def _rndPtUniformPolyXY(poly: poly) -> pt:
    # Get list of triangles ===================================================
    # TODO: tripy.earclip() to be replaced
    lstTriangle = tripy.earclip(poly)

    # Weight them and make draws ==============================================
    lstWeight = []
    for i in range(len(lstTriangle)):
        lstWeight.append(calTriangleAreaXY(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2]))

    # Select a triangle and randomize a point in the triangle =================
    idx = rndPick(lstWeight)
    (x, y) = _rndPtUniformTriangleXY(lstTriangle[idx])

    return (x, y)

def _rndPtUniformPolyXYs(polys: polys) -> pt:
    # Get all triangulated triangles ==========================================
    # TODO: tripy.earclip() to be replaced
    lstTriangle = []
    for p in polys:
        lstTriangle.extend(tripy.earclip(p))

    # Weight them and make draws ==============================================
    lstWeight = []
    for i in range(len(lstTriangle)):
        lstWeight.append(calTriangleAreaXY(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2]))

    # Select a triangle and randomize a point in the triangle =================
    idx = rndPick(lstWeight)
    (x, y) = _rndPtUniformTriangleXY(lstTriangle[idx])

    return (x, y)

def _rndPtUniformCircleXY(radius: float, center: pt) -> pt:
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, radius ** 2))
    x = center[0] + r * math.cos(theta)
    y = center[1] + r * math.sin(theta)

    return (x, y)

def _rndPtUniformCircleLatLon(radius: float, center: pt) -> pt:
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, radius ** 2))
    (lat, lon) = ptInDistLatLon(center, theta, r)

    return (lat, lon)

def _rndPtRoadNetworkPolyLatLon(N: int, road: dict, poly: poly, roadClass: str | list[str]) -> list[pt]:
    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for rID in road:
        roadLength = 0
        includedFlag = False
        if ('class' in road[rID] and road[rID]['class'] in roadClass):
            if (poly == None):
                includedFlag = True
            else:
                for i in range(len(road[rID]['shape'])):
                    if (isPtInPoly(road[rID]['shape'][i], poly)):
                        includedFlag = True
                        break
        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(road[rID]['shape']) - 1):
                roadLength += distLatLon(road[rID]['shape'][i], road[rID]['shape'][i + 1])['dist']
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(rID)

    # Check if there are roads included =======================================
    if (sum(lengths) == 0):
        raise EmptyError("ERROR: No road is found.")

    # Use accept-denial to test if the node is within poly ====================
    # FIXME: Inefficient approach, will need to be rewritten
    # TODO: Truncate the roads that partially inside polygon
    nodeLocs = []
    for i in range(N):
        lat = None
        lon = None
        if (poly == None):
            idx = rndPick(lengths)
            edgeLength = lengths[idx]
            edgeDist = random.uniform(0, 1) * edgeLength
            (lat, lon) = locInSeq(road[roadIDs[idx]]['shape'], edgeDist, 'LatLon')
        else:
            insideFlag = False
            while (not insideFlag):
                idx = rndPick(lengths)
                edgeLength = lengths[idx]
                edgeDist = random.uniform(0, 1) * edgeLength
                (lat, lon) = locInSeq(road[roadIDs[idx]]['shape'], edgeDist, 'LatLon')
                if (isPtInPoly([lat, lon], poly)):
                    insideFlag = True
        nodeLocs.append((lat, lon))

    return nodeLocs

def _rndPtRoadNetworkCircleLatLon(N: int, road: dict, radius: float, center: pt, roadClass: str | list[str]) -> list[pt]:
    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for rID in road:
        roadLength = 0
        includedFlag = False
        for i in range(len(road[rID]['shape'])):
            if (road[rID]['class'] in roadClass and distLatLon(road[rID]['shape'][i], center) <= radius)['dist']:
                includedFlag = True
                break

        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(road[rID]['shape']) - 1):
                roadLength += distLatLon(road[rID]['shape'][i], road[rID]['shape'][i + 1])['dist']
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(rID)


    # Check if there are roads included =======================================
    if (sum(lengths) == 0):
        return []

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
            (lat, lon) = locInSeq(road[roadIDs[idx]]['shape'], edgeDist, 'LatLon')
            if (distLatLon([lat, lon], center)['dist'] <= radius):
                insideFlag = True
        nodeLocs.append((lat, lon))

    return nodeLocs

def rndPlainArcs(
    A: int|None = None,
    arcIDs: list[int|str] = [],
    distr: dict = {
            'method': 'UniformLengthInSquareXY',
            'xRange': (0, 100),
            'yRange': (0, 100),
            'minLen': 0,
            'maxLen': 10
        }
    ) -> dict:

    """Randomly create a set of arcs (to be visited) 

    Parameters
    ----------

    A: integer, optional, default None
        Number of arcs to be visited
    arcIDs: list, optional, default None
        Alternative input parameter of `A`. A list of arc IDs, `A` will be overwritten if `arcIDs` is given
    distr: dictionary, optional, default {'method': 'UniformLengthInSquareXY', 'xRange': (0, 100), 'yRange': (0, 100), 'minLen': 0, 'maxLen': 10}
        Spatial distribution of arcs, optional are as following:
        1) (default) Uniformly sample from a square on the Euclidean space, with uniformly selected length
        >>> distr = {
        ...     'method': 'UniformLengthInSquareXY',
        ...     'xRange': (0, 100),
        ...     'yRange': (0, 100),
        ...     'minLen': 0,
        ...     'maxLen': 10
        ... }

    """

    # Sanity check ============================================================
    if (distr == None or 'method' not in distr):
        raise MissingParameterError(ERROR_MISSING_ARCS_DISTR)

    arcs = {}
    if (arcIDs == [] and A == None):
        raise MissingParameterError(ERROR_MISSING_N)
    elif (arcIDs == [] and A != None):
        arcIDs = [i for i in range(A)]

    # Generate instance =======================================================
    if (distr['method'] == 'UniformLengthInSquareXY'):
        if ('minLen' not in distr or 'maxLen' not in distr):
            raise MissingParameterError("ERROR: Missing required field 'minLen' and/or 'maxLen' in `distr`")
        xRange = None
        yRange = None
        if ('xRange' not in distr or 'yRange' not in distr):
            xRange = [0, 100]
            yRange = [0, 100]
            warnings.warn("WARNING: Set sample area to be default as a (0, 100) x (0, 100) square")
        else:
            xRange = [float(distr['xRange'][0]), float(distr['xRange'][1])]
            yRange = [float(distr['yRange'][0]), float(distr['yRange'][1])]
        for n in arcIDs:
            arcs[n] = {
                'arc': _rndArcUniformSquareXY(xRange, yRange, distr['minLen'], distr['maxLen'])
            }
    else:
        raise UnsupportedInputError(ERROR_MISSING_ARCS_DISTR)

    return arcs

def _rndArcUniformSquareXY(xRange: list[int]|list[float], yRange: list[int]|list[float], minLen: int|float, maxLen: int|float) -> tuple[pt, pt]:
    length = random.uniform(minLen, maxLen)
    direction = random.uniform(0, 360)
    xStart = random.uniform(xRange[0], xRange[1])
    yStart = random.uniform(yRange[0], yRange[1])
    (xEnd, yEnd) = ptInDistXY((xStart, yStart), direction, length)
    return ((xStart, yStart), (xEnd, yEnd))

