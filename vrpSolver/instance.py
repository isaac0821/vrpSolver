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
    method: dict = {
            'distr': 'UniformSquareXY', 
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
    method: dictionary, optional, default {'distr': 'UniformSquareXY', 'xRange': (0, 100), 'yRange': (0, 100)}
        Spatial distribution of nodes, options are as following:
            1) (default) Uniformly sample from a square on the Euclidean space
            >>> method = {
            ...     'distr': 'UniformSquareXY', 
            ...     'xRange': (0, 100), # A 2-tuple with minimum/maximum range of x, default as (0, 100), 
            ...     'yRange': (0, 100), # A 2-tuple with minimum/maximum range of y, default as (0, 100), 
            ... }
            2) Uniformly sample from a given polygon on the Euclidean space
            >>> method = {
            ...     'distr': 'UniformPolyXY', 
            ...     'polyXY': poly, # polygon of the area, (no holes)
            ...     'polyXYs': polys, # alternative option for 'polyXY', as a list of polygons 
            ... }
            3) Uniformly sample from a circle on the Euclidean space
            >>> method = {
            ...     'distr': 'UniformCircleXY',
            ...     'centerXY': (0, 0), # centering location, default as (0, 0), 
            ...     'radius': 100, # radius of the circle , default as 100
            ... }
            4) Uniformly sample from a given polygon by lat/lon
            >>> method = {
            ...     'distr': 'UniformPolyLatLon', 
            ...     'polyLatLon': polygon of the area, (no holes)
            ...     'polyLatLons': alternative option for 'polyLatLon', as a list of polygons 
            ... }
            5) Uniformly sample from a given circle by lat/lon,
            >>> method = {
            ...     'distr': 'UniformCircleLatLon', 
            ...     'centerLatLon': required, centering location in lat/lon, 
            ...     'radiusInMeters': radius of the circle in meters 
            ... }
            6) Uniformly generate from a given polygon on a road network
            >>> method = {
            ...     'distr': 'RoadNetworkPolyLatLon'
            ...     'RoadNetwork': list of arcs that can be sampled 
            ...     'polyLatLon': nodes should generated within the polygon, if not provided, will consider the entire network, 
            ...     'roadClass': list of classes that can be sampled 
            ... }
            7) Uniformly generate from a given circle on a road network
            >>> method = {
            ...     'distr': 'RoadNetworkCircleLatLon', 
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
        Missing `method` field or values in `method`.
    UnsupportedInputError
        Option is not supported for `method['distr']`
    NotAvailableError
        Functions/options that are not ready yet.
    EmptyError
        The sample area is empty.
    """

    # Sanity checks ===========================================================
    if (method == None or 'distr' not in method):
        raise MissingParameterError(ERROR_MISSING_NODES_DISTR)

    nodes = {}
    if (nodeIDs == [] and N == None):
        raise MissingParameterError(ERROR_MISSING_N)
    elif (nodeIDs == [] and N != None):
        nodeIDs = [i for i in range(N)]

    # Generate instance =======================================================
    # Uniformly sample from a square on the Euclidean space
    if (method['distr'] == 'UniformSquareXY'):
        xRange = None
        yRange = None
        if ('xRange' not in method or 'yRange' not in method):
            xRange = [0, 100]
            yRange = [0, 100]
            warnings.warn("WARNING: Set sampled area to be default as a (0, 100) x (0, 100) square")
        else:
            xRange = [float(method['xRange'][0]), float(method['xRange'][1])]
            yRange = [float(method['yRange'][0]), float(method['yRange'][1])]
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformSquareXY(xRange, yRange)
            }

    # Uniformly sample from a polygon/a list of polygons on the Euclidean space
    elif (method['distr'] == 'UniformPolyXY'):
        if ('polyXY' not in method and 'polyXYs' not in method):
            raise MissingParameterError("ERROR: Missing required key 'polyXY' or 'polyXYs' in field `method`, which indicates a polygon / a list of polygons in the Euclidean space")
        if ('polyXY' in method):
            for n in nodeIDs:
                nodes[n] = {
                    'loc': _rndPtUniformPolyXY(method['polyXY'])
                }
        elif ('polyXYs' in method):
            for n in nodeIDs:
                nodes[n] = {
                    'loc': _rndPtUniformPolyXYs(method['polyXYs'])
                }

    # Uniformly sample from a circle on the Euclidean space
    elif (method['distr'] == 'UniformCircleXY'):
        centerXY = None
        radius = None
        if ('centerXY' not in method or 'radius' not in method):
            centerXY = (0, 0)
            radius = 100
            warnings.warn("WARNING: Set sample area to be default as a circle with radius of 100 centering at (0, 0)")
        else:
            centerXY = method['centerXY']
            radius = method['radius']
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformCircleXY(radius, centerXY)
            }

    # Uniformly sample from a polygon by lat/lon
    elif (method['distr'] == 'UniformPolyLatLon'):
        if ('polyLatLon' not in method and 'polyLatLons' not in method):
            raise MissingParameterError("ERROR: Missing required key 'polyXY' or 'polyXYs' in field `method`, which indicates a polygon / a list of polygons in the Euclidean space")
        # TODO: Mercator projection
        raise VrpSolverNotAvailableError("ERROR: 'UniformPolyLatLon' is not available yet, please stay tune.")

    # Uniformly sample from a circle by lat/lon
    elif (method['distr'] == 'UniformCircleLatLon'):
        if ('centerLatLon' not in method or 'radiusInMeters' not in method):
            raise MissingParameterError("ERROR: Missing required key 'centerLatLon' or 'radiusInMeters' in field `method`.")
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformCircleLatLon(method['radiusInMeters'], method['centerLatLon'])
            }

    # Uniformly sample from the roads/streets within a polygon/a list of polygons from given road networks
    elif (method['distr'] == 'RoadNetworkPolyLatLon'):
        if ('polyLatLon' not in method):
            raise MissingParameterError("ERROR: Missing required key 'polyXY' or 'polyXYs' in field `method`, which indicates a polygon / a list of polygons in the Euclidean space")
        elif ('roadNetwork' not in method):
            raise MissingParameterError("ERROR: Missing required key 'RoadNetwork' in field `method`. Need to provide the road network where the nodes are generated.")
        elif ('roadClass' not in method):
            warnings.warn("WARNING: Set 'roadClass' to be default as ['residential']")
        nodeLocs = _rndPtRoadNetworkPolyLatLon(
            N if N != None else len(nodeIDs),
            method['roadNetwork'], 
            method['polyLatLon'],
            method['roadClass'] if 'roadClass' in method else ['residential'])
        for n in range(len(nodeIDs)):
            nodes[nodeIDs[n]] = {
                'loc': nodeLocs[n]
            }

    # Uniformly sample from the roads/streets within a circle from given road network
    elif (method['distr'] == 'RoadNetworkCircleLatLon'):
        if ('centerLatLon' not in method or 'radiusInMeters' not in method):
            raise MissingParameterError("ERROR: Missing required key 'centerLatLon' or 'radiusInMeters' in field `method`.")
        elif ('roadNetwork' not in method):
            raise MissingParameterError("ERROR: Missing required key 'RoadNetwork' in field `method`. Need to provide the road network where the nodes are generated.")
        elif ('roadClass' not in method):
            warnings.warn("WARNING: Set 'roadClass' to be default as ['residential']")
        nodeLocs = _rndPtRoadNetworkCircleLatLon(
            N if N != None else len(nodeIDs),
            method['roadNetwork'], 
            method['radiusInMeters'],
            method['centerLatLon'],
            method['roadClass'] if 'roadClass' in method else ['residential'])
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

def rndTimeWindowsNodes(
    nodes: dict,
    nodeIDs: list[int|str]|str = 'All',
    method: dict = {
            'mode': 'RandomStartInt',
            'startTW': [0, 100],
            'endLimit': float('inf')
        }
    ) -> dict:

    """Add 'startTime', 'endTime' fields for given nodes

    """

    # Sanity check ============================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
        else:
            for i in nodeIDs:
                if (i not in nodes):
                    raise OutOfRangeError("ERROR: Node %s is not in `nodes`." % i)

    if (method == None or 'mode' not in method):
        raise MissingParameterError("ERROR: Missing required field `method`, or missing required key 'mode' in `method`.")
    
    if (method['mode'] == 'RandomStartInt'):
        if ('startTW' not in method or type(method['startTW']) not in [list, tuple] or len(method['startTW']) != 2):
            raise MissingParameterError("ERROR: Missing required key 'startTW' in `method` or incorrect format.")
        startTW = list(method['startTW'])

        endLimit = None
        if ('endLimit' not in method):
            endLimit = float('inf')
        else:
            endLimit = method['endLimit']
        if (startTW[1] > endLimit):
            warnings.warn("WARNING: Start time cannot be later than end time, auto adjust last posible start time to be 'endLimie'")
            startTW[1] = endLimit

        for n in nodeIDs:
            start = random.randrange(int(startTW[0]), int(startTW[1]))
            nodes[n]['timeWindow'] = [start, endLimit]

    elif (method['mode'] == 'RandomStart'):
        startTW = None
        if ('startTW' not in method or type(method['startTW']) not in [list, tuple] or len(method['startTW']) != 2):
            raise MissingParameterError("ERROR: Missing required key 'startTW' in `method` or incorrect format.")
        startTW = list(method['startTW'])

        endLimit = None
        if ('endLimit' not in method):
            endLimit = float('inf')
        else:
            endLimit = method['endLimit']
        if (startTW[1] > endLimit):
            warnings.warn("WARNING: Start time cannot be later than end time, auto adjust last posible start time to be 'endLimie'")
            startTW[1] = endLimit

        for n in nodeIDs:
            start = startTW[0] + (startTW[1] - startTW[0]) * random.random()
            nodes[n]['timeWindow'] = [start, endLimit]

    elif (method['mode'] == 'Random'):
        if ('timeWindow' not in method or type(method['timeWindow']) not in [list, tuple] or len(method['timeWindow']) != 2):
            raise MissingParameterError("ERROR: Missing required key 'timeWindow' in `method` or incorrect format.")
        timeWindow = method['timeWindow']

        if (timeWindow[0] >= timeWindow[1]):
            raise UnsupportedInputError("ERROR: End time should be later than start time.")

        for n in nodeIDs:
            rnd1 = timeWindow[0] + (timeWindow[1] - timeWindow[0]) * random.random()
            rnd2 = timeWindow[0] + (timeWindow[1] - timeWindow[0]) * random.random()
            startTime = min(rnd1, rnd2)
            endTime = max(rnd1, rnd2)
            nodes[n]['timeWindow'] = [startTime, endTime]

    elif (method['mode'] == 'RandomInt'):
        if ('timeWindow' not in method or type(method['timeWindow']) not in [list, tuple] or len(method['timeWindow']) != 2):
            raise MissingParameterError("ERROR: Missing required key 'timeWindow' in `method` or incorrect format.")
        timeWindow = method['timeWindow']

        if (timeWindow[0] >= timeWindow[1]):
            raise UnsupportedInputError("ERROR: End time should be later than start time.")

        for n in nodeIDs:
            rnd1 = random.randrange(timeWindow[0], timeWindow[1])
            rnd2 = random.randrange(timeWindow[0], timeWindow[1])
            startTime = min(rnd1, rnd2)
            endTime = max(rnd1, rnd2)
            nodes[n]['timeWindow'] = [startTime, endTime]

    else:
        raise UnsupportedInputError("ERROR: Unsupported option for `method`. Supported 'mode' includes: 'Random', 'RandomInt', 'RandomStart', 'RandomStartInt'.")

    return nodes

def rndDemandNodes(
    nodes: dict,
    nodeIDs: list[int|str]|str = 'All',
    method: dict = {
            'mode': 'RandomInt',
            'range': [0, 100]
        }
    ) -> dict:

    """Assign 'demand' for given nodes

    """

    # Sanity check ============================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
        else:
            for i in nodeIDs:
                if (i not in nodes):
                    raise OutOfRangeError("ERROR: Node %s is not in `nodes`." % i)

    if (method == None or 'mode' not in method):
        raise MissingParameterError("ERROR: Missing required field `method`, or missing required key 'mode' in `method`.")

    if (method['mode'] == 'Random'):
        if ('range' not in method or type(method['range']) not in [list, tuple] or len(method['range']) != 2):
            raise MissingParameterError("ERROR: Missing required key 'range' in `method` or incorrect format.")
        Range = method['range']

        if (Range[0] >= Range[1]):
            raise UnsupportedInputError("ERROR: Upper range should be greater than lower range.")

        for n in nodeIDs:
            rnd = Range[0] + (Range[1] - Range[0]) * random.random()
            nodes[n]['demand'] = rnd

    elif (method['mode'] == 'RandomInt'):
        if ('range' not in method or type(method['range']) not in [list, tuple] or len(method['range']) != 2):
            raise MissingParameterError("ERROR: Missing required key 'range' in `method` or incorrect format.")
        Range = method['range']

        if (Range[0] >= Range[1]):
            raise UnsupportedInputError("ERROR: Upper range should be greater than lower range.")

        for n in nodeIDs:
            rnd = random.randrange(Range[0], Range[1])
            nodes[n]['demand'] = rnd

    else:
        raise UnsupportedInputError("ERROR: Unsupported option for `method`. Supported 'mode' includes: 'Random', and 'RandomInt'.")

    return nodes

def rndPlainArcs(
    A: int|None = None,
    arcIDs: list[int|str] = [],
    method: dict = {
            'distr': 'UniformLengthInSquareXY',
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
    method: dictionary, optional, default {'distr': 'UniformLengthInSquareXY', 'xRange': (0, 100), 'yRange': (0, 100), 'minLen': 0, 'maxLen': 10}
        Spatial distribution of arcs, optional are as following:
        1) (default) Uniformly sample from a square on the Euclidean space, with uniformly selected length
        >>> method = {
        ...     'distr': 'UniformLengthInSquareXY',
        ...     'xRange': (0, 100),
        ...     'yRange': (0, 100),
        ...     'minLen': 0,
        ...     'maxLen': 10
        ... }

    """

    # Sanity check ============================================================
    if (method == None or 'distr' not in method):
        raise MissingParameterError(ERROR_MISSING_ARCS_DISTR)

    arcs = {}
    if (arcIDs == [] and A == None):
        raise MissingParameterError(ERROR_MISSING_N)
    elif (arcIDs == [] and A != None):
        arcIDs = [i for i in range(A)]

    # Generate instance =======================================================
    if (method['distr'] == 'UniformLengthInSquareXY'):
        if ('minLen' not in method or 'maxLen' not in method):
            raise MissingParameterError("ERROR: Missing required field 'minLen' and/or 'maxLen' in `method`")
        xRange = None
        yRange = None
        if ('xRange' not in method or 'yRange' not in method):
            xRange = [0, 100]
            yRange = [0, 100]
            warnings.warn("WARNING: Set sample area to be default as a (0, 100) x (0, 100) square")
        else:
            xRange = [float(method['xRange'][0]), float(method['xRange'][1])]
            yRange = [float(method['yRange'][0]), float(method['yRange'][1])]
        for n in arcIDs:
            arcs[n] = {
                'arc': _rndArcUniformSquareXY(xRange, yRange, method['minLen'], method['maxLen'])
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

