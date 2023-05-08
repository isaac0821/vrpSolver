import geopy
import math
import random
import tripy
import warnings

from .common import *
from .const import *
from .geometry import *
from .msg import *
from .error import *
from .relation import *

def rndPlainNodes(
    N:          "Number of vertices" = None,
    nodeIDs:    "(alternative of `N`) A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
    distr:      "Spatial distribution of nodes" = None
    ) -> "Randomly create a set of locations":

    """Randomly create a set of locations

    Parameters
    ----------

    N: integer, optional, default None
        Number of locations/vertices/customers to be randomly created
    nodeIDs: list, optional, default None
        Alternative input parameter of `N`. A list of node IDs, `N` will be overwritten if `nodeIDs` is given
    distr: dictionary, required, default {'method': 'uniformSquareXY', 'xRange': [0, 100], 'yRange': [0, 100]}
        Spatial distribution of nodes, options are as following:
            1) (default) Uniformly sample from a square on the Euclidean space
            >>> distr = {
            ...     'method': 'uniformSquareXY', 
            ...     'xRange': (0, 100), # A 2-tuple with minimum/maximum range of x, default as (0, 100), 
            ...     'yRange': (0, 100), # A 2-tuple with minimum/maximum range of y, default as (0, 100), 
            ... }
            2) Uniformly sample from a given polygon on the Euclidean space
            >>> distr = {
            ...     'method': 'uniformPolyXY', 
            ...     'polyXY': poly, # polygon of the area, (no holes)
            ...     'polyXYs': polys, # alternative option for 'polyXY', as a list of polygons 
            ... }
            3) Uniformly sample from a circle on the Euclidean space
            >>> distr = {
            ...     'method': 'uniformCircleXY',
            ...     'centerXY': (0, 0), # centering location, default as (0, 0), 
            ...     'radius': 100, # radius of the circle , default as 100
            ... }
            4) Uniformly sample from a given polygon by lat/lon
            >>> distr = {
            ...     'method': 'uniformPolyLatLon', 
            ...     'polyLatLon': polygon of the area, (no holes)
            ...     'polyLatLons': alternative option for 'polyLatLon', as a list of polygons 
            ... }
            5) Uniformly sample from a given circle by lat/lon,
            >>> distr = {
            ...     'method': 'uniformCircleLatLon', 
            ...     'centerLatLon': required, centering location in lat/lon, 
            ...     'radiusInMeters': radius of the circle in meters 
            ... }
            6) Uniformly generate from a given polygon on a road network
            >>> distr = {
            ...     'method': 'roadNetworkPolyLatLon'
            ...     'roadNetwork': list of arcs that can be sampled 
            ...     'polyLatLon': nodes should generated within the polygon, if not provided, will consider the entire network, 
            ...     'roadClass': list of classes that can be sampled 
            ... }
            7) Uniformly generate from a given circle on a road network
            >>> distr = {
            ...     'method': 'roadNetworkCircleLatLon', 
            ...     'roadNetwork': list of arcs that can be sampled 
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
    NotAvailableError
        Functions/options that are not ready yet.
    EmptyError
        The sample area is empty.
    """

    # Sanity checks ===========================================================
    if (N == None and nodeIDs == None):
        raise MissingParameterError(ERROR_MISSING_N)
    if (distr == None or 'method' not in distr):
        raise MissingParameterError(ERROR_MISSING_NODES_DISTR)

    # Initialize ==============================================================
    nodes = {}
    if (nodeIDs == None):
        nodeIDs = [i for i in range(N)]

    # Generate instance =======================================================
    # Uniformly sample from a square on the Euclidean space
    if (distr['method'] == 'uniformSquareXY'):
        xRange = None
        yRange = None
        if ('xRange' not in distr or 'yRange' not in distr):
            xRange = [0, 100]
            yRange = [0, 100]
            warnings.warn("Set sampled area to be default as a (0, 100) x (0, 100) square")
        else:
            xRange = distr['xRange']
            yRange = distr['yRange']
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformSquareXY(xRange, yRange)
            }

    # Uniformly sample from a polygon/a list of polygons on the Euclidean space
    elif (distr['method'] == 'uniformPolyXY'):
        if ('polyXY' not in distr and 'polyXYs' not in distr):
            raise MissingParameterError(ERROR_MISSING_NODES_DISTR_POLYXY)
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
    elif (distr['method'] == 'uniformCircleXY'):
        centerXY = None
        radius = None
        if ('centerXY' not in distr or 'radius' not in distr):
            centerXY = (0, 0),
            radius = 100
            warnings.warn("Set sample area to be default as a circle with radius of 100")
        else:
            centerXY = distr['centerXY']
            radius = dist['radius']
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformCircleXY(radius, centerXY)
            }

    # Uniformly sample from a polygon by lat/lon
    elif (distr['method'] == 'uniformPolyLatLon'):
        if ('polyLatLon' not in distr and 'polyLatLons' not in distr):
            raise MissingParameterError(ERROR_MISSING_NODES_DISTR_POLYLATLON)
        # TODO: Mercator projection
        raise NotAvailableError("`uniformPolyLatLon` is not available yet, please stay tune.")

    # Uniformly sample from a circle by lat/lon
    elif (distr['method'] == 'uniformCircleLatLon'):
        if ('centerLatLon' not in distr or 'radiusInMeters' not in distr):
            raise MissingParameterError(ERROR_MISSING_NODES_DISTR_CIRCLELATLON)
        for n in nodeIDs:
            nodes[n] = {
                'loc': _rndPtUniformCircleLatLon(distr['radiusInMeters'], distr['centerLatLon'])
            }

    # Uniformly sample from the roads/streets within a polygon/a list of polygons from given road networks
    elif (distr['method'] == 'roadNetworkPolyLatLon'):
        if ('polyLatLon' not in distr):
            raise MissingParameterError(ERROR_MISSING_NODES_DISTR_POLYLATLON)
        elif ('roadNetwork' not in distr):
            raise MissingParameterError(ERROR_MISSING_NODES_DISTR_ROADNETWORK)
        elif ('roadClass' not in distr):
            warnings.warn("WARNING: Set `roadClass` to be default as ['residential']")
        nodeLocs = _rndPtRoadNetworkPolyLatLon(
            N if N != None else len(nodeIDs),
            distr['roadNetwork'], 
            distr['polyLatLon'] if 'polyLatLon' in distr else None,
            distr['roadClass'] if 'roadClass' in distr else ['residential'])
        for n in range(len(nodeIDs)):
            nodes[nodeIDs[n]] = {
                'loc': nodeLocs[n]
            }

    # Uniformly sample from the roads/streets within a circle from given road network
    elif (distr['method'] == 'roadNetworkCircleLatLon'):
        if ('centerLatLon' not in distr or 'radiusInMeters' not in distr):
            raise MissingParameterError(ERROR_MISSING_NODES_DISTR_CIRCLELATLON)
        elif ('roadNetwork' not in distr):
            raise MissingParameterError(ERROR_MISSING_NODES_DISTR_ROADNETWORK)
        elif ('roadClass' not in distr):
            warnings.warn("WARNING: Set `roadClass` to be default as ['residential']")
        nodeLocs = _rndPtRoadNetworkCircleLatLon(
            N if N != None else len(nodeIDs),
            distr['roadNetwork'], 
            distr['centerLatLon'],
            distr['radiusInMeters'],
            distr['roadClass'] if 'roadClass' in distr else ['residential'])
        for n in range(len(nodeIDs)):
            nodes[nodeIDs[n]] = {
                'loc': nodeLocs[n]
            }
    else:
        raise MissingParameterError(ERROR_MISSING_NODES_DISTR)

    return nodes

def _rndPtUniformSquareXY(
    xRange:    "The range of x coordinates",
    yRange:    "The range of y coordinates"
    ) -> "Given the range of x, y, returns a random point in the square defined by the ranges":
    x = random.randrange(xRange[0], xRange[1])
    y = random.randrange(yRange[0], yRange[1])
    return (x, y)

def _rndPtUniformTriangleXY(
    triangle:   "The triangle for generating random points"
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

def _rndPtUniformPolyXY(
    poly:       "The polygon for generating random points"
    ) -> "Given a polygon, generate a random point in the polygons uniformly":

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

def _rndPtUniformPolyXYs(
    polys:       "A list of polygons for generating random points"
    ) -> "Given a list of polygons, generate a random point in the polygons uniformly":

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

def _rndPtUniformCircleXY(
    radius:     "Radius of the circle",
    center:  "Center location of the circle"
    ):
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, radius ** 2))
    x = center[0] + r * math.cos(theta)
    y = center[1] + r * math.sin(theta)
    return (x, y)

def _rndPtUniformCircleLatLon(
    radius:     "Radius of the circle",
    center:  "Center location of the circle"
    ):
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, radius ** 2))
    (lat, lon) = ptInDistLatLon(center, theta, r)
    return (lat, lon)

def _rndPtRoadNetworkPolyLatLon(
    N:          "Number of nodes",
    road:       "Dictionary of road network in the format of \
                {\
                    roadID: {\
                        'shape': [[lat, lon], [lat, lon], ...],\
                        'class': 'residential', etc. as categorized in OpenStreetMap,\
                    }\
                }",
    poly:       "Nodes should also within this polygon",
    roadClass:  "List of classes that can be sampled"
    ) -> "Given a road network, generate customers that locates on the road network":
    
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
                    if (isPtOnPoly(road[rID]['shape'][i], poly)):
                        includedFlag = True
                        break
        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(road[rID]['shape']) - 1):
                roadLength += distLatLon(road[rID]['shape'][i], road[rID]['shape'][i + 1])
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(rID)

    # Check if there are roads included =======================================
    if (sum(lengths) == 0):
        raise EmptyError("No road is found.")

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
            (lat, lon) = mileageInPathLatLon(road[roadIDs[idx]]['shape'], edgeDist)
        else:
            insideFlag = False
            while (not insideFlag):
                idx = rndPick(lengths)
                edgeLength = lengths[idx]
                edgeDist = random.uniform(0, 1) * edgeLength
                (lat, lon) = mileageInPathLatLon(road[roadIDs[idx]]['shape'], edgeDist)
                if (isPtOnPoly([lat, lon], poly)):
                    insideFlag = True
        nodeLocs.append((lat, lon))
    return nodeLocs

def _rndPtRoadNetworkCircleLatLon(
    N:          "Number of nodes",
    road:       "Dictionary of road network in the format of \
                {\
                    roadID: {\
                        'shape': [[lat, lon], [lat, lon], ...]\
                    }\
                }",\
    center:  "Center location",
    radius:     "Radius in [m]",
    roadClass:  "List of classes that can be sampled"
    ) -> "Given a road network, generate customers that locates on the road network":
    
    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for rID in road:
        roadLength = 0
        includedFlag = False
        for i in range(len(road[rID]['shape'])):
            if (road[rID]['class'] in roadClass and distLatLon(road[rID]['shape'][i], center) <= radius):
                includedFlag = True
                break

        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(road[rID]['shape']) - 1):
                roadLength += distLatLon(road[rID]['shape'][i], road[rID]['shape'][i + 1])
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(rID)


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
            (lat, lon) = mileageInPathLatLon(road[roadIDs[idx]]['shape'], edgeDist)
            if (distLatLon([lat, lon], center) <= radius):
                insideFlag = True
        nodeLocs.append((lat, lon))
    return nodeLocs

