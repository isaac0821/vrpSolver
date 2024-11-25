import math
import random
import warnings

from .common import *
from .geometry import *
from .msg import *
from .road import *

def rndLocs(
    N: int, 
    distr = 'UniformSquareXY', 
    **kwargs) -> list:

    """Randomly create a list of N locations

    Parameters
    ----------

    N: integer, required
        Number of locations/vertices/customers to be randomly created
    distr: string, optional, default as 'UniformSquareXY'
        Spatial distribution of locations, options and required additional inputs are as follows:

        1) (default) 'UniformSquareXY', uniformly sample from a square on the Euclidean space
            - xRange: 2-tuple, with minimum/maximum range of x, default as (0, 100)
            - yRange: 2-tuple, with minimum/maximum range of y, default as (0, 100)
        2) 'UniformPolyXY', uniformly sample from a given polygon
            - polyXY: poly, the polygon of the area, (no holes)
            - polyXYs: list of polys, alternative option for `polyXY`
        3) 'UniformAvoidPolyXY', uniformly sample from a square avoiding some polygons
            - xRange: 2-tuple, with minimum/maximum range of x, default as (0, 100)
            - yRange: 2-tuple, with minimum/maximum range of y, default as (0, 100)
            - polyXY: poly, the polygon of the area, (no holes)
            - polyXYs: list of polys, alternative option for `polyXY`
        4) 'UniformCircleXY', uniformly sample from a circle on the Euclidean space
            - centerXY: 2-tuple, the center of circle
            - radius: float, the radius of the circle
        5) 'UniformPolyLatLon', uniformly sample from a polygon by lat/lon
            - polyLatLon: poly, the polygon of the area, (no holes)
            - polyLatLons: list of polys, alternative option for `polyLatLon`
        6) 'UniformCircleLatLon', uniformly sample from a circle by lat/lon
            - centerLatLon: 2-tuple, the (lat, lon) for the center
            - radiusInMeters: float, the radius of the circle in meters
        7) 'RoadNetworkPolyLatLon', uniformly generate within a given polygon on a road network
            - roads: dict, the road network dictionary
            - polyLatLon: poly, optional, the polygon on the map to sample from
            - polyLatLons: list of polys, optional, alternative for `polyLatLon`
            - roadClass: list[str], the road classes that allows to sample from
        8) 'RoadNetworkCircleLatLon', uniformly generate within a circle on a road network
            - roads: dict, the road network dictionary
            - centerLatLon: 2-tuple, the (lat, lon) for the center
            - radiusInMeters: float, the radius of the circle in meters
            - roadClass: list[str], the road classes that allows to sample from
    **kwargs: optional
        Provide additional inputs for different `distr` options

    Returns
    -------
    list
        A list of randomly created locations

    Raises
    ------
    MissingParameterError
        Missing required inputs in `**kwargs`.
    UnsupportedInputError
        Option is not supported for `distr`
    NotAvailableError
        Functions/options that are not ready yet.
    EmptyError
        The sample area is empty.
    """

    nodeLocs = []
    # Uniformly sample from a square on the Euclidean space
    if (distr == 'UniformSquareXY'):
        xRange = None
        yRange = None
        if ('xRange' not in kwargs or 'yRange' not in kwargs):
            xRange = [0, 100]
            yRange = [0, 100]
            warnings.warn("WARNING: Set sampled area to be default as a (0, 100) x (0, 100) square")
        else:
            xRange = [float(kwargs['xRange'][0]), float(kwargs['xRange'][1])]
            yRange = [float(kwargs['yRange'][0]), float(kwargs['yRange'][1])]
        for n in range(N):
            nodeLocs.append(_rndPtUniformSquareXY(xRange, yRange))

    elif (distr == 'UniformCubeXYZ'):
        xRange = None
        yRange = None
        zRange = None
        if ('xRange' not in kwargs or 'yRange' not in kwargs or 'zRange' not in kwargs):
            xRange = [0, 100]
            yRange = [0, 100]
            zRange = [0, 100]
            warnings.warn("WARNING: Set sampled area to be default as a (0, 100) x (0, 100) x (0, 100) cube")
        else:
            xRange = [float(kwargs['xRange'][0]), float(kwargs['xRange'][1])]
            yRange = [float(kwargs['yRange'][0]), float(kwargs['yRange'][1])]
            zRange = [float(kwargs['zRange'][0]), float(kwargs['zRange'][1])]
        for n in range(N):
            nodeLocs.append(_rndPtUniformCubeXYZ(xRange, yRange, zRange))

    # Uniformly sample from a polygon/a list of polygons on the Euclidean space
    elif (distr == 'UniformPolyXY'):
        if ('polyXY' not in kwargs and 'polyXYs' not in kwargs):
            raise MissingParameterError("ERROR: Missing required args 'polyXY' or 'polyXYs', which indicates a polygon / a list of polygons in the Euclidean space")
        if ('polyXY' in kwargs):
            for n in range(N):
                nodeLocs.append(_rndPtUniformPolyXY(kwargs['polyXY']))
        elif ('polyXYs' in kwargs):
            for n in range(N):
                nodeLocs.append(_rndPtUniformPolyXYs(kwargs['polyXY']))

    # Uniformly sample from the Euclidean space avoiding polygons
    elif (distr == 'UniformAvoidPolyXY'):
        if ('polyXY' not in kwargs and 'polyXYs' not in kwargs):
            raise MissingParameterError("ERROR: Missing required args 'polyXY' or 'polyXYs', which indicates a polygon / a list of polygons in the Euclidean space")
        xRange = None
        yRange = None
        if ('xRange' not in kwargs or 'yRange' not in kwargs):
            xRange = [0, 100]
            yRange = [0, 100]
            warnings.warn("WARNING: Set sampled area to be default as a (0, 100) x (0, 100) square")
        else:
            xRange = [float(kwargs['xRange'][0]), float(kwargs['xRange'][1])]
            yRange = [float(kwargs['yRange'][0]), float(kwargs['yRange'][1])]
        if ('polyXY' in kwargs):
            for n in range(N):
                nodeLocs.append(_rndPtUniformAvoidPolyXY(kwargs['polyXY'], xRange, yRange))
        elif ('polyXYs' in kwargs):
            for n in range(N):
                nodeLocs.append(_rndPtUniformAvoidPolyXYs(kwargs['polyXYs'], xRange, yRange))

    # Uniformly sample from a circle on the Euclidean space
    elif (distr == 'UniformCircleXY'):
        centerXY = None
        radius = None
        if ('centerXY' not in kwargs or 'radius' not in kwargs):
            centerXY = (0, 0)
            radius = 100
            warnings.warn("WARNING: Set sample area to be default as a circle with radius of 100 centering at (0, 0)")
        else:
            centerXY = kwargs['centerXY']
            radius = kwargs['radius']
        for n in range(N):
            nodeLocs.append(_rndPtUniformCircleXY(radius, centerXY))

    # Uniformly sample from a polygon by lat/lon
    elif (distr == 'UniformPolyLatLon'):
        if ('polyLatLon' not in kwargs and 'polyLatLons' not in kwargs):
            raise MissingParameterError("ERROR: Missing required args 'polyXY' or 'polyXYs', which indicates a polygon / a list of polygons in the Euclidean space")
        # TODO: Mercator projection
        raise VrpSolverNotAvailableError("ERROR: 'UniformPolyLatLon' is not available yet, please stay tune.")

    # Uniformly sample from a circle by lat/lon
    elif (distr == 'UniformCircleLatLon'):
        if ('centerLatLon' not in kwargs or 'radiusInMeters' not in kwargs):
            raise MissingParameterError("ERROR: Missing required args 'centerLatLon' or 'radiusInMeters'.")
        for n in range(N):
            nodeLocs.append(_rndPtUniformCircleLatLon(kwargs['radiusInMeters'], kwargs['centerLatLon']))

    # Uniformly sample from the roads/streets within a polygon/a list of polygons from given road networks
    elif (distr == 'RoadNetworkPolyLatLon'):
        if ('polyLatLon' not in kwargs and 'polyLatLons' not in kwargs):
            raise MissingParameterError("ERROR: Missing required args 'polyLatLon' or 'polyLatLons', which indicates a polygon / a list of polygons in the Euclidean space")
        elif ('roads' not in kwargs):
            raise MissingParameterError("ERROR: Missing required args 'RoadNetwork'. Need to provide the road network where the nodes are generated.")
        elif ('roadClass' not in kwargs):
            warnings.warn("WARNING: Set 'roadClass' to be default as ['residential']")
        if ('polyLatLon' in kwargs):
            nodeLocs = _rndPtRoadNetworkPolyLatLon(
                N if N != None else len(nodeIDs),
                kwargs['roads'], 
                kwargs['polyLatLon'],
                kwargs['roadClass'] if 'roadClass' in kwargs else ['residential'])
        elif ('polyLatLons' in kwargs):
            nodeLocs = _rndPtRoadNetworkPolyLatLons(
                N if N != None else len(nodeIDs),
                kwargs['roads'], 
                kwargs['polyLatLons'],
                kwargs['roadClass'] if 'roadClass' in kwargs else ['residential'])

    # Uniformly sample from the roads/streets within a circle from given road network
    elif (distr == 'RoadNetworkCircleLatLon'):
        if ('centerLatLon' not in kwargs or 'radiusInMeters' not in kwargs):
            raise MissingParameterError("ERROR: Missing required args 'centerLatLon' or 'radiusInMeters'.")
        elif ('roads' not in kwargs):
            raise MissingParameterError("ERROR: Missing required args 'RoadNetwork'. Need to provide the road network where the nodes are generated.")
        elif ('roadClass' not in kwargs):
            warnings.warn("WARNING: Set 'roadClass' to be default as ['residential']")
        nodeLocs = _rndPtRoadNetworkCircleLatLon(
            N if N != None else len(nodeIDs),
            kwargs['roads'], 
            kwargs['radiusInMeters'],
            kwargs['centerLatLon'],
            kwargs['roadClass'] if 'roadClass' in kwargs else ['residential'])
    
    else:
        raise UnsupportedInputError(ERROR_MISSING_NODES_DISTR)

    return nodeLocs

def rndNodes(
    N: int|None = None, 
    nodeIDs: list[int|str] = [], 
    nodes: dict|None = None,
    distr = 'UniformSquareXY',
    locFieldName = 'loc',
    **kwargs
    ) -> dict:

    """Randomly create a nodes dictionary

    Parameters
    ----------

    N: integer, optional
        Number of locations/vertices/customers to be randomly created
    nodeIDs: list of int|str, optional
        A list of ids for the locations to be created, an alternative option for `N`
    nodes: dict, optional
        A nodes dictionary, if given, new locations will be append into this dictionary
    distr: string, optional, default as 'UniformSquareXY'
        See `distr` docstring in :func:`~vrpSolver.instance.rndLocs()`
    locFieldName: str, optional, default as 'loc'
        The key in nodes dictionary to indicate the locations
    **kwargs: optional
        Provide additional inputs for different `distr` options

    Returns
    -------
    list
        A list of randomly created locations

    Raises
    ------
    MissingParameterError
        Missing required inputs in **kwargs.
    UnsupportedInputError
        Option is not supported for `distr`
    NotAvailableError
        Functions/options that are not ready yet.
    EmptyError
        The sample area is empty.
    """

    # Sanity checks ===========================================================
    if (nodes == None):
        nodes = {}
    
    if (nodeIDs == [] and N == None):
        raise MissingParameterError(ERROR_MISSING_N)
    elif (nodeIDs == [] and N != None):
        nodeIDs = [i for i in range(N)]

    # Generate instance =======================================================
    nodeLocs = rndLocs(
        N = len(nodeIDs), 
        distr = distr,
        **kwargs
        )
    for n in range(len(nodeIDs)):
        if (nodeIDs[n] in nodes):
            warnings.warn("WARNING: %s already exists, will be replaced" % n)
        nodes[nodeIDs[n]] = {
            locFieldName: nodeLocs[n]
        }

    return nodes

def _rndPtUniformSquareXY(xRange: list[int]|list[float], yRange: list[int]|list[float]) -> pt:
    x = random.uniform(xRange[0], xRange[1])
    y = random.uniform(yRange[0], yRange[1])
    return (x, y)

def _rndPtUniformCubeXYZ(xRange: list[int]|list[float], yRange: list[int]|list[float], zRange: list[int]|list[float]) -> pt3D:
    x = random.uniform(xRange[0], xRange[1])
    y = random.uniform(yRange[0], yRange[1])
    z = random.uniform(zRange[0], zRange[1])
    return (x, y, z)

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
    # TODO: polyTriangulation() to be replaced
    lstTriangle = polyTriangulation(poly)

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
    # TODO: polyTriangulation() to be replaced
    lstTriangle = []
    for p in polys:
        lstTriangle.extend(polyTriangulation(p))

    # Weight them and make draws ==============================================
    lstWeight = []
    for i in range(len(lstTriangle)):
        lstWeight.append(calTriangleAreaXY(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2]))

    # Select a triangle and randomize a point in the triangle =================
    idx = rndPick(lstWeight)
    (x, y) = _rndPtUniformTriangleXY(lstTriangle[idx])

    return (x, y)

def _rndPtUniformAvoidPolyXY(poly: poly, xRange: list[int]|list[float], yRange: list[int]|list[float]) -> pt:
    # Use the low efficient accept-denial approach
    while (True):
        x = random.uniform(xRange[0], xRange[1])
        y = random.uniform(yRange[0], yRange[1])
        if (isPtInPoly((x, y), poly)):
            return (x, y)

def _rndPtUniformAvoidPolyXYs(polys: polys, xRange: list[int]|list[float], yRange: list[int]|list[float]) -> pt:
    while (True):
        x = random.uniform(xRange[0], xRange[1])
        y = random.uniform(yRange[0], yRange[1])
        notInPolys = True
        for p in polys:
            if (isPtInPoly((x, y), p)):
                notInPolys = False
                break
        if (notInPolys):
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

def _rndPtRoadNetworkPolyLatLon(N: int, roads: dict, poly: poly, roadClass: str | list[str]) -> list[pt]:
    # If poly is given, clip road networks by poly ============================
    clipRoad = {}
    if (poly != None):
        clipRoad = clipRoadsByPoly(roads, poly)
    else:
        clipRoad = roads

    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for rID in clipRoad:
        roadLength = 0
        includedFlag = False
        if ('class' in clipRoad[rID] and clipRoad[rID]['class'] in roadClass):
            for i in range(len(clipRoad[rID]['shape']) - 1):
                roadLength += distLatLon(clipRoad[rID]['shape'][i], clipRoad[rID]['shape'][i + 1])['dist']
            lengths.append(roadLength)
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
        idx = rndPick(lengths)
        edgeLength = lengths[idx]
        edgeDist = random.uniform(0, 1) * edgeLength
        (lat, lon) = ptInSeqMileage(clipRoad[roadIDs[idx]]['shape'], edgeDist, 'LatLon')
        nodeLocs.append((lat, lon))

    return nodeLocs

def _rndPtRoadNetworkPolyLatLons(N: int, roads: dict, polys: polys, roadClass: str | list[str]) -> list[pt]:
    # If poly is given, clip road networks by poly ============================
    clipRoad = {}
    if (poly != None):
        clipRoad = clipRoadsByPolys(roads, polys)
    else:
        clipRoad = roads

    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for rID in clipRoad:
        roadLength = 0
        includedFlag = False
        if ('class' in clipRoad[rID] and clipRoad[rID]['class'] in roadClass):
            for i in range(len(clipRoad[rID]['shape']) - 1):
                roadLength += distLatLon(clipRoad[rID]['shape'][i], clipRoad[rID]['shape'][i + 1])['dist']
            lengths.append(roadLength)
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
        idx = rndPick(lengths)
        edgeLength = lengths[idx]
        edgeDist = random.uniform(0, 1) * edgeLength
        (lat, lon) = ptInSeqMileage(clipRoad[roadIDs[idx]]['shape'], edgeDist, 'LatLon')
        nodeLocs.append((lat, lon))

    return nodeLocs

def _rndPtRoadNetworkCircleLatLon(N: int, roads: dict, radius: float, center: pt, roadClass: str | list[str]) -> list[pt]:
    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for rID in roads:
        roadLength = 0
        includedFlag = False
        for i in range(len(roads[rID]['shape'])):
            if ('class' in roads[rID] and roads[rID]['class'] in roadClass and distLatLon(roads[rID]['shape'][i], center)['dist'] <= radius):
                includedFlag = True
                break

        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(roads[rID]['shape']) - 1):
                roadLength += distLatLon(roads[rID]['shape'][i], roads[rID]['shape'][i + 1])['dist']
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(rID)


    # Check if there are roads included =======================================
    if (sum(lengths) == 0):
        raise EmptyError("ERROR: No road is found.")

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
            (lat, lon) = ptInSeqMileage(roads[roadIDs[idx]]['shape'], edgeDist, 'LatLon')
            if (distLatLon([lat, lon], center)['dist'] <= radius):
                insideFlag = True
        nodeLocs.append((lat, lon))

    return nodeLocs

def rndNodeNeighbors(
    nodes: dict,
    nodeIDs: list[int|str]|str = 'All', 
    shape: str = 'Circle',
    locFieldName = 'loc',
    neighborFieldName = 'neighbor',
    **kwargs
    ) -> dict:

    """Given a node dictionary, create neighborhood to selected nodes

    WARNING
    -------    
    This function will modify the input dictionary `nodes`

    Parameters
    ----------
    nodes: dictionary, required
        A plain nodes dictionary to add neighborhoods.
    nodeIDs: string|list[int|str], optional, default 'All'
        A list of node IDs to add neighborhood, leave it as 'All' to indicate adding such information to all nodes.
    shape: str, optional, default as 'Circle'
        The shape of neighborhoods, options and required addtional inputs are as follows:

        1) (default) 'Cirlce', add circle surrounding nodes
            - 'radius': The radius, default as 1
            - 'lod': The level of details, circle will be approximated as a x-gon polygon, default as 30
        2) 'Poly', add polygon surrounding nodes
            - 'poly': In relative axis where node locates in [0, 0]
        3) 'Egg', add egg shape to nodes. The curve function: :math:`\\frac{x^2}{(a - b)x + ab} + \\frac{y^2}{c^2} = 1`
            - 'a': required
            - 'b': required
            - 'c': required
            - 'direction': default as 0
            - 'lod': default as 30
        4) 'RndSquare', add random size squares around nodes
            - 'minLen': required, minimum length
            - 'maxLen': required, maximum length
        4) 'RndCurvy', add random curvy shapes around nodes
            - 'maxRadius': default as 1.2
            - 'minRadius': default as 0.8
            - 'N': default as 5
            - 'w': default as 3
            - 'lod': default as 30
        5) 'RndConvexPoly', add convex polygons with random size arond nodes
            - 'maxNumSide': maximum number of sides
            - 'maxDiag': maximum length of the diagonal
            - 'minDiag': minimum length of the diagonal
    **kwargs: optional
        Provide additional inputs for different `distr` options

    Returns
    -------
    dict
        Changes will apply to the original `nodes` dictionary

    """

    # Sanity check ============================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
        else:
            for i in nodeIDs:
                if (i not in nodes):
                    raise OutOfRangeError("ERROR: Node %s is not in `nodes`." % i)
    
    # Add neighborhood by 'shape' =============================================
    if (shape == 'Poly'):
        for n in nodeIDs:
            if ('poly' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'poly'")
            poly = [[i[0] + nodes[n][locFieldName][0], i[1] + nodes[n][locFieldName][1]] for i in kwargs['poly']]
            nodes[n][neighborFieldName] = [poly[i] for i in range(len(poly)) if distEuclideanXY(poly[i], poly[i - 1])['dist'] > CONST_EPSILON]
            
    elif (shape == 'Circle'):
        for n in nodeIDs:
            if ('radius' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'radius'")
            # By default, a circle is plotted by a 30-gon
            lod = 30
            if ('lod' in kwargs and type(kwargs['lod']) == int):
                lod = kwargs['lod']
            poly = [[
                nodes[n][locFieldName][0] + kwargs['radius'] * math.sin(2 * d * math.pi / lod),
                nodes[n][locFieldName][1] + kwargs['radius'] * math.cos(2 * d * math.pi / lod),
            ] for d in range(lod + 1)]
            nodes[n][neighborFieldName] = [poly[i] for i in range(len(poly)) if distEuclideanXY(poly[i], poly[i - 1])['dist'] > CONST_EPSILON]
        
    elif (shape == 'Egg'):
        for n in nodeIDs:
            if ('a' not in kwargs or 'b' not in kwargs or 'c' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'a', 'b', and/or 'c'.")
            direction = 0
            if ('direction' in kwargs):
                direction = kwargs['direction']
            lod = 30
            if ('lod' in kwargs and type(kwargs['lod']) == int):
                lod = kwargs['lod']
            # Formulation:
            # \frac{x^2}{(a - b)x + ab} + \frac{y^2}{c^2} = 1
            a = kwargs['a']
            b = kwargs['b']
            c = kwargs['c']
            
            vHLod = math.ceil(lod * 2 / 9)
            vTLod = math.ceil(lod / 9)
            hLod = math.ceil(lod * 2 / 3)

            polyL = []
            polyM = []
            polyR = []
            for d in range(vHLod + 1):
                y = c * 0.75 * d / vHLod
                A = 1
                B = (y ** 2 / c ** 2 - 1) * (a - b)
                C = (y ** 2 / c ** 2 - 1) * a * b
                X = (-B - math.sqrt(B ** 2 - 4 * A * C)) / (2 * A)
                polyL.append((X, y))
                xStart = X
            for d in range(vTLod + 1):
                y = c * 0.4 * d / vHLod
                A = 1
                B = (y ** 2 / c ** 2 - 1) * (a - b)
                C = (y ** 2 / c ** 2 - 1) * a * b
                X = (-B + math.sqrt(B ** 2 - 4 * A * C)) / (2 * A)
                polyR.insert(0, (X, y))
                xEnd = X
            for d in range(hLod + 1):
                x = xStart + (xEnd - xStart) * d / hLod
                Y = math.sqrt(c * c * (1 - (x * x) / ((a - b) * x + a * b)))
                polyM.append((x, Y))
            polyHf = []
            polyHf.extend(polyL)
            polyHf.extend(polyM)
            polyHf.extend(polyR)

            polyB4Rot = [i for i in polyHf]
            polyB4Rot.extend([(polyHf[len(polyHf) - 1 - k][0], - polyHf[len(polyHf) - 1 - k][1]) for k in range(len(polyHf))])
            
            poly = []
            for d in range(len(polyB4Rot)):
                di = headingXY((0, 0), polyB4Rot[d])
                r = distEuclideanXY((0, 0), polyB4Rot[d])['dist']
                pt = ptInDistXY(nodes[n][locFieldName], di + direction, r)
                poly.append(pt)

            nodes[n][neighborFieldName] = [poly[i] for i in range(len(poly)) if distEuclideanXY(poly[i], poly[i - 1])['dist'] > CONST_EPSILON]
            
    elif (shape == 'RndSquare'):
        for n in nodeIDs:
            if ('maxLen' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'maxLen'")
            if ('minLen' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'minLen'")
            if (kwargs['minLen'] > kwargs['maxLen']):
                warnings.warn("WARNING: 'minLen' is greater than 'maxLen', will be swapped")
                kwargs['maxLen'], kwargs['minLen'] = kwargs['minLen'], kwargs['maxLen']
            
            width = random.uniform(kwargs['minLen'], kwargs['maxLen'])
            height = random.uniform(kwargs['minLen'], kwargs['maxLen'])

            nodes[n][neighborFieldName] = [
                [nodes[n][loc][0] - width / 2, nodes[n][loc][1] - height / 2], 
                [nodes[n][loc][0] + width / 2, nodes[n][loc][1] - height / 2], 
                [nodes[n][loc][0] + width / 2, nodes[n][loc][1] + height / 2], 
                [nodes[n][loc][0] - width / 2, nodes[n][loc][1] + height / 2]
            ]

    elif (shape == 'RndCurvy'):
        for n in nodeIDs:
            if ('maxRadius' not in kwargs or 'minRadius' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'maxRadius' or 'minRadius'")
            lod = 30
            if ('lod' in kwargs and type(kwargs['lod']) == int):
                lod = kwargs['lod']

            r = []
            for i in range(lod + 1):
                r.append(kwargs['minRadius'])
            N = 4
            if ('N' in kwargs and type(kwargs['N']) == int):
                N = kwargs['N']
            w = 3
            if ('w' in kwargs and type(kwargs['w']) == int):
                w = kwargs['w']

            for k in range(N):
                a = random.uniform(0, 1)
                b = random.randint(1, w)
                c = random.uniform(0, 2)
                for i in range(lod + 1):
                    r[i] += a * math.sin(b * 2 * i * math.pi / lod + math.pi * c)

            maxRI = max(r)
            for i in range(len(r)):
                r[i] = r[i] * (kwargs['maxRadius'] - kwargs['minRadius']) / maxRI

            poly = [[
                nodes[n][locFieldName][0] + (r[d] + kwargs['minRadius']) * math.sin(2 * d * math.pi / lod),
                nodes[n][locFieldName][1] + (r[d] + kwargs['minRadius']) * math.cos(2 * d * math.pi / lod),
            ] for d in range(lod + 1)]
            nodes[n][neighborFieldName] = [poly[i] for i in range(len(poly)) if distEuclideanXY(poly[i], poly[i - 1])['dist'] > CONST_EPSILON]
    
    elif (shape == 'RndConvexPoly'):
        for n in nodeIDs:
            if ('maxNumSide' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'maxNumSide'")
            if ('maxDiag' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'maxDiag'")
            if ('minDiag' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'minDiag'")
            polyPts = []
            for i in range(kwargs['maxNumSide']):
                deg = random.uniform(0, 1) * 360
                r = kwargs['minDiag'] / 2 + random.uniform(0, 1) * (kwargs['maxDiag'] - kwargs['minDiag']) / 2
                polyPts.append(ptInDistXY(
                    pt = nodes[n][locFieldName], direction = deg, dist = r))
            polyShapely = shapely.convex_hull(shapely.MultiPoint(points = polyPts))
            poly = [i for i in mapping(polyShapely)['coordinates'][0]]
            nodes[n][neighborFieldName] = [poly[i] for i in range(len(poly)) if distEuclideanXY(poly[i], poly[i - 1])['dist'] > CONST_EPSILON]

    elif (shape == 'RndStar'):
        for n in nodeIDs:
            if ('maxNumSide' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'maxNumSide'")
            if ('maxDiag' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'maxDiag'")
            if ('minDiag' not in kwargs):
                raise MissingParameterError("ERROR: Missing required args 'minDiag'")

            degs = []
            for i in range(kwargs['maxNumSide']):
                degs.append(random.uniform(0, 1) * 360)
            degs.sort()

            polyPts = []
            for i in range(kwargs['maxNumSide']):
                r = kwargs['minDiag'] / 2 + random.uniform(0, 1) * (kwargs['maxDiag'] - kwargs['minDiag']) / 2
                polyPts.append(ptInDistXY(
                    pt = nodes[n][locFieldName], direction = degs[i], dist = r))
            nodes[n][neighborFieldName] = [polyPts[i] for i in range(len(polyPts)) if distEuclideanXY(polyPts[i], polyPts[i - 1])['dist'] > CONST_EPSILON]

    else:
        raise UnsupportedInputError("ERROR: Unsupported option for `kwargs`. Supported 'shape' includes: 'Poly', 'Circle', 'Egg', 'RndSquare', 'RndConvexPoly' and 'RndCurvy'.")

    return nodes

def rndArcs(
    A: int|None = None,
    arcIDs: list[int|str] = [],
    distr = 'UniformLengthInSquareXY',
    arcFieldName: str = 'arc',
    **kwargs
    ) -> dict:

    """Randomly create a set of arcs 

    Parameters
    ----------

    A: integer, optional, default as None
        Number of arcs to be visited
    arcIDs: list, optional, default as None
        Alternative input parameter of `A`. A list of arc IDs, `A` will be overwritten if `arcIDs` is given
    distr: str, optional, default as 'UniformLengthInSquareXY'
        The distribution of arcs. Options and required additional inputs are as follows:

        1) (default) 'UniformLengthInSquareXY', uniformly sample from a square on the Euclidean space, with uniformly selected length
            - xRange: 2-tuple, with minimum/maximum range of x, default as (0, 100)
            - yRange: 2-tuple, with minimum/maximum range of y, default as (0, 100)
            - minLen: float, minimum length of the arcs
            - maxLen: float, maximum length of the arcs
    **kwargs: optional
        Provide additional inputs for different `distr` options

    Returns
    -------
    dict
        A dictionary of randomly created arcs.

    """

    # Sanity check ============================================================
    arcs = {}
    if (arcIDs == [] and A == None):
        raise MissingParameterError(ERROR_MISSING_N)
    elif (arcIDs == [] and A != None):
        arcIDs = [i for i in range(A)]

    # Generate instance =======================================================
    if (distr == 'UniformLengthInSquareXY'):
        if ('minLen' not in kwargs or 'maxLen' not in kwargs):
            raise MissingParameterError("ERROR: Missing required field 'minLen' and/or 'maxLen'")
        if ('minDeg' not in kwargs):
            kwargs['minDeg'] = 0
        if ('maxDeg' not in kwargs):
            kwargs['maxDeg'] = 360
        xRange = None
        yRange = None
        if ('xRange' not in kwargs or 'yRange' not in kwargs):
            xRange = [0, 100]
            yRange = [0, 100]
            warnings.warn("WARNING: Set sample area to be default as a (0, 100) x (0, 100) square")
        else:
            xRange = [float(kwargs['xRange'][0]), float(kwargs['xRange'][1])]
            yRange = [float(kwargs['yRange'][0]), float(kwargs['yRange'][1])]
        for n in arcIDs:
            arcs[n] = {
                arcFieldName : _rndArcUniformSquareXY(xRange, yRange, kwargs['minLen'], kwargs['maxLen'], kwargs['minDeg'], kwargs['maxDeg'])
            }
    else:
        raise UnsupportedInputError(ERROR_MISSING_ARCS_DISTR)

    return arcs

def _rndArcUniformSquareXY(xRange, yRange, minLen, maxLen, minDeg, maxDeg) -> tuple[pt, pt]:
    length = random.uniform(minLen, maxLen)
    direction = random.uniform(minDeg, maxDeg)
    xStart = random.uniform(xRange[0], xRange[1])
    yStart = random.uniform(yRange[0], yRange[1])
    (xEnd, yEnd) = ptInDistXY((xStart, yStart), direction, length)
    return ((xStart, yStart), (xEnd, yEnd))

def rndPolys(
    P: int|None = None,
    polyIDs: list[int|str]|None = None,
    distr = 'UniformSquareXY',
    shape = 'Circle',
    anchorFieldName = 'anchor',
    polyFieldName = 'poly',    
    allowOverlapFlag = True,
    returnAsListFlag = True,
    **kwargs
    ) -> dict:

    """
    Randomly create polygons

    Parameters
    ----------

    P: int|str, optional, default as None
        Number of polygons to create
    polyIDs: list[int|str]|None, optional, default as None
        A list of ids for the polygons to be created, an alternative option for `P`
    distr: str, optional, default as 'UniformSquareXY'
        Anchor locations of each polygon. Options and required additional information are referred to :func:`~vrpSolver.instance.rndLocs()`.
    shape: str, optional, default as 'Circle',
        Shape of the polygons. Options and required additional information are referred to :func:`~vrpSolver.instance.rndNodeNeighbors()`.
    anchorFieldName: str, optional, default as 'anchor'
        The key value of the anchor location
    polyFieldName: str, optional, default as 'poly',
        The key value of the polygons
    allowOverlapFlag: bool, optional, default as True
        True if allows the polygons to overlap
    returnAsListFlag: bool, optional, default as True
        True if returns a list of polygons instead of a dictionary

    Returns
    -------
    dict
        A dictionary with polygon information
    """

    # Sanity check ============================================================
    if (polyIDs == None and P == None):
        raise MissingParameterError("ERROR: Missing required field `P` and `polyIDs`.")
    
    elif (polyIDs == None and P != None):
        polyIDs = [i for i in range(P)]

    # If overlapping is allowed ===============================================
    if (allowOverlapFlag):
        polygons = rndNodes(
            N = P,
            nodeIDs = polyIDs,
            distr = distr,
            locFieldName = anchorFieldName,
            **kwargs)

        # Next, create P polygons relative to anchor points ===================
        polygons = rndNodeNeighbors(
            nodes = polygons,
            shape = shape,
            locFieldName = anchorFieldName,
            neighborFieldName = polyFieldName,
            **kwargs)

        if (not returnAsListFlag):
            return polygons
        else:
            return [polygons[i][polyFieldName] for i in polygons]
    
    # If overlapping is not allowed ===========================================
    maxNumOfFailedTrial = 20
    anchor = []
    polys = []

    numOfFailedTrial = 0
    while (len(polys) < len(polyIDs)):
        addedFlag = True
        p = rndNodes(
            N = 1,
            distr = distr,
            **kwargs)
        p = rndNodeNeighbors(
            nodes = p,
            shape = shape,
            **kwargs)
        newPoly = p[0]['neighbor']
        for poly in polys:
            if (isPolyIntPoly(poly, newPoly) == True):
                numOfFailedTrial += 1
                addedFlag = False
                break

        if (addedFlag == True):
            polys.append(newPoly)
            anchor.append(p[0]['loc'])
            numOfFailedTrial = 0
        else:
            if (numOfFailedTrial >= maxNumOfFailedTrial):
                break
    if (len(polys) < len(polyIDs)):
        warnings.warn("WARNING: Space is too limited, only %s of polygons are created." % len(polys))

    if (returnAsListFlag):        
        return polys
    else:
        polygons = {}
        for p in range(len(polys)):
            polygons[polyIDs[p]] = {
                anchorFieldName: anchor[p],
                polyFieldName: polys[p]
            }
        return polygons