import requests

from .common import *
from .geometry import *

def matrixDist(nodes: dict, locFieldName: str = 'loc', nodeIDs: list|str = 'All', edges: str = 'Euclidean', detailFlag: bool = False, **kwargs) -> dict:
    """
    Given a `nodes` dictionary, returns the traveling matrix between nodes

    Parameters
    ----------

    nodes: dict, required
        A `nodes`dictionary with location information
    locFieldName: str, optional, default as 'loc'
        The key in nodes dictionary to indicate the locations
    nodeIDs: list of int|str, or 'All', optional, default as 'All'
        A list of nodes in `nodes` that needs to be considered, other nodes will be ignored
    edges: str, optional, default as 'Euclidean'
        The methods for the calculation of distances between nodes. Options and required additional information are as follows:

        1) (default) 'Euclidean', using Euclidean distance, no additional information needed
        2) 'EuclideanBarrier', using Euclidean distance, if `polys` is provided, the path between nodes will consider them as barriers and by pass those areas.
            - polys: list of poly, the polygons to be considered as barriers
        3) 'LatLon', calculate distances by lat/lon, no additional information needed
            - distUnit: str, the unit of distance, default as 'meter'
        4) 'ManhattenXY', calculate distance by Manhatten distance            
        5) 'Dictionary', directly provide the travel matrix
            - tau: the traveling matrix
        6) 'Grid', traveling on a grid with barriers, usually used in warehouses
            - column: number of columns
            - row: number of rows
            - barrier: a list of coordinates on the grid indicating no-entrance
    **kwargs: optional
        Provide additional inputs for different `edges` options

    Returns
    -------

    tuple
        Two dictionaries, the first one is the travel matrix, index by (nodeID1, nodeID2), the second one is the dictionary for path between start and end locations (useful for 'EuclideanBarrier').

    """

    # Define tau
    tau = {}
    pathLoc = {}

    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    if (edges == 'Euclidean'):
        tau, pathLoc = _matrixDistEuclideanXY(
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            locFieldName = locFieldName)
    elif (edges == 'EuclideanBarrier'):
        if ('polys' not in kwargs or kwargs['polys'] == None):
            warnings.warning("WARNING: No barrier provided.")
            tau, pathLoc = _matrixDistEuclideanXY(
                nodes = nodes, 
                nodeIDs = nodeIDs, 
                locFieldName = locFieldName)
        else:
            tau, pathLoc = _matrixDistBtwPolysXY(
                nodes = nodes, 
                nodeIDs = nodeIDs, 
                polys = kwargs['polys'], 
                locFieldName = locFieldName)
    elif (edges == 'LatLon'):
        distUnit = 'meter' if 'distUnit' not in kwargs else kwargs['distUnit']
        tau, pathLoc = _matrixDistLatLon(
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            distUnit = distUnit,
            locFieldName = locFieldName)
    elif (edges == 'Manhatten'):
        tau, pathLoc = _matrixDistManhattenXY(
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            locFieldName = locFieldName)
    elif (edges == 'Dictionary'):
        if ('tau' not in kwargs or kwargs['tau'] == None):
            raise MissingParameterError("ERROR: 'tau' is not specified")
        for p in kwargs['tau']:
            tau[p] = kwargs['tau'][p]
            pathLoc[p] = kwargs['path'][p] if 'path' in kwargs else [nodes[p[0]][locFieldName], nodes[p[1]][locFieldName]]
    elif (edges == 'Grid'):
        if ('grid' not in kwargs or kwargs['grid'] == None):
            raise MissingParameterError("ERROR: 'grid' is not specified")
        if ('column' not in kwargs['grid'] or 'row' not in kwargs['grid']):
            raise MissingParameterError("'column' and 'row' need to be specified in 'grid'")
        tau, pathLoc = _matrixDistGrid(
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            grids = kwargs['grid'], 
            locFieldName = locFieldName)
    elif (edges == 'RoadNetwork'):
        if ('source' not in kwargs or kwargs['source'] == None):
            raise MissingParameterError("ERROR: 'source' is not specified")
        if ('APIkey' not in kwargs or kwargs['APIkey'] == None):
            raise MissingParameterError("ERROR: 'APIkey' is not specified")

        if (kwargs['source'] == 'Baidu'):
            if (detailFlag):
                raise UnsupportedInputError("ERROR: Stay tune")
            else:
                try:
                    tau = _matrixBaiduSimple(nodes, nodeIDs, API, locFieldName)
                except:
                    raise UnsupportedInputError("ERROR: Failed to fetch data, check network connection and API key")
        else:
            raise UnsupportedInputError("ERROR: Right now we support 'Baidu'")
    else:
        raise UnsupportedInputError(ERROR_MISSING_EDGES)        

    return tau, pathLoc

def _matrixBaiduSimple(nodes: dict, nodeIDs: list, API, locFieldName = 'loc'):
    locStr = ""
    for i in range(len(nodeIDs)):
        locStr += str(nodes[nodeIDs[i]][locFieldName][0]) + "," + str(nodes[nodeIDs[i]][locFieldName][1]) + "|"
    locStr = locStr[:-1]

    url = "https://api.map.baidu.com/routematrix/v2/driving"
    params = {
        "origins": locStr,
        "destinations": locStr,
        "ak": API,
    }
    response = requests.get(url=url, params=params)
    
    tau = {}
    if (response):
        k = 0
        for i in range(len(nodeIDs)):
            for j in range(len(nodeIDs)):
                tau[nodeIDs[i], nodeIDs[j]] = response.json()['result'][k]['distance']['value']
                k += 1

    return tau

def _matrixDistEuclideanXY(nodes: dict, nodeIDs: list, locFieldName = 'loc'):
    tau = {}
    pathLoc = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                d = distEuclideanXY(nodes[i][locFieldName], nodes[j][locFieldName])
                tau[i, j] = d['dist']
                tau[j, i] = d['dist']
                pathLoc[i, j] = [nodes[i][locFieldName], nodes[j][locFieldName]]
                pathLoc[j, i] = [nodes[j][locFieldName], nodes[i][locFieldName]]
            else:
                tau[i, j] = CONST_EPSILON
                tau[j, i] = CONST_EPSILON
                pathLoc[i, j] = []
                pathLoc[j, i] = []
    return tau, pathLoc

def _matrixDistManhattenXY(nodes: dict, nodeIDs: list, locFieldName = 'loc'):
    tau = {}
    pathLoc = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                d = distManhattenXY(nodes[i][locFieldName], nodes[j][locFieldName])
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

def _matrixDistLatLon(nodes: dict, nodeIDs: list, distUnit = 'meter', locFieldName = 'loc'):
    tau = {}
    pathLoc = {}
    for i in nodeIDs:
        for j in nodeIDs:
            if (i != j):
                d = distLatLon(nodes[i][locFieldName], nodes[j][locFieldName], distUnit)
                tau[i, j] = d['dist']
                tau[j, i] = d['dist']
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

def vectorDist(loc: pt, nodes: dict, locFieldName: str = 'loc', nodeIDs: list|str = 'All', edges: str = 'Euclidean', **kwargs) -> dict:
    """
    Given a location and a `nodes` dictionary, returns the traveling distance and path between the location to each node.

    Parameters
    ----------

    loc: pt, required
        Origin/destination location.
    nodes: dict, required
        A `nodes`dictionary with location information. See :ref:`nodes` for reference.
    locFieldName: str, optional, default as 'loc'
        The key in nodes dictionary to indicate the locations
    nodeIDs: list of int|str, or 'All', optional, default as 'All'
        A list of nodes in `nodes` that needs to be considered, other nodes will be ignored
    edges: str, optional, default as 'Euclidean'
        The methods for the calculation of distances between nodes. Options and required additional information are referred to :func:`~vrpSolver.geometry.matrixDist()`.
    **kwargs: optional
        Provide additional inputs for different `edges` options

    Returns
    -------

    tuple
        tau, revTau, pathLoc, revPathLoc. Four dictionaries, the first one is the travel distance from loc to each node index by nodeID, the second the travel distance from each node back to loc. The third and fourth dictionaries are the corresponded path.

    """

    # Define tau
    tau = {}
    revTau = {}
    pathLoc = {}
    revPathLoc = {}

    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    if (edges == 'Euclidean'):
        tau, revTau, pathLoc, revPath = _vectorDistEuclideanXY(
            loc = loc,
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            locFieldName = locFieldName)
    elif (edges == 'EuclideanBarrier'):
        if ('polys' not in kwargs or kwargs['polys'] == None):
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
                polys = kwargs['polys'], 
                locFieldName = locFieldName)
    elif (edges == 'LatLon'):
        tau, revTau, pathLoc, revPath = _vectorDistLatLon(
            loc = loc,
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            locFieldName = locFieldName)
    elif (edges == 'Manhatten'):
        tau, revTau, pathLoc, revPath = _vectorDistManhattenXY(
            loc = loc,
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            locFieldName = locFieldName)
    elif (edges == 'Grid'):
        if ('grid' not in kwargs or kwargs['grid'] == None):
            raise MissingParameterError("'grid' is not specified")
        if ('column' not in kwargs['grid'] or 'row' not in kwargs['grid']):
            raise MissingParameterError("'column' and 'row' need to be specified in 'grid'")
        tau, revTau, pathLoc, revPath = _vectorDistGrid(
            loc = loc,
            nodes = nodes, 
            nodeIDs = nodeIDs, 
            grids = kwargs['grid'], 
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

def scaleDist(loc1: pt, loc2: pt, edges: str = 'Euclidean', **kwargs) -> dict:
    """
    Given a two locations, returns the traveling distance and path between locations.

    Parameters
    ----------

    loc1: pt, required
        The first location.
    loc2: pt, required
        The second location.
    edges: str, optional, default as 'Euclidean'
        The methods for the calculation of distances between nodes. Options and required additional information are referred to :func:`~vrpSolver.geometry.matrixDist()`.
    **kwargs: optional
        Provide additional inputs for different `edges` options

    Returns
    -------

    dict
        tau, revTau, pathLoc, revPathLoc. Four keys, the first one is the travel distance from loc to each node index by nodeID, the second the travel distance from each node back to loc. The third and fourth dictionaries are the corresponded path.

    """

    # Define tau
    dist = None
    revDist = None
    pathLoc = []
    revPathLoc = []

    if (edges == 'Euclidean'):
        dist = distEuclideanXY(loc1, loc2)['dist']
        revDist = dist
        pathLoc = [loc1, loc2]
        revPathLoc = [loc2, loc1]
    elif (edges == 'EuclideanBarrier'):
        if ('polys' not in kwargs or kwargs['polys'] == None):
            warnings.warning("WARNING: No barrier provided.")
            dist = distEuclideanXY(loc1, loc2)['dist']
            revDist = dist
            pathLoc = [loc1, loc2]
            revPathLoc = [loc2, loc1]
        else:
            res = distBtwPolysXY(pt1, pt2, kwargs['polys'])
            dist = res['dist']
            revDist = dist
            pathLoc = [i for i in res['path']]
            revPathLoc = [pathLoc[len(pathLoc) - i - 1] for i in range(len(pathLoc))]
    elif (edges == 'LatLon'):
        dist = distLatLon(loc1, loc2)['dist']
        revDist = dist
        pathLoc = [loc1, loc2]
        revPathLoc = [loc2, loc1]
    elif (edges == 'Manhatten'):
        dist = distManhattenXY(loc1, loc2)['dist']
        revDist = dist
        pathLoc = [loc1, (pt1[0], pt2[1]), loc2]
        revPathLoc = [loc2, (pt1[0], pt2[1]), loc1]
    elif (edges == 'Grid'):
        if ('grid' not in kwargs or kwargs['grid'] == None):
            raise MissingParameterError("'grid' is not specified")
        if ('column' not in kwargs['grid'] or 'row' not in kwargs['grid']):
            raise MissingParameterError("'column' and 'row' need to be specified in 'grid'")
        res = distOnGrid(loc1, loc2, kwargs['grid'])
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
