import math
import shapely
from shapely.geometry import mapping
import coptpy as cp

from .geometry import *

def createNeighbor(
    nodes: dict,
    nodeIDs: list[int | str] | str = 'All',
    neighbor: dict = {'shape': 'circle', 'radius': 1, 'lod': 30}
    ) -> dict:

    """Given a node dictionary, add neighborhood to selected nodes

    Parameters
    ----------
    node: dictionary, required
        A node dictionary indicating the locations in 'loc'
    nodeIDs: string | list[int|str], optional, default as 'All'
        A list of node IDs to add neighborhood, leave it as 'All' to indication adding to all.
    neighbor: dictionary, optional, default {'shape': 'circle', 'radius': 1, 'lod': 30}
        The shape of dictionary. Options includes
        1) Adding polygon surrounding nodes
            >>> neighbor = {
            ...     'shape': 'poly',
            ...     'poly': poly, # In relative axis where node locates in [0, 0]
            ... }
        2) Adding disk circle surrounding nodes
            >>> neighbor = {
            ...     'shape': 'circle',
            ...     'radius': 1,
            ...     'lod': 30 # Optional, 'lod' = 'level of detail', default to use a 30-gon representing cirlce
            }
        3) Add egg shape surrounding nodes
            >>> neighbor = {
            ... 
            ... }

    Returns
    -------

    dictionary
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

    if (neighbor == None or 'shape' not in neighbor):
        raise MissingParameterError("ERROR: Missing required field `neighbor`")
    
    for n in nodeIDs:
        if (neighbor['shape'] == 'poly'):
            if ('poly' not in neighbor):
                raise MissingParameterError("ERROR: Missing required key 'poly' in `neighbor`")
            poly = [[i[0] + nodes[n]['loc'][0], i[1] + nodes[n]['loc'][1]] for i in neighbor['poly']]
            nodes[n]['neighbor'] = poly
        if (neighbor['shape'] == 'circle'):
            if ('radius' not in neighbor):
                raise MissingParameterError("ERROR: Missing required key 'poly' in `neighbor`")
        # By default, a circle is plotted by a 30-gon
        lod = 30
        if ('lod' in neighbor and type(neighbor['lod']) == int):
            lod = neighbor['lod']
        poly = [[
            nodes[n]['loc'][0] + neighbor['radius'] * math.sin(2 * d * math.pi / lod),
            nodes[n]['loc'][1] + neighbor['radius'] * math.cos(2 * d * math.pi / lod),
        ] for d in range(lod)]
        nodes[n]['neighbor'] = poly
    return nodes

def cutNodesNeighbor(
    nodes: dict) -> dict:

    """Given a node dictionary with neighbors, reduce the neiborhood to the convex hull

    WARNING
    -------
    This function will modify the input dictionary `nodes`


    Parameter
    ---------

    nodes: dictionary, required
        The nodes dictionary with neighborhood.

    Returns
    -------

    dictionary
        node dictionary that neighborhood has been reduced

    """

    # First, create the convex hull for all customers
    lstNodeLoc = []
    for n in nodes:
        lstNodeLoc.append(shapely.Point(nodes[n]['loc'][0], nodes[n]['loc'][1]))
    ch = shapely.convex_hull(shapely.MultiPoint(points = lstNodeLoc))

    # Next, for all nodes with neighborhood, update it with intersection
    for n in nodes:
        if ('neighbor' in nodes[n]):
            neiPoly = shapely.Polygon([[p[0], p[1]] for p in nodes[n]['neighbor']])
            intSectPoly = shapely.intersection(ch, neiPoly)
            newNei = [i for i in mapping(intSectPoly)['coordinates'][0]]
            nodes[n]['neighbor'] = newNei
    
    return nodes

def createSteinerZone(
    nodes: dict, 
    order: int|None = None
    ) -> list[dict]:
    
    """Given a node dictionary, returns a list of Steiner zones

    Parameters
    ----------

    nodes: dictionary, required
        The nodes dictionary with neighborhood.
    order: int, optional, 5
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
    for i in nodes:
        for j in nodes:
            if (i < j):
                neiI = None
                if ('neighbor' in nodes[i]):
                    neiI = shapely.Polygon([[p[0], p[1]] for p in nodes[i]['neighbor']])
                else:
                    neiI = shapely.Point([nodes[i]['loc'][0], nodes[i]['loc'][1]])
                neiJ = None
                if ('neighbor' in nodes[j]):
                    neiJ = shapely.Polygon([[p[0], p[1]] for p in nodes[j]['neighbor']])
                else:
                    neiJ = shapely.Point([nodes[j]['loc'][0], nodes[j]['loc'][1]])

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
    if (sum(overlapMatrix.values()) == 0):
        return [
            {
                'poly': nodes[n]['neighbor'] if 'neighbor' in nodes[n] else [nodes[n]['loc']],
                'repPt': list(shapely.Polygon([[p[0], p[1]] for p in nodes[n]['neighbor']]).centroid.coords[0]) if 'neighbor' in nodes[n] else [nodes[n]['loc']],
                'nodeIDs': [n]
            } for n in nodes]

    pointer = 0
    checkOrder = 2
    while (checkOrder <= (order if order != None else len(nodes))):
        # Q: How many SZ needs to be checked? A: From `pointer` to `endPointer`
        endPointer = len(lstSteinerZoneShape)

        # Check each SZ in this order, to see if it can be increase by order 1
        for p in range(pointer, endPointer):
            # Get a SZ, see if there is a node n intersect with this SZ
            SZShape = lstSteinerZoneShape[p]
            for n in nodes:
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
                            neiN = shapely.Polygon([[p[0], p[1]] for p in nodes[n]['neighbor']]) if 'neighbor' in nodes[n] else shapely.Point(nodes[n]['loc'])
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
        'poly': nodes[n]['neighbor'] if 'neighbor' in nodes[n] else [nodes[n]['loc']],
        'repPt': list(shapely.Polygon([[p[0], p[1]] for p in nodes[n]['neighbor']]).centroid.coords[0]) if 'neighbor' in nodes[n] else nodes[n]['loc'],
        'nodeIDs': [n]
    } for n in nodes]

    for n in lstSteinerZoneShape:
        # print(n['polyShape'])
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


