import math
import shapely

from .geometry import *

def createNeighbor(
    nodes: dict,
    nodeIDs: list[int | str] | str = 'All',
    neighbor: dict = {'shape': 'circle', 'radius': 1, 'lod': 30}
    ) -> dict:

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

def neighborIntersection(
    poly1: poly,
    poly2: poly):
    poly = None
    return poly

def shortestPt2Poly2Pt(
    pt1: pt,
    poly: poly,
    pt2: pt) -> pt:

    return pt

