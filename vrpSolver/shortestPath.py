import heapq

from .common import *
from .msg import *

def graphShortestPath(
    weightArcs: "A list of 3-tuples, (ID1, ID2, weight), indexes of vertices must start from 0" = None,
    oID:        "Original node id" = None,
    dID:        "Destination node id" = None,
    algo:       "1) String, (default) 'Dijkstra' for single-source non-negative edge weight or, \
                 2) String, 'BellmanFord' for single-source or, \
                 2) String, 'FloydWarshall' for all pairs or, \
                 3) String, 'AStar' or, \
                 4) String, 'Johnson'" = 'Dijkstra'    
    ) -> "Return a shortest path":

    arcType = ""


    return res

def _graphShortestPathDijkstra(weightArcs, oID, dID):
    sp = {}

    return sp