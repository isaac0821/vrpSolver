import heapq
import math
import re
import shapely
from shapely.geometry import mapping

from .const import *
from .common import *
from .neighbor import *
from .graph import *
from .geometry import *
from .msg import *

# History =====================================================================
# 20230519 - Initialize
# =============================================================================

def heuCETSP(
    nodes: dict,
    algo: dict,
    tolerance: float = 0.01,
    ) -> dict | None:

    """Use heuristic method to find suboptimal CETSP solution

    Parameters
    ----------

    nodes: dictionary, required, default None
        The coordinates of given nodes, in the following format::
            >>> nodes = {
            ...     'nodeID1': {'loc': (x, y), 'neighbor': poly}, 
            ...     'nodeID2': {'loc': (x, y), 'neighbor': poly}, # ...
            ... }

    Returns
    -------

    dictionary

    """

    # Create convex hull of the nodes, truncate neighborhoods
    nodes = cutNodesNeighbor(nodes)

    # Create a list of Steiner zones
    sz = createSteinerZone(nodes)

    # Constructive phase
    # Step 1: Start from SZ that has the highest order
    remain = [i for i in range(12)]
    szIDs = []
    for i in range(len(sz)):
        szID = len(sz) - i - 1

        allInRemain = True
        for j in sz[szID]['nodeIDs']:
            if (j not in remain):
                allInRemain = False
        if (allInRemain):
            szIDs.append(szID)
            for k in sz[szID]['nodeIDs']:
                remain.remove(k)

    # Step 2: Get centroids of selected SZs as representative point
    repNodes = []
    for szID in szIDs:
        repNodes[szID] = {
            'loc': sz[szID]['centroid']
        }

    # Step 3: Transform into shortest path problem
    repSeq = []