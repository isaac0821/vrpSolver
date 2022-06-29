from .geometry import *
from .vector import *

def getTau(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...} or \
                 4) String 'Grid', will need to add arguments using `edgeArgs`" = "Euclidean",
    edgeArgs:   "Dictionary, to help defining traveling cost matrix \
                 1) for 'Euclidean', 'LatLon', and the travel cost dictionary,\
                    {\
                        'scale': 1, \
                    }\
                 2) for 'Grid'\
                    {\
                        'colRow': (numCol, numRow),\
                        'barriers': [(coordX, coordY), ...], \
                    }" = None,
    depotID:    "DepotID, default to be 0" = 0,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0):
    # Define tau ==============================================================
    tau = {}
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            if (edgeArgs == None):
                tau = _getTauEuclidean(nodes, nodeIDs)
            elif ('scale' not in edgeArgs):
                msgWarning("WARNING: Missing 'scale' in `edgeArgs`. Set to be default value as 1")
                tau = _getTauEuclidean(nodes, nodeIDs)
            else:
                tau = _getTauEuclidean(nodes, nodeIDs, edgeArgs['scale'])
        elif (edges == 'LatLon'):
            if (edgeArgs == None):
                tau = _getTauLatLon(nodes, nodeIDs)
            elif ('scale' not in edgeArgs):
                msgWarning("WARNING: Missing 'scale' in `edgeArgs`. Set to be default value as 1")
                tau = _getTauLatLon(nodes, nodeIDs)
            else:
                tau = _getTauLatLon(nodes, nodeIDs, edgeArgs['scale'])
        elif (edges == 'Grid'):
            if (edgeArgs == None or 'colRow' not in edgeArgs or 'barriers' not in edgeArgs):
                msgError("ERROR: Need more information to define the grid.")
                return None
            tau = _getTauGrid(nodes, nodeIDs, edgeArgs)
        else:
            msgError(ERROR_INCOR_TAU)
            return None
    else:
        if (edgeArgs == None):
            for p in edges:
                tau[p] = edges[p]
        elif ('scale' not in edgeArgs):
            msgWarning("WARNING: Missing 'scale' in `edgeArgs`. Set to be default value as 1")
            for p in edges:
                tau[p] = edges[p]
        else:
            for p in edges:
                tau[p] = edges[p] * edgeArgs['scale']

    # Service time ============================================================
    if (serviceTime != None and serviceTime > 0):
        for (i, j) in tau:
            if (i != depotID and j != depotID and i != j):
                tau[i, j] += serviceTime
            elif (i == depotID or j == depotID and i != j):
                tau[i, j] += serviceTime / 2 

    return tau

def _getTauEuclidean(nodes, nodeIDs, speed = 1):
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

def _getTauLatLon(nodes, nodeIDs, distUnit = 'meter', speed = 1):
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

def _getTauGrid(nodes, nodeIDs, grid):
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
                t = gridPathFinding(
                    gridColRow = grid['colRow'],
                    barriers = grid['barriers'],
                    startGridCoord = nodes[i]['loc'],
                    endGridCoord = nodes[j]['loc'])['dist']
                tau[i, j] = t
                tau[j, i] = t
            else:
                tau[i, j] = 0
                tau[j, i] = 0
    return tau

def _getTauRoadNetwork(
    nodes:      "Dictionary, returns the coordinate of given nodeIDs, \
                    {\
                        nodeIDs1: {'loc': (x, y)}, \
                        nodeIDs2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All', 
    network:    "List of dictionary of networks in the format of \
                    {\
                        'start': (lat, lon),\
                        'end': (lat, lon),\
                    }" = None,
    ) -> "Given a network, e.g., road network, and a set of locations on network (could be on arc), returns distance matrix":

    return tau

def getNeighborCluster(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon'" = "Euclidean",
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    diameter:   "The diameter of the neighborhood" = None,
    maxSize:    "Maximum number of nodes in each cluster" = 5
    ) -> "Given a dictionary of locations, an a radius, return the sets that any two locations are within the distance":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define edges ============================================================
    if (edges == 'Euclidean'):
        edges = getTauEuclidean(nodes)
    elif (edges == 'LatLon'):
        edges = getTauLatLon(nodes)
    else:
        print("Error: Incorrect type `edges`")
        return None

    # Initialize ==============================================================
    seedClique = []
    # For each node get the neighbors within diameter
    neighbor = {}
    for n1 in nodeIDs:
        for n2 in nodeIDs:
            if (n1 != n2 and edges[n1, n2] <= diameter):
                if (n1 not in neighbor):
                    neighbor[n1] = []
                neighbor[n1].append(n2)

    # Find seed cliques in the neighbor graph =================================
    # FIXME: Now using stupid method, will be replaced by clique searching algorithm, or should I?
    # 3-clique
    for n1 in neighbor:
        for n2 in neighbor[n1]:
            for n3 in neighbor[n2]:
                if (n1 < n2 < n3 and n1 in neighbor[n3]):
                    seedClique.append([n1, n2, n3])

    # Try to union existing seed-clique to find larger clique =================
    canUnionFlag = True
    while (canUnionFlag):
        canUnionFlag = False
        for c1 in range(len(seedClique) - 1):
            for c2 in range(c1 + 1, len(seedClique)):
                # Try to merge two cliques
                clique1 = [i for i in seedClique[c1]]
                clique2 = [i for i in seedClique[c2]]
                intersect = listSetIntersect(clique1, clique2)                
                # Two cliques can be merged iff they have intersection
                if (len(intersect) > 0):
                    diff1 = listSetMinus(clique1, intersect)
                    diff2 = listSetMinus(clique2, intersect)                    
                    # Try to see if the nodes that are not in the intersection are close
                    mergeFlag = True
                    for n1 in diff1:
                        for n2 in diff2:
                            if (n1 not in neighbor[n2]):
                                mergeFlag = False
                                break
                        if (not mergeFlag):
                            break
                    if (mergeFlag):
                        newClique = listSetUnion(clique1, clique2)
                        if (maxSize == None or len(newClique) <= maxSize):
                            newClique.sort()
                            seedClique.remove(clique1)
                            seedClique.remove(clique2)
                            seedClique.append(newClique)
                            canUnionFlag = True
                            break
            if (canUnionFlag):
                break
    return seedClique

def getSortedNodesByDist(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon'" = "Euclidean",
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    refNodeID:  "List, [x, y], the reference location to calculate distance" = None    
    ) -> "Given a set of locations, and a referencing node, sort the nodes by distance to this referencing node":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'LatLon'):
            edges = getTauLatLon(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    # Sort distance ===========================================================
    sortedSeq = []
    sortedSeqHeap = []
    for n in nodeIDs:
        dist = edges[refNodeID, n]
        heapq.heappush(sortedSeqHeap, (dist, n))
    while (len(sortedSeqHeap) > 0):
        sortedSeq.append(heapq.heappop(sortedSeqHeap)[1])  

    return sortedSeq

def getSweepSeq(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    centerLoc:  "1) (Default) String, 'Centroid' the centroid of nodes, or\
                 2) List, [x, y], the center point" = 'Centroid',
    isClockwise: "True if the sweeping direction is clock-wise, False otherwise" = True,
    initDeg:    "Starting direction of the sweeping, 0 as North" = 0
    ) -> "Given a set of locations, and a center point, gets the sequence from sweeping":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Initialize centroid =====================================================
    if (centerLoc == 'Centroid'):
        centerLoc = getCentroid(nodes = nodes, nodeIDs = nodeIDs)

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

def getCentroid(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    algo:       "1) String, 'Weiszfeld'" = 'Weiszfeld'
    ) -> "Get centroid location for given nodes":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)

    # Call subroutines ========================================================
    centroid = None
    if (algo == "Weiszfeld"):
        centroid = _getCentroidWeiszfeld(nodes, nodeIDs)
    else:
        return None

    return centroid

def _getCentroidWeiszfeld(nodes, nodeIDs):
    # Initialize ==============================================================
    q = [1 for i in range(len(nodeIDs))]
    a = [nodes[nodeIDs[i]]['loc'][0] for i in range(len(nodeIDs))]
    b = [nodes[nodeIDs[i]]['loc'][1] for i in range(len(nodeIDs))]
    x = sum(a) / len(a)
    y = sum(b) / len(b)
    f = 0
    for i in range(len(nodeIDs)):
        f += math.sqrt((x - a[i])**2 + (y - b[i])**2)

    # Iterations ==============================================================
    canGoFlag = True
    while (canGoFlag):
        # update q
        q = []
        for i in range(len(nodeIDs)):
            q.append(1 / math.sqrt((x - a[i])**2 + (y - b[i])**2))

        # update x, y
        x = 0
        y = 0
        for i in range(len(nodeIDs)):
            x += q[i] * a[i]
            y += q[i] * b[i]
        x /= sum(q)
        y /= sum(q)

        # update f
        newF = 0
        for i in range(len(nodeIDs)):
            newF += math.sqrt((x - a[i])**2 + (y - b[i])**2)
        if (abs(newF - f) < CONST_EPSILON):
            canGoFlag = False
        f = newF

    # Output ==================================================================
    centroid = (x, y)

    return centroid

def getScan(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    direction:  "Direction of scanning, 0 as North" = 0
    ) -> "Scan nodes from one direction":

    # Initialize ==============================================================
    baseline = []
    maxDist = 0
    centroid = getCentroid(nodes)
    for n in nodes:
        d = distEuclidean2D(nodes[n]['loc'], centroid)
        if (maxDist == None or d > maxDist):
            maxDist = 1.2 * d
    basePt = ptInDistXY(centroid, direction, maxDist)
    baseline = getPerpendicularLine(basePt, vecPolar2XY([10, direction]))

    # Distance to the baseline ================================================
    distHeap = []
    for n in nodes:
        dist2Baseline = distPt2Line(nodes[n]['loc'], baseline)
        heapq.heappush(distHeap, (dist2Baseline, n))

    # Scan ====================================================================
    scanSeq = []
    while (len(distHeap)):
        scanSeq.append(heapq.heappop(distHeap)[1])

    return scanSeq

def getConvexHull(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    algo:       "1) String, 'Jarvis', O(nH) (current implementation is O(nHlog n)) or\
                 2) String, (not available) 'DNC', O(nlog n) or\
                 3) String, (not available) 'Graham', O(nlog n) or\
                 4) String, (not available) 'Melkman'" = "Jarvis"
    ) -> "Given a set of node locations, return a list of nodeID which construct the convex hull":

    # Initialize ==============================================================
    chSeq = None

    # Some extreme cases ======================================================
    if (len(nodes) == 0):
        return None
    elif (len(nodes) <= 3):
        chSeq = []
        for n in nodes:
            chSeq.append(n)
        return chSeq

    # Call subroutines ========================================================
    if (algo == "Jarvis"):
        chSeq = _getConvexHullJavis(nodes)
    else:
        return None
    
    return chSeq

def _getConvexHullJavis(nodes):
    # References ==============================================================
    # 1. https://blog.csdn.net/Bone_ACE/article/details/46239187
    # 2. chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/viewer.html?pdfurl=http%3A%2F%2Fwww.ams.sunysb.edu%2F~jsbm%2Fcourses%2F345%2F13%2Fmelkman.pdf&clen=46562&chunk=true
    # 3. https://en.wikipedia.org/wiki/Convex_hull_algorithms

    # Initialize ==============================================================
    chSeq = []

    # Get the location of the left-most nodeID ================================
    # Find an initial point which guaranteed to be in convex hull
    leftMostID = None
    leftMostX = None
    for n in nodes:
        if (leftMostID == None or nodes[n]['loc'][0] < leftMostX):
            leftMostID = n
            leftMostX = nodes[n]['loc'][0]

    # Jarvis march ============================================================
    curNodeID = leftMostID
    curDirection = 0
    marchFlag = True
    while (marchFlag):
        sweepSeq = getSweepSeq(
            nodes = nodes,
            nodeIDs = [i for i in nodes if i not in chSeq],
            centerLoc = nodes[curNodeID]['loc'],
            initDeg = curDirection)
        if (sweepSeq[0] == leftMostID):
            marchFlag = False
        chSeq.append(sweepSeq[0])
        curDirection = headingXY(nodes[curNodeID]['loc'], nodes[sweepSeq[0]]['loc'])    
        curNodeID = sweepSeq[0]
    return chSeq
