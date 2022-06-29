import heapq
import gurobipy as grb
import math

from .common import *

def gridPathFinding(
    grid:       "Dictionary, includes the grid info\
                {\
                    'gridColRow': gridColRow,\
                    'barriers': barriers,\
                    'type': 'Square' # or 'Hexagon'\
                }" = None,
    startCoord: "Start coordinate" = (0, 0),
    endCoord:   "End coordinate" = (0, 0),
    algo:       "1) String, (default) 'A*', or\
                 2) String, (not available) 'Dijkstra', or\
                 3) String, (not available) 'BFS', or\
                 4) String, (not available) 'Greed-BFS', or\
                 5) String, (not available) 'B*'" = 'A*',
    algoArgs:   "Dictionary, args for path finding algorithm \
                 1) For 'A*', \
                    {\
                        'distMeasure': (default) 'Manhatten', other options are 'Euclidean';\
                    }" = None
    ) -> "Given two coordinates on the grid, finds the 'shortest' path to travel":

    # Decode ==================================================================
    gridColRow = grid['gridColRow']
    barriers = grid['barriers']

    # Call path finding =======================================================
    if (algo == 'A*'):
        if (algoArgs == None or 'distMeasure' not in algoArgs):
            msgWarning("Warning: Missing `algoArgs` or missing 'distMeasure' in `algoArgs`. Set 'distMeasure' to be 'Manhatten'")
            algoArgs = {'distMeasure': 'Manhatten'}
        res = _gridPathFindingAStar(gridColRow, barriers, startCoord, endCoord, algoArgs['distMeasure'])
    else:
        print("Error: Incorrect or not available grid path finding option!")
    return res

def _gridPathFindingAStar(gridColRow, barriers, startCoord, endCoord, distMeasure):
    # Heuristic measure ==================================================-
    def _calManhattenDist(coord1, coord2):
        return abs(coord1[0] - coord2[0]) + abs(coord1[1] - coord2[1])
    def _calEuclideanDist(coord1, coord2):
        return math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2)

    # Initialize grid ====================================================-
    # Evaluate value f(n) = g(n) + h(n)
    gridStatus = {}
    for col in range(gridColRow[0]):
        for row in range(gridColRow[1]):
            if ((col, row) not in barriers):
                # Content in the dictionary (g(n), h(n), fromCoord)
                # At this stage, no need to calculate h(n) 
                gridStatus[(col, row)] = (None, None, None)
            else:
                gridStatus[(col, row)] = 'block'
    if (distMeasure == 'Manhatten'):
        gridStatus[startCoord] = (0, _calManhattenDist(startCoord, endCoord), None)
    elif (distMeasure == 'Euclidean'):
        gridStatus[startCoord] = (0, _calEuclideanDist(startCoord, endCoord), None)
    gridStatus[endCoord] = (None, 0, None)

    # Open/close set ======================================================
    openList = [startCoord]
    closeList = [i for i in barriers]

    # Find smallest Fn ====================================================
    def _findSmallestFnGrid():
        bestFn = None
        bestCoord = None
        for coord in openList:
            if (gridStatus[coord] != None
                and gridStatus[coord] != 'block' 
                and (bestFn == None or gridStatus[coord][0] + gridStatus[coord][1] < bestFn)):
                bestFn = gridStatus[coord][0] + gridStatus[coord][1]
                bestCoord = coord
        return bestCoord

    # For each grid in open set, update g(n) ==============================
    while (len(openList) > 0):
        tmpOpenList = []
        coord = _findSmallestFnGrid()
        # Up
        upCoord = (coord[0], coord[1] + 1)
        if (coord[1] + 1 < gridColRow[1] and gridStatus[upCoord] != None and gridStatus[upCoord] != 'block' and upCoord not in closeList):
            if (gridStatus[upCoord][0] == None or gridStatus[upCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[upCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(upCoord, endCoord), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[upCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(upCoord, endCoord), coord)
                if (upCoord == endCoord):
                    break
                else:
                    tmpOpenList.append(upCoord)
        # Down
        downCoord = (coord[0], coord[1] - 1)
        if (coord[1] - 1 >= 0 and gridStatus[downCoord] != None and gridStatus[downCoord] != 'block' and downCoord not in closeList):
            if (gridStatus[downCoord][0] == None or gridStatus[downCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[downCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(downCoord, endCoord), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[downCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(downCoord, endCoord), coord)
                if (downCoord == endCoord):
                    break
                else:
                    tmpOpenList.append(downCoord)
        # Left
        leftCoord = (coord[0] - 1, coord[1])
        if (coord[0] - 1 >= 0 and gridStatus[leftCoord] != None and gridStatus[leftCoord] != 'block' and leftCoord not in closeList):
            if (gridStatus[leftCoord][0] == None or gridStatus[leftCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[leftCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(leftCoord, endCoord), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[leftCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(leftCoord, endCoord), coord)
                if (leftCoord == endCoord):
                    break
                else:
                    tmpOpenList.append(leftCoord)
        # Right
        rightCoord = (coord[0] + 1, coord[1])
        if (coord[0] + 1 < gridColRow[0] and gridStatus[rightCoord] != None and gridStatus[rightCoord] != 'block' and rightCoord not in closeList):
            if (gridStatus[rightCoord][0] == None or gridStatus[rightCoord][0] > gridStatus[coord][0] + 1):
                if (distMeasure == 'Manhatten'):
                    gridStatus[rightCoord] = (gridStatus[coord][0] + 1, _calManhattenDist(rightCoord, endCoord), coord)
                if (distMeasure == 'Euclidean'):
                    gridStatus[rightCoord] = (gridStatus[coord][0] + 1, _calEuclideanDist(rightCoord, endCoord), coord)
                if (rightCoord == endCoord):
                    break
                else:
                    tmpOpenList.append(rightCoord)
        openList.remove(coord)
        openList.extend(tmpOpenList)
        closeList.append(coord)

    # Recover path ========================================================
    path = []
    curCoord = endCoord
    finishReconstructFlag = True
    while (finishReconstructFlag):
        finishReconstructFlag = False
        path.insert(0, curCoord)
        curCoord = gridStatus[curCoord][2]
        if (curCoord != None):
            finishReconstructFlag = True
    return {
        'dist': len(path) - 1,
        'path': path
    }

def treeTraversal(
    tree:       "Dictionary, returns the children of given nodeID, None if no child, \
                {\
                    nodeID1: {child1: weight, child2: weight, ...}, \
                    nodeWithNoChild: None, # We can ignore these items in the dictionary\
                    ... \
                }" = None, 
    oID:        "String/Integer, nodeID of the root" = None,
    algo:       "1) String, (default) 'DepthFirst' or, \
                 2) String, 'BreadthFirst'" = 'DepthFirst'
    ) -> "Return a sequence of node ids that traverses the tree":

    # Solve by different algorithms ===========================================
    res = None
    if (algo == 'DepthFirst'):
        res = _treeTraversalDepthFirst(tree, oID)
    elif (algo == 'BreadthFirst'):
        res = _treeTraversalBreadthFirst(tree, oID)
    return res

def graphTraversal(
    arcs:       "1) A list of 3-tuple (nodeID1, nodeID2, weight) or, \
                 2) A list of 2-tuple (nodeID1, nodeID2)",
    oID:        "1) String/Integer, nodeID of the root or, \
                 2) None, (default) the first nodeID in `tree`" = None,
    algo:       "1) String, (default) 'DepthFirst' or, \
                 2) String, 'BreadthFirst'" = 'DepthFirst'
    ) -> "Return a sequence of node ids that traverses the tree":
    # Convert arcs into adjList ===============================================
    weightArcs = []
    for arc in arcs:
        if (len(arc) == 3):
            weightArcs.append(arc)
        elif (len(arc) == 2):
            weightArcs.append((arc[0], arc[1], 1))
        else:
            msgError("ERROR: Incorrect `arcs` format")
            return

    mst = graphMST(weightArcs = weightArcs, algo = 'Krusal', exportAs = 'Tree')['mst']

    # Set Default oID =========================================================
    if (oID == None):
        msgWarning("Warning: Missing `oID` as the root of traversing, default set as the smallest ID")
        oID = min(mst)

    # Solve by different algorithms ===========================================
    res = treeTraversal(mst, oID, algo)
    return res

def _treeTraversalDepthFirst(tree, oID):
    visited = []

    # Visit children recursively ==============================================
    def _visitNode(nodeID):
        if (nodeID not in visited):
            visited.append(nodeID)
            if (nodeID in tree and tree[nodeID] != None):
                for child in tree[nodeID]:
                    _visitNode(child)

    # Start search from root ==================================================
    _visitNode(oID)

    return {
        'seq': visited
    }

def _treeTraversalBreadthFirst(tree, oID):
    visited = [oID]
    pointer = 0

    # Main iterations =========================================================
    while (pointer < len(visited)):
        # Scan though visited, add the children to the end of visited
        if (visited[pointer] in tree and tree[visited[pointer]] != None):             
            for v in tree[visited[pointer]]:
                if (v not in visited):
                    visited.append(v)
        pointer += 1

    return {
        'seq': visited
    }
    
def graphComponents(
    arcs:       "1) A list of 2-tuples, (ID1, ID2), or \
                 2) A list of 3-tuples, (ID1, ID2, weight)" = None
    ) -> "A list of components":
    # Create adj list, each vertex start with an empty list ===================
    adjList = graphArcs2AdjList(arcs)

    # Initialize ==============================================================
    found = {}
    for node in adjList:
        found[node] = 0
    components = []

    # Main algorithm, mark neighbors ==========================================
    for i in adjList:
        comp = []
        q = []
        if (found[i] == 0):
            found[i] = 1
            comp.append(i)
            q.append(i)
            while (q):
                v = q.pop(0)
                if (v in adjList and adjList[v] != None):
                    for u in adjList[v]:
                        if (found[u] == 0):
                            found[u] = 1
                            comp.append(u)
                            q.append(u)
            components.append(comp)
    return components

# [Constructing]
def graphShortestPath(
    weightArcs: "A list of 3-tuples, (ID1, ID2, weight)" = None,
    algo:       "1) String, (not available) 'Dijkstra' with list, or, \
                 2) String, (not available) 'Johnson', Dijkstra algorithm with binary heap, or, \
                 3) String, (not available) 'FredmanTarjan', Dijkstra algorithm with Fibonacci heap, or, \
                 4) String, (not available) 'BellmanFold', or, \
                 5) String, (not available) 'FlyodWarshall'" = 'A*',
    oID:        "ID of origin vertex" = None,
    sID:        "1) Integer, ID of target vertex, or, \
                 2) List of integers, IDs of target vertices, or\
                 3) None, calculates distance to all other vertices" = None
    ) -> "Given a graph and the source vertex, returns a dictionary that includes the value and shapepoints of \
          1) shortest path to target vertex if target vertex is given, \
          2) shortest path to all other vertices if target vertex is not given":

    # Number of vertices ======================================================
    vertices = []
    for i in range(len(weightArcs)):
        if (weightArcs[i][0] not in vertices):
            vertices.append(weightArcs[i][0])
        if (weightArcs[i][1] not in vertices):
            vertices.append(weightArcs[i][1])

    # Check if the path weights are all nonnegative ===========================
    allNonnegativeFlag = True
    for arc in weightArcs:
        if (arc[2] < 0):
            allNonnegativeFlag = False
            break
    if (not allNonnegativeFlag):
        if (algo in ['Dijkstra', 'Johnson', 'FredmanTarjan']):
            print(ERROR_OPTS_SHORTESTPATH_ALGO)

    return res

def graphMST(
    weightArcs: "A list of 3-tuples, (ID1, ID2, weight), indexes of vertices must start from 0" = None,
    algo:       "1) String, (default) 'Krusal' or, \
                 2) String, (not available) 'Prim_AdjList' or, \
                 3) String, (not available) 'Prim_AdjMat' or, \
                 4) String, (not available) 'Boruvka' or, \
                 5) String, (not available) 'ReverseDelete'" = 'Krusal',
    exportAs:   "1) String, 'Arcs' or,\
                 2) String, 'Tree'" = 'Arcs'
    ) -> "A list of weightArcs/tree which forms a minimal spanning tree":

    # Number of vertices ======================================================
    vertices = []
    for i in range(len(weightArcs)):
        if (weightArcs[i][0] not in vertices):
            vertices.append(weightArcs[i][0])
        if (weightArcs[i][1] not in vertices):
            vertices.append(weightArcs[i][1])
    numVertices = len(vertices)

    # Call MST ================================================================
    if (algo == 'Krusal'):
        mst = _graphMSTKrusal(weightArcs, numVertices)
    else:
        print("Error: Incorrect or not available MST option!")

    # Decide export format ====================================================
    if (exportAs == 'Arcs'):
        return mst
    elif (exportAs == 'Tree'):
        mstAsTree = graphArcs2AdjList(mst['mst'])
        return {
            'mst': mstAsTree,
            'value': mst['value']
        }
    else:
        return

def _graphMSTKrusal(weightArcs, numVertices):
    # Initialize ==============================================================
    mst = []
    val = 0
    compList = []

    # Arc sorting =============================================================
    sortedWeightArcs = []
    for i in range(len(weightArcs)):
        heapq.heappush(sortedWeightArcs, (weightArcs[i][2], weightArcs[i]))
    
    # Krusal algorithm, add weightArcs between components =====================
    # NOTICE: Will return incorrect result if there are multiple components in the graph
    while(len(mst) < numVertices - 1 and len(sortedWeightArcs) > 0):
        # Uninserted arc with minimal weight
        currArc = heapq.heappop(sortedWeightArcs)

        # Mark two nodes
        nodeID1 = currArc[1][0]
        nodeID2 = currArc[1][1]
        weight = currArc[1][2]
        compID1 = None
        compID2 = None
        findNodeFlag1 = False
        findNodeFlag2 = False

        # Find component that nodes belong to
        for i in range(len(compList)):
            if (nodeID1 in compList[i]):
                findNodeFlag1 = True
                compID1 = i
            if (nodeID2 in compList[i]):
                findNodeFlag2 = True
                compID2 = i
            if (findNodeFlag1 and findNodeFlag2):
                break

        # If two nodes are not in the same component, merge components
        if ((not findNodeFlag1) and (not findNodeFlag2)):
            mst.append(currArc[1])
            val += weight
            compList.append([nodeID1, nodeID2])
        elif (findNodeFlag1 and (not findNodeFlag2)):
            mst.append(currArc[1])
            val += weight
            compList[compID1].append(nodeID2)
        elif ((not findNodeFlag1) and findNodeFlag2):
            mst.append(currArc[1])
            val += weight
            compList[compID2].append(nodeID1)
        elif (findNodeFlag1 and findNodeFlag2):
            if (compID1 != compID2):
                mst.append(currArc[1])
                val += weight
                compList[compID1].extend(compList[compID2].copy())
                compList.remove(compList[compID2])

    return {
        'mst': mst,
        'value': val
    }

def graphCheckBipartite(
    arcs:       "1) A list of 3-tuple (nodeID1, nodeID2, weight) or, \
                 2) A list of 2-tuple (nodeID1, nodeID2)" = None
    ) -> "Given a graph, check if the graph is bipartite":
    
    # FIXME: This function can be simplified, but not urgent
    adjList = graphArcs2AdjList(arcs = arcs)

    # Initialize setA and setB
    minID = min(adjList)
    setA = [minID]
    setB = [i for i in adjList[minID]]

    # Use BFS to color the graph
    tvs = graphTraversal(arcs = arcs, oID = minID, algo = 'BreadthFirst')['seq']
    for v in tvs:
        childrenV = [i for i in adjList[v]]
        if (v in setA):
            for c in childrenV:
                if (c in setA):
                    return {
                        'bipartiteFlag': False,
                        'setA': None,
                        'setB': None
                    }
                elif (c not in setA and c not in setB):
                    setB.append(c)
        elif (v in setB):
            for c in childrenV:
                if (c in setB):
                    return {
                        'bipartiteFlag': False,
                        'setA': None,
                        'setB': None
                    }
                elif (c not in setA and c not in setB):
                    setA.append(c)
    return {
        'bipartiteFlag': True,
        'setA': setA,
        'setB': setB
    }

# [Constructing]
def graphMaximumFlow():
    return

# [Constructing]
def graphMatching(
    weightArcs: "A list of 3-tuples, (ID1, ID2, weight), indexes of vertices must start from 0" = None,
    bipartiteFlag: "Bypass bipartite graph checking, leave it as None if we do not know, preferably set it to be True/False if we know" = None,
    mType:      "1) String, 'Minimum' if to find minimum matching, or\
                 2) String, 'Maximum' if to find maximum matching" = 'Minimum',
    algo:       "1) String, (not available) 'Blossom' or, \
                 2) String, (default if bipartite) 'Hungarian' or 'Kuhn_Munkres', O(|V|^3), for bipartite graph or, \
                 3) String, (default if not bipartite) 'IP', NPC" = None,
    algoArgs:   "1) If `algo` == 'Hungarian' or 'Kuhn_Munkres'\
                {\
                    'setA': setA, \
                    'setB': setB\
                }" = None
    ) -> "Return a set of vertices that forms a Minimum/Maximum Matching": 

    # Check bipartite =========================================================
    checkBipartiteFlag = False
    # Case 1: algo is not designated, bipartiteFlag is None, check bipartite to determine default
    if (algo == None and bipartiteFlag == None):
        checkBipartiteFlag = True
    # Case 2: algo is not designated, bipartiteFlag is True, but missing setA and/or setB info
    elif (algo == None and bipartiteFlag == True 
        and (algoArgs == None or 'setA' not in algoArgs or 'setB' not in algoArgs)):
        checkBipartiteFlag = True
    # Case 2: algo is designated as bipartite method, bipartiteFlag is None
    elif (algo in ['Hungarian', 'Kuhn_Munkres'] and bipartiteFlag == None):
        checkBipartiteFlag = True
    # Case 3: algo is designated as bipartite method, bipartiteFlag is True, but missing setA and/or setB info
    elif (algo in ['Hungarian', 'Kuhn_Munkres'] and bipartiteFlag == True 
        and (algoArgs == None or 'setA' not in algoArgs or 'setB' not in algoArgs)):
        checkBipartiteFlag = True

    if (checkBipartiteFlag):
        bpt = graphCheckBipartite(arcs = weightArcs)
        bipartiteFlag = bpt['bipartiteFlag']
        algoArgs['setA'] = bpt['setA']
        algoArgs['setB'] = bpt['setB']

    if (algo == None):
        if (bipartiteFlag):
            algo = 'Hungarian'
        else:
            algo = 'IP'

    # Convert arcs based on mType =============================================
    convertedArcs = []
    if (mType == 'Maximum'):
        convertedArcs = [i for i in weightArcs]
    elif (mType == 'Minimum'):
        maxWeight = max([arc[2] for arc in weightArcs]) + 1
        convertedArcs = [(i[0], i[1], maxWeight - i[2]) for i in weightArcs]
    else:
        msgError("ERROR: Please select `mType` from ['Maximum', 'Minimum']")

    # Calculate matching using different algorithms ===========================
    res = None
    if (algo == 'IP'):
        res = _graphMaxMatchingIP(convertedArcs)
    elif (algo in ['Hungarian', 'Kuhn_Munkres']):
        if (bipartiteFlag != True):
            msgError("ERROR: 'Hungarian' option only applies for bipartite graph")
            return
        res = _graphMaxMatchingHungarian(convertedArcs, algoArgs['setA'], algoArgs['setB'])

    return res

def _graphMaxMatchingIP(weightArcs):
    matching = []
    M = grb.Model('Matching')

    # Decision variables ==================================================
    x = {}
    for e in range(len(weightArcs)):
        x[e] = M.addVar(vtype = grb.GRB.BINARY, obj = weightArcs[e][2])

    # Matching objective function ========================================-
    M.modelSense = grb.GRB.MAXIMUM
    M.update()

    # Perfect matching ====================================================
    # First find neighborhoods
    neighborhoods = graphArcs2AdjList(weightArcs)
    for node in neighborhoods:
        neis = neighborhoods[node]
        neiArcs = []
        for nei in neis:
            for i in range(len(weightArcs)):
                if ((node == weightArcs[i][0] and nei == weightArcs[i][1]) 
                    or (node == weightArcs[i][1] and nei == weightArcs[i][0])):
                    neiArcs.append(i)
        M.addConstr(grb.quicksum(x[e] for e in neiArcs) == 1)

    # Matching ============================================================
    M.optimize()

    # Construct solution ==================================================
    ofv = None
    if (M.status == grb.GRB.status.OPTIMAL):
        ofv = M.getObjective().getValue()
        for e in x:
            if (x[e].x > 0.8):
                matching.append(weightArcs[e])

    return {
        'ofv': ofv, 
        'matching': matching
    }

def _graphMaxMatchingHungarian(weightArcs, setA, setB):
    # Ref: https://cse.hkust.edu.hk/~golin/COMP572/Notes/Matching.pdf

    # Initialize ==============================================================
    adjList = graphArcs2AdjList(arcs = weightArcs)
    # Initialize labeling, E_l and M
    # NOTE: Assume all arcs in E_l is from B -> A, i.e., from X -> Y
    labelA = {}
    labelB = {}
    curE = []
    for v in setA:
        labelA[v] = 0
    for v in setB:
        labelB[v] = max(adjList[v].values())
        curE.append((v, max(adjList[v]), max(adjList[v].values())))
    # Find initial match
    curE = sorted(curE, key = lambda x: x[2], reverse = True)
    coveredY = []
    coveredX = []
    curM = []
    for arc in curE:
        if (arc[0] not in coveredX and arc[1] not in coveredY):
            curM.append(arc)
            coveredX.append(arc[0])
            coveredY.append(arc[1])
    # Find initial set S and set T
    S = [v for v in setB if v not in coveredX]  # S \subseteq X(or setB)
    T = []                                      # T <= N_l(S) \neq Y(or setA)

    # Subroutine ==============================================================
    def getNlS(Sl, curEl):
        Tl = []
        for arc in curEl:
            if (arc[0] in Sl and arc[1] not in Tl):
                Tl.append(arc[1])
        return Tl

    # Improve matching ========================================================
    canImproveFlag = True
    while (canImproveFlag):
        canImproveFlag = False

        # Get N_l(S)
        NlS = getNls(S, curE)

        # If N_l(S) = T update labels
        if (Nls == T):
            # Find alpha_l = \min_{x \in S, y \notin T} \{l(x) + l(y) - w(x, y)\}
            alphal = None
            for x in S:
                for y in setB:
                    if (y not in T):
                        if (alphal == None or labelB[x] + label[y] - adjList[x][y] < alphal):
                            alphal = labelB[x] + label[y] - adjList[x][y]

        else:
            pass

    return {
        'ofv': ofv,
        'matching': matching 
    }

def graphArcs2AdjList(
    arcs:       "1) A list of 3-tuple (nodeID1, nodeID2, weight) or, \
                 2) A list of 2-tuple (nodeID1, nodeID2)"
    ) -> "Dictionary of neighbors of each node":

    adjList = {}
    for i in range(len(arcs)):
        if (arcs[i][0] not in adjList):
            adjList[arcs[i][0]] = {}
        if (arcs[i][1] not in adjList):
            adjList[arcs[i][1]] = {} 
        adjList[arcs[i][0]][arcs[i][1]] = arcs[i][2] if (len(arcs[i]) == 3) else 1
        adjList[arcs[i][1]][arcs[i][0]] = arcs[i][2] if (len(arcs[i]) == 3) else 1

    return adjList

def graphAdjList2Arcs(
    adjList:    "A Dictionary that has the neighbors info" = None,
    directedFlag: "True if the graph is directed, False otherwise" = False
    ) -> "adjList convert to arcs":

    arcs = []
    for n in adjList:
        if (adjList[n] != None):
            for v in adjList[n]:
                if (directedFlag or n < v):
                    arcs.append((n, v, adjList[n][v]))
    return arcs