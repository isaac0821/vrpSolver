import heapq
from gurobipy import *

from .common import *

def findComponentsUndirected(
    arcs:       "A list of 2-tuples", 
    ) -> "A list of components":
    # Create adj list, each vertex start with an empty list
    adjList = arcs2AdjList(arcs)

    # Initialize
    found = {}
    for node in adjList:
        found[node] = 0
    components = []

    # Main algorithm, mark neighbors
    for i in adjList:
        comp = []
        q = []
        if (found[i] == 0):
            found[i] = 1
            comp.append(i)
            q.append(i)
            while (q):
                v = q.pop(0)
                for u in adjList[v]:
                    if (found[u] == 0):
                        found[u] = 1
                        comp.append(u)
                        q.append(u)
            components.append(comp)
    return components

def traversalTree(
    tree:   "Dictionary, returns the children of given nodeID, tuple if in order, list otherwise, \
                {\
                    nodeID1: (child1, child2, ...), \
                    nodeID2: (child1, child2, ...), \
                    nodeWithNoChild: None, \
                    ... \
                }" = None, 
    oID:    "1) String/Integer, nodeID of the root or, \
             2) None, (default) the first nodeID in `tree`" = None,
    algo:   "1) String, (default) 'DepthFirst' or, \
             2) String, 'BreadthFirst'" = 'DepthFirst'
    ) -> "Return a sequence of node ids that traverses the tree":

    # Subroutines =============================================================
    def traversalTreeDepthFirst():
        visited = []

        # Visit children recursively ------------------------------------------
        def visitNode(nodeID):
            visited.append(nodeID)
            children = tree[nodeID]
            if (children != None and children not in visited):
                if (type(children) == int or type(children) == str):
                    visitNode(children)
                else:
                    for child in children:
                        visitNode(child)

        # Start search from root ----------------------------------------------
        # FIXME! Incorrect for dictionary that root is not the first element
        if (oID == None):
            oID = list(tree.keys())[0]
        visitNode(oID)

        return {
            'seq': visited
        }

    # Solve by different algorithms ===========================================
    res = None
    if (algo == 'DepthFirst'):
        res = traversalTreeDepthFirst()

    return res

def traversalGraph(
    arcs:   "1) A list of 3-tuple (nodeID1, nodeID2, weight) or, \
             2) A list of 2-tuple (nodeID1, nodeID2)",
    oID:    "1) String/Integer, nodeID of the root or, \
             2) None, (default) the first nodeID in `tree`" = None,
    algo:   "1) String, (default) 'DepthFirst' or, \
             2) String, 'BreadthFirst'" = 'DepthFirst'
    ) -> "Return a sequence of node ids that traverses the tree":

    # Convert arcs into adjList ===============================================
    neighbors = arcs2AdjList(arcs)

    # Subroutines =============================================================
    def traversalGraphDepthFirst():
        visited = []
        
        # Visit neighbors that has not been visited ===============================
        def visitNode(nodeID):
            visited.append(nodeID)
            neis = neighbors[nodeID]
            for nei in neis:
                if (nei not in visited):
                    visitNode(nei)

        # Start search from root ==================================================
        if (oID == None):
            oID = list(neighbors.keys())[0]
        visitNode(oID)

        return {
            'seq': visited,
            'oID': oID
        }

    # Solve by different algorithms ===========================================
    res = None
    if (algo == 'DepthFirst'):
        res = traversalGraphDepthFirst()

    return res

def graphMST(
    weightArcs: "A list of 3-tuples, (ID1, ID2, weight), indexes of vertices must start from 0" = None,
    numVertices:"1) Integer, number of vertices, or \
                 2) (default) None, assuming all vertices are mentioned in `weightArcs`" = None,
    algo:       "1) String, (default) 'Krusal' or, \
                 2) String, (not available) 'Prim_AdjList' or, \
                 3) String, (not available) 'Prim_AdjMat' or, \
                 4) String, (not available) 'Boruvka' or, \
                 5) String, (not available) 'ReverseDelete'" = 'Krusal'
    ) -> "A list of weightArcs which forms a minimal spanning tree":

    # Number of vertices ======================================================
    if (numVertices == None):
        vertices = []
        for i in range(len(weightArcs)):
            if (weightArcs[i][0] not in vertices):
                vertices.append(weightArcs[i][0])
            if (weightArcs[i][1] not in vertices):
                vertices.append(weightArcs[i][1])
        numVertices = len(vertices)

    # Subroutines =============================================================
    def mstKrusal(weightArcs, numVertices):
        # Initialize ----------------------------------------------------------
        mst = []
        val = 0
        compList = []

        # Arc ranking ---------------------------------------------------------
        sortedWeightArcs = []
        for i in range(len(weightArcs)):
            heapq.heappush(sortedWeightArcs, (weightArcs[i][2], weightArcs[i]))
        
        # Krusal algorithm, add weightArcs between components -----------------
        while(len(mst) < numVertices - 1 and sortedWeightArcs):    
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

    # Call MST ================================================================
    if (algo == 'Krusal'):
        res = mstKrusal(weightArcs, numVertices)
    else:
        print("Error: Incorrect or not available MST option!")
    return res

def graphMatching(
    weightArcs: "A list of 3-tuples, (ID1, ID2, weight), indexes of vertices must start from 0" = None,
    mType:      "1) String, (default) 'Maximum' or \
                 2) String, 'Minimum'" = 'Maximum',
    algo:       "1) String, (default) 'Blossom' or, \
                 2) String, (not available) 'MAlterning' or, \
                 3) String, (not available) 'IP'" = 'Blossom'
    ) -> "Return a set of vertices that forms a Maximum/Minimum Matching": 

    # Subroutine ==============================================================
    def matchingIP(weightArcs, mType):
        matching = []
        M = Model('Matching')

        # Decision variables --------------------------------------------------
        x = {}
        for e in range(len(weightArcs)):
            x[e] = M.addVar(vtype = GRB.BINARY, obj = weightArcs[e][2])

        # Matching objective function -----------------------------------------
        if (mType == 'Maximum'):
            M.modelSense = GRB.MAXIMIZE
        elif (mType == 'Minimum'):
            M.modelSense = GRB.MINIMIZE
        M.update()

        # Perfect matching ----------------------------------------------------
        # First find neighborhoods
        neighborhoods = arcs2AdjList(weightArcs)
        for node in neighborhoods:
            neis = neighborhoods[node]
            neiArcs = []
            for nei in neis:
                for i in range(len(weightArcs)):
                    if ((node == weightArcs[i][0] and nei == weightArcs[i][1]) 
                        or (node == weightArcs[i][1] and nei == weightArcs[i][0])):
                        neiArcs.append(i)
            M.addConstr(quicksum(x[e] for e in neiArcs) == 1)

        # Matching ------------------------------------------------------------
        M.optimize()

        # Construct solution --------------------------------------------------
        ofv = None
        if (M.status == GRB.status.OPTIMAL):
            ofv = M.getObjective().getValue()
            for e in x:
                if (x[e].x > 0.8):
                    matching.append(weightArcs[e])

        return {
            'ofv': ofv, 
            'matching': matching
        }

    # Calculate matching using different algorithms ===========================
    res = None
    if (mType == 'Minimum'):
        if (algo == 'IP'):
            res = matchingIP(weightArcs, mType)
    elif (mType == 'Maximum'):
        if (algo == 'IP'):
            res = matchingIP(weightArcs, mType)

    return res

