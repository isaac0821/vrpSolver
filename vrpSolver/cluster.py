import networkx as nx

from .common import *

# 【Constructing]
def getClusterByRadius(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'SphereEuclidean' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    critera:    "The shortest distances between nodes that two nodes still considered to be in the cluster" = None,
    ) -> "Given a set of plain nodes, a radius, returns a list of clusters that can be covered by a circle of that radius":

    # Initialize ==============================================================
    cluster = []

    # Get distance from one node to all the other nodes =======================
    # Define edges ============================================================
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            edges = getTauEuclidean(nodes)
        elif (edges == 'SphereEuclidean'):
            edges = getTauSphereEuclidean(nodes)
        else:
            print("Error: Incorrect type `edges`")
            return None

    dist = {}
    for n1 in nodes:
        for n2 in nodes:
            if (n2 > n1 and edges[n1, n2] <= criteria['closeEnough']):
                dist[n1].append(n2)
                dist[n2].append(n1)

    # Create the graph ========================================================
    G = nx.Graph()
    for n1 in nodes:
        G.add_node(n1)
        for n2 in dist[n1]:
            G.add_edge(n1, n2)

    # Find all cliques ========================================================
    cluster = nx.enumerate_all_cliques(G)


    return cluster