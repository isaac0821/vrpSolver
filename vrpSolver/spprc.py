import heapq
import math
import warnings
import networkx as nx

from .common import *
from .geometry import *
from .msg import *

def initialLabelFunc(oID):
    return {
        'path': [oID],
        'dist': 0,
        'load': 0,
        'time': 0
    }

def feasiblityFunc(curPath, nextNode):
    # Check if node has been covered
    if (nextNode in curPath['path']):
        return False
    # Check time windows
    lastNode = curPath['path'][-1]
    arrTime = max(
        curPath['time'] + g.edges[lastNode, nextNode]['travelDist'],
        g.nodes[nextNode]['timeWindow'][0])

    availTW = []
    if (nextNode != dupDepotID):
        availTW = g.nodes[nextNode]['timeWindow']
    else:         `                                                                                                                                                                                                                                                                     
        availTW = g.nodes[depotID]['timeWindow']

    if (arrTime < availTW[0] or arrTime > availTW[1]):
        return False
    # Check loads
    accLoad = curPath['load'] + g.nodes[nextNode]['weight']
    if (accLoad > vehCap):
        return False
    return True

def solveSPPRC(
    nodes: dict,
    tau: dict,
    oID: int|str,
    dID: int|str,
    initialLabelFunc,
    feasiblityFunc,
    extendFunc,
    dominateFunc 
    ) -> dict: 

    # Create graph ============================================================
    g = nx.DiGraph()
    for n in nodes:
        g.add_node(n)
    for (i, j) in tau:
        g.add_edge(i, j)
