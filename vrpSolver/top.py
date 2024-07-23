import heapq
import math
import warnings
import networkx as nx

from .common import *
from .geometry import *
from .msg import *

def ipTOP(
    nodes: dict, 
    maxBudget: float,
    locFieldName: str = 'loc',
    priceFieldName: str = 'price',
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    vehicles: dict = {
        0: {'speed': 1, 'depotID': 0, 'startID': 0, 'endID': 0}
    },
    vehicleID: int|str = 0,
    edges: dict = {
        'method': "Euclidean", 
        'ratio': 1
    },
    method: dict = {
        'fml': 'MTZ',
        'solver': 'Gurobi',
        'timeLimit': None,
        'outputFlag': False,
        'env': None
    },
    detailsFlag: bool = False,
    metaFlag: bool = False
    ) -> dict|None:

    return