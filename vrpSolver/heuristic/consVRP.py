###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 12/04/2020                             #
# consVRP() - Constructive Heuristic for VRP, and typical variants            #
# - 11/09/2020 DepthFirst Method                                              #
# - 11/09/2020 Christofides Algorithm                                         #
# - 11/11/2020 Nearest Neighborhood Method                                    #
# - 11/11/2020 Farthest Neighborhood Method                                   #
# - 11/11/2020 Random Sequence                                                #
###############################################################################

from vrpSolver.common import *
from vrpSolver.graph.mst import *
from vrpSolver.graph.basic import *
from vrpSolver.graph.search import *
from vrpSolver.graph.matching import *

def consVRP(
	nodes: 	"Dictionary, returns the coordinate of given nodeID, \
				{\
					nodeID1: {'loc': (x, y), 'demand': d, 'tStart': t1, 'tEnd': t2}, \
					nodeID2: {'loc': (x, y), 'demand': d, 'tStart': t1, 'tEnd': t2}, \
					... \
				}" = None,
	tau:	"1) String 'Euclidean' or \
			 2) String (default) 'SphereEuclidean' or \
			 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "SphereEuclidean",
	nodeIDs:"1) String (default) 'All', or \
			 2) A list of node IDs" = 'All',
	algo:	"1) String 'ClarkeWright' or \
			 2) String 'Sweep' or \
			 3) String 'SortFirstSplitSecond' or \
			 4) String 'ClusterFirstRouteSecond'" = 'ClarkeWright'
	) -> "Heuristic solution for VRP":

	return

def _consVRPClarkeWright(nodes):

	return