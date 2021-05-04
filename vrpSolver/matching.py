from gurobipy import *

from .common import *

def graphMatching(
	weightArcs:	"A list of 3-tuples, (ID1, ID2, weight), indexes of vertices must start from 0" = None,
	numVertices:"1) Integer, number of vertices, or \
				 2) (default) None, assuming all vertices are mentioned in `weightArcs`" = None,
	mType:		"1) String, (default) 'Maximum' or \
				 2) String, 'Minimum'" = 'Maximum',
	algo:		"1) String, (default) 'Blossom' or, \
				 2) String, (not available) 'M_Alterning' or, \
				 2) String, (not available) 'IP'" = 'Blossom'
	) -> "Return a set of vertices that forms a Maximum/Minimum Matching": 

	# Calculate matching using different algorithms ===========================
	res = None
	if (algo == 'IP'):
		res = _matchingIP(weightArcs, mType)

	return res

def _matchingIP(weightArcs, mType):
	matching = []
	M = Model('Matching')

	# Decision variables ======================================================
	x = {}
	for e in range(len(weightArcs)):
		x[e] = M.addVar(vtype = GRB.BINARY, obj = weightArcs[e][2])

	# Matching objective function =============================================
	if (mType == 'Maximum'):
		M.modelSense = GRB.MAXIMIZE
	elif (mType == 'Minimum'):
		M.modelSense = GRB.MINIMIZE
	M.update()

	# Perfect matching ========================================================
	# First find neighborhoods
	neighborhoods = convertArcs2Neighbor(weightArcs)
	for node in neighborhoods:
		neis = neighborhoods[node]
		neiArcs = []
		for nei in neis:
			for i in range(len(weightArcs)):
				if ((node == weightArcs[i][0] and nei == weightArcs[i][1]) 
					or (node == weightArcs[i][1] and nei == weightArcs[i][0])):
					neiArcs.append(i)
		M.addConstr(quicksum(x[e] for e in neiArcs) == 1)

	# Matching ================================================================
	M.optimize()

	# Construct solution ======================================================
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