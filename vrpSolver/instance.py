import random

def rndTSPInstance(
	N: "Number of vertices" = None,
	nodeIDs: "A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
	xRange: "A 2-tuple with minimum/maximum range of x" = (0, 100),
	yRange: "A 2-tuple with minimum/maximum range of y" = (0, 100),
	) -> "A set of nodes with id start from 0 to N":

	# The Dictionary
	nodeLoc = {}

	# Node IDs
	if (nodeIDs == None):
		nodeIDs = [i for i in range(N)]

	# Generate instance
	for i in nodeIDs:
		x = random.randrange(xRange[0], xRange[1])
		y = random.randrange(yRange[0], yRange[1])
		nodeLoc[i] = {'loc': (x, y)}
	return nodeLoc

def rndVRPInstance(

	):

	return {
		'nodes': nodes,
		'edges': edges,
		'vehicles': vehicles
	}