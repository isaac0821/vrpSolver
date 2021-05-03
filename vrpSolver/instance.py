import random

def rndPlainNodes(
	N:			 "Number of vertices" = None,
	nodeIDs: 	"A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
	xRange: 	"A 2-tuple with minimum/maximum range of x" = (0, 100),
	yRange: 	"A 2-tuple with minimum/maximum range of y" = (0, 100),
	) -> "A set of nodes with id start from 0 to N":

	# The Dictionary
	nodes = {}

	# Node IDs
	if (nodeIDs == None):
		nodeIDs = [i for i in range(N)]

	# Generate instance
	for i in nodeIDs:
		x = random.randrange(xRange[0], xRange[1])
		y = random.randrange(yRange[0], yRange[1])
		nodes[i] = {'loc': (x, y)}
		
	return nodes

def rndTimeWindowsNodes(
	N:			"Number of vertices" = None,
	nodeIDs: 	"A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
	xRange: 	"A 2-tuple with minimum/maximum range of x" = (0, 100),
	yRange: 	"A 2-tuple with minimum/maximum range of y" = (0, 100),
	twType:		"Type of time windows\
				 1) String, 'Fixed' \
				 2) String, 'Random'" = 'Random',
	twFixArgs: 	"Dictionary, \
				{\
					'startTime': startTime,\
					'endTime': endTime,\
				}" = None,
	twRndArgs:	"Dictionary, \
				{\
					'avgIntervalDuration': average interval time between time windows,\
					'avgTimeWindowDuration': average length of time windows\
				}" = None
	):

	return nodes