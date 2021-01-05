import matplotlib.pyplot as plt
import random

def randomColor():
	color = "#%06x" % random.randint(0, 0xFFFFFF)
	return color

def plotTSP(
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	seq:	"List of visiting sequences" = None,
	color:	"1) String 'Random', or\
			 2) List of color of seq i" = 'Random'
	) -> "Draw TSP solution":

	# Draw nodes ==============================================================
	x = []
	y = []
	for node in nodes:
		x.append(nodes[node]['loc'][0])
		y.append(nodes[node]['loc'][1])
		plt.annotate(node, (nodes[node]['loc'][0], nodes[node]['loc'][1]))
	plt.plot(x, y, 'ro')

	# Draw Seq ================================================================
	if (seq != None):
		for i in range(len(seq)):
			lx = []
			ly = []
			newColor = randomColor()
			for s in range(len(seq[i]) - 1):
				lx.append(nodes[seq[i][s]]['loc'][0])
				ly.append(nodes[seq[i][s]]['loc'][1])
			lx.append(nodes[seq[i][len(seq[i]) - 1]]['loc'][0])
			ly.append(nodes[seq[i][len(seq[i]) - 1]]['loc'][1])
			if (color == 'Random'):
				plt.plot(lx, ly, color=newColor)
			else:
				plt.plot(lx, ly, color=color[i])

	return

def plotSeq(
	nodeLoc:"Dictionary, returns the coordinate of given nodeID, \
			{\
				nodeID1: (lat, lon), \
				nodeID2: (lat, lon), \
				... \
			}" = None, 
	seq:	"List of visiting sequences" = None,
	color:	"1) String 'Random', or\
			 2) List of color of seq i" = 'Random'
	) -> "Draw sequences of visiting":

	# Define nodeIDs ==========================================================
	if (type(nodeIDs) is not list):
		if (nodeIDs == 'All'):
			nodeIDs = []
			for i in nodeLoc:
				nodeIDs.append(i)

	# Draw nodes ==============================================================
	x = []
	y = []
	for node in nodeIDs:
		x.append(nodeLoc[node][0])
		y.append(nodeLoc[node][1])
		plt.annotate(node, (nodeLoc[node][0], nodeLoc[node][1]))
	plt.plot(x, y, 'ro')

	# Draw Seq ================================================================
	for i in range(len(seq)):
		lx = []
		ly = []
		newColor = randomColor()
		for s in range(len(seq[i]) - 1):
			lx.append(nodeLoc[seq[i][s]][0])
			ly.append(nodeLoc[seq[i][s]][1])
		lx.append(nodeLoc[seq[i][len(seq[i]) - 1]][0])
		ly.append(nodeLoc[seq[i][len(seq[i]) - 1]][1])
		if (color == 'Random'):
			plt.plot(lx, ly, color=newColor)
		else:
			plt.plot(lx, ly, color=color[i])

	return

def plotArcs(
	nodeLoc:"Dictionary, returns the coordinate of given nodeID, \
			{\
				nodeID1: (lat, lon), \
				nodeID2: (lat, lon), \
				... \
			}" = None, 
	nodeIDs:"1) String (default) 'All', or \
			 2) A list of node IDs" = 'All',
	arcs:	"List of arcs" = None,
	color:	"Color of seq" = 'green'
	) -> "Draw sequences of visiting":
	
	# Define nodeIDs ==========================================================
	if (type(nodeIDs) is not list):
		if (nodeIDs == 'All'):
			nodeIDs = []
			for i in nodeLoc:
				nodeIDs.append(i)

	# Draw nodes ==============================================================
	x = []
	y = []
	for node in nodeIDs:
		x.append(nodeLoc[node][0])
		y.append(nodeLoc[node][1])
		plt.annotate(node, (nodeLoc[node][0], nodeLoc[node][1]))
	plt.plot(x, y, 'ro')

	# Draw arcs ===============================================================
	for arc in arcs:
		lx = [nodeLoc[arc[0]][0], nodeLoc[arc[1]][0]]
		ly = [nodeLoc[arc[0]][1], nodeLoc[arc[1]][1]]
		plt.plot(lx, ly, color=color)

	return