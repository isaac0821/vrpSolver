import matplotlib.pyplot as plt

def plotSeq(
	nodeLoc:"Dictionary, returns the coordinate of given nodeID, \
			{\
				nodeID1: (lat, lon), \
				nodeID2: (lat, lon), \
				... \
			}" = None, 
	nodeIDs:"1) String (default) 'All', or \
			 2) A list of node IDs" = 'All',
	seq:	"List of visiting sequences" = None,
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

	# Draw Seq ================================================================
	lx = []
	ly = []
	for s in range(len(seq) - 1):
		lx.append(nodeLoc[seq[s]][0])
		ly.append(nodeLoc[seq[s]][1])
	lx.append(nodeLoc[seq[len(seq) - 1]][0])
	ly.append(nodeLoc[seq[len(seq) - 1]][1])
	plt.plot(lx, ly, color=color)

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