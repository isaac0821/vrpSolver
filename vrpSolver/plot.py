import matplotlib.pyplot as plt
import random

def randomColor():
	color = "#%06x" % random.randint(0, 0xFFFFFF)
	return color

def plotGantt(
	fig:		"Based matplotlib figure object" = None,
	ax:			"Based matplotlib ax object" = None,
	gantt:		"List of dictionaries, in the following format\
					[{\
						'entityID': entityID, \
						'timeWindow': [startTime, endTime], \
						'color': color, \
						'style': 'solid' \
					}, ... , \
					{\
						'entityID': entityID, \
						'timeStamps': [timeStamp1, timeStamp2, ..., timeStampN], \
						'color': color, \
						'style': 'solid' \
					}]\
				" = None,
	entities:	"The Gantt chart will be drawn in this order, if None, takes the entities in `gantt`" = None,
	startTime:	"Start time of Gantt, default to be 0, if None, use the earliest time in `gantt`" = 0,
	endTime:	"End time of Gantt, default to be None, if None, use the latest time in `gantt`" = None,
	width:		"Width of the figure" = 12,
	height:		"Height of the figure" = 5,
	saveFigPath:"1) None, if not exporting image, or \
				 2) String, the path for exporting image" = None
	) -> "Given a Gantt dictionary, plot Gantt":

	# Precalculate ============================================================
	s = None
	e = None
	entList = []
	for g in gantt:
		if (g['entityID'] not in entList):
			entList.append(g['entityID'])
		if ('timeWindow' in g):
			if (s == None or s > g['timeWindow'][0]):
				s = g['timeWindow'][0]
			if (e == None or e < g['timeWindow'][1]):
				e = g['timeWindow'][1]
		elif ('timeStamps' in g):
			if (s == None or s > g['timeStamps'][0]):
				s = g['timeStamps'][0]
			if (e == None or e < g['timeStamps'][-1]):
				e = g['timeStamps'][-1]

	if (entities != None):
		for g in entList:
			if (g not in entities):
				print("ERROR: Missing entity in `entities`")
				return
	else:
		entities = [i for i in entList]
	if (startTime != None):
		if (startTime > s):
			print("WARNING: `startTime` later than earliest time in `gantt`, auto-fixed")
			startTime = s
	else:
		startTime = s
	if (endTime != None):
		if (endTime < e):
			print("WARNING: `endTime` earlier than latest time in `gantt`, auto-fixed")
			endTime = e
	else:
		endTime = e

	# If no based matplotlib figure, define fig size ==========================
	if (fig == None or ax == None):
		fig, ax = plt.subplots()
		fig.set_figheight(height)
		fig.set_figwidth(width)
		ax.set_xlim(startTime, endTime)
		ax.set_ylim(0, len(entities))

	# Set axis ================================================================
	ax.set_yticks([i + 0.5 for i in range(len(entities))])
	entities.reverse()
	ax.set_yticklabels(entities)
	entities.reverse()
	ax.set_xlabel("Time")

	# Loop through `gantt` and draw gantt =====================================
	for g in gantt:
		bottom = len(entities) - entities.index(g['entityID']) - 1
		top = len(entities) - entities.index(g['entityID'])
		if ('timeWindow' in g):
			s = g['timeWindow'][0]
			e = g['timeWindow'][1]
			x = [s, s, e, e, s]
			y = [bottom, top, top, bottom, bottom]
			ax.plot(x, y, color = 'black', linewidth = 2)
			ax.fill(x, y, color = g['color'])
			if (g['style'] == 'shadow'):
				ax.fill(x, y, hatch = "/", fill=False)
		elif ('timeStamps' in g):
			for i in range(len(g['timeStamps']) - 1):
				s = g['timeStamps'][i]
				e = g['timeStamps'][i + 1]
				x = [s, s, e, e, s]
				y = [bottom, top, top, bottom, bottom]
				ax.plot(x, y, color = 'black', linewidth = 2)
				ax.fill(x, y, color = g['color'])
				if (g['style'] == 'shadow'):
					ax.fill(x, y, hatch = "/", fill=False)

	# Save figure =============================================================
	if (saveFigPath != None):
		fig.savefig(saveFigPath)

	return fig, ax

def plotNodes(
	fig:		"Based matplotlib figure object" = None,
	ax:			"Based matplotlib ax object" = None,
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	color:		"1) String 'Random', or\
				 2) String, color" = 'Random',
	saveFigPath:"1) None, if not exporting image, or \
				 2) String, the path for exporting image" = None
	) -> "Draw nodes":

	# If no based matplotlib figure, define boundary ==========================
	if (fig == None or ax == None):
		fig, ax = plt.subplots()
		allX = []
		allY = []
		for i in nodes:
			allX.append(nodes[i]['loc'][0])
			allY.append(nodes[i]['loc'][1])
		xMin = min(allX) - 0.25
		xMax = max(allX) + 0.25
		yMin = min(allY) - 0.25
		yMax = max(allY) + 0.25
		xSpan = None
		ySpan = None
		if (xMax - xMin > yMax - yMin):
			xSpan = 20
			ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
		else:
			xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
			ySpan = 20
		fig.set_figheight(ySpan)
		fig.set_figwidth(xSpan)
		ax.set_xlim(xMin, xMax)
		ax.set_ylim(yMin, yMax)

	# Draw nodes ==============================================================
	x = []
	y = []
	for node in nodes:
		x.append(nodes[node]['loc'][0])
		y.append(nodes[node]['loc'][1])
		ax.annotate(node, (nodes[node]['loc'][0], nodes[node]['loc'][1]))
	ax.plot(x, y, 'ro')

	# Save figure =============================================================
	if (saveFigPath != None):
		fig.savefig(saveFigPath)

	return fig, ax

def plotArcs(
	fig:		"Based matplotlib figure object" = None,
	ax:			"Based matplotlib ax object" = None,
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	arcs:		"List of 2-tuples, arcs in format of [(nodeID1, nodeID2), ...]" = None,
	arrowFlag: 	"Boolean, whether or not add arrows to route" = True,
	color:		"1) String 'Random', or\
				 2) String, color" = 'Random',
	saveFigPath:"1) None, if not exporting image, or \
				 2) String, the path for exporting image" = None
	) -> "Draw arcs":

	# If no based matplotlib figure, define boundary ==========================
	if (fig == None or ax == None):
		fig, ax = plt.subplots()
		allX = []
		allY = []
		for i in arcs:
			allX.append(i[0][0])
			allY.append(i[0][1])
			allX.append(i[1][0])
			allY.append(i[1][1])
		xMin = min(allX) - 0.25
		xMax = max(allX) + 0.25
		yMin = min(allY) - 0.25
		yMax = max(allY) + 0.25
		xSpan = None
		ySpan = None
		if (xMax - xMin > yMax - yMin):
			xSpan = 20
			ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
		else:
			xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
			ySpan = 20
		fig.set_figheight(ySpan)
		fig.set_figwidth(xSpan)
		ax.set_xlim(xMin, xMax)
		ax.set_ylim(yMin, yMax)

	# Draw arcs ===============================================================
	for arc in arcs:
		lx = [nodeLoc[arc[0]][0], nodeLoc[arc[1]][0]]
		ly = [nodeLoc[arc[0]][1], nodeLoc[arc[1]][1]]
		if (color == 'Random'):
			rndColor = randomColor()
			plt.plot(lx, ly, color = rndColor)
		else:
			plt.plot(lx, ly, color = color)

	# Save figure =============================================================
	if (saveFigPath != None):
		fig.savefig(saveFigPath)

	return fig, ax

def plotRoutes(
	fig:		"Based matplotlib figure object" = None,
	ax:			"Based matplotlib ax object" = None,
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	seq:		"List of nodeIDs, visiting sequences" = None,
	arrowFlag:	"Boolean, whether or not add arrows to route" = True,
	color:		"1) String 'Random', or\
				 2) String, color" = 'Random',
	saveFigPath:"1) None, if not exporting image, or \
				 2) String, the path for exporting image" = None
	) -> "Draw a route, e.g., TSP solution":

	# If no based matplotlib figure, define boundary ==========================
	if (fig == None or ax == None):
		fig, ax = plt.subplots()
		allX = []
		allY = []
		for i in seq:
			allX.append(nodes[i]['loc'][0])
			allY.append(nodes[i]['loc'][1])
		xMin = min(allX) - 0.25
		xMax = max(allX) + 0.25
		yMin = min(allY) - 0.25
		yMax = max(allY) + 0.25
		xSpan = None
		ySpan = None
		if (xMax - xMin > yMax - yMin):
			xSpan = 20
			ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
		else:
			xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
			ySpan = 20
		fig.set_figheight(ySpan)
		fig.set_figwidth(xSpan)
		ax.set_xlim(xMin, xMax)
		ax.set_ylim(yMin, yMax)

	# Draw Seq ================================================================
	rndColor = randomColor()
	if (seq != None and len(seq) > 0):
		lx = []
		ly = []
		for s in range(len(seq) - 1):
			lx.append(nodes[seq[s]]['loc'][0])
			ly.append(nodes[seq[s]]['loc'][1])
		lx.append(nodes[seq[len(seq) - 1]]['loc'][0])
		ly.append(nodes[seq[len(seq) - 1]]['loc'][1])
		if (color == 'Random'):
			ax.plot(lx, ly, color=rndColor)
		else:
			ax.plot(lx, ly, color=color)

	# Draw arrows =============================================================
	if (seq != None and len(seq) > 0):
		for s in range(len(seq) - 1):
			x1 = nodes[seq[s]]['loc'][0]
			y1 = nodes[seq[s]]['loc'][1]
			x2 = nodes[seq[s + 1]]['loc'][0]
			y2 = nodes[seq[s + 1]]['loc'][1]
			dx = x2 - x1
			dy = y2 - y1
			if (color == 'Random'):
				ax.arrow(x=x1, y=y1, dx=dx, dy=dy, linewidth=1, color=rndColor)
				ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=1, head_width=0.1, head_length=0.2, color=rndColor)
			else:
				ax.arrow(x=x1, y=y1, dx=dx, dy=dy, linewidth=1, color=color)
				ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=1, head_width=0.1, head_length=0.2, color=color)

	# Save figure =============================================================
	if (saveFigPath != None):
		fig.savefig(saveFigPath)

	return fig, ax

def plotReqs(
	fig:		"Based matplotlib figure object" = None,
	ax:			"Based matplotlib ax object" = None,
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	reqs:		"Dictionary, includes all pickup and delivery requests, in the format of \
				{\
					'reqID': {\
						'pickup': pickupZipCode,\
						'delivery': deliveryZipCode,\
						'size': itemSize, \
					}, \
					'reqID': ... \
				}" = None,
	color:		"1) String 'Random', or \
				 2) String of color" = 'Random',
	saveFigPath:"1) None, if not exporting image, or \
				 2) String, the path for exporting image" = None
	) -> "Plot delivery requests":

	# If no based matplotlib figure, define boundary ==========================
	if (fig == None or ax == None):
		fig, ax = plt.subplots()
		allX = []
		allY = []
		for r in reqs:
			allX.append(nodes[reqs[r]['pickup']]['loc'][0])
			allX.append(nodes[reqs[r]['delivery']]['loc'][0])
			allY.append(nodes[reqs[r]['pickup']]['loc'][1])
			allY.append(nodes[reqs[r]['delivery']]['loc'][1])
		xMin = min(allX) - 0.25
		xMax = max(allX) + 0.25
		yMin = min(allY) - 0.25
		yMax = max(allY) + 0.25
		xSpan = None
		ySpan = None
		if (xMax - xMin > yMax - yMin):
			xSpan = 20
			ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
		else:
			xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
			ySpan = 20
		fig.set_figheight(ySpan)
		fig.set_figwidth(xSpan)
		ax.set_xlim(xMin, xMax)
		ax.set_ylim(yMin, yMax)

	# Each request will be plot as arrow ======================================
	for r in reqs:
		x1 = nodes[reqs[r]['pickup']]['loc'][0]
		y1 = nodes[reqs[r]['pickup']]['loc'][1]
		x2 = nodes[reqs[r]['delivery']]['loc'][0]
		y2 = nodes[reqs[r]['delivery']]['loc'][1]
		dx = x2 - x1
		dy = y2 - y1

		# Define color
		if (color == 'Random'):
			rndColor = randomColor()
			ax.arrow(x=x1, y=y1, dx=dx, dy=dy, linewidth=1, color=rndColor)
			ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=1, head_width=0.1, head_length=0.2, color=rndColor)
		else:
			ax.arrow(x=x1, y=y1, dx=dx, dy=dy, linewidth=1, color=color)
			ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=1, head_width=0.1, head_length=0.2, color=color)	
		ax.annotate(r, (x1 + dx / 2, y1 + dy / 2))
	
	# Save figure =============================================================
	if (saveFigPath != None):
		fig.savefig(saveFigPath)

	return fig, ax

def plotActions(
	fig:		"Based matplotlib figure object" = None,
	ax:			"Based matplotlib ax object" = None,
	nodes:		"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: {'loc': (x, y)}, \
						nodeID2: {'loc': (x, y)}, \
						... \
					}" = None, 
	actions:	"List of 3-tuples, in format of [('pickup', reqID, nodeID), ...]" = [],
	color:		"1) String 'Random', or \
				 2) String of color" = 'Random',
	saveFigPath:"1) None, if not exporting image, or \
				 2) String, the path for exporting image" = None
	) -> "Plot delivery requests":

	# If no based matplotlib figure, define boundary ==========================
	if (fig == None or ax == None):
		fig, ax = plt.subplots()
		allX = []
		allY = []
		for r in actions:
			allX.append(nodes[r[2]]['loc'][0])
			allY.append(nodes[r[2]]['loc'][1])
		xMin = min(allX) - 0.25
		xMax = max(allX) + 0.25
		yMin = min(allY) - 0.25
		yMax = max(allY) + 0.25
		xSpan = None
		ySpan = None
		if (xMax - xMin > yMax - yMin):
			xSpan = 20
			ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
		else:
			xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
			ySpan = 20
		fig.set_figheight(ySpan)
		fig.set_figwidth(xSpan)
		ax.set_xlim(xMin, xMax)
		ax.set_ylim(yMin, yMax)

	# Draw actions ============================================================
	route = []
	for act in actions:
		route.append(act[2])
	fig, ax = plotRoutes(
		fig = fig,
		ax = ax,
		nodes = nodes,
		seq = route,
		color = color,
		arrowFlag = True,
		saveFigPath = saveFigPath)

	return fig, ax
