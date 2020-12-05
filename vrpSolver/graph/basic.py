import random

def findComponentsUndirected(
	arcs:		"A list of 2-tuples", 
	) -> "A list of components":
	# Create adj list, each vertex start with an empty list
	adjList = {}
	for e in arcs:
		if (e[0] not in adjList):
			adjList[e[0]] = [e[1]]
		else:
			adjList[e[0]].append(e[1])
		if (e[1] not in adjList):
			adjList[e[1]] = [e[0]]
		else:
			adjList[e[1]].append(e[0])

	# Initialize
	found = {}
	for node in adjList:
		found[node] = 0
	components = []

	# Main algorithm, mark neighbors
	for i in adjList:
		comp = []
		q = []
		if (found[i] == 0):
			found[i] = 1
			comp.append(i)
			q.append(i)
			while (q):
				v = q.pop(0)
				for u in adjList[v]:
					if (found[u] == 0):
						found[u] = 1
						comp.append(u)
						q.append(u)
			components.append(comp)

	return components

def calSeqCostArcs(
	weightArcs:	"A list of 3-tuple (nodeID1, nodeID2, weight)",
	seq:		"List, sequence of visiting node ids"
	) -> "Return the cost on the graph given a list of arcs weights":

	# Accumulate costs ========================================================
	cost = 0
	for i in range(len(seq) - 1):
		c = None
		for j in range(len(weightArcs)):
			if (seq[i] == weightArcs[j][0] and seq[i + 1] == weightArcs[j][1]):
				c = weightArcs[j][2]
				break
			elif (seq[i] == weightArcs[j][1] and seq[i + 1] == weightArcs[j][0]):
				c = weightArcs[j][2]
				break
		if (c == None):
			print("Error: Missing arc (%s, %s) in `weightArcs`" % (seq[i], seq[i + 1]))
			return
		else:
			cost += c

	return cost

def calSeqCostMatrix(
	tau: "Dictionary {(nodeID1, nodeID2): dist, ...}", 
	seq: "List, sequence of visiting node ids"
	) -> "Return the cost on the graph given cost matrix/dictionary tau":

	# Accumulate costs ========================================================
	cost = 0
	for i in range(len(seq) - 1):
		cost += tau[seq[i], seq[i + 1]]

	return cost

def convertArcs2Neighbor(
	arcs:	"1) A list of 3-tuple (nodeID1, nodeID2, weight) or, \
			 2) A list of 2-tuple (nodeID1, nodeID2)"
	) -> "Dictionary of neighbors of each node":

	neighbors = {}
	for i in range(len(arcs)):
		if (arcs[i][0] not in neighbors):
			neighbors[arcs[i][0]] = [arcs[i][1]]
		else:
			neighbors[arcs[i][0]].append(arcs[i][1])
		if (arcs[i][1] not in neighbors):
			neighbors[arcs[i][1]] = [arcs[i][0]]
		else:
			neighbors[arcs[i][1]].append(arcs[i][0])

	return neighbors