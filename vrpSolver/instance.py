import random
import math

from .msg import *

def rndPlainNodes(
	N:			"Number of vertices" = None,
	nodeIDs: 	"A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
	xRange: 	"A 2-tuple with minimum/maximum range of x" = (0, 100),
	yRange: 	"A 2-tuple with minimum/maximum range of y" = (0, 100),
	) -> "A set of nodes with id start from 0 to N":

	# Check for required fields ===============================================
	if (N == None):
		print(ERROR_MISSING_N)
		return

	# Initialize ==============================================================
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
	nodes:		"If nodes are provided, time windows will be applied on those nodes" = None,
	nodeIDs: 	"A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
	xRange: 	"A 2-tuple with minimum/maximum range of x" = (0, 100),
	yRange: 	"A 2-tuple with minimum/maximum range of y" = (0, 100),
	timeSpan:	"Time span for generating time windows, starting from 0" = None,
	twType:		"Type of time windows\
				 1) String 'Known', all time windows are given priorly, or \
				 2) String 'Periodic', time windows are generated periodically, or \
				 3) String 'Random', time windows are generated randomly " = 'Random',
	twKnownArgs: "Dictionary, \
				{\
					'timeWindows': A list of 2-tuples, e.g., [(s1, e1), (s2, e2), ...]\
				}" = None,
	twPeriodicArgs: "Dictionary, \
				{\
					'cycleTime': cycleTime, \
					'startTimeInCycle': startTimeInCycle, \
					'endTimeInCycle': endTimeInCycle \
				}" = None,
	twRandomArgs: "Dictionary, \
				{\
					'indFlag': True if different nodes have different time windows, False otherwise,\
					'avgBtwDur': average interval time between time windows, exponential distributed,\
					'avgTWDur': average length of time windows, exponential distributed\
				}" = None
	) -> "A set of nodes with time windows":

	# Check for required fields ===============================================
	if (N == None and nodes == None):
		print(ERROR_MISSING_N)
		return
	if (timeSpan == None):
		print(ERROR_MISSING_TIMESPAN)
		return
	validTWOpts = ['Known', 'Periodic', 'Random']
	if (twType not in validTWOpts):
		print(ERROR_OPTS_TWTYPE % validTWOpts)
		return
	elif (twType == 'Known'):
		if (twKnownArgs == None):
			print(ERROR_MISSING_TWKNOWNARGS)
			return
		elif ('timeWindows' not in twKnownArgs):
			print (ERROR_MISSING_TWKNOWNARGS_TW)
			return
	elif (twType == 'Periodic'):
		if (twPeriodicArgs == None):
			print(ERROR_MISSING_TWPERIODICARGS)
			return
		elif ('cycleTime' not in twPeriodicArgs):
			print(ERROR_MISSING_TWPERIODICARGS_CT)
			return
		elif ('startTimeInCycle' not in twPeriodicArgs):
			print(ERROR_MISSING_TWPERIODICARGS_ST)
			return
		elif ('endTimeInCycle' not in twPeriodicArgs):
			print(ERROR_MISSING_TWPERIODICARGS_ET)
			return
	elif (twType == 'Random'):
		if (twRandomArgs == None):
			print(ERROR_MISSING_TWRANDOMARGS)
			return
		elif ('indFlag' not in twRandomArgs):
			print(ERROR_MISSING_TWRANDOMARGS_FG)
			return
		elif ('avgBtwDur' not in twRandomArgs):
			print(ERROR_MISSING_TWRANDOMARGS_BTW)
			return
		elif ('avgTWDur' not in twRandomArgs):
			print(ERROR_MISSING_TWRANDOMARGS_TW)
			return

	# Initialize nodes ========================================================
	if (nodes == None):
		nodes = rndPlainNodes(N, nodeIDs, xRange, yRange)

	# Generate time windows for different types ===============================
	for n in nodes:
		nodes[n]['timeWindows'] = []
	if (twType == "Known"):
		for n in nodes:
			nodes[n]['timeWindows'] = twKnownArgs['timeWindows']
	elif (twType == 'Periodic'):
		repeatNum = math.ceil(timeSpan / twPeriodicArgs['cycleTime'])
		for tw in range(repeatNum):
			for n in nodes:
				nodes[n]['timeWindows'].append((
					tw * twPeriodicArgs['cycleTime'] + twPeriodicArgs['startTimeInCycle'], 
					tw * twPeriodicArgs['cycleTime'] + twPeriodicArgs['endTimeInCycle']))
	elif (twType == 'Random'):
		avgBtwDur = twRandomArgs['avgBtwDur']
		avgTWDur = twRandomArgs['avgTWDur']
		# If all nodes have different time windows
		if (twRandomArgs['indFlag'] == False):
			now = 0
			pre = 0
			availFlag = True if (random.random() < (avgTWDur / (avgTWDur + avgBtwDur))) else False
			while (now < timeSpan):			
				interval = 0
				if (availFlag):
					interval = random.expovariate(1 / avgTWDur)
					now += interval
					if (interval > 0.0001):
						for n in nodes:
							nodes[n]['timeWindows'].append((pre, now))
				else:
					interval = random.expovariate(1 / avgBtwDur)
					now += interval
				availFlag = not availFlag
				pre = now
		# If all nodes have the same time windows
		else:
			for n in nodes:
				now = 0
				pre = 0
				availFlag = True if (random.random() < (avgTWDur / (avgTWDur + avgBtwDur))) else False
				while (now < timeSpan):
					interval = 0
					if (availFlag):
						interval = random.expovariate(1 / avgTWDur)
						now += interval
						if (interval > 0.0001):
							nodes[n]['timeWindows'].append((pre, now))
					else:
						interval = random.expovariate(1 / avgBtwDur)
						now += interval
					availFlag = not availFlag
					pre = now

	# Truncate last time window to fit time span ==============================
	for n in nodes:
		while (nodes[n]['timeWindows'][-1][0] >= timeSpan or nodes[n]['timeWindows'][-1][1] > timeSpan):
			lastStart = nodes[n]['timeWindows'][-1][0]
			if (lastStart >= timeSpan):
				for n in nodes:
					nodes[n]['timeWindows'] = nodes[n]['timeWindows'][:-1]
			lastEnd = nodes[n]['timeWindows'][-1][1]
			if (lastEnd >= timeSpan):
				for n in nodes:
					nodes[n]['timeWindows'][-1] = (
						nodes[n]['timeWindows'][-1][0], 
						timeSpan)

	return nodes