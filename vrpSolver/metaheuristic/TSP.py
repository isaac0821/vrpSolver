###############################################################################
#                       Written by Isaac0821 (Lan Peng)                       #
#                             Version: 11/11/2020                             #
# saTSP - Use Simulated Annealing to find TSP solution                        #
###############################################################################

import datetime

from vrpSolver.common import *
from vrpSolver.graph.basic import *
from vrpSolver.heuristic.TSP import *

def saTSP(
	nodeLoc:	"Dictionary, returns the coordinate of given nodeID, \
					{\
						nodeID1: (lat, lon), \
						nodeID2: (lat, lon), \
						... \
					}" = None, 
	tau:		"1) String 'Euclidean' or \
				 2) String (default) 'SphereEuclidean' or \
				 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "SphereEuclidean",
	nodeIDs:	"1) String (default) 'All', or \
				 2) A list of node IDs" = 'All',
	initSol:	"1) String, 'NearestNeighbor' or, \
				 2) String, 'Random' or, \
				 3) String, 'FarthestNeighbor'" = 'Random',
	initTemp:	"Float, Initial temperature" = None, 
	lengTemp:	"Integer, Temperature length" = None,
	neighRatio:	"A 4-tuple, sum to 1, (swap, exchange, rotate)" = (0.2, 0.2, 0.3),
	coolRate:	"Float, Temperature drop rate, (0, 1)" = None,
	stopType:	"List, with options as follows: \
				 1) String, 'Final_Temperature' or \
				 2) String, 'Num_Iterations_Without_Improving' or \
				 3) String, 'Percent_of_Accepted_Move' or \
				 4) String, 'Num_Iterations' or \
				 5) String, 'Executed_Time' " = ['Final_Temperature'],
	stopTemp:	"Float, stopping temperature, if stopType is 'Final_Temperature'" = None,
	stopNoImp:	"Integer, number of iteration without improving, if stopType is 'Num_Iterations_Without_Improving" = None,
	stopAptRate:"Float, ratio of acceptance, (0, 1), if stopType is 'Percent_of_Accepted_Move'" = None,
	stopIter:	"Integer, number of iteration, if stopType is 'Num_Iterations'" = None, 
	stopTime:	"Float, number of seconds before it stops, if stopType is 'Executed_Time'" = None
	) -> "TSP tour":

	# Initial clock ===========================================================
	startTime = datetime.datetime.now()

	# Define nodeIDs ==========================================================
	if (type(nodeIDs) is not list):
		if (nodeIDs == 'All'):
			nodeIDs = []
			for i in nodeLoc:
				nodeIDs.append(i)

	# Define tau ==============================================================
	if (type(tau) is not dict):
		lstNodeID = nodeIDs.copy()
		if (tau == 'Euclidean'):
			tau = getTauEuclidean(nodeLoc, lstNodeID)
		elif (tau == 'SphereEuclidean'):
			tau = getTauSphereEuclidean(nodeLoc, lstNodeID)
		else:
			print("Error: Incorrect type `tau`")
			return None

	# Subroutines to generate neighborhoods ===================================
	# Swap two nearby vertices
	def swap(seq):
		N = len(seq)
		i = random.randint(0, N - 1)

		# newSeq
		newSeq = [k for k in seq]
		t = newSeq[i]
		j = iterSeq(N, i, 'next')
		newSeq[i] = newSeq[j]
		newSeq[j] = t

		# deltaC = newC - preC
		deltaC = ((tau[seq[iterSeq(N, i, 'prev')], seq[j]]
				 + tau[seq[i], seq[iterSeq(N, j, 'next')]])
			    - (tau[seq[iterSeq(N, i, 'prev')], seq[i]]
				 + tau[seq[j], seq[iterSeq(N, j, 'next')]]))

		return {
			'seq': newSeq,
			'deltaC': deltaC
		}
	# Randomly exchange two vertices
	def exchange(seq):
		# Randomly choose i, j
		N = len(seq)
		i = None
		j = None
		while (i == None or j == None or abs(i - j) <= 2 or (i == 0 and j == len(seq) - 1) or (i == len(seq) - 1 and j == 0)):
			i = random.randint(0, N - 1)
			j = random.randint(0, N - 1)

		# new seq
		newSeq = [k for k in seq]
		t = newSeq[i]
		newSeq[i] = newSeq[j]
		newSeq[j] = t

		# deltaC = newC - preC	
		deltaC = ((tau[seq[iterSeq(N, i, 'prev')], seq[j]] 
				 + tau[seq[j], seq[iterSeq(N, i, 'next')]] 
				 + tau[seq[iterSeq(N, j, 'prev')], seq[i]]
				 + tau[seq[i], seq[iterSeq(N, j, 'next')]])
			    - (tau[seq[iterSeq(N, i, 'prev')], seq[i]] 
				 + tau[seq[i], seq[iterSeq(N, i, 'next')]] 
				 + tau[seq[iterSeq(N, j, 'prev')], seq[j]]
				 + tau[seq[j], seq[iterSeq(N, j, 'next')]]))

		return {
			'seq': newSeq,
			'deltaC': deltaC
		}
	# Randomly rotate part of seq
	def rotate(seq):
		# randomize i, j
		N = len(seq)
		i = None
		j = None
		while (i == None or j == None or j - i <= 2 or (i == 0 and j == len(seq) - 1)):
			i = random.randint(0, N - 1)
			j = random.randint(0, N - 1)

		# new seq
		newSeq = [seq[k] for k in range(i)]
		newSeq.append(seq[j])
		newSeq.extend([seq[j - k - 1] for k in range(j - i - 1)])
		newSeq.append(seq[i])
		newSeq.extend([seq[k] for k in range(j + 1, N)])

		# deltaC = newC - preC
		deltaC = ((tau[seq[iterSeq(N, i, 'prev')], seq[j]]
				 + tau[seq[i], seq[iterSeq(N, j, 'next')]])
			    - (tau[seq[iterSeq(N, i, 'prev')], seq[i]]
				 + tau[seq[j], seq[iterSeq(N, j, 'next')]]))

		return {
			'seq': newSeq,
			'deltaC': deltaC
		}

	# Initialize ==============================================================
	# Initial temperature
	T = initTemp
	# Temperature length (maximum temperature iteration)
	L = lengTemp
	# Initial Solution
	curSeq = None
	ofv = None
	if (initSol in ['NearestNeighbor', 'Random', 'FarthestNeighbor']):
		res = consTSP(nodeLoc, tau, nodeIDs, initSol)
		curSeq = res['seq'][:-1] # To avoid all kind of trouble, seq here is not closed
		ofv = res['ofv']
	else:
		return None

	# Main cooling ============================================================
	contFlag = True
	iterTotal = 0
	iterNoImp = 0
	iterAcc = 0
	apRate = 1
	reportTime = 0
	ofvCurve = []
	while (contFlag):
		# Repeat in the same temperature
		for l in range(L):
			# Export time
			currTime = (datetime.datetime.now() - startTime).total_seconds()
			if (currTime > reportTime):
				print('t: %ss \t T: %s \t iter %s \t noImp %s \t apRate %s \t ofv: %s' % (
					round(currTime, 2), 
					round(T, 2), 
					iterTotal,
					iterNoImp,
					round(apRate, 3),
					ofv))
				reportTime += 10
				ofvCurve.append(ofv)

			# Increment iterator
			iterTotal += 1

			# Generate a neighbor using different type
			typeOfNeigh = randomPick(list(neighRatio))
			newSeq = None
			deltaC = None
			res = None
			if (typeOfNeigh == 0):
				res = swap(curSeq)
			elif (typeOfNeigh == 1):
				res = exchange(curSeq)
			elif (typeOfNeigh == 2):
				res = rotate(curSeq)
			elif (typeOfNeigh == 3):
				res = insert(curSeq)
			newSeq = res['seq']
			deltaC = res['deltaC']

			# If this new neighbor is good, accept it, 
			#     otherwise accept it with probability
			if (deltaC <= 0): # deltaC = newC - preC, <0 means improve
				curSeq = [i for i in newSeq]				
				ofv += deltaC
				iterAcc += 1
				iterNoImp = 0
			else:
				sample = random.random()
				if (sample < math.exp(- deltaC / T)):
					curSeq = [i for i in newSeq]
					ofv += deltaC
					iterAcc += 1
				else:
					iterNoImp += 1
			apRate = iterAcc / iterTotal

			# Check stopping criteria
			endCriteria = None
			if ('Final_Temperature' in stopType):
				if (T < stopTemp):
					contFlag = False
					endCriteria = 'Final_Temperature'
					break
			if ('Num_Iterations_Without_Improving' in stopType):
				if (iterNoImp > stopNoImp):
					contFlag = False
					endCriteria = 'Num_Iterations_Without_Improving'
					break
			if ('Percent_of_Accepted_Move' in stopType):
				if (iterTotal > 0 and apRate < stopAptRate):
					contFlag = False
					endCriteria = 'Percent_of_Accepted_Move'
					break
			if ('Num_Iterations' in stopType):
				if (iterTotal > stopIter):
					contFlag = False
					endCriteria = 'Num_Iterations'
					break
			if ('Executed_Time' in stopType):
				if ((datetime.datetime.now() - startTime).total_seconds() > stopTime):
					contFlag = False
					endCriteria = 'Executed_Time'
					break
		
		# Cool down
		T = coolRate * T

	curSeq.append(curSeq[0])
	runtime = (datetime.datetime.now() - startTime).total_seconds()
	return {
		'seq': curSeq,
		'ofv': ofv,
		'runtime': runtime,
		'endCriteria': endCriteria,
		'ofvCurve': ofvCurve
	}

