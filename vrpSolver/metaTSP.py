import datetime
import math

from .const import *
from .common import *
from .graph import *
from .geometry import *
from .msg import *
from .heuTSP import *

def metaTSP(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...} or \
                 4) String 'Grid', will need to add arguments using `edgeArgs`" = "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                    {\
                        'colRow': (numCol, numRow),\
                        'barriers': [(coordX, coordY), ...], \
                    }" = None,
    depotID:    "DepotID, default to be 0" = 0,
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    metaAlgo:   "1) String (default) 'SimulatedAnnealing' or, \
                 2) String (not available) 'GenaticAlgorithm', or\
                 3) String (not available) 'TabuSearch'" = 'SimulatedAnnealing',
    metaAlgoArgs: "Dictionary, args for metaheuristic\
                 1) for 'SimulatedAnnealing'\
                    {\
                        'initAlgo': acceptable constructive heuristic option for heuTSP(),\
                        'initAlgoArgs': usually not needed, use when 'initAlgo' needs additional args in heuTSP()\
                        'initTemp': initial temperature,\
                        'iterTemp': number of iterations per temperature,\
                        'optRatio': ratio of randomize local search operations,\
                        'coolRate': cooling rate\
                        'stopCriteria': list of 2-tuple, options of the first element in tuple are as follows: \
                            1) 'Final_Temperature', \
                            2) 'Num_Iterations_Without_Improving', \
                            3) 'Percent_of_Accepted_Move', \
                            4) 'Num_Iterations', \
                            5) 'Executed_Time'  \
                            second element is the correspond value of stop criteria, e.g. [('Executed_Time', 100), ...]\
                    }" = None
    ) -> "Using metaheuristic for solving TSP tour":

    # Define nodeIDs ==========================================================
    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = []
            for i in nodes:
                nodeIDs.append(i)
        else:
            msgError(ERROR_INCOR_NODEIDS)
            return

    # Define tau ==============================================================
    tau = {}
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            tau = getTauEuclidean(nodes, nodeIDs)
        elif (edges == 'LatLon'):
            tau = getTauLatLon(nodes, nodeIDs)
        elif (edges == 'Grid'):
            if (edgeArgs == None or 'colRow' not in edgeArgs or 'barriers' not in edgeArgs):
                msgError("ERROR: Need more information to define the grid.")
                return None
            tau = getTauGrid(nodes, nodeIDs, edgeArgs['colRow'], edgeArgs['barriers'])
        else:
            msgError(ERROR_INCOR_TAU)
            return None
    else:
        tau = dict(edges)

    # Service time ============================================================
    if (serviceTime != None and serviceTime > 0):
        for (i, j) in tau:
            if (i != depotID and j != depotID and i != j):
                tau[i, j] += serviceTime
            elif (i == depotID or j == depotID and i != j):
                tau[i, j] += serviceTime / 2 
    else:
        serviceTime = 0

    # Configuration ===========================================================
    if (metaAlgo == 'SimulatedAnnealing'):
        if (metaAlgoArgs == None):
            metaAlgoArgs = {
                'initAlgo': 'Sweep', 
                'initTemp': 1000, 
                'iterTemp': 20,
                'optRatio': (0.2, 0.2, 0.3),
                'coolRate': 0.99,
                'stopCriteria': [('Executed_Time', 10)]
            }
            msgWarning("WARNING: Missing `metaAlgoArgs`. Using default setting as follows, results may not be satisfactory.")
            msgWarning(metaAlgoArgs)
        else:
            canGo = True
            warningGo = False
            if ('initAlgo' not in metaAlgoArgs):
                msgWarning("WARNING: Missing 'initAlgo' in `metaAlgoArgs`. Initial solution set to be 'Sweep'")
                metaAlgoArgs['initAlgo'] = 'Sweep'
                warningGo = True
            if ('initTemp' not in metaAlgoArgs):
                msgError("ERROR: Missing 'initTemp' in `metaAlgoArgs`. Require initial temperature")
                canGo = False
            if ('iterTemp' not in metaAlgoArgs):
                msgError("ERROR: Missing 'iterTemp' in `metaAlgoArgs`. Require number of iterations per temperature")
                canGo = False
            if ('optRatio' not in metaAlgoArgs): 
                msgWarning("WARNING: Missing 'optRatio' in `metaAlgoArgs`. 'optRatio' is the percentage of taking (swap, exchange, 2Opt) operators")
                metaAlgoArgs['optRatio'] = (0.3, 0.3, 0.4)
                warningGo = True
            if ('coolRate' not in metaAlgoArgs):
                msgError("ERROR: Missing 'coolRate' in `metaAlgoArgs`. Require cooling rate, e.g., 0.99")
                canGo = False
            if ('stopCriteria' not in metaAlgoArgs):
                msgWarning("WARNING: Missing 'stopCriteria' in `metaheuristic`. Set to be default as 'Executed_Time' for 10 seconds")
                metaAlgoArgs['stopCriteria'] = [('Executed_Time', 10)]
                warningGo = True
            if (not canGo):
                return None
            if (warningGo):
                msgWarning("WARNING: Adjusted `metaAlgoArgs` is as follows:")
                msgWarning(metaAlgoArgs)

    # Subroutines to generate neighborhoods ===================================
    # Swap two nearby vertices
    def _swap(seq):
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
    def _exchange(seq):
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
    def _2Opt(seq):
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

    # Subroutines for different metaheuristic =================================
    def _metaTSPSimulatedAnnealing(initAlgo, initAlgoArgs, initTemp, iterTemp, optRatio, coolRate, stopCriteria):
        # Initialize
        # Initial temperature
        T = initTemp
        # Temperature length (maximum temperature iteration)
        L = iterTemp
        # Initial Solution
        initSol = heuTSP(
            nodes = nodes, 
            edges = edges, 
            edgeArgs = edgeArgs, 
            depotID = depotID, 
            nodeIDs = nodeIDs, 
            serviceTime = serviceTime, 
            consAlgo = initAlgo, 
            consAlgoArgs = initAlgoArgs)
        # To avoid all kind of trouble, seq here is not closed
        curSeq = initSol['seq'][:-1] 
        ofv = initSol['ofv'] - tau[curSeq[-1], depotID]

        # Main cooling --------------------------------------------------------
        startTime = datetime.datetime.now()
        contFlag = True
        iterTotal = 0
        iterNoImp = 0
        iterAcc = 0
        apRate = 1
        while (contFlag):
            # Repeat in the same temperature
            for l in range(L):
                # Increment iterator
                iterTotal += 1

                # Generate a neighbor using different type
                typeOfNeigh = rndPick(list(optRatio))
                newSeq = None
                deltaC = None
                res = None
                if (typeOfNeigh == 0):
                    res = _swap(curSeq)
                elif (typeOfNeigh == 1):
                    res = _exchange(curSeq)
                elif (typeOfNeigh == 2):
                    res = _2Opt(curSeq)
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
                    if (sample < math.exp(-deltaC / T)):
                        curSeq = [i for i in newSeq]
                        ofv += deltaC
                        iterAcc += 1
                    else:
                        iterNoImp += 1
                apRate = iterAcc / iterTotal

                # Check stopping criteria
                endCriteria = None
                for end in stopCriteria:
                    if (end[0] == 'Final_Temperature'):
                        if (T < end[1]):
                            contFlag = False
                            endCriteria = 'Final_Temperature'
                            break
                    if (end[0] == 'Num_Iterations_Without_Improving'):
                        if (iterNoImp > end[1]):
                            contFlag = False
                            endCriteria = 'Num_Iterations_Without_Improving'
                            break
                    if (end[0] == 'Percent_of_Accepted_Move'):
                        if (iterTotal > 0 and apRate < end[1]):
                            contFlag = False
                            endCriteria = 'Percent_of_Accepted_Move'
                            break
                    if (end[0] == 'Num_Iterations'):
                        if (iterTotal > end[1]):
                            contFlag = False
                            endCriteria = 'Num_Iterations'
                            break
                    if (end[0] == 'Executed_Time'):
                        if ((datetime.datetime.now() - startTime).total_seconds() > end[1]):
                            contFlag = False
                            endCriteria = 'Executed_Time'
                            break
            
            # Cool down
            T = coolRate * T
        curSeq.append(curSeq[0])
        metaStat = {
            'endCriteria': endCriteria,
            'temperature': T,
            'iterTotal': iterTotal,
            'iterNoImp': iterNoImp,
            'acceptRate': apRate,
            'runtime': (datetime.datetime.now() - startTime).total_seconds()
        }

        return {
            'seq': curSeq,
            'metaStat': metaStat
        }

    meta = None
    if (metaAlgo == 'SimulatedAnnealing'):
        meta = _metaTSPSimulatedAnnealing(
            initAlgo = metaAlgoArgs['initAlgo'],
            initAlgoArgs = metaAlgoArgs['initAlgoArgs'] if 'initAlgoArgs' in metaAlgoArgs else None,
            initTemp = metaAlgoArgs['initTemp'],
            iterTemp = metaAlgoArgs['iterTemp'],
            optRatio = metaAlgoArgs['optRatio'],
            coolRate = metaAlgoArgs['coolRate'],
            stopCriteria = metaAlgoArgs['stopCriteria'])

    # Fix the sequence to make it start from and end with the depot ===========
    # NOTE: nodeID gets duplicated, if nodeID == 0, the sequence starts and ends with a 0
    startIndex = None
    truckSeq = []
    seq = meta['seq']
    for k in range(len(seq)):
        if (seq[k] == depotID):
            startIndex = k
    if (startIndex <= len(seq) - 1):
        for k in range(startIndex, len(seq)):
            truckSeq.append(seq[k])
    if (startIndex >= 0):
        for k in range(0, startIndex):
            truckSeq.append(seq[k])
    truckSeq.append(depotID)

    ofv = calSeqCostMatrix(tau, truckSeq)

    return {
        'ofv': ofv,
        'seq': truckSeq,
        'metaStat': meta['metaStat'],
        'serviceTime': serviceTime
    }