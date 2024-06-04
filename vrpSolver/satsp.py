import heapq
import math
import warnings
import networkx as nx

from .common import *
from .geometry import *
from .msg import *
from .ds import *

def metaTSP(
    nodes: dict, 
    locFieldName: str = 'loc',
    depotID: int|str = 0,
    nodeIDs: list[int|str]|str = 'All',
    serviceTime: float = 0,
    vehicles: dict = {
        0: {'speed': 1}
    },
    vehicleID: int|str = 0,
    edges: dict = {
        'method': "Euclidean", 
        'ratio': 1
    },
    method: dict = {
        'algo': 'SimulatedAnnealing',
        'cons': 'Insertion',
        'initTemp': None,
        'lengTemp': None,
        'neighRatio': {
            '2Opt': 0.6,
            'swap': 0.2,
            'removeInsert': 0.2
        },
        'coolRate': None,
        'stop': {
            'stopTemp': None,
            'stopNoImp': None,
            'stopAptRate': None,
            'stopIter': None,
            'stopRuntime': None
        }
    },
    detailsFlag: bool = False,
    metaFlag: bool = False
    ) -> dict|None:
    
    # Sanity check ============================================================
    if (nodes == None or type(nodes) != dict):
        raise MissingParameterError(ERROR_MISSING_NODES)

    if (type(nodeIDs) is not list):
        if (nodeIDs == 'All'):
            nodeIDs = [i for i in nodes]
        else:
            for i in nodeIDs:
                if (i not in nodes):
                    raise OutOfRangeError("ERROR: Node %s is not in `nodes`." % i)
    if ((type(nodeIDs) == list and depotID != None and depotID not in nodeIDs)
        or (nodeIDs == 'All' and depotID != None and depotID not in nodes)):
        raise OutOfRangeError("ERROR: Cannot find `depotID` in given `nodes`/`nodeIDs`")

    if (vehicles == None):
        raise MissingParameterError("ERROR: Missing required field `vehicles`.")
    if (vehicleID not in vehicles):
        raise MissingParameterError("ERROR: Cannot find `vehicleID` in `vehicles`.")

    # Define tau ==============================================================
    tau = None
    path = None
    if (detailsFlag):
        tau, path = matrixDist(nodes, edges, depotID, nodeIDs, locFieldName)
    else:
        tau, _ = matrixDist(nodes, edges, depotID, nodeIDs, locFieldName)

    # Check symmetric =========================================================
    asymFlag = False
    for (i, j) in tau:
        if (tau[i, j] != tau[j, i]):
            asymFlag = True
            break

    return

def _metaTSPSimulatedAnnealing(nodeLoc, tau, nodeIDs, initSol, initTemp, lengTemp, neighRatio, coolRate, stop) -> dict:

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
        res = consTSP(nodeLoc, tau, initSol)
        curSeq = res['seq'][:-1] # To avoid all kind of trouble, seq here is not closed
        ofv = res['ofv']
    else:
        return None    

    # Main cooling ============================================================
    contFlag = True
    iterTotal = 0
    iterNoImp = 0
    iterAcc = 0
    reportTime = 0
    ofvCurve = []
    while (contFlag):
        # Repeat in the same temperature
        for l in range(L):
            # Export time
            currTime = (datetime.datetime.now() - startTime).total_seconds()
            if (currTime > reportTime):
                print('t: %ss \t T: %s \t iter %s \t noImp %s \t ofv: %s' % (
                    round(currTime, 2), 
                    round(T, 2), 
                    iterTotal,
                    iterNoImp,
                    ofv))
                reportTime += 10
                ofvCurve.append(ofv)

            # Increment iterator
            iterTotal += 1

            # Generate a neighbor using different type
            typeOfNeigh = rndPick(list(neighRatio))
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