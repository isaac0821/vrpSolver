from .const import *
from .geometry import *

def VRP2Gantt(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    nodeIDs:    "1) String (default) 'All', or \
                 2) A list of node IDs" = 'All',
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                {\
                    'colRow': (numCol, numRow),\
                    'barriers': [(coordX, coordY), ...], \
                }" = None,

    vrpSol:     "Dictionary, solution of VRP, by heuVRP(), heuVRPTW(), heuVRPPD(), heuVRPD(), etc." = None,
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    ) -> "Given a VRP result, returns the gantt dictionary for plotGantt()":

    # Define tau ==============================================================
    tau = {}
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            tau = getTauEuclidean(nodes, nodeIDs)
        elif (edges == 'LatLon'):
            tau = getTauLatLon(nodes, nodeIDs)
        elif (edges == 'Grid'):
            tau = getTauGrid(nodes, nodeIDs, edgeArgs['colRow'], edgeArgs['barriers'])
        else:
            msgError(ERROR_INCOR_TAU)
            return None
    else:
        tau = dict(edges)

    # Create gantt ============================================================
    gantt = []
    for veh in vrpSol['route']:
        gantt.extend(visitSeq2Gantt(
            entityID = 'Truck_%s' % veh,
            tau = tau,
            visitSeq = vrpSol['route'][veh],
            serviceTime = serviceTime))

    return gantt

def TSP2Gantt(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                {\
                    'colRow': (numCol, numRow),\
                    'barriers': [(coordX, coordY), ...], \
                }" = None,

    tspSol:     "Dictionary, solution of TSP, by ipTSP() or heuTSP()" = None,
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    ) -> "Given a TSP result, returns the gantt dictionary for plotGantt()":

    # Define tau ==============================================================
    tau = {}
    if (type(edges) is not dict):
        if (edges == 'Euclidean'):
            tau = getTauEuclidean(nodes, nodeIDs)
        elif (edges == 'LatLon'):
            tau = getTauLatLon(nodes, nodeIDs)
        elif (edges == 'Grid'):
            tau = getTauGrid(nodes, nodeIDs, edgeArgs['colRow'], edgeArgs['barriers'])
        else:
            msgError(ERROR_INCOR_TAU)
            return None
    else:
        tau = dict(edges)

    # Create gantt ============================================================
    gantt = visitSeq2Gantt(
        entityID = 'Truck',
        tau = tau,
        visitSeq = tspSol['seq'],
        serviceTime = serviceTime)

    return gantt

def visitSeq2Gantt(
    entityID:   "Truck" = 'Truck',
    tau:        "Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    visitSeq:    "Travel sequence" = None,
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    ) -> "Given a TSP result, returns the gantt dictionary for plotGantt()":

    # Process TSP sequence ====================================================
    gantt = []
    acc = 0
    for i in range(len(visitSeq) - 2):
        dt = tau[visitSeq[i], visitSeq[i + 1]]
        gantt.append({
            'entityID': entityID,
            'timeWindow': [acc, acc + dt],
            'desc': 'cus_' + str(visitSeq[i + 1]),
            'color': 'lightgreen',
            'style': '///'
        })
        gantt.append({
            'entityID': entityID,
            'timeWindow': [acc + dt, acc + dt + serviceTime],
            'desc': '',
            'color': 'lightgray',
            'style': '////'
        })
        acc += dt + serviceTime
    dt = tau[visitSeq[-2], visitSeq[-1]]
    gantt.append({
        'entityID': entityID,
        'timeWindow': [acc, acc + dt],
        'desc': 'Return',
        'color': 'lightgreen',
        'style': '///'
    })

    return gantt

def getStepFromTWs(
    res:        "Name of the resource" = None,
    tws:        "A list of time windows that one of the resource is occupied, could (and highly likely) be overlapped" = None,
    ) -> "Given a set of occupied time windows of given resource, returns a step sequence represent utilization": 
    
    # Initialize ==============================================================
    timeStamp = []
    useLevel = []

    # Get the timeStamps ======================================================
    for tw in tws:
        if (tw[0] not in timeStamp):
            timeStamp.append(tw[0])
        if (tw[1] not in timeStamp):
            timeStamp.append(tw[1])
    timeStamp.sort()
    for i in range(len(timeStamp)):
        useLevel.append(0)

    # Update useLevel =========================================================
    for tw in tws:
        # Use level between start time and end time will increase by 1
        for i in range(len(timeStamp)):
            if (timeStamp[i] >= tw[0] and timeStamp[i] < tw[1]):
                useLevel[i] += 1

    return {
        'resID': res,
        'timeStamp': timeStamp,
        'useLevel': useLevel
    }
