from .const import *
from .geometry import *

def stepFromTWs(
    res:        "Name of the resource" = None,
    tws:        "A list of time windows that one of the resource is occupied, could (and highly likely) be overlapped" = None,
    share:      "number of user sharing one item" = 1,
    stepInt:    "Interval of stat" = 1,
    ) -> "Given a set of occupied time windows of given resource, returns a step sequence represent utilization": 
    
    # Initialize ==============================================================
    timeStamp = []
    useLevel = []
    occLevel = []

    # Get the timeStamps ======================================================
    for tw in tws:
        if (tw[0] not in timeStamp):
            timeStamp.append(tw[0])
        if (tw[1] not in timeStamp):
            timeStamp.append(tw[1])
    timeStamp.sort()
    for i in range(len(timeStamp)):
        useLevel.append(0)
        occLevel.append(0)

    # Update useLevel =========================================================
    for tw in tws:
        # Use level between start time and end time will increase by 1
        for i in range(len(timeStamp)):
            if (timeStamp[i] >= tw[0] and timeStamp[i] < tw[1]):
                useLevel[i] += 1
                occLevel[i] = math.ceil(useLevel[i] / share)

    # Get stat ================================================================
    rangeByEntity = [min(useLevel), max(useLevel)]
    stat = []
    for r in range(0, math.ceil(rangeByEntity[1] - rangeByEntity[0]), stepInt):
        rangeStat = [rangeByEntity[0] + r * stepInt, rangeByEntity[0] + (r + 1) * stepInt]
        stat.append({
            'rangeStat': rangeStat,
            'totalLength': 0
        })
    for i in range(len(timeStamp) - 1):
        lengthOfLevel = timeStamp[i + 1] - timeStamp[i]
        for j in range(len(stat)):
            if (useLevel[i] > stat[j]['rangeStat'][0] and useLevel[i] <= stat[j]['rangeStat'][1]):
                stat[j]['totalLength'] += lengthOfLevel

    return {
        'resID': res,
        'timeStamp': timeStamp,
        'useLevel': occLevel,
        'stat': stat
    }

def ganttFromVRPSol(
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

    vrpSol:     "Dictionary, solution of VRP, by heuVRP(), heuVRPTW(), heuVRPPD(), heuVRPD(), etc." = None
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
        gantt.extend(ganttFromVisitSeq(
            entityID = 'Truck_%s' % veh,
            edges = tau,
            visitSeq = vrpSol['route'][veh]['route'],
            serviceTime = vrpSol['serviceTime'])['gantt'])

    return {
        'gantt': gantt
    }

def ganttFromTSPSol(
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

    tspSol:     "Dictionary, solution of TSP, by ipTSP() or heuTSP()" = None
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
    gantt = ganttFromVisitSeq(
        entityID = 'Truck',
        edges = tau,
        visitSeq = tspSol['seq'],
        serviceTime = tspSol['serviceTime'])['gantt']

    return {
        'gantt': gantt
    }

def ganttFromVisitSeq(
    entityID:   "Truck" = 'Truck',
    edges:        "Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    visitSeq:    "Travel sequence" = None,
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    ) -> "Given a TSP result, returns the gantt dictionary for plotGantt()":

    # Process TSP sequence ====================================================
    gantt = []
    acc = 0
    for i in range(len(visitSeq) - 2):
        dt = edges[visitSeq[i], visitSeq[i + 1]]
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
    dt = edges[visitSeq[-2], visitSeq[-1]]
    gantt.append({
        'entityID': entityID,
        'timeWindow': [acc, acc + dt],
        'desc': 'Return',
        'color': 'lightgreen',
        'style': '///'
    })

    return {
        'gantt': gantt
    }
