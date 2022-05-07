
from .geometry import *
from .ipTSP import *
from .heuTSP import *

def VRP2Gantt(

    ) -> "Given a VRP result, returns the gantt dictionary for plotGantt()":

    return

def TSP2Gantt(

    )

def visitSeq2Gantt(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    depotID:    "Depot ID (serviceTime will not apply)" = 0,
    entityID:   "Truck" = 'Truck',
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...}" = "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                {\
                    'colRow': (numCol, numRow),\
                    'barriers': [(coordX, coordY), ...], \
                }" = None,
    visitSeq:    "Travel sequence" = None,
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    ) -> "Given a TSP result, returns the gantt dictionary for plotGantt()":

    # Define tau ==============================================================
    tau = {}
    if (type(edges) is not dict):
        if (tau == 'Euclidean'):
            tau = getTauEuclidean(nodes, nodeIDs)
        elif (tau == 'LatLon'):
            tau = getTauLatLon(nodes, nodeIDs)
        elif (tau == 'Grid'):
            tau = getTauGrid(nodes, nodeIDs, edgeArgs['colRow'], edgeArgs['barriers'])
        else:
            msgError(ERROR_INCOR_TAU)
            return None
    else:
        tau = dict(edges)

    # Process TSP sequence ====================================================
    ganttTSP = []
    acc = 0
    for i in range(len(visitSeq) - 2):
        dt = tau[visitSeq[i], visitSeq[i + 1]]
        ganttTSP.append({
            'entityID': entityID,
            'timeWindow': [acc, acc + dt],
            'desc': 'cus_' + str(visitSeq[i + 1]),
            'color': 'lightgreen',
            'style': '///'
        })
        ganttTSP.append({
            'entityID': entityID,
            'timeWindow': [acc + dt, acc + dt + serviceTime],
            'desc': '',
            'color': 'lightgray',
            'style': '////'
        })
        acc += dt + serviceTime
    dt = tau[visitSeq[-2], visitSeq[-1]]
    ganttTSP.append({
        'entityID': entityID,
        'timeWindow': [acc, acc + dt],
        'desc': 'Return',
        'color': 'lightgreen',
        'style': '///'
    })

    return ganttTSP
