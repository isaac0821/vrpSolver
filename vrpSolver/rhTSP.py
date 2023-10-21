import datetime

from .heuTSP import *
from .geometry import *

# History =====================================================================
# 20231022 - Start writing script
# =============================================================================
def rhTSP(
    rhTime: tuple[int, int],
    # Environment settings ----------------------------------------------------
    nodes: dict = None, 
    locFieldName: str = 'loc',
    depotID: int|str = 0, 
    nodeIDs: list[int|str]|str = 'All', 
    serviceTime: float = 0,
    vehicles: dict = {
        0: {'speed': 1}
    },
    vehicleID: int|str = 0,
    edges: dict = {
        'method': 'Euclidean', 
        'ratio': 1
    }, 
    # Solver engine settigs ---------------------------------------------------
    method: dict = {
        'cons': 'Insertion', 
        'impv': '2Opt'
    },
    reoptPolicy: dict = {
        'method': 'FixedTime',
        'timeInterval': 100
    },
    # Output configs ----------------------------------------------------------
    detailsFlag: bool = False,
    metaFlag: bool = False,
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

    if (reoptPolicy == None or 'method' not in reoptPolicy):
        raise MissingParameterError("ERROR: Missing required field `reoptPolicy`.")

    # Initialize phase ========================================================
    # 初始时钟
    clock = 0
    # 初始状态
    curVehTimedSeq = [(nodes[depotID][locFieldName], 0), (nodes[depotID][locFieldName], rhTime)]
    # 初始待访问集合
    notVisitIDs = [n for n in nodes if n not in visit and nodes[n]['timeWindow'][0] <= clock]
    # 所有可能节点间的距离+路径
    tau, path = matrixDist(
        nodes = nodes, 
        edges = edges, 
        depotID = depotID, 
        nodeIDs = nodeIDs)
    nodes['snap'] = {'loc': nodes[depotID][locFieldName]}

    # Get time stamps =========================================================
    arrivalTimes = [nodes[i]['timeWindow'][0] for i in nodes]
    fixedTimeList = []
    if (reoptPolicy['method'] == 'FixedTime'):
        fixedTimeList = [i * reoptPolicy['timeInterval'] for i in range(int(rhTime / reoptPolicy['timeInterval']))]
        arrivalTimes.extend(fixedTimeList)
    timeList = list(set(arrivalTimes))
    timeList.sort()
    timeStampInList = 0

    # Main loop ===============================================================
    reoptFlag = True
    reoptPrev = 0
    while (clock <= rhTime and timeStampInList < len(timeList) - 1):
        # Initialize: current status ------------------------------------------
        # 在当前时刻vehicle的位置
        snap = snapInTimedSeq(curVehTimedSeq, clock)
        nodes['snap'] = {'loc': snap['loc']}
        # 在当前时刻尚未被访问的位置
        unvisitIDs = [i for i in nodes if nodes[i]['timeWindow'][0] <= clock and nodes[i]['status'] != 'Visited']
        
        # update tau and path
        scaleTau, _, scalePathLoc, _ = scaleDist(
            loc = snap['loc'], 
            nodes = nodes, 
            nodeIDs = unvisitIDs,
            edges = edges,
            locFieldName = locFieldName)
        for n in unvisitIDs:
            tau['snap', n] = scaleTau[n]
            tau[n, 'snap'] = CONST_EPSILON
            path['snap', n] = [snap['loc'], nodes[n][locFieldName]]
            path[n, 'snap'] = []

        # Check if reoptimize at this time ------------------------------------
        reoptFlag = False
        if (reoptPolicy['method'] == 'FixedTime' and clock - reoptPrev >= reoptPolicy['timeInterval']):
            reoptFlag = True
        elif (reoptPolicy['method'] == 'AccuArr' and unvisitIDs >= reoptPolicy['accuLevel']):
            reoptFlag = True

        # Find TSP route for unvisited nodes ----------------------------------
        if (reoptFlag):
            if (len(unvisitIDs) > 0):
                reoptPrev = clock
                reoptTsp = heuTSP(
                    nodes = nodes,
                    edges = {'method': 'Dictionary', 'tau': tau, 'path': path},
                    depotID = 'snap',
                    nodeIDs = unvisitIDs,
                    detailsFlag = True)
                curVehTimedSeq = tsp['vehicles']['timedSeq']
            else:
                # Stay at the same location
                curVehTimedSeq = [(snap['loc'], clock), (snap['loc'], rhTime)]

        # Post optimization ---------------------------------------------------
        

        # Clock advance -------------------------------------------------------
        timeStampInList += 1
        clock = timeList[timeStampInList]

    return {
        'ofv': ofv,
        'nodes': nodes,
        'vehicles': vehicles,
        'reoptPolicy': reoptPolicy
    }