import geopy
import math
import random
import tripy

from .common import *
from .const import *
from .geometry import *
from .msg import *
from .relation import *

def rndPlainNodes(
    N:          "Number of vertices" = None,
    nodeIDs:    "A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
    distr:      "Spatial distribution of nodes, the options are \
                 1) String, (default) 'uniformSquare', or \
                 2) String, 'uniformCircleXY', or \
                 3) String, 'uniformCircleLatLon', or \
                 4) String, 'uniformPoly', or\
                 5) String, 'uniformRoadNetwork', or\
                 6) String, 'clustered'" = 'uniformSquare',
    distrArgs:  "Dictionary that describes the distribution of the nodes\
                 1) for 'uniformSquare'\
                    {\
                        'xRange': A 2-tuple with minimum/maximum range of x, default as (0, 100), \
                        'yRange': A 2-tuple with minimum/maximum range of y, default as (0, 100), \
                    }\
                 2) for 'uniformCircle'\
                    {\
                        'centerLoc': centering location \
                        'radius': radius of the circle \
                    }\
                 3) for 'uniformPoly'\
                    {\
                        'poly': polygon of the area, (no holes)\
                        (or 'polys': list of polygons) \
                    }\
                 4) for 'uniformRoadNetwork'\
                    {\
                        'network': list of arcs that can be sampled \
                        'poly': nodes should generated within the polygon, if not provided, will consider the entire network \
                    }\
                 5) for 'clustered\
                    {\
                        'numCluster': number of cluster centers\
                        'xRange': xRange of cluster centroid\
                        'yRange': yRange of cluster centroid\
                        'poly': polygon for customers\
                        'centroidLocs': list of cluster center locations\
                        'clusterDiameter': the spread of nodes, in diameter\
                    }\
                " = {
                    'xRange': (0, 100),
                    'yRange': (0, 100)
                 }
    ) -> "A set of nodes with id start from 0 to N":

    # Check for required fields ===============================================
    if (N == None and nodeIDs == None):
        print(ERROR_MISSING_N)
        return

    # Initialize ==============================================================
    nodes = {}
    if (nodeIDs == None):
        nodeIDs = [i for i in range(N)]

    # Generate instance =======================================================
    if (distr == 'uniformSquare'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            distrArgs = {
                'xRange': (0, 100),
                'yRange': (0, 100)
            }
        if ('xRange' not in distrArgs):
            distrArgs['xRange'] = (0, 100)
        if ('yRange' not in distrArgs):
            distrArgs['yRange'] = (0, 100)
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformSquare(distrArgs['xRange'], distrArgs['yRange'])
            }

    elif (distr == 'uniformCircleXY'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            distrArgs = {
                'centerLoc': [0, 0],
                'radius': 100
            }
        if ('centerLoc' not in distrArgs or 'radius' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_UNICC)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformCircleXY(distrArgs['radius'], distrArgs['centerLoc'])
            }

    elif (distr == 'uniformCircleLatLon'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformCircleLatLon(distrArgs['radius'], distrArgs['centerLoc'])
            }

    elif (distr == 'uniformPoly'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return
        if ('poly' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_UNIPOLY)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformPoly(distrArgs['poly'])
            }

    elif (distr == 'uniformPolys'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return
        if ('polys' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_UNIPOLY)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformPolys(distrArgs['polys'])
            }

    elif (distr == 'uniformRoadNetwork'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return
        if ('roadNetwork' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_UNIROADNETWORK)
            return

        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtRoadNetwork(
                    distrArgs['roadNetwork'], 
                    distrArgs['poly'] if 'poly' in distrArgs else None)
            }

    elif (distr == 'clusterXY'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return
        if ('centroidLocs' not in distrArgs or 'clusterDiameter' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_CLUSTER)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtClusterXY(distrArgs['centroidLocs'], distrArgs['clusterDiameter'])
            }

    elif (distr == 'clusterLatLon'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return
        if ('centroidLocs' not in distrArgs or 'clusterDiameterInMeters' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_CLUSTER)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtClusterLatLon(distrArgs['centroidLocs'], distrArgs['clusterDiameterInMeters'])
            }
    else:
        print(ERROR_INCOR_DISTARG)
        return

    return nodes

def rndTimeWindowsNodes(
    plainNodes: "Needs to provide plain nodes" = None,
    timeSpan:   "Time span for generating time windows, starting from 0" = None,
    twType:     "Type of time windows\
                 1) String 'Known', all time windows are given priorly, or \
                 2) String 'Periodic', time windows are generated periodically, or \
                 3) String 'Random', time windows are generated randomly " = 'Random',
    twArgs:     "Dictionary, \
                 1) For 'Known' twType:\
                    {\
                        'timeWindows': A list of 2-tuples, e.g., [(s1, e1), (s2, e2), ...]\
                    }\
                 2) For 'Periodic' twType:\
                    {\
                        'cycleTime': cycleTime, \
                        'startTimeInCycle': startTimeInCycle, \
                        'endTimeInCycle': endTimeInCycle \
                    }\
                 3) For 'Random' twType:\
                    {\
                        'indFlag': True if different nodes have different time windows, False otherwise,\
                        'avgBtwDur': average interval time between time windows, exponential distributed,\
                        'avgTWDur': average length of time windows, exponential distributed\
                    }" = None
    ) -> "A set of nodes with time windows":
    # Initialize ==============================================================
    nodes = copy.deepcopy(plainNodes)

    # Check for required fields ===============================================
    if (plainNodes == None):
        print(ERROR_MISSING_PLAINNODES)
        return
    if (timeSpan == None):
        print(ERROR_MISSING_TIMESPAN)
        return
    validTWOpts = ['Known', 'Periodic', 'Random']
    if (twType not in validTWOpts):
        print(ERROR_OPTS_TWTYPE % validTWOpts)
        return
    elif (twType == 'Known'):
        if (twArgs == None):
            print(ERROR_MISSING_TWKNOWNARGS)
            return
        elif ('timeWindows' not in twArgs):
            print (ERROR_MISSING_TWKNOWNARGS_TW)
            return
    elif (twType == 'Periodic'):
        if (twArgs == None):
            print(ERROR_MISSING_TWPERIODICARGS)
            return
        elif ('cycleTime' not in twArgs):
            print(ERROR_MISSING_TWPERIODICARGS_CT)
            return
        elif ('startTimeInCycle' not in twArgs):
            print(ERROR_MISSING_TWPERIODICARGS_ST)
            return
        elif ('endTimeInCycle' not in twArgs):
            print(ERROR_MISSING_TWPERIODICARGS_ET)
            return
    elif (twType == 'Random'):
        if (twArgs == None):
            print(ERROR_MISSING_TWRANDOMARGS)
            return
        elif ('indFlag' not in twArgs):
            print(ERROR_MISSING_TWRANDOMARGS_FG)
            return
        elif ('avgBtwDur' not in twArgs):
            print(ERROR_MISSING_TWRANDOMARGS_BTW)
            return
        elif ('avgTWDur' not in twArgs):
            print(ERROR_MISSING_TWRANDOMARGS_TW)
            return

    # Generate time windows for different types ===============================
    for n in nodes:
        nodes[n]['timeWindows'] = []
    if (twType == "Known"):
        for n in nodes:
            nodes[n]['timeWindows'] = twArgs['timeWindows']
    elif (twType == 'Periodic'):
        repeatNum = math.ceil(timeSpan / twArgs['cycleTime'])
        for tw in range(repeatNum):
            for n in nodes:
                nodes[n]['timeWindows'].append((
                    tw * twArgs['cycleTime'] + twArgs['startTimeInCycle'], 
                    tw * twArgs['cycleTime'] + twArgs['endTimeInCycle']))
    elif (twType == 'Random'):
        avgBtwDur = twArgs['avgBtwDur']
        avgTWDur = twArgs['avgTWDur']
        # If all nodes have different time windows
        if (twArgs['indFlag'] == False):
            now = 0
            pre = 0
            availFlag = True if (random.random() < (avgTWDur / (avgTWDur + avgBtwDur))) else False
            while (now < timeSpan):            
                interval = 0
                if (availFlag):
                    if (avgTWDur > 0):
                        interval = random.expovariate(1 / avgTWDur)
                        now += interval
                        if (interval > 0.0001):
                            for n in nodes:
                                nodes[n]['timeWindows'].append((pre, now))
                else:
                    if (avgBtwDur > 0):
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
                        if (avgTWDur > 0):
                            interval = random.expovariate(1 / avgTWDur)
                            now += interval
                            if (interval > 0.0001):
                                nodes[n]['timeWindows'].append((pre, now))
                    else:
                        if (avgBtwDur > 0):
                            interval = random.expovariate(1 / avgBtwDur)
                            now += interval
                    availFlag = not availFlag
                    pre = now

    # Truncate last time window to fit time span ==============================
    for n in nodes:
        while (len(nodes[n]['timeWindows']) > 0 and (nodes[n]['timeWindows'][-1][0] >= timeSpan or nodes[n]['timeWindows'][-1][1] > timeSpan)):
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

def rndWinds(
    degAvg:     "Prevailing direction of wind" = None,
    degSpread:  "Wind direction difference in 10th to 90th percentile bonds" = None,
    degPattern: "1) String, 'Consistent' or, \
                 2) String, 'Normal' or, \
                 3) String, 'Swap'" = None,
    spdAvg:     "Average wind speed" = None,
    spdSpread:  "Wind speed difference in 10th to 90th percentile bonds" = None,
    spdPattern: "1) String, 'Consistent' or, \
                 2) String, 'Increasing'  or, \
                 3) String, 'Decreasing' or, \
                 4) String, 'Normal' or" = None,
    duration:   "Time window of creating wind dictionaries" = [0, None],
    sampleInt:  "Interval between sampling in [sec]" = 3600
    ) -> "Generate a list of wind dictionaries":

    # Initialize ==============================================================
    winds = []
    now = duration[0]
    while (now <= duration[1]):
        winds.append({
            'startTime': now,
            'endTime': now + sampleInt,
            'windSpd': None,
            'windDeg': None
        })
        now += sampleInt

    # Direction ===============================================================
    for w in range(len(winds)):
        deg = None
        if (degPattern == 'Consistent'):
            deg = degAvg
        elif (degPattern == 'Normal'):
            deg = random.normalvariate(degAvg, degSpread / 1.28) # 10th to 90th percentile bonds
        elif (degPattern == 'Swap'):
            deg = degAvg - degSpread + w * (2 * degSpread / len(winds))
        winds[w]['windDeg'] = deg
    
    # Speed ===================================================================
    for w in range(len(winds)):
        spd = None
        if (spdPattern == 'Consistent'):
            spd = spdAvg
        elif (spdPattern == 'Increasing'):
            spd = spdAvg - spdSpread + w * (2 * spdSpread / len(winds))

        elif (spdPattern == 'Decreasing'):
            spd = spdAvg + spdSpread - (w + 1) * (2 * spdSpread / len(winds))

        elif (spdPattern == 'Normal'):
            spd = random.normalvariate(spdAvg, spdSpread / 1.28)
            while (spd <= 0):
                spd = random.normalvariate(spdAvg, spdSpread / 1.28)
        winds[w]['windSpd'] = spd

    return winds

def rndClouds(
    polyLatLon: "The polygon that clouds will be covering" = None,
    cloudMode:  "1) String, 'Cumulus', or\
                 2) String, 'Cumulonimbus'" = 'Cumulus',
    cloudDeg:   "The direction of where the clouds are entering" = [CLOUD_SPD_DEG_RANGE[0], CLOUD_SPD_DEG_RANGE[1]],
    numOfClouds: "Total number of clouds" = None,
    duration:   "Duration of the instance, in [sec]" = [0, 36000],
    ) -> "Given a polygon, returns a list of clouds that move through the polygon":

    # Time stamps =============================================================
    appearTime = []
    for i in range(numOfClouds):
        appearTime.append(random.uniform(duration[0], duration[1]))
    appearTime.sort()

    # Randomize the centroid of the cloud =====================================
    cloudCentroid = rndPlainNodes(
        N = len(appearTime),
        nodeIDs = [i for i in range(len(appearTime))],
        distr = 'uniformPoly',
        distrArgs = {'poly': polyLatLon})
    
    clouds = []
    for t in range(len(appearTime)):
        # Entering locations --------------------------------------------------
        enteringDist = None
        if (cloudMode == 'Cumulus'):
            enteringDist = random.randrange(CLOUD_CUMULUS_SIZE[0], CLOUD_CUMULUS_SIZE[1])
        elif (cloudMode == 'Cumulonimbus'):
            enteringDist = random.randrange(CLOUD_CUMULONIMBUS_SIZE[0], CLOUD_CUMULONIMBUS_SIZE[1])

        # Initial location of the cloud ---------------------------------------
        direction = random.randrange(cloudDeg[0], cloudDeg[1])
        enterLoc = pointInDistLatLon(
            pt = cloudCentroid[t]['loc'],
            direction = direction - 180,
            distMeters = enteringDist + CONST_EPSILON)
        
        # Size of cloud -------------------------------------------------------
        widthInMeter = None
        lengthInMeter = None
        if (cloudMode == 'Cumulus'):
            widthInMeter = random.randrange(CLOUD_CUMULUS_SIZE[0], CLOUD_CUMULUS_SIZE[1])
            lengthInMeter = random.randrange(CLOUD_CUMULUS_SIZE[0], CLOUD_CUMULUS_SIZE[1])
            existDur = random.randrange(CLOUD_CUMULUS_DURATION[0], CLOUD_CUMULUS_DURATION[1])
        elif (cloudMode == 'Cumulonimbus'):
            widthInMeter = random.randrange(CLOUD_CUMULONIMBUS_SIZE[0], CLOUD_CUMULONIMBUS_SIZE[1])
            lengthInMeter = random.randrange(CLOUD_CUMULONIMBUS_SIZE[0], CLOUD_CUMULONIMBUS_SIZE[1])
            existDur = random.randrange(CLOUD_CUMULONIMBUS_DURATION[0], CLOUD_CUMULONIMBUS_DURATION[1])

        # Cloud generation ----------------------------------------------------
        initPolyLatLon = rectInWidthLengthOrientationLatLon(
            centroidLatLon = enterLoc,
            widthInMeter = widthInMeter,
            lengthInMeter = lengthInMeter,
            oriDeg = random.randrange(-45, 45))        
        clouds.append({
            'appearTime': appearTime[t],
            'initPolyLatLon': initPolyLatLon,
            'movingVecPolar': (random.randrange(CLOUD_SPD_ABS_RANGE[0], CLOUD_SPD_ABS_RANGE[1]), direction),
            'existDur': existDur
        })

    return clouds
