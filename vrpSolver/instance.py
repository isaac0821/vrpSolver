import random
import math
import tripy
import geopy

from .msg import *
from .geometry import *
from .common import *
from .matrices import *

def nodesEuc2Geo(
    nodesInEuc: "Node dictionary in euclidean space" = None,
    ratio:      "Distance ratio between euclidean traveling matrix (`edgeInEuc`) and geological distance (as in meters)" = 1,
    refNodeID:  "Reference nodeID" = None,
    refNodeLoc:  "Lat/lon of reference node" = None
    ) -> "Given nodes in euclidean space, transform it into lat/lon coordinates":

    # Initialize ==============================================================
    # NOTICE: For lat/lon, the x, y coordinate is in different order
    nodesInGeo = {}

    # Get the relative distance to each node ==================================
    for n in nodesInEuc:
        deg2RefNode = getHeadingXY(nodesInEuc[refNodeID]['loc'], nodesInEuc[n]['loc'])
        dist2RefNode = distEuclidean2D(nodesInEuc[refNodeID]['loc'], nodesInEuc[n]['loc'])
        geoLoc = pointInDistLatLon(refNodeLoc, deg2RefNode, dist2RefNode)
        (lat, lon) = (geoLoc[0], geoLoc[1])
        nodesInGeo[n] = {'loc': (lat, lon)}

    return nodesInGeo

def rndPlainNodes(
    N:          "Number of vertices" = None,
    nodeIDs:    "A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
    distr:      "Spatial distribution of nodes, the options are \
                 1) String, (default) 'uniformSquare', or \
                 2) String, 'uniformCircle', or \
                 2) String, 'uniformPoly', or\
                 3) String, 'uniformOnNetwork', or\
                 2) String, 'clustered, or \
                 3) String, 'ring', or \
                 4) String, 'normal2D'" = 'uniform',
    distrArgs:  "Dictionary that describes the distribution of the nodes\
                 1) for 'uniformSquare'\
                    {\
                        'xRange': A 2-tuple with minimum/maximum range of x, default as (0, 100), \
                        'yRange': A 2-tuple with minimum/maximum range of y, default as (0, 100), \
                    }\
                 2) for 'uniformCircle'\
                 2) for 'uniformPoly'\
                    {\
                        'poly': polygon of the area, (no holes)\
                    }\
                 3) for 'uniformOnNetwork'\
                 2) for 'clustered\
                 3) for 'ring'\
                    {\
                        'radius': radius of the ring,\
                        'degOffset': clock-wise rotate the nodes, default as 0, which is pointing north\
                        'centerLoc': centering location of the ring, default as (0, 0)\
                    }\
                 4) for 'normal2D'\
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
    if (distr == "uniformSquare"):
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
            x = random.randrange(distrArgs['xRange'][0], distrArgs['xRange'][1])
            y = random.randrange(distrArgs['yRange'][0], distrArgs['yRange'][1])
            nodes[n] = {'loc': (x, y)}
    elif (distr == "uniformCircle"):
        print("Stay tune")
        return
    elif (distr == "uniformPoly"):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return
        if ('poly' not in distrArgs and 'polys' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_UNIPOLY)
            return
        # Create nodes --------------------------------------------------------
        polys = []
        if ('poly' in distrArgs and 'polys' not in distrArgs):
            polys = [distrArgs['poly']]
        elif ('poly' in distrArgs and 'polys' in distrArgs):
            polys = [distrArgs['poly']]
            polys.extend(distrArgs['polys'])
        elif ('polys' in distrArgs):
            polys = distrArgs['polys']
        # Get all triangulated triangles
        lstTriangle = []
        
        for p in polys:
            lstTriangle.extend(tripy.earclip(p))
        # Weight them and make draws
        lstWeight = []
        for i in range(len(lstTriangle)):
            lstWeight.append(calTriangleAreaByCoords(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2]))
        for n in nodeIDs:
            idx = rndPick(lstWeight)
            [x1, y1] = lstTriangle[idx][0]
            [x2, y2] = lstTriangle[idx][1]
            [x3, y3] = lstTriangle[idx][2]
            rndR1 = np.random.uniform(0, 1)
            rndR2 = np.random.uniform(0, 1)
            rndX = (1 - math.sqrt(rndR1)) * x1 + math.sqrt(rndR1) * (1 - rndR2) * x2 + math.sqrt(rndR1) * rndR2 * x3
            rndY = (1 - math.sqrt(rndR1)) * y1 + math.sqrt(rndR1) * (1 - rndR2) * y2 + math.sqrt(rndR1) * rndR2 * y3
            nodes[n] = {'loc': (rndX, rndY)}
    elif (distr == "uniformOnNetwork"):
        print("Stay tune")
        return
    elif (distr == "clustered"):
        print("Stay tune")
        return
    elif (distr == "ring"):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return        
        if ('radius' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_RING)
            return
        if ('degOffset' not in distrArgs):
            distrArgs['degOffset'] = 0
        if ('centerLoc' not in distrArgs):
            distrArgs['centerLoc'] = (0, 0)
        # Create nodes --------------------------------------------------------
        initDeg = distrArgs['degOffset']
        deltaDeg = 360 / N
        for i in range(N):
            deg = initDeg + i * deltaDeg
            (x, y) = pointInDistXY(distrArgs['centerLoc'], deg, distrArgs['radius'])
            nodes[nodeIDs[i]] = {'loc': (x, y)}
    elif (distr == "normal2D"):
        print("Stay tune")
        return
    else:
        print(ERROR_INCOR_DISTARG)
        return

    return nodes

def rndTimeWindowsNodes(
    N:          "Number of vertices" = None,
    nodes:      "If nodes are provided, will overwrite `nodeIDs`, `xRange`, `yRange`, time windows will be applied on those nodes" = None,
    nodeIDs:    "A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
    xRange:     "A 2-tuple with minimum/maximum range of x" = (0, 100),
    yRange:     "A 2-tuple with minimum/maximum range of y" = (0, 100),
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

    # Initialize nodes ========================================================
    if (nodes == None):
        nodes = rndPlainNodes(N, nodeIDs, xRange, yRange)

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