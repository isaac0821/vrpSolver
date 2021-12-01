import random
import math
import tripy
import geopy
import copy

from .msg import *
from .geometry import *
from .common import *

def rndPlainNodes(
    N:          "Number of vertices" = None,
    nodeIDs:    "A list of node IDs, `N` will be overwritten if `nodeIDs` is given" = None,
    distr:      "Spatial distribution of nodes, the options are \
                 1) String, (default) 'uniformSquare', or \
                 2) String, 'uniformCircle', or \
                 3) String, 'uniformPoly', or\
                 4) String, (not available) 'uniformOnNetwork', or\
                 5) String, 'clustered, or \
                 6) String, (not available) 'normal2D'" = 'uniformSquare',
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
                 4) for 'uniformOnNetwork'\
                    {\
                        'network': list of arcs that can be sampled \
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
                 6) for 'normal2D'\
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
            theta = random.uniform(0, 2 * math.pi)
            r = math.sqrt(random.uniform(0, distrArgs['radius'] ** 2))
            x = distrArgs['centerLoc'][0] + r * math.cos(theta)
            y = distrArgs['centerLoc'][1] + r * math.sin(theta)
            nodes[n] = {'loc': (x, y)}
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
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            print(ERROR_MISSING_DISTRARGS)
            return
        if (('numCluster' not in distrArgs or 'xRange' not in distrArgs or 'yRange' not in distrArgs) and 'centroidLocs' not in distrArgs):
            print(ERROR_MISSING_DISTRARGS_CLUSTER)
            return
        if ('clusterDiameter' not in distrArgs):
            distrArgs['clusterDiameter'] = 20
        # Create nodes --------------------------------------------------------
        centroidLocs = []
        if ('centroidLocs' in distrArgs):
            centroidLocs = distrArgs['centroidLocs']
        else:
            for cl in range(distrArgs['numCluster']):
                x = random.randrange(distrArgs['xRange'][0], distrArgs['xRange'][1])
                y = random.randrange(distrArgs['yRange'][0], distrArgs['yRange'][1])
                centroidLocs.append([x, y])
        for n in nodeIDs:
            # First pick a centroidLoc
            idx = random.randint(0, len(centroidLocs) - 1)
            ctrLoc = centroidLocs[idx]
            theta = random.uniform(0, 2 * math.pi)
            r = math.sqrt(random.uniform(0, distrArgs['clusterDiameter']))
            x = ctrLoc[0] + r * math.cos(theta)
            y = ctrLoc[1] + r * math.sin(theta)
            nodes[n] = {'loc': (x, y)}
    elif (distr == "normal2D"):
        print("Stay tune")
        return
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
