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
                 4) for 'uniformRoadNetworkPoly'\
                    {\
                        'network': list of arcs that can be sampled \
                        'poly': nodes should generated within the polygon, if not provided, will consider the entire network \
                    }\
                 5) for 'uniformRoadNetworkCircle'\
                    {\
                        'network': list of arcs that can be sampled \
                        'circle': {\
                            'centerLoc': [lat, lon], \
                            'radius': radius in [m] \
                        }\
                    }\
                 6) for 'clustered\
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
        msgError(ERROR_MISSING_N)
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
            msgError(ERROR_MISSING_DISTRARGS_UNICC)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformCircleXY(distrArgs['radius'], distrArgs['centerLoc'])
            }

    elif (distr == 'uniformCircleLatLon'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            msgError(ERROR_MISSING_DISTRARGS)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformCircleLatLon(distrArgs['radius'], distrArgs['centerLoc'])
            }

    elif (distr == 'uniformPoly'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            msgError(ERROR_MISSING_DISTRARGS)
            return
        if ('poly' not in distrArgs):
            msgError(ERROR_MISSING_DISTRARGS_UNIPOLY)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformPoly(distrArgs['poly'])
            }

    elif (distr == 'uniformPolys'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            msgError(ERROR_MISSING_DISTRARGS)
            return
        if ('polys' not in distrArgs):
            msgError(ERROR_MISSING_DISTRARGS_UNIPOLY)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtUniformPolys(distrArgs['polys'])
            }

    elif (distr == 'uniformRoadNetworkPoly'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            msgError(ERROR_MISSING_DISTRARGS)
            return
        if ('roadNetwork' not in distrArgs):
            msgError(ERROR_MISSING_DISTRARGS_UNIROADNETWORK)
            return

        # Create nodes --------------------------------------------------------
        nodeLocs = getRndPtRoadNetworkPoly(
            N if N != None else len(nodeIDs),
            distrArgs['roadNetwork'], 
            distrArgs['poly'] if 'poly' in distrArgs else None)
        for n in range(len(nodeIDs)):
            nodes[nodeIDs[n]] = {
                'loc': nodeLocs[n]
            }

    elif (distr == 'uniformRoadNetworkCircle'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            msgError(ERROR_MISSING_DISTRARGS)
            return
        if ('roadNetwork' not in distrArgs):
            msgError(ERROR_MISSING_DISTRARGS_UNIROADNETWORK)
            return
        if ('centerLoc' not in distrArgs or 'radius' not in distrArgs):
            print("Missing circle definition")
            return

        # Create nodes --------------------------------------------------------
        nodeLocs = getRndPtRoadNetworkCircle(
            N if N != None else len(nodeIDs),
            distrArgs['roadNetwork'], 
            distrArgs['centerLoc'],
            distrArgs['radius'])
        for n in range(len(nodeIDs)):
            nodes[nodeIDs[n]] = {
                'loc': nodeLocs[n]
            }

    elif (distr == 'clusterXY'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            msgError(ERROR_MISSING_DISTRARGS)
            return
        if ('centroidLocs' not in distrArgs or 'clusterDiameter' not in distrArgs):
            msgError(ERROR_MISSING_DISTRARGS_CLUSTER)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtClusterXY(distrArgs['centroidLocs'], distrArgs['clusterDiameter'])
            }

    elif (distr == 'clusterLatLon'):
        # Sanity check --------------------------------------------------------
        if (distrArgs == None):
            msgError(ERROR_MISSING_DISTRARGS)
            return
        if ('centroidLocs' not in distrArgs or 'clusterDiameterInMeters' not in distrArgs):
            msgError(ERROR_MISSING_DISTRARGS_CLUSTER)
            return
        # Create nodes --------------------------------------------------------
        for n in nodeIDs:
            nodes[n] = {
                'loc': getRndPtClusterLatLon(distrArgs['centroidLocs'], distrArgs['clusterDiameterInMeters'])
            }
    else:
        msgError(ERROR_INCOR_DISTARG)
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
    nodes = dict(plainNodes)

    # Check for required fields ===============================================
    if (plainNodes == None):
        msgError(ERROR_MISSING_PLAINNODES)
        return
    if (timeSpan == None):
        msgError(ERROR_MISSING_TIMESPAN)
        return
    validTWOpts = ['Known', 'Periodic', 'Random']
    if (twType not in validTWOpts):
        msgError(ERROR_OPTS_TWTYPE % validTWOpts)
        return
    elif (twType == 'Known'):
        if (twArgs == None):
            msgError(ERROR_MISSING_TWKNOWNARGS)
            return
        elif ('timeWindows' not in twArgs):
            print (ERROR_MISSING_TWKNOWNARGS_TW)
            return
    elif (twType == 'Periodic'):
        if (twArgs == None):
            msgError(ERROR_MISSING_TWPERIODICARGS)
            return
        elif ('cycleTime' not in twArgs):
            msgError(ERROR_MISSING_TWPERIODICARGS_CT)
            return
        elif ('startTimeInCycle' not in twArgs):
            msgError(ERROR_MISSING_TWPERIODICARGS_ST)
            return
        elif ('endTimeInCycle' not in twArgs):
            msgError(ERROR_MISSING_TWPERIODICARGS_ET)
            return
    elif (twType == 'Random'):
        if (twArgs == None):
            msgError(ERROR_MISSING_TWRANDOMARGS)
            return
        elif ('indFlag' not in twArgs):
            msgError(ERROR_MISSING_TWRANDOMARGS_FG)
            return
        elif ('avgBtwDur' not in twArgs):
            msgError(ERROR_MISSING_TWRANDOMARGS_BTW)
            return
        elif ('avgTWDur' not in twArgs):
            msgError(ERROR_MISSING_TWRANDOMARGS_TW)
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

def getRndPtRoadNetworkPoly(
    N:          "Number of nodes" = 1,
    roadNetwork: "Dictionary of road network in the format of \
                {\
                    roadID: {\
                        'line': [[lat, lon], [lat, lon], ...]\
                    }\
                }" = None,
    poly:       "Nodes should also within this polygon" = None,
    ) -> "Given a road network, generate customers that locates on the road network":
    
    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for road in roadNetwork:
        roadLength = 0
        includedFlag = False
        if (poly == None):
            includedFlag = True
        else:
            for i in range(len(roadNetwork[road]['line'])):
                if (isPtOnPoly(roadNetwork[road]['line'][i], poly)):
                    includedFlag = True
                    break

        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(roadNetwork[road]['line']) - 1):
                roadLength += distLatLon(roadNetwork[road]['line'][i], roadNetwork[road]['line'][i + 1])
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(road)

    # Check if there are roads included =======================================
    if (sum(lengths) == 0):
        return None

    # Use accept-denial to test if the node is within poly ====================
    # FIXME: Inefficient approach, will need to be rewritten
    nodeLocs = []
    for i in range(N):
        lat = None
        lon = None
        if (poly == None):
            idx = rndPick(lengths)
            edgeLength = lengths[idx]
            edgeDist = random.uniform(0, 1) * edgeLength
            (lat, lon) = getMileageInPathLatLon(roadNetwork[roadIDs[idx]]['line'], edgeDist)
        else:
            insideFlag = False
            while (not insideFlag):
                idx = rndPick(lengths)
                edgeLength = lengths[idx]
                edgeDist = random.uniform(0, 1) * edgeLength
                (lat, lon) = getMileageInPathLatLon(roadNetwork[roadIDs[idx]]['line'], edgeDist)
                if (isPtOnPoly([lat, lon], poly)):
                    insideFlag = True
        nodeLocs.append((lat, lon))
    return nodeLocs

def getRndPtRoadNetworkCircle(
    N:          "Number of nodes" = 1,
    roadNetwork: "Dictionary of road network in the format of \
                {\
                    roadID: {\
                        'line': [[lat, lon], [lat, lon], ...]\
                    }\
                }" = None,\
    centerLoc:  "Center location" = None,
    radius:     "Radius in [m]" = None
    ) -> "Given a road network, generate customers that locates on the road network":
    
    # Calculate the length of each edge =======================================
    lengths = []
    roadIDs = []
    for road in roadNetwork:
        roadLength = 0
        includedFlag = False
        for i in range(len(roadNetwork[road]['line'])):
            if (distLatLon(roadNetwork[road]['line'][i], centerLoc) <= radius):
                includedFlag = True
                break

        # Check if this road is inside polygon
        if (includedFlag):
            for i in range(len(roadNetwork[road]['line']) - 1):
                roadLength += distLatLon(roadNetwork[road]['line'][i], roadNetwork[road]['line'][i + 1])
            lengths.append(roadLength)            
        else:
            lengths.append(0)

        roadIDs.append(road)

    # Check if there are roads included =======================================
    if (sum(lengths) == 0):
        return None

    # Use accept-denial to test if the node is within poly ====================
    # FIXME: Inefficient approach, will need to be rewritten
    nodeLocs = []
    for i in range(N):
        lat = None
        lon = None
        insideFlag = False
        while (not insideFlag):
            idx = rndPick(lengths)
            edgeLength = lengths[idx]
            edgeDist = random.uniform(0, 1) * edgeLength
            (lat, lon) = getMileageInPathLatLon(roadNetwork[roadIDs[idx]]['line'], edgeDist)
            if (distLatLon([lat, lon], centerLoc) <= radius):
                insideFlag = True
        nodeLocs.append((lat, lon))
    return nodeLocs

def getRndPtUniformSquare(
    xRange:    "The range of x coordinates" = (0, 100),
    yRange:    "The range of y coordinates" = (0, 100)
    ) -> "Given the range of x, y, returns a random point in the square defined by the ranges":
    x = random.randrange(xRange[0], xRange[1])
    y = random.randrange(yRange[0], yRange[1])
    return (x, y)

def getRndPtUniformTriangle(
    triangle:   "The triangle for generating random points" = None
    ) -> "Given a triangle, generate a random point in the triangle uniformly":
    
    # Get three extreme points ================================================
    [x1, y1] = triangle[0]
    [x2, y2] = triangle[1]
    [x3, y3] = triangle[2]

    # Generate random points ==================================================
    rndR1 = random.uniform(0, 1)
    rndR2 = random.uniform(0, 1)
    x = (1 - math.sqrt(rndR1)) * x1 + math.sqrt(rndR1) * (1 - rndR2) * x2 + math.sqrt(rndR1) * rndR2 * x3
    y = (1 - math.sqrt(rndR1)) * y1 + math.sqrt(rndR1) * (1 - rndR2) * y2 + math.sqrt(rndR1) * rndR2 * y3

    return (x, y)

def getRndPtUniformPoly(
    poly:       "The polygon for generating random points" = None
    ) -> "Given a polygon, generate a random point in the polygons uniformly":

    # Get list of triangles ===================================================
    lstTriangle = tripy.earclip(poly)

    # Weight them and make draws ==============================================
    lstWeight = []
    for i in range(len(lstTriangle)):
        lstWeight.append(calTriangleAreaByCoords(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2]))

    # Select a triangle and randomize a point in the triangle =================
    idx = rndPick(lstWeight)
    (x, y) = getRndPtUniformTriangle(lstTriangle[idx])

    return (x, y)

def getRndPtUniformPolys(
    polys:       "A list of polygons for generating random points" = None
    ) -> "Given a list of polygons, generate a random point in the polygons uniformly":

    # Get all triangulated triangles ==========================================
    lstTriangle = []
    for p in polys:
        lstTriangle.extend(tripy.earclip(p))

    # Weight them and make draws ==============================================
    lstWeight = []
    for i in range(len(lstTriangle)):
        lstWeight.append(calTriangleAreaByCoords(lstTriangle[i][0], lstTriangle[i][1], lstTriangle[i][2]))

    # Select a triangle and randomize a point in the triangle =================
    idx = rndPick(lstWeight)
    (x, y) = getRndPtUniformTriangle(lstTriangle[idx])

    return (x, y)

def getRndPtClusterXY(
    centroidLocs: "A list of center locs of clusters" = None,
    clusterDiameter: "Diameter of cluster" = None
    ):
    idx = random.randint(0, len(centroidLocs) - 1)
    ctrLoc = centroidLocs[idx]
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, clusterDiameter))
    x = ctrLoc[0] + r * math.cos(theta)
    y = ctrLoc[1] + r * math.sin(theta)
    return (x, y)

def getRndPtClusterLatLon(
    centroidLocs: "A list of center locs of clusters" = None,
    clusterDiameterInMeters: "Diameter of cluster" = None
    ):
    idx = random.randint(0, len(centroidLocs) - 1)
    ctrLoc = centroidLocs[idx]
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, clusterDiameterInMeters))
    (lat, lon) = pointInDistLatLon(ctrLoc, theta, r)
    return (lat, lon)

def getRndPtUniformCircleXY(
    radius:     "Radius of the circle" = None,
    centerLoc:  "Center location of the circle" = None
    ):
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, radius ** 2))
    x = centerLoc[0] + r * math.cos(theta)
    y = centerLoc[1] + r * math.sin(theta)
    return (x, y)

def getRndPtUniformCircleLatLon(
    radius:     "Radius of the circle" = None,
    centerLoc:  "Center location of the circle" = None
    ):
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(0, radius ** 2))
    (lat, lon) = pointInDistLatLon(centerLoc, theta, r)
    return (lat, lon)
