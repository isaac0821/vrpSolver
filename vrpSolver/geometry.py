import math
import heapq

from .const import *
from .common import *

def valDeg2Vec(
    vVal:       "Norm of the vector",
    vDeg:       "Degree of the vector, 0 as North, in [0, 360)"
    ) -> "Given vector's norm and its degree to North, convert it into a 2-tuple vector":

    vX = None
    vY = None

    while(vDeg < 0):
        vDeg = vDeg + 360

    while(vDeg >= 360):
        vDeg = vDeg - 360

    vX = vVal * math.sin(math.radians(vDeg))
    vY = vVal * math.cos(math.radians(vDeg))

    return (vX, vY)

def vec2ValDeg(
    vec:        "2-tuple, as the vector"
    ) -> "Given a 2-tuple, convert it into a norm and a direction in degree":
    
    (vX, vY) = vec
    
    vDeg = None
    vVal = None
    if (abs(vX) <= 0.0001):
        if (vY >= 0):
            vDeg = 0
            vVal = vY
        elif (vY < 0):
            vDeg = 180
            vVal = -vY
    elif (abs(vY) <= 0.0001):
        if (vX >= 0):
            vVal = vX
            vDeg = 90
        elif (vX < 0):
            vVal = -vX
            vDeg = 270
    else:
        vVal = math.sqrt(vX**2 + vY**2)
        # 1st quad
        if (vX > 0 and vY >= 0):
            vDeg = math.degrees(math.atan(vX / vY))
        # 2nd quad
        elif (vX > 0 and vY < 0):
            vDeg = 180 + math.degrees(math.atan(vX / vY))
        # 3rd quad
        elif (vX < 0 and vY < 0):
            vDeg = 180 + math.degrees(math.atan(vX / vY))
        # 4th quad
        elif (vX < 0 and vY >= 0):
            vDeg = 360 + math.degrees(math.atan(vX / vY))

    return vVal, vDeg

def calVecAddition(
    v1Val:      "Norm of vector 1", 
    v1Deg:      "Degree of the vector 1, 0 as North, in [0, 360)", 
    v2Val:      "Norm of vector 2", 
    v2Deg:      "Degree of the vector 2, 0 as North, in [0, 360)",
    ) -> "Given two vectors' norm and their meteorological degrees, get v3 that v3 = v1 + v2":
    # Change to 2-tuple vector ================================================
    (v1X, v1Y) = valDeg2Vec(v1Val, v1Deg)
    (v2X, v2Y) = valDeg2Vec(v2Val, v2Deg)

    # Get v3 ==================================================================
    (v3X, v3Y) = (v1X + v2X, v1Y + v2Y)

    # Get vector norm and direction ===========================================
    v3Val, v3Deg = vec2ValDeg((v3X, v3Y))

    return v3Val, v3Deg

def calVecSubtract(
    v1Val:      "Norm of vector 1", 
    v1Deg:      "Degree of the vector 1, 0 as North, in [0, 360)", 
    v2Val:      "Norm of vector 2", 
    v2Deg:      "Degree of the vector 2, 0 as North, in [0, 360)", 
    ) -> "Given two vectors' norm and their meteorological degrees, get v3 that v1 = v2 + v3":
    # Change to 2-tuple vector ================================================
    (v1X, v1Y) = valDeg2Vec(v1Val, v1Deg)
    (v2X, v2Y) = valDeg2Vec(v2Val, v2Deg)

    # Get v3 ==================================================================
    (v3X, v3Y) = (v1X - v2X, v1Y - v2Y)

    # Get vector norm and direction ===========================================
    v3Val, v3Deg = vec2ValDeg((v3X, v3Y))

    return v3Val, v3Deg

def getSweepSeq(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None,
    centerLoc:  "List, [x, y], the center point" = None,
    isClockwise: "True if the sweeping direction is clock-wise, False otherwise" = True,
    initDeg:    "Starting direction of the sweeping, 0 as North" = 0
    ) -> "Given a set of locations, and a center point, gets the sequence from sweeping":

    # Initialize heap =========================================================
    degHeap = []
    centerLocNodes = []
    
    # Build heap ==============================================================
    for nodeID in nodes:
        dist = distEuclidean2D(nodes[nodeID]['loc'], centerLoc)
        # If the nodes are too close, separate it/them
        if (dist <= CONST_EPSILON):
            centerLocNodes.append(nodeID)
        else:
            dx = nodes[nodeID]['loc'][0] - centerLoc[0]
            dy = nodes[nodeID]['loc'][1] - centerLoc[1]
            (_, deg) = vec2ValDeg([dx, dy])
            # Calculate angles
            evalDeg = None
            if (isClockwise):
                evalDeg = deg - initDeg
            else:
                evalDeg = initDeg - deg
            while(evalDeg > 360):
                evalDeg -= 360
            while(evalDeg < 0):
                evalDeg += 360
            heapq.heappush(degHeap, (evalDeg, nodeID))

    # Sweep ===================================================================
    sweepSeq = []
    while (len(degHeap)):
        sweepSeq.append(heapq.heappop(degHeap)[1])
    sweepSeq.extend(centerLocNodes)

    return sweepSeq
                