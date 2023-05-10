import math



def vecPolar2XY(
    vecPolar:   "2-tuple (vVal, vDeg), `vVal` is the norm and `vDeg` is the direction, 0 as North, clockwise, in [0, 360)" = None
    ) -> "Given vector's norm and its degree to North, convert it into a 2-tuple vector":

    # Initialize ==============================================================
    (vVal, vDeg) = vecPolar

    vX = None
    vY = None

    while(vDeg < 0):
        vDeg = vDeg + 360

    while(vDeg >= 360):
        vDeg = vDeg - 360

    vX = vVal * math.sin(math.radians(vDeg))
    vY = vVal * math.cos(math.radians(vDeg))

    return (vX, vY)

def vecXY2Polar(
    vecXY:      "2-tuple (vX, vY), the coordinate of vector" = None
    ) -> "Given a 2-tuple, convert it into a norm and a direction in degree":
    
    (vX, vY) = vecXY    
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

    return (vVal, vDeg)

def calPolarVecAddition(
    vecPolar1:  "First vector in polar system" = None,
    vecPolar2:  "Second vector in polar system" = None
    ) -> "Given two vectors' norm and their meteorological degrees, get v3 that v3 = v1 + v2":
    # Change to 2-tuple vector ================================================
    (v1X, v1Y) = vecPolar2XY(vecPolar1)
    (v2X, v2Y) = vecPolar2XY(vecPolar2)

    # Get v3 ==================================================================
    (v3X, v3Y) = (v1X + v2X, v1Y + v2Y)

    # Get vector norm and direction ===========================================
    (v3Val, v3Deg) = vecXY2Polar((v3X, v3Y))

    return (v3Val, v3Deg)

def calPolarVecSubtract(
    vecPolar1:  "First vector in polar system" = None,
    vecPolar2:  "Second vector in polar system" = None
    ) -> "Given two vectors' norm and their meteorological degrees, get v3 that v1 = v2 + v3":
    # Change to 2-tuple vector ================================================
    (v1X, v1Y) = vecPolar2XY(vecPolar1)
    (v2X, v2Y) = vecPolar2XY(vecPolar2)

    # Get v3 ==================================================================
    (v3X, v3Y) = (v1X - v2X, v1Y - v2Y)

    # Get vector norm and direction ===========================================
    v3Val, v3Deg = vecXY2Polar((v3X, v3Y))

    return (v3Val, v3Deg)

def calXYVecAddition(
    vecXY1:     "First vector in xy system" = None,
    vecXY2:     "Second vector in xy system" = None
    ) -> "Given two vectors in XY, return the addition of two vectors":
    return (vecXY1[0] + vecXY2[0], vecXY1[1] + vecXY2[1])

def calXYVecSubtract(
    vecXY1:     "First vector in xy system" = None,
    vecXY2:     "Second vector in xy system" = None
    ) -> "Given two vectors in XY, return the subtract of two vectors":
    return (vecXY1[0] - vecXY2[0], vecXY1[1] - vecXY2[1])
