import math
import random
import geopy.distance

from .geometry import *
from .timeWindows import *
from .const import *

'''
# Description =================================================================
# This script is to simulate the movement of clouds
# Assumption:
# 1. The rain and snow happens only when there are clouds above
# 2. Wind speed is not related to the speed of clouds
# 3. Spinning of clouds are not modeled

# Reference ===================================================================
# Taylor's Hypothesis ---------------------------------------------------------
#     A series of changes in time at a fixed place is due to the passage of an unchanging 
#     spatial pattern over that locale. That is, when observing a cloud or thunderstorm 
#     passing overhead, the clouds of the thunderstorm floating by are effectively 
#     "unchanging spatial patterns" (big blocks of cloud). Put another way, the lateral 
#     change in the cloud conditions across regions of the meteorological event (e.g., 
#     cloud-sky-cloud-sky-cloud) can be directly connected locally with a variation in 
#     irradiance measurement over time (e.g., a periodic measure of dark-bright-dark-bright-dark).
# Ref: https://www.e-education.psu.edu/eme810/node/583

# Key parameters --------------------------------------------------------------
# Average speed of meteorological phenomena: 17 m/s
# Four types of meteorological phenomena and their scale:
# - Cumulus: 2-5 km, 10-100 minutes
# - Cumulonimbus: 10-200 km, 1-5 hours
# - Cumulonimbus cluster: 50-1000 km, 3-36 hours
# - Synoptic: 1000-400000 km, 2-15+ days
# We only consider: 1) Cumulus, 2) Cumulonimbus, and 3) Cumulus + Cumulonimbus (?) situations, 
#     the other two are too large
'''

def getCloudCurrentPosition(
    cloud:      "A cloud dictionary" = None,
    timestamp:  "Absolute time stamp" = None
    ) -> "Given a cloud dictionary and needed time stamps, returns the current position of cloud polygon":

    # Time elapsed ============================================================
    elapsed = timestamp - cloud['appearTime']
    if (elapsed > cloud['existDur'] + CONST_EPSILON or elapsed < 0):
        return None

    # Calculate current location ==============================================
    currentPoly = []
    dist = elapsed * cloud['movingVecPolar'][0]
    for pt in cloud['initPolyLatLon']:
        currentPoly.append(pointInDistLatLon(
            pt = pt, 
            direction = cloud['movingVecPolar'][1],
            distMeters = dist))

    return currentPoly

def getLocCoverByCloudsTW(
    clouds:     "A list of cloud dictionaries" = None,
    loc:        "A location, in [lat, lon]" = None
    ) -> "Given a cloud and a location, returns the time window that the location is covered by cloud":
    
    # Initialize ==============================================================
    coverTW = []

    # Calculate for each cloud ================================================
    for cloud in clouds:
        tws = twMovingPtInsidePolyLatLon(
            ptLatLon = loc,
            polyLatLon = cloud['initPolyLatLon'],
            vecPolarPt = (0, 0),
            vecPolarPoly = cloud['movingVecPolar'])
        if (len(tws) > 0):
            for tw in tws:
                if (tw[0] != None and tw[1] != None and tw[0] < cloud['existDur']):
                    coverTW.append([cloud['appearTime'] + tw[0], cloud['appearTime'] + min(cloud['existDur'], tw[1])])
    
    # Merge time windows of clouds ============================================
    coverTW = mergeTimeWindows(coverTW)

    return coverTW

def getSegCoverByCloudsTW(
    clouds:     "A list of cloud dictionaries" = None,
    loc1:       "Lat/lon of the first location, in [lat, lon]" = None,
    loc2:       "Lat/lon of the first location, in [lat, lon]" = None
    ) -> "Given a list of clouds and two locations, returns the time windows that anywhere between two locations is covered by cloud":

    # FIXME: Will need to be rewritten, now using stupid approach, but it works

    # Initialize ==============================================================
    coverTWs = []

    # For each cloud, it should generate at most one cover time window ========
    for cloud in clouds:
        ts = cloud['appearTime']
        te = cloud['appearTime'] + cloud['existDur']
        distInterval = 100 # Try every 100 m along loc1 to loc2 to see if it will be covered, only for now
        heading = getHeadingLatLon(loc1, loc2)
        accDist = 0
        remainDist = distLatLon(loc1, loc2)

        # Try to see if any of the location between loc1 and loc2 will be covered by cloud
        curLoc = loc1
        intersectFlag = False
        while (remainDist > 0 and not intersectFlag):
            tmpTW = None
            tmLeft = None
            tmRight = None
            # First, if this loc is already inside, search right directions
            if (isPtInsidePoly(curLoc, cloud['initPolyLatLon'])):
                twpTW = [ts, ts] 
                intersectFlag = True
            else:
                # Then, if this loc is not inside, try to see if it can be covered            
                tryCoverTW = getLocCoverByCloudsTW([cloud], curLoc)
                
                if (len(tryCoverTW) > 0):
                    tmpTW = tryCoverTW[0] # There should be only one time window
                    intersectFlag = True
                else:
                    if (remainDist > distInterval):
                        remainDist -= distInterval
                        accDist += distInterval
                        curLoc = pointInDistLatLon(curLoc, heading, distInterval)
                    else:
                        remainDist = 0
                        curLoc = loc2

            # If we find a loc will be covered, search in both direction to see when will be covered
            if (tmpTW != None):
                # Search left
                # print("Left")
                tsLeft = cloud['appearTime']         # Suppose to be not covered
                teLeft = tmpTW[0]   # Suppose to be covered
                while (teLeft - tsLeft > 2): # Search until the 1 sec of precision
                    tmLeft = (teLeft - tsLeft) / 2 + tsLeft
                    curCloud = getCloudCurrentPosition(cloud, tmLeft)
                    segIntCloud = isSegCrossPoly([loc1, loc2], curCloud)
                    if (segIntCloud):
                        teLeft = tmLeft
                    else:
                        tsLeft = tmLeft
                    # print(segIntCloud, tsLeft, tmLeft, teLeft, curCloud)
                tmLeft = (teLeft - tsLeft) / 2 + tsLeft
                # Search right
                # print("Right")
                tsRight = tmpTW[1]
                teRight = cloud['appearTime'] + cloud['existDur']
                while (teRight - tsRight > 2):
                    tmRight = (teRight - tsRight) / 2 + tsRight
                    curCloud = getCloudCurrentPosition(cloud, tmRight)
                    segIntCloud = isSegCrossPoly([loc1, loc2], curCloud)
                    # print([loc1, loc2], curCloud)
                    if (segIntCloud):
                        tsRight = tmRight
                    else:
                        teRight = tmRight
                    # print(segIntCloud, tsRight, tmRight, teRight)
                tmRight = (teRight - tsRight) / 2 + tsRight

            if (intersectFlag):
                if (tmLeft != None and tmRight != None):
                    coverTWs.append([tmLeft, tmRight])

    # Merge time windows of clouds ============================================
    coverTWs = mergeTimeWindows(coverTWs)

    return coverTWs
