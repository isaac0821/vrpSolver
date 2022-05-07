import math
import random
import geopy.distance

from .geometry import *
from .timeWindows import *
from .const import *

'''
# Description =================================================================
# This script is to simulate the movement of clouds
# FIXME: This script may need to be rewritten for more detailed clouds
# Assumption:
# 1. The rain and snow happens only when there are clouds above
# 2. Wind speed is not related to the movement of clouds, in both speed and direction
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
'''
# Parameters for clouds
CLOUD_SPD_ABS_RANGE = [16.5, 17.5] # 17 m/s
CLOUD_SPD_DEG_RANGE = [30, 60] # [deg, deg]

CLOUD_CUMULUS_SIZE = [2000, 5000] # 2-5 km
CLOUD_CUMULONIMBUS_SIZE = [10000, 200000] # 10-200 km

CLOUD_CUMULUS_DURATION = [10 * 60, 100 * 60] # 10 - 100 mins
CLOUD_CUMULONIMBUS_DURATION = [60 * 60, 5 * 60 * 60] # 1-5 hours

def rndWinds(
    degAvg:     "Prevailing direction of wind - in the direction of wind vector, not meteorology term" = None,
    degSpread:  "Wind direction difference in 10th to 90th percentile bonds" = None,
    degPattern: "1) String, 'Consistent' or, \
                 2) String, 'Normal' or, \
                 3) String, 'Swap' or, \
                 4) String, 'Uniform'" = None,
    spdAvg:     "Average wind speed" = None,
    spdSpread:  "Wind speed difference in 10th to 90th percentile bonds" = None,
    spdPattern: "1) String, 'Consistent' or, \
                 2) String, 'Increasing'  or, \
                 3) String, 'Decreasing' or, \
                 4) String, 'Normal' or, \
                 5) String, 'Uniform'" = None,
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
        elif (degPattern == 'Uniform'):
            deg = random.uniform(degAvg - degSpread, degAvg + degSpread)
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

        elif (spdPattern == 'Uniform'):
            spd = random.uniform(spdAvg - spdSpread, spdAvg + spdSpread)
        winds[w]['windSpd'] = spd

    return winds

def rndClouds(
    polyLatLon: "The polygon that clouds will be covering" = None,
    cloudSize:  "Size of the clouds, e.g., [3000, 12000], in [m, m]" = None,
    cloudShape: "Overwrite cloudSize, if provided, it will be the shape of the clouds" = None,
    cloudDur:   "Exist duration of clouds, e.g., [1200, 1800], in [sec, sec]" = None,
    cloudDeg:   "The direction of where the clouds are entering" = [CLOUD_SPD_DEG_RANGE[0], CLOUD_SPD_DEG_RANGE[1]],
    cloudOri:   "Orientation of the cloud, if provided as None, the orientation will be randomized" = None,
    numOfClouds: "Total number of clouds" = None,
    duration:   "Duration of the instance, in [sec]" = [0, 36000],
    ) -> "Given a polygon, returns a list of clouds that move through the polygon":

    # Key parameters ==========================================================
    # Average speed of meteorological phenomena: 17 m/s
    # Four types of meteorological phenomena and their scale:
    # - Cumulus: 2-5 km, 10-100 minutes
    # - Cumulonimbus: 10-200 km, 1-5 hours
    # - Cumulonimbus cluster: 50-1000 km, 3-36 hours
    # - Synoptic: 1000-400000 km, 2-15+ days

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
        # Size of cloud -------------------------------------------------------
        widthInMeter = None
        lengthInMeter = None
        if (cloudShape != None and type(cloudShape) == list):
            widthInMeter = cloudShape[0]
            lengthInMeter = cloudShape[1]
        elif (cloudSize != None and type(cloudSize) == list):
            widthInMeter = random.randrange(cloudSize[0], cloudSize[1])
            lengthInMeter = random.randrange(cloudSize[0], cloudSize[1])

        # Moving direction ----------------------------------------------------
        direction = None
        if (cloudDeg != None and type(cloudDeg) == list):
            direction = random.randrange(cloudDeg[0], cloudDeg[1])
        elif (type(cloudDeg) == int or type(cloudDeg) == float):
            direction = cloudDeg

        # Existing duration ---------------------------------------------------
        existDur = None
        if (cloudDur != None and type(cloudDur) == list):
            existDur = random.randrange(cloudDur[0], cloudDur[1])
        elif (type(cloudDur) == int or type(cloudDur) == float):
            existDur = cloudDur

        # Initial location of the cloud ---------------------------------------        
        enterLoc = pointInDistLatLon(
            pt = cloudCentroid[t]['loc'],
            direction = direction - 180,
            distMeters = math.sqrt(widthInMeter**2 + lengthInMeter**2) + CONST_EPSILON)
        
        # Cloud generation ----------------------------------------------------
        initPolyLatLon = rectInWidthLengthOrientationLatLon(
            centroidLatLon = enterLoc,
            widthInMeter = widthInMeter,
            lengthInMeter = lengthInMeter,
            oriDeg = random.randrange(-45, 45) if cloudOri == None else cloudOri)        
        clouds.append({
            'appearTime': appearTime[t],
            'initPolyLatLon': initPolyLatLon,
            'widthInMeter': widthInMeter,
            'lengthInMeter': lengthInMeter,
            'movingVecPolar': (17, direction), # cloud moving speed: 17 m/s
            'existDur': existDur
        })

    return clouds

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
                    # print(curLoc)
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
                while (teLeft - tsLeft > 1): # Search until the 1 sec of precision
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
                while (teRight - tsRight > 1):
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
