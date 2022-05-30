import math

from .common import *
from .const import *

'''
# Description =================================================================
# This script is for processing time windows
#     including: TW (time window) and TWs (time windows)
# Notice:
# 1. Time windows can be unbounded, the unbounded side will be None, e.g., [None, 9] represents t <= 9
# 2. All the time windows include the edge (close set)
'''

def sortTWs(
    unsortedTWs: "List of unsorted time windows" = None,
    ) -> "Given a list of time windows, returns None if the TWs are overlapping returns a sorted list of time windows":

    # Example =================================================================
    # In: -> unsortedTWs = [[17.8, 20], [12, 14], [21, None], [17, 17.5], [None, 10]]
    #     -> sortTWs(unsortedTWs)
    # Out:-> [[None, 10], [12, 14], [17, 17.5], [17.8, 20], [21, None]]

    # Initialize ==============================================================
    sortedTWs = []

    # Separate infinite time windows ==========================================
    leftInfTW = None
    rightInfTW = None
    bothInfTW = None
    startTimeList = []
    endTimeList = []
    for tw in unsortedTWs:
        if (tw[0] == None and tw[1] == None):
            if (bothInfTW != None):
                msgError(ERROR_TW_OVERLAP)
                return None
            else:
                bothInfTW = tw # tw = [None, None]
        if (tw[0] == None and tw[1] != None):
            if (leftInfTW != None):
                msgError(ERROR_TW_OVERLAP)
                return None
            else:
                leftInfTW = tw # tw = [None, time]
        if (tw[0] != None and tw[1] == None):
            if (rightInfTW != None):
                msgError(ERROR_TW_OVERLAP)
                return None
            else:
                rightInfTW = tw # tw = [time, None]
        if (tw[0] != None and tw[1] != None):
            startTimeList.append(tw[0])
            endTimeList.append(tw[1])

    # Quick filter for bothInfTW ==============================================
    if (bothInfTW != None):
        if (len(unsortedTWs) > 1):
            msgError(ERROR_TW_OVERLAP)
            return None
        else:
            return [bothInfTW]

    # Sort start/end time of each finite time windows =========================
    # Get the list of start times and end times and all timestamps
    sortedStartIndexList = sorted(range(len(startTimeList)), key=lambda k: startTimeList[k])
    sortedEndIndexList = sorted(range(len(endTimeList)), key=lambda k: endTimeList[k])

    # Check overlapping =======================================================
    if (sortedStartIndexList != sortedEndIndexList):
        msgError(ERROR_TW_OVERLAP)
        return None

    # Sort time windows =======================================================
    sortedTWs = [[startTimeList[sortedStartIndexList[i]], endTimeList[sortedStartIndexList[i]]] for i in range(len(sortedStartIndexList))]

    # Check overlapping =======================================================
    for i in range(len(sortedTWs) - 1):
        if (sortedTWs[i][1] > sortedTWs[i + 1][0]):
            msgError(ERROR_TW_OVERLAP)
            return None 

    # Add back infinite time windows ==========================================
    if (leftInfTW != None):
        if (len(sortedTWs) > 0):
            if (leftInfTW[1] > sortedTWs[0][0]):
                msgError(ERROR_TW_OVERLAP)
                return None
            else:
                sortedTWs.insert(0, leftInfTW)
        else:
            sortedTWs.insert(0, leftInfTW)
    if (rightInfTW != None):
        if (len(sortedTWs) > 0):
            if (rightInfTW[0] < sortedTWs[-1][1]):
                msgError(ERROR_TW_OVERLAP)
                return None
            else:
                sortedTWs.append(rightInfTW)
        else:
            sortedTWs.append(rightInfTW)

    return sortedTWs

def unionTWs(
    tws:        "List of time windows, might be overlapping" = None
    ) -> "Given a list of time windows, which might be overlapping, merge the time windows that are overlapping":

    # Examples ================================================================
    # Example 1:
    # In: -> tws1 = [[0, 10], [12, 14], [18, 20]]
    #     -> tws2 = [[9, 13], [19, 23]]
    #     -> twsBeforeMerging = []
    #     -> twsBeforeMerging.extend(tws1)
    #     -> twsBeforeMerging.extend(tws2)
    #     -> twsAfterMerging = unionTWs(twsBeforeMerging)
    # Out:-> [[0, 14], [18, 23]]
    # Example 2:
    # In: -> tws = [[31, None], [17.5, 17.8], [15.6, 17.6], [2, 14], [None, 3]]
    #     -> unionTimeWindows(tws)
    # Out:-> [[None, 14], [15.6, 17.8], [31, None]]

    # Initialize ==============================================================
    unionedTWs = []

    # Separate infinite time windows ==========================================
    leftInfTW = None
    rightInfTW = None
    bothInfTW = None
    for tw in tws:
        if (tw[0] == None and tw[1] == None):
            bothInfTW = tw # tw = [None, None]
            return bothInfTW
        if (tw[0] == None and tw[1] != None):
            if (leftInfTW == None or tw[1] > leftInfTW[1]):
                leftInfTW = tw
        if (tw[0] != None and tw[1] == None):
            if (rightInfTW == None or tw[0] < rightInfTW[0]):
                rightInfTW = tw

    # Merge finite tws ========================================================
    for tw in tws:
        # Find the index of tw that involves the s/e of tw, if there is any
        sIndex = None
        eIndex = None
        for bgTwIndex in range(len(unionedTWs)):
            if (insideInterval(tw[0], unionedTWs[bgTwIndex])):
                sIndex = bgTwIndex
            if (insideInterval(tw[1], unionedTWs[bgTwIndex])):
                eIndex = bgTwIndex
            if (sIndex != None and eIndex != None):
                break
        # Find the s/e of merged tw
        ts = None
        te = None
        if (sIndex != None):
            ts = unionedTWs[sIndex][0]
        else:
            ts = tw[0]
        if (eIndex != None):
            te = unionedTWs[eIndex][1]
        else:
            te = tw[1]
        unionedTWs.append([ts, te])

        # Remove the tws that partially covered by [ts, te], remove them
        tmpMergedTWs = [i for i in unionedTWs]
        unionedTWs = [[ts, te]] # merged tw should be included
        for newTW in tmpMergedTWs:
            # For existed tws, both s and e should not inside the [ts, te]
            if (not insideInterval(newTW[0], [ts, te]) and not insideInterval(newTW[1], [ts, te])):
                unionedTWs.append(newTW)
        
        # Sort tws
        unionedTWs = sortTWs(unionedTWs)

    return unionedTWs

def intersectTWs(
    tws1:       "List of time windows, must not be overlapping" = None,
    tws2:       "List of time windows, must not be overlapping" = None
    ) -> "Given two lists of time windows, return a new list of time windows that are the intersection of those two lists":

    intTWs = []
    for tw1 in tws1:
        for tw2 in tws2:
            intTW = None
            # Case by case ============================================================
            if (not insideInterval(tw2[0], tw1) and not insideInterval(tw2[1], tw1)):
                # Case 1: tw1 is entirely inside tw2
                if (tw2[0] < tw1[0] and tw1[1] < tw2[1]):
                    intTW = tw1
                # Case 2: No overlapping
                else:
                    intTW = None
            elif (insideInterval(tw2[0], tw1) and not insideInterval(tw2[1], tw1)):
                # Case 3: First half of tw2 is inside tw1
                intTW = [tw2[0], tw1[1]]
            elif (not insideInterval(tw2[0], tw1) and insideInterval(tw2[1], tw1)):
                # Case 4: Last half of tw2 is inside tw1
                intTW = [tw1[0], tw2[1]]
            else:
                # Case 5: tw2 is entirely inside tw1
                intTW = [tw2[0], tw2[1]]
            if (intTW != None):
                intTWs.append(intTW)
    intTWs = sortTWs(intTWs)

    return intTWs

def lenTWOverlap(
    tw1:        "The first time window" = None,
    tw2:        "The second time window" = None
    ) -> "Given two time windows, returns the length of time that these two time windows are overlapping":

    # Initialize ==============================================================
    overlapTime = 0

    # Case by case ============================================================
    if (not insideInterval(tw2[0], tw1) and not insideInterval(tw2[1], tw1)):
        # Case 1: tw1 is entirely inside tw2
        if (tw2[0] < tw1[0] and tw1[1] < tw2[1]):
            overlapTime = tw1[1] - tw1[0]
        # Case 2: No overlapping
        else:
            overlapTime = 0
    elif (insideInterval(tw2[0], tw1) and not insideInterval(tw2[1], tw1)):
        # Case 3: First half of tw2 is inside tw1
        overlapTime = tw1[1] - tw2[0]
    elif (not insideInterval(tw2[0], tw1) and insideInterval(tw2[1], tw1)):
        # Case 4: Last half of tw2 is inside tw1
        overlapTime = tw2[1] - tw1[0]
    else:
        # Case 5: tw2 is entirely inside tw1
        overlapTime = tw2[1] - tw2[0]

    return overlapTime

def reverseTWs(
    tws:        "List of existing time windows" = None,
    epoch:      "Another time window that we wanna find opposite time windows of `tws` inside `epoch`" = None
    ) -> "Given a list of time windows `tws`, given another time window `epoch`, returns the opposite time windows of `tws` in `epoch":

    # Initialize ==============================================================
    oppTws = [[i for i in epoch]]

    # Reverse time windows ====================================================
    # Start with epoch, using `blockTWs` to reverse time windows
    for tw in tws:
        oppTws = blockTWs(oppTws, tw)
    return oppTws

def blockTWs(
    tws:        "List of existing available time windows" = None,
    block:      "A piece of time make part of tws unavailable" = None,
    ) -> "Given a list of non-overlapping time windows, use a time window to block (make the period unavailable) the existing `tws`":
    newTWs = []
    for tw in tws:
        # Use '<' in this section
        # Case 1: Block entirely inside the time window
        if (tw[0] < block[0] and block[1] < tw[1]):
            # Split into two time windows
            newTWs.append([tw[0], block[0]])
            newTWs.append([block[1], tw[1]])
        # Case 2: No overlap
        elif (block[1] < tw[0] or tw[1] < block[0]):
            # Keep origin time window
            newTWs.append(tw)
        # Case 3: Partially overlapping on the left of the time window
        elif (insideInterval(tw[0], block) and block[1] < tw[1]):
            # From the end of block to the end of tw
            newTWs.append([block[1], tw[1]])
        # Case 4: Partially overlapping on the right of the time window
        elif (tw[0] < block[0] and insideInterval(tw[1], block)):
            # From the start of tw to the start of block
            newTWs.append([tw[0], block[0]])
        # Case 5: Block entirely cover the time window
        elif (insideInterval(tw[0], block) and insideInterval(tw[1], block)):
            # No window added
            pass
    return newTWs

def getIndexInTimeWindow(
    bgTws:      "List of non-overlapping time windows as background" = None,
    t:          "Time stamp" = None
    ) -> "Given a time stamp and a list of time windows, check if the time stamp is in time windows":
    for i in range(len(bgTws)):
        if (insideInterval(t, bgTws[i])):
            return i
    return None

def getEarliestNextAvailTime(
    bgTws:      "List of non-overlapping time windows as background" = None,
    t:          "Time stamp" = None
    ) -> "Given a time stamp, get next launchable time":
    earliestNext = None
    for i in range(len(bgTws)):
        if (not insideInterval(t, bgTws[i])):
            if (t < bgTws[i][0]):
                if (earliestNext is None or earliestNext > bgTws[i][0]):
                    earliestNext = bgTws[i][0]
        else:
            return t
    return earliestNext

def getAvailStartTW(
    bgTws:      "List of non-overlapping time windows as background (CANNOT be covered by inserting time window)" = None,
    epoch:      "The outer epoch of `bgTws` time windows" = None,
    twLen:      "Length of time window that trying to add into bgTws (within epoch)" = None
    ) -> "Given a list of background time windows and an outer epoch, given the length of a new time window, \
          return a list of time windows that can start the tw with the length of `twLen`":

    # Initialize ==============================================================
    startTws = []

    # Get the list of time windows that can insert new time window ============
    availTws = reverseTWs(bgTws, epoch)

    # For each available time windows, try to see if it is long enough ========
    for tw in availTws:
        bgTwLen = tw[1] - tw[0]
        if (bgTwLen >= twLen):
            startTws.append([tw[0], tw[0] + (bgTwLen - twLen)])

    return startTws

def getTWInsertFlexibility(
    bgTws:      "List of non-overlapping time windows as background (CANNOT be covered by inserting time window)" = None,
    epoch:      "The outer epoch of `bgTws` time windows" = None,
    tw:         "The time windows to be inserted" = None,
    tdRange:    "This is for time-dependent tw, which is the time that the tw can move advance/delay without \
                 changing length. Also, tdRange[0] <= 0 and tdRange[1] >= 0" = [None, None]
    ) -> "Given a list of non-overlapping time windows (and the outer epoch), \
          for a new time window to be inserted, find \
          1) Can it be inserted? \
          2) if it can be inserted, the maximum time the tw can be advanced/delayed; \
          3) if it can not be inserted, the minimum time advanced/delayed to make the tw can be inserted":

    # Check if the tw is overlapping with bgTws ===============================
    canInsertFlag = True
    insertedTW = [[epoch[0] - CONST_EPSILON, epoch[0]]]
    insertedTW.extend([i for i in bgTws])
    insertedTW.append(tw)
    insertedTW.append([epoch[1], epoch[1] + CONST_EPSILON])
    insertedTW = sortTWs(insertedTW)
    if (insertedTW is None):
        canInsertFlag = False

    # Calculate flexibility ===================================================
    # If insert-able, find the gap that can fit the entire tw    
    if (canInsertFlag):
        # Get reversion of bgTWs, return the gap ------------------------------
        reverseTW = reverseTWs(bgTws, epoch)
        gapInserted = None
        for gap in reverseTW:
            if (gap[0] <= tw[0] and tw[1] <= gap[1]):
                gapInserted = gap
        if (gapInserted is None):
            print("ERROR: Cannot find inserting gap.")
            return None
        # Calculate allowed advancing/delaying time ---------------------------
        maxAdvanceAllowed = 0
        if (tdRange[0] is None):
            maxAdvanceAllowed = tw[0] - gapInserted[0]
        else:
            # NOTICE: tdRange[0] <= 0
            maxAdvanceAllowed = min(-tdRange[0], tw[0] - gapInserted[0])
        maxDelayAllowed = 0
        if (tdRange[1] is None):
            maxDelayAllowed = gapInserted[1] - tw[1]
        else:
            maxDelayAllowed = min(tdRange[1], gapInserted[1] - tw[1])

        return {
            'canInsertFlag': True,
            'maxAdvanceAllowed': maxAdvanceAllowed,
            'maxDelayAllowed': maxDelayAllowed
        }

    # If not insert-able, find if there is a gap that can fit in
    else:
        # Get the reverse of bgTWs, return all candidate gaps -----------------
        # NOTICE: tdRange[0] <= 0
        candidateEpoch = None
        if (tdRange[0] is None and tdRange[1] is None):
            candidateEpoch = [i for i in epoch]
        elif (tdRange[0] is not None and tdRange[1] is None):
            candidateEpoch = [tw[0] + tdRange[0], epoch[1]]
        elif (tdRange[0] is None and tdRange[1] is not None):
            candidateEpoch = [epoch[0], tw[1] + tdRange[1]]
        elif (tdRange[0] is not None and tdRange[1] is not None):
            candidateEpoch = [tw[0] + tdRange[0], tw[1] + tdRange[1]]
        reverseTW = reverseTWs(bgTws, candidateEpoch)

        # Find gaps within the candidateEpoch, and calculate the minimum time needed for inserting the tw by advancing/delaying
        minDelayNeeded = None
        minAdvanceNeeded = None
        for gap in reverseTW:
            # First, the gap should be large enough for `tw`
            if ((gap[1] - gap[0]) > (tw[1] - tw[0])):
                # If the tw is (entirely or partially) on the left side of the gap, the tw can fit in by delaying
                if (tw[0] < gap[0]):
                    if (minDelayNeeded == None or minDelayNeeded > (gap[0] - tw[0])):
                        minDelayNeeded = gap[0] - tw[0]
                # Else if the tw is (entirely or partially) on the right side of the gap, the tw can fit in by advancing
                elif (gap[1] < tw[1]):
                    if (minAdvanceNeeded == None or minAdvanceNeeded > (tw[1] - gap[1])):
                        minAdvanceNeeded = tw[1] - gap[1]

        return {
            'canInsertFlag': False,
            'minAdvanceNeeded': minAdvanceNeeded,
            'minDelayNeeded': minDelayNeeded
        }
