import math

from .common import *
from .const import *

def sortTimeWindows(
    unsortedTW: "List of unsorted time windows" = None,
    ) -> "Given a list of time windows, returns None if the TWs are overlapping returns a sorted list of time windows":

    # Example =================================================================
    # In: -> unsortedTW = [[17.8, 20], [12, 14], [21, None], [17, 17.5], [None, 10]]
    #     -> sortTimeWindows(unsortedTW)
    # Out:-> [[None, 10], [12, 14], [17, 17.5], [17.8, 20], [21, None]]

    # Initialize ==============================================================
    sortedTW = []

    # Separate infinite time windows ==========================================
    leftInfTW = None
    rightInfTW = None
    bothInfTW = None
    for tw in unsortedTW:
        if (tw[0] == None and tw[1] == None):
            if (bothInfTW != None):
                return None
            else:
                bothInfTW = tw # tw = [None, None]
        if (tw[0] == None and tw[1] != None):
            if (leftInfTW != None):
                return None
            else:
                leftInfTW = tw
        if (tw[0] != None and tw[1] == None):
            if (rightInfTW != None):
                return None
            else:
                rightInfTW = tw

    # Quick filter for bothInfTW ==============================================
    if (bothInfTW != None):
        if (len(unsortedTW) > 1):
            return None
        else:
            return [bothInfTW]

    # Sort start/end time of each finite time windows =========================
    # Get the list of start times and end times and all timestamps
    startTimeList = []
    endTimeList = []
    for tw in unsortedTW:
        if (tw[0] != None and tw[1] != None):
            startTimeList.append(tw[0])
            endTimeList.append(tw[1])
    sortedStartIndexList = sorted(range(len(startTimeList)), key=lambda k: startTimeList[k])
    sortedEndIndexList = sorted(range(len(endTimeList)), key=lambda k: endTimeList[k])

    # Check overlapping =======================================================
    if (sortedStartIndexList != sortedEndIndexList):
        return None

    # Sort time windows =======================================================
    for i in range(len(sortedStartIndexList)):
        sortedTW.append([startTimeList[sortedStartIndexList[i]], endTimeList[sortedStartIndexList[i]]])

    # Check overlapping =======================================================
    for i in range(len(sortedTW) - 1):
        if (sortedTW[i][1] > sortedTW[i + 1][0]):
            return None 

    # Add back infinite time windows ==========================================
    if (leftInfTW != None):
        if (len(sortedTW) > 0):
            if (leftInfTW[1] > sortedTW[0][0]):
                return None
            else:
                sortedTW.insert(0, leftInfTW)
        else:
            sortedTW.insert(0, leftInfTW)
    if (rightInfTW != None):
        if (len(sortedTW) > 0):
            if (rightInfTW[0] < sortedTW[-1][1]):
                return None
            else:
                sortedTW.append(rightInfTW)
        else:
            sortedTW.append(rightInfTW)

    return sortedTW

def mergeTimeWindows(
    tws:        "List of time windows, might be overlapping" = None,
    ) -> "Given a list of time windows, which might be overlapping, merge the time windows that are overlapping":

    # Example =================================================================
    # In: -> tws = [[31, None], [17.5, 17.8], [15.6, 17.6], [2, 14], [None, 3]]
    #     -> mergeTimeWindows(tws)
    # Out:-> [[None, 14], [15.6, 17.8], [31, None]]

    # Initialize ==============================================================
    mergedTWs = []

    # Separate infinite time windows ==========================================
    leftInfTW = None
    rightInfTW = None
    bothInfTW = None
    for tw in tws:
        if (tw[0] == None and tw[1] == None):
            bothInfTW = tw # tw = [None, None]
        if (tw[0] == None and tw[1] != None):
            if (leftInfTW != None):
                return None
            else:
                leftInfTW = tw
        if (tw[0] != None and tw[1] == None):
            if (rightInfTW != None):
                return None
            else:
                rightInfTW = tw

    # Quick filter for bothInfTW ==============================================
    if (bothInfTW != None):
        return [bothInfTW]

    # Merge finite tws ========================================================
    for tw in tws:
        # If the tw is finite
        if (tw[0] != None and tw[1] != None):
            # Find the index of tw that involves the s/e of tw, if there is any
            sIndex = None
            eIndex = None
            for bgTwIndex in range(len(mergedTWs)):
                if (insideInterval(tw[0], mergedTWs[bgTwIndex])):
                    sIndex = bgTwIndex
                if (insideInterval(tw[1], mergedTWs[bgTwIndex])):
                    eIndex = bgTwIndex
                if (sIndex != None and eIndex != None):
                    break
            # Find the s/e of merged tw
            ts = None
            te = None
            if (sIndex != None):
                ts = mergedTWs[sIndex][0]
            else:
                ts = tw[0]
            if (eIndex != None):
                te = mergedTWs[eIndex][1]
            else:
                te = tw[1]
            mergedTWs.append([ts, te])

            # Remove the tws that partially covered by [ts, te], remove them
            tmpMergedTWs = [i for i in mergedTWs]
            mergedTWs = [[ts, te]] # merged tw should be included
            for newTW in tmpMergedTWs:
                # For existed tws, both s and e should not inside the [ts, te]
                if (not insideInterval(newTW[0], [ts, te]) and not insideInterval(newTW[1], [ts, te])):
                    mergedTWs.append(newTW)
            
    # Process the infinite tws ================================================
    if (leftInfTW != None):
        updatedRight = leftInfTW[1]
        for tw in mergedTWs:
            if (insideInterval(updatedRight, tw)):
                updatedRight = tw[1]
        if (updatedRight != None):
            mergedTWs = [i for i in mergedTWs if i[0] >= updatedRight]
            mergedTWs.insert(0, [None, updatedRight])
        else:
            return [[None, None]]

    if (rightInfTW != None):
        updatedLeft = rightInfTW[0]
        for tw in mergedTWs:
            if (insideInterval(updatedLeft, tw)):
                updatedLeft = tw[0]
        if (updatedLeft != None):
            mergedTWs = [i for i in mergedTWs if i[1] <= updatedLeft]
            mergedTWs.append([updatedLeft, None])
        else:
            return [[None, None]]

    # Sort tws ================================================================
    mergedTWs = sortTimeWindows(mergedTWs)

    return mergedTWs

def lenTWOverlap(
    tw1:        "The first time window" = None,
    tw2:        "The second time window" = None
    ) -> "Given two time windows, returns the length of time that these two time windows are overlapping":

    # Initialize ==============================================================
    overlapTime = 0

    return overlapTime

def reverseTimeWindows(
    tws:        "List of existing time windows" = None,
    epoch:      "Another time window that we wanna find opposite time windows of `tws` inside `epoch`" = None
    ) -> "Given a list of time windows `tws`, given another time window `epoch`, returns the opposite time windows of `tws` in `epoch":

    # Initialize ==============================================================
    oppTws = [[epoch[0], epoch[1]]]

    # Reverse time windows ====================================================
    # Start with epoch, using `blockTimeWindows` to reverse time windows
    for tw in tws:
        oppTws = blockTimeWindows(oppTws, tw)
    return oppTws

def blockTimeWindows(
    tws:        "List of existing available time windows" = None,
    block:      "A piece of time make part of tws unavailable" = None,
    ) -> "Given a list of non-overlapping time windows, use a time window to block the existing `tws`":
    
    # Initialize ==============================================================
    newTWs = []

    # For finite block ========================================================
    if (block[0] != None and block[1] != None):
        for tw in tws:
            # Case 1: tw is finite
            if (tw[0] != None and tw[1] != None):
                # Case 1.1: block entirely inside tw
                if (insideInterval(block[0], tw) and insideInterval(block[1], tw)):
                    newTWs.append([tw[0], block[0]])
                    newTWs.append([block[1], tw[1]])
                # Case 1.2: Partially overlapping on the left of tw
                elif (not insideInterval(block[0], tw) and insideInterval(block[1], tw)):
                    newTWs.append([tw[1], block[1]])
                # Case 1.3: Partially overlapping on the right of tw
                elif (insideInterval(block[0], tw) and not insideInterval(block[1], tw)):
                    newTWs.append([tw[0], block[0]])
                # Case 1.4: No overlapping
                elif (not insideInterval(block[0], tw) and not insideInterval(block[1], tw)
                    and not insideInterval(tw[0], block) and not insideInterval(tw[1], block)):
                    newTWs.append(tw)
                # Case 1.5: block entirely cover tw
                elif (insideInterval(tw[0], block) and insideInterval(tw[1], block)):
                    pass
            # Case 2: tw is infinite from left
            elif (tw[0] == None and tw[1] != None):
                # Case 2.1: block entirely inside tw
                if (insideInterval(block[0], tw) and insideInterval(block[1], tw)):
                    newTWs.append([None, block[0]])
                    newTWs.append([block[1], tw[1]])
                # Case 2.2: partially overlapping
                elif (insideInterval(block[0], tw) and not insideInterval(block[1], tw)):
                    newTWs.append([None, block[0]])
                # Case 2.3: No overlapping
                elif (not insideInterval(block[0], tw) and not insideInterval(block[1], tw)):
                    pass
            # Case 3: tw is infinite from right
            elif (tw[0] != None and tw[1] == None):
                # Case 3.1: block entirely inside tw
                if (insideInterval(block[0], tw) and insideInterval(block[1], tw)):
                    newTWs.append([tw[0], block[0]])
                    newTWs.append([block[1], None])
                # Case 3.2: Partially overlapping on the right of tw
                elif (not insideInterval(block[0], tw) and insideInterval(block[1], tw)):
                    newTWs.append([block[1], None])
                # Case 3.3: No overlapping
                elif (not insideInterval(block[0], tw) and not insideInterval(block[1], tw)):
                    pass
            # Case 4: tw is [-inf, inf]
            elif (tw[0] == None and tw[1] == None):
                newTWs.append([None, block[0]])
                newTWs.append([block[1], None])

    # Truncate left side ======================================================
    elif (block[0] == None and block[1] != None):
        for tw in tws:
            # Case 1: tw is finite
            if (tw[0] != None and tw[1] != None):
                # Case 1.1: tw is entirely clear
                if (not insideInterval(tw[0], block) and not insideInterval(tw[1], block)):
                    newTWs.append(tw)
                # Case 1.2: tw is half blocked
                elif (insideInterval(tw[0], block) and not insideInterval(tw[1], block)):
                    newTWs.append([block[1], tw[1]])
                # Case 1.3: tw is entirely blocked
                elif (insideInterval(tw[0], block) and insideInterval(tw[1], block)):
                    pass
            # Case 2: tw is infinite from left
            elif (tw[0] == None and tw[1] != None):
                # Case 2.1: tw is half blocked
                if (not insideInterval(tw[1], block)):
                    newTWs.append([block[1], tw[1]])
                # Case 2.2: tw is entirely blocked
                elif (insideInterval(tw[1], block)):
                    pass
            # Case 3: tw is infinite from right
            elif (tw[0] != None and tw[1] == None):
                # Case 3.1: tw is entirely clear
                if (not insideInterval(tw[0], block)):
                    newTWs.append(tw)
                # Case 3.2: tw is partially blocked
                elif (insideInterval(tw[1], block)):
                    newTWs.append([block[1], None])
            # Case 4: tw is infinite from both
            elif (tw[0] == None and tw[1] == None):
                newTWs.append([block[1], None])

    # Truncate right side =====================================================
    elif (block[0] != None and block[1] == None):
        for tw in tws:
            # Case 1: tw is finite
            if (tw[0] != None and tw[1] != None):
                # Case 1.1: tw is entirely clear
                if (not insideInterval(tw[0], block) and not insideInterval(tw[1], block)):
                    newTWs.append(tw)
                # Case 1.2: tw is half blocked
                elif (not insideInterval(tw[0], block)):
                    newTWs.append([tw[0], block[0]])
                # Case 1.3: tw is entirely blocked
                elif (insideInterval([tw[0], block])):
                    pass
            # Case 2: tw is infinite from left side
            elif (tw[0] == None and tw[1] != None):
                # Case 2.1: tw is half blocked
                if (insideInterval(tw[1], block)):
                    newTWs.append([None, block[0]])
                elif (not insideInterval(tw[1], block)):
                    newTWs.append(tw)
            # Case 3: tw is infinite from right side
            elif (tw[0] != None and tw[1] == None):
                # Case 1: tw is half blocked
                if (not insideInterval(tw[0], block)):
                    newTWs.append([tw[0], block[0]])
                # Case 2: tw is entirely blocked
                if (insideInterval(tw[1], block)):
                    pass
            elif (tw[0] == None and tw[1] == None):
                newTWs.append([None, block[0]])

    # Entirely blocked ========================================================
    else:
        return [[]]

    return newTWs

def getIndexInTimeWindow(
    bgTws:      "List of non-overlapping time windows as background (CANNOT be covered by inserting time window)" = None,
    t:          "Time stamp" = None
    ) -> "Given a time stamp and a list of time windows, check if the time stamp is in time windows":
    for i in range(len(bgTws)):
        if (insideInterval(t, bgTws[i])):
            return i
    return None

def getEarliestNextAvailTime(
    bgTws:      "List of non-overlapping time windows as background (CANNOT be covered by inserting time window)" = None,
    t:          "Time stamp" = None
    ) -> "Given a time stamp, get next launchable time":
    earliestNext = None
    for i in range(len(bgTws)):
        if (not insideInterval(t, bgTws[i])):
            if (bgTws[i][0] != None and t < bgTws[i][0]):
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
    availTws = reverseTimeWindows(bgTws, epoch)

    # For each available time windows, try to see if it is long enough ========
    for tw in availTws:
        # For finite tws
        if (tw[0] != None and tw[1] != None):
            bgTwLen = tw[1] - tw[0]
            if (bgTwLen >= twLen):
                startTws.append([tw[0], tw[0] + (bgTwLen - twLen)])
        # For left-infinite tws
        elif (tw[0] == None and tw[1] != None):
            startTws.append([None, tw[1] - twLen])
        # For right-infinite tws
        elif (tw[0] != None and tw[1] == None):
            startTws.append([tw[0], None])
        # For both-end-infinite tws
        elif (tw[0] == None and tw[1] == None):
            startTws.append([None, None])

    return startTws

def getTWInsertFlexibility(
    bgTws:      "List of non-overlapping time windows as background (CANNOT be covered by inserting time window)" = None,
    epoch:      "The outer epoch of `bgTws` time windows" = [0, None],
    tw:         "The time windows to be inserted" = None,
    tdRange:    "This is for time-dependent tw, which is the time that the tw can move advance/delay without \
                 changing length. Also, tdRange[0] <= 0 and tdRange[1] >= 0" = [None, None]
    ) -> "Given a list of non-overlapping time windows (and the outer epoch), \
          for a new time window to be inserted, find \
          1) if it can be inserted; \
          2) if it can be inserted, the maximum time the tw can be advanced/delayed; \
          3) if it can not be inserted, the minimum time advanced/delayed to make the tw can be inserted":

    # Check if the tw is overlapping with bgTws ===============================
    canInsertFlag = True

    insertedTW.extend([i for i in bgTws])
    insertedTW.append(tw)
    insertedTW.append([epoch[1], epoch[1] + EPSILON])
    insertedTW = sortTimeWindows(insertedTW)
    if (insertedTW is None):
        canInsertFlag = False

    # Calculate flexibility ===================================================
    # If insert-able, find the gap that can fit the entire tw    
    if (canInsertFlag):
        # Get reversion of bgTWs, return the gap ------------------------------
        reverseTW = reverseTimeWindows(bgTws, epoch)
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
        reverseTW = reverseTimeWindows(bgTws, candidateEpoch)

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

