import math

from .common import *
from .const import *

# [Done]
def sortTimeWindows(
    unsortedTW: "List of unsorted time windows" = None,
    ) -> "Given a list of time windows, returns None if the TWs are overlapping returns a sorted list of time windows":

    # Initialize ==============================================================
    sortedTW = []

    # Sort start/end time of each time windows ================================
    # Get the list of start times and end times and all timestamps
    startTimeList = []
    endTimeList = []
    for tw in unsortedTW:
        startTimeList.append(tw[0])
        endTimeList.append(tw[1])
    sortedStartIndexList = sorted(range(len(startTimeList)), key=lambda k: startTimeList[k])
    sortedEndIndexList = sorted(range(len(endTimeList)), key=lambda k: endTimeList[k])

    # Check overlapping =======================================================
    if (sortedStartIndexList != sortedEndIndexList):
        return None

    # Sort time windows =======================================================
    for i in range(len(sortedStartIndexList)):
        sortedTW.append(unsortedTW[sortedStartIndexList[i]])

    # Check overlapping again =================================================
    for i in range(len(sortedTW) - 1):
        if (sortedTW[i][1] > sortedTW[i + 1][0]):
            return None 

    return sortedTW

# [Done]
def reverseTimeWindows(
    tws:        "List of existing time windows" = None,
    epoch:      "Another time window that we wanna find opposite time windows of `tws` inside `epoch`" = None
    ) -> "Given a list of time windows `tws`, given another time window `epoch`, returns the opposite time windows of `tws` in `epoch":

    # Initialize ==============================================================
    oppTws = [[i for i in epoch]]

    # Reverse time windows ====================================================
    # Start with epoch, using `blockTimeWindows` to reverse time windows
    for tw in tws:
        oppTws = blockTimeWindows(oppTws, tw)
    return oppTws

# [Done]
def blockTimeWindows(
    tws:        "List of existing time windows" = None,
    block:      "A piece of time make part of tws unavailable" = None,
    ) -> "Given a list of non-overlapping time windows, use a time window to block the existing `tws`":
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

# [Done]
def getIndexInTimeWindow(
    bgTws:      "List of non-overlapping time windows as background (CANNOT be covered by inserting time window)" = None,
    t:          "Time stamp" = None
    ) -> "Given a time stamp and a list of time windows, check if the time stamp is in time windows":
    for i in range(len(bgTws)):
        if (insideInterval(t, bgTws[i])):
            return i
    return None

# [Done]
def getEarliestNextAvailTime(
    bgTws:      "List of non-overlapping time windows as background (CANNOT be covered by inserting time window)" = None,
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

# [Done] NOTICE: this is for sorted and non-overlapping time windows
def getTWInsertFlexibility(
    bgTws:      "List of non-overlapping time windows as background (CANNOT be covered by inserting time window)" = None,
    epoch:      "The outer epoch of `bgTws` time windows" = None,
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
    insertedTW = [[epoch[0] - EPSILON, epoch[0]]]
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
