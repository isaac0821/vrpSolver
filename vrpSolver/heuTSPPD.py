import heapq
import math
import datetime
import warnings

from .const import *
from .common import *
from .graph import *
from .geometry import *
from .msg import *
from .operator import *
from .calculate import *
from .plot import *

def heuTSPPD():
    return

def consTSPPD(
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    reqs:       "Pickup and delivery pairs\
                    {\
                        reqID1: {\
                            'pickup': nodeID, \
                            'delivery': nodeID, \
                            'size': size, \
                            'deadline': deadline}, \
                        ... \
                    }" = None,
    edges:      "1) String (default) 'Euclidean' or \
                 2) String 'LatLon' or \
                 3) Dictionary {(nodeID1, nodeID2): dist, ...} or \
                 4) String 'Grid', will need to add arguments using `edgeArgs`" = "Euclidean",
    edgeArgs:   "If choose 'Grid' as tau option, we need to provide the following dictionary \
                    {\
                        'colRow': (numCol, numRow),\
                        'barriers': [(coordX, coordY), ...], \
                    }" = None,
    depotID:    "DepotID, default to be 0" = 0,
    capacity:   "Capacity of loads" = 3,
    serviceTime: "Service time spent on each customer (will be added into travel matrix)" = 0,
    ) -> "TSP with PD, DDL, capacity":

    # Algorithm outline =======================================================
    # Structures:
    # - `backwardStack`: Stack, store "next feasible action" ordered by distance 
    #    to previous action. Each element is in format of (preAction, availAction).
    #    for the same `preAction`, `availAction` on top of stack is closer in distance
    # - `forwardStack`: Stack, store the best feasible sequences actions so far.
    # Algorithm:
    # Step 0: Initialize
    #     Step 0.1: backwardStack = [], forwardStack = []
    #     Step 0.2: Get an initial location, e.g., depot, and find the closest pickup 
    #               action as starting point, and push it into `forwardStack`
    # Step 1: A = forwardStack.peek()
    #     Step 1.1: If (len(forwardStack)) == 2 * requests, Return `forwardStack`
    # Step 2: Find all available action for A: lstNextA = availNextActions()
    #     Step 2.1: If (len(lstNextA)) == 0, dead end, try rolling back
    #         Step 2.1.1: If (len(backwardStack) > 0): go to Step 3
    #         Step 2.2.1: If (len(backwardStack) == 0): Return infeasible
    #     Step 2.2: If (len(lstNextA)) > 0, push the closest action to `forwardStack`
    #               push all the rest of actions to `backwardStack`, A is `preAction`
    #               go to Step 1
    # Step 3: Rolling back
    #     Step 3.1: ' while (forwardStack.peek() != backwardStack.peek().preAction):
    #               '     forwardStack.pop()
    #     Step 3.2: B = backwardStack.pop(), forwardStack.push(B), go to 1.

    # Initialize ==============================================================
    # In backward stack, it has all candidate next steps
    # backwardStack[0]  - stack bottom
    # backwardStack[-1] - stack top
    backwardStack = []

    # In forward stack, it has the best actions so far    
    # forwardStack[0]  - stack bottom
    # forwardStack[-1] - stack top
    forwardStack = []
    traveled = []

    # Given a current location find all available next steps ==================
    # The "current" location is based on `forwardStack`
    def availNextActions():
        # Initialize, available action list -----------------------------------
        availActions = []

        # Initialize, find locations of items/requests ------------------------
        # lstPicked + lstNotYet = All
        # lstOnboard + lstDelivered = lstPicked
        # lstDelivered + lstUnfulfilled = All
        lstPicked = []
        lstNotYet = []
        lstOnboard = []
        lstDelivered = []
        lstUnfulfilled = []
        for act in forwardStack:
            if (act[0] == 'pickup'):
                lstPicked.append(act[1])
            elif (act[0] == 'deliver'):
                lstDelivered.append(act[1])
        for r in reqs:
            if (r not in lstPicked):
                lstNotYet.append(r)
            if (r not in lstDelivered):
                lstUnfulfilled.append(r)
        for r in lstPicked:
            if (r not in lstDelivered):
                lstOnboard.append(r)

        # Check if `req` will violate DDL when `act` is taken -----------------
        # Notice that `act` has not been taken yet, the last action "as of this moment"
        #     is the last action in "forwardStack"
        def checkBestCaseDDL(req, act):
            # True if the DDL will not be violated after `act`
            ddlFlag = False

            # Calculate the best case accumulated distances for request
            lastNodeID = forwardStack[-1][2]
            leastAccDist = traveled[-1]
            
            # Case 1: If delivered in `act`, leastAccDist is the distance from first
            #    pickup in the bundle to this delivery            
            if (act[1] == req and act[0] == 'deliver'):
                leastAccDist += tau[lastNodeID, reqs[req]['delivery']]
            # Case 2: If picked up in `act`, leastAccDist is the distance between pickup
            #    and delivery location of the request
            elif (act[1] == req and act[0] == 'pickup'):
                leastAccDist += tau[lastNodeID, reqs[req]['pickup']]
                leastAccDist += reqs[req]['travelDist']
            # Case 3: If picked up before `act`, not delivered in `act`, leastAccDist
            #    is the distance from first pickup in the bundle, to current loc, then
            #    to deliver loc
            elif (req in lstOnboard):
                leastAccDist += tau[lastNodeID,act[2]]
                leastAccDist += tau[act[2], reqs[req]['delivery']]
            # Case 4: If not picked up until `act`, leastAccDist is the distance from 
            #    first pickup in the bundle, to current loc, to pickup loc, to deliver loc
            elif (req in lstNotYet):
                leastAccDist += tau[lastNodeID,act[2]]
                leastAccDist += tau[act[2], reqs[req]['pickup']]
                leastAccDist += reqs[req]['travelDist']
            # Error: e.g. req has been delivered
            else:
                print("ERROR: Failed to check DDL constraints. Most likely trying to pickup/deliver some delivered request.")
                return

            # Check DDL
            ddlFlag = (leastAccDist < req['deadline'])
            return ddlFlag

        # Find current location and current load ------------------------------
        if (len(forwardStack) == 0):
            curNodeID = initNodeID
            curLoad = 0
        else:
            curNodeID = forwardStack[-1][2]
            curLoad = 0
            for i in range(len(forwardStack)):
                if (forwardStack[i][0] == 'pickup'):
                    curLoad += reqs[forwardStack[i][1]]['size']
                elif (forwardStack[i][0] == 'deliver'):
                    curLoad -= reqs[forwardStack[i][1]]['size']

        # Find nearest feasible next steps ------------------------------------
        # A heap, each item in format of (dist, ('pickup/delivery', reqID, nodeID))
        heapAvailActions = []
        # 1. Find pickup that will not jeopardize DDL, will not violate capacity
        for r in lstNotYet:
            tryAct = ('pickup', r, reqs[r]['pickup'])

            # Check capacity
            capacityFlag = (curLoad + reqs[r]['size'] <= capacity)

            # Check DDL - for every unfulfilled requests, cannot violate ddl
            ddlFlag = False
            if (capacityFlag):
                ddlFlag = True
                for unfulfilled in lstUnfulfilled:
                    if (not checkBestCaseDDL(unfulfilled, tryAct)):
                        ddlFlag = False
                        break

            # If this pickup action is doable, add it
            if (ddlFlag):
                d = tau[curNodeID, tryAct[2]]
                heapq.heappush(heapAvailActions, (d, tryAct))

        # 2. Find delivery that will not jeopardize DDL
        for r in lstOnboard:
            tryAct = ('deliver', r, reqs[r]['delivery'])

            # Check DDL - for every unfulfilled requests
            ddlFlag = True
            for unfulfilled in lstUnfulfilled:
                if (not checkBestCaseDDL(unfulfilled, tryAct)):
                    ddlFlag = False
                    break

            # If this deliver action is doable, add it
            if (ddlFlag):
                d = tau[curNodeID, tryAct[2]]
                heapq.heappush(heapAvailActions, (d, tryAct))

        # Sorted actions (from furthest to nearest) ---------------------------
        while (len(heapAvailActions) > 0):
            availActions.insert(0, heapq.heappop(heapAvailActions)[1])

        return {
            'availActions': availActions
        }

    # Initial action ==========================================================
    if (depot != None):
        initialClosest = []
        for r in reqs:
            d = tau[depot, reqs[r]['pickup']]
            heapq.heappush(initialClosest, (d, ('pickup', r, reqs[r]['pickup'])))
        forwardStack.append(heapq.heappop(initialClosest)[1])
        # Clear memory
        initialClosest = []
    else:
        if (firstReqID == None):
            return {
                'dist': None,
                'actions': None
            }
        forwardStack.append(('pickup', firstReqID, reqs[firstReqID]['pickup']))

    # Depth-search + nearest neighborhood =====================================
    while (len(forwardStack) < 2 * len(reqs)):
        # push available next actions into backward stack, the top of stack will
        #     will be the nearest action that is feasible
        lstNextActions = availNextActions()['availActions']
        # If at current location, there still feasible next steps, continue search child tree
        if (len(lstNextActions) > 0):
            for act in lstNextActions:
                backwardStack.append((forwardStack[-1], act))
            nextAction = backwardStack.pop()[1]
            forwardStack.append(nextAction)
        # If at current location, there is no where to go, trace back to previous location
        elif (len (lstNextActions) == 0):
            # If no where to trace back, return infeasible
            if (len(backwardStack) == 0):
                return {
                    'dist': None,
                    'actions': None
                }
            else:
                # If forwardStack.peek() is not the `preAction` of backwardStack.peek()
                # Pop until we find correct node that has feasible next step
                while (forwardStack[-1][0] != backwardStack[-1][0][0] 
                    or forwardStack[-1][1] != backwardStack[-1][0][1]):
                    forwardStack.pop()
                # Trace back to where there are feasible next steps
                newAction = backwardStack.pop()[1]
                forwardStack.append(newAction)
    dist = evalAccDist(reqs, forwardStack, curTime)['accDist']

    return {
        'dist': dist,
        'actions': forwardStack
    }

def localTSPPD(
    reqs:       "Dictionary, returns the pairs of pickup and delivery node IDs, \
                    {\
                        reqID1: {\
                            'pickup': nodeID, \
                            'delivery': nodeID, \
                            'size': size, \
                            'arrivalTime': arrivalTime}, \
                        ... \
                    }" = None,
    initActions:"Sequence of actions that comes from constructive heuristic, needs to be feasible" = [],
    curTime:    "Current time, unit: hour" = None,
    ) -> "TSPPD with DDL, capacity":

    # Initialize ==============================================================
    actions = [i for i in initActions]
    newDist = evalAccDist(reqs, actions, curTime)['accDist']
    if (newDist == None):
        return {
            'dist': None,
            'actions': None
        }

    # Swapping operation ======================================================
    # For an action (`act`), find the best actions, if any, that can swap with it to get improvement
    def swap2Locs(actions):
        # Initialize ----------------------------------------------------------
        currentBestMove = None
        currentBestDist = evalAccDist(reqs, actions, curTime)['accDist']
        if (currentBestDist == None):
            print("Swap2Loc failed", reqs, actions, curTime)

        # O(n^2) swapping -----------------------------------------------------
        # Try to swap action i with action j to see if there is improvement
        for i in range(len(actions) - 1):
            # If action i is 'pickup', it cannot be swapped with action after 'delivery'
            # If the swapping cause infeasibility of other actions, it will be detected in `evalAccDist()`
            if (actions[i][0] == 'pickup'):
                for j in range(i, len(actions)):
                    if (actions[j][0] == 'deliver' and actions[j][1] == actions[i][1]):
                        break
                    else:
                        swappedActions = [i for i in actions]
                        swappedActions[i], swappedActions[j] = swappedActions[j], swappedActions[i]
                        swappedDist = evalAccDist(reqs, swappedActions, curTime)['accDist']
                        # If swapped is feasible
                        if (swappedDist != None and swappedDist < currentBestDist):
                            currentBestDist = swappedDist
                            currentBestMove = (i, j)
            # If action i is 'delivery', it cannot be swapped with action before 'pickup'
            # However, we can safely try to swap with the rest of seq because 'pickup' should be 
            #     in front of this action.
            elif (actions[i][0] == 'deliver'):
                for j in range(i, len(actions)):
                    swappedActions = [i for i in actions]
                    swappedActions[i], swappedActions[j] = swappedActions[j], swappedActions[i]
                    swappedDist = evalAccDist(reqs, swappedActions, curTime)['accDist']
                    # If swapped is feasible
                    if (swappedDist != None and swappedDist < currentBestDist):
                        currentBestDist = swappedDist
                        currentBestMove = (i, j)        

        return {
            'newSwap': currentBestMove,
            'newDist': currentBestDist
        }

    # 2-opt operation =========================================================
    # NOTICE: Need to have at least 4 nodes/stops/actions to try 2-opts
    def revert(i, j, actions):
        revertAction = []
        revertAction.extend([actions[k] for k in range(i)])
        revertAction.extend([actions[j - k] for k in range(j - i + 1)])
        revertAction.extend([actions[k] for k in range(j + 1, len(actions))])
        return revertAction

    def linKernighan(actions):
        # Initialize ----------------------------------------------------------
        currentBestMove = None
        currentBestDist = evalAccDist(reqs, actions, curTime)['accDist']

        # O(n^2) 2-opt --------------------------------------------------------
        # First arc: (i, i + 1)
        # Second arc: (j, j + 1)
        # New actions: (start -> i) -> (j -> i + 1) -> j + 1 -> end 
        for i in range(len(actions) - 3):
            for j in range(i + 2, len(actions) - 1):
                optedActions = revert(i, j, actions)
                optedDist = evalAccDist(reqs, optedActions, curTime)['accDist']
                # If 2-opt is feasible
                if (optedDist != None and optedDist < currentBestDist):
                    currentBestDist = optedDist
                    currentBestMove = (i, j)

        return {
            'new2Opt': currentBestMove,
            'newDist': currentBestDist
        }

    # Local improvement with different methods ================================
    improveFlag = True
    while (improveFlag):
        improveFlag = False

        # Try swap improvement
        newImprove = swap2Locs(actions)
        if (newImprove['newSwap'] != None):
            improveFlag = True
            (i, j) = newImprove['newSwap']
            newDist = newImprove['newDist']
            actions[i], actions[j] = actions[j], actions[i]

        # Try Lin-Kernighan  2-opt
        newImprove = linKernighan(actions)
        if (newImprove['new2Opt'] != None):
            improveFlag = True
            (i, j) = newImprove['new2Opt']
            newDist = newImprove['newDist']
            actions = revert(i, j, actions)

    return {
        'dist': newDist,
        'actions': actions
    }

def insertTSPPD(
    reqs:       "Dictionary, returns the pairs of pickup and delivery node IDs, \
                    {\
                        reqID1: {\
                            'pickup': nodeID, \
                            'delivery': nodeID, \
                            'size': size, \
                            'arrivalTime': arrivalTime}, \
                        ... \
                    }" = None,
    insertReq:        "Dictionary, returns the pairs of pickup and delivery node IDs, only has one request\
                    {\
                        reqID1: {\
                            'pickup': nodeID, \
                            'delivery': nodeID, \
                            'size': size, \
                            'arrivalTime': arrivalTime}, \
                        ... \
                    }" = None,
    actions:    "Sequence of actions" = [],
    curTime:    "Current time, the time of bundle/singleton generated (to be posted)" = None,
    ) -> "Insert `insertReq` into `actions` which is the route for `reqs`":

    insertReqID = list(insertReq.keys())[0]

    # First try to insert to existing route ===================================
    # (i, j), positions for inserted pickup and delivery
    # NOTE: position j is BEFORE inserting i. So when inserting delivery position, it should be j + 1
    currentBestMove = None
    currentBestDist = None
    for i in range(len(actions) - 1):
        for j in range(i + 1, len(actions)):
            tryActions = [i for i in actions]
            tryActions.insert(i, ('pickup', insertReqID, insertReq[insertReqID]['pickup']))
            tryActions.insert(j + 1, ('deliver', insertReqID, insertReq[insertReqID]['delivery']))
            tryActionDist = evalAccDist({**reqs, **insertReq}, tryActions, curTime)['accDist']
            if (tryActionDist != None and (currentBestDist == None or tryActionDist < currentBestDist)):
                currentBestDist = tryActionDist
                currentBestMove = (i, j)

    # If insertion route is feasible, return inserted route ===================
    if (currentBestMove == None):
        return {
            'dist': None,
            'actions': None,
            'potential': None,
            'reliableLevel': None,
            'expireTime': None,
            'deliverTime': None
        }

    # Local improve this new route ============================================
    insertedActions = [i for i in actions]
    (i, j) = currentBestMove
    insertedActions.insert(i, ('pickup', insertReqID, insertReq[insertReqID]['pickup']))
    insertedActions.insert(j + 1, ('deliver', insertReqID, insertReq[insertReqID]['delivery']))    
    localReq = localTSPPD(
        reqs = {**reqs, **insertReq},
        initActions = insertedActions,
        curTime = curTime)

    # Calculate expire time ===================================================
    remainHours = {}
    deliverTime = {}
    if (localReq['dist'] != None):
        newReqs = {**reqs, **insertReq}
        travelDist = evalAccDist(newReqs, localReq['actions'], curTime)['travelDist']
        for r in newReqs:
            deliverTime[r] = curTime + int(24 * travelDist[r] / config["SETTING_MILESPERDAY"])
            remainHours[r] = slaCheckRemainTime(
                arrivalTime = newReqs[r]['arrivalTime'], 
                curTime = curTime, 
                reqDist = newReqs[r]['travelDist'],
                travelDist = travelDist[r])['remainHours']
    expireTime = curTime + min(remainHours.values())

    if (localReq['dist'] == None):
        return {
            'dist': None,
            'actions': None,
            'potential': None,
            'reliableLevel': None,
            'expireTime': None,
            'deliverTime': None
        }
    else:
        potential = forecastBundlePotential(insertedActions)['potential']
        return {
            'dist': currentBestDist,
            'actions': insertedActions,
            'potential': potential,
            'reliableLevel': 'Ins',
            'expireTime': expireTime,
            'deliverTime': deliverTime
        }

def removeTSPPD(
    reqs:       "Dictionary, returns the pairs of pickup and delivery node IDs, \
                    {\
                        reqID1: {\
                            'pickup': nodeID, \
                            'delivery': nodeID, \
                            'size': size, \
                            'arrivalTime': arrivalTime}, \
                        ... \
                    }" = None,
    removeReqID:"RequestID of the request gets removed" = None,
    actions:    "Sequence of actions" = [],
    curTime:    "Current time, the time of bundle/singleton generated (to be posted)" = None,
    ) -> "Remove a request from a bundle and update actions":

    # Validation ==============================================================
    if (removeReqID not in reqs):
        print("ERROR: removeReqID not in reqs")
        return

    # If there are only two requests, there will be one singleton left ========
    if (len(reqs) == 2):
        restReqID = [r for r in reqs if r != removeReqID][0]
        dist = reqs[restReqID]['travelDist']
        actions = [
            ('pickup', restReqID, reqs[restReqID]['pickup']),
            ('deliver', restReqID, reqs[restReqID]['delivery'])
        ]
        remainHours = slaCheckRemainTime(
            arrivalTime = reqs[restReqID]['arrivalTime'], 
            curTime = curTime, 
            reqDist = reqs[restReqID]['travelDist'],
            travelDist = reqs[restReqID]['travelDist'])['remainHours']
        return {
            'dist': dist,
            'actions': actions,
            'reliableLevel': 'Ins',
            'expireTime': curTime + remainHours
        }
    elif (len(reqs) <= 1):
        print("ERROR: For `removeTSPPD`, number of reqs should >= 2")
        return {
            'dist': None,
            'actions': None,
            'potential': None,
            'reliableLevel': None,
            'expireTime': None,
            'deliverTime': None
        }

    # First remove this request, it will give another feasible route
    removedActions = [i for i in actions if (i[1] != removeReqID)]
    reqsWithoutR = {r: reqs[r] for r in reqs if r != removeReqID}
    removedDist = evalAccDist(reqsWithoutR, removedActions, curTime)['accDist']
    if (removedDist == None):
        return {
            'dist': None,
            'actions': None,
            'potential': None,
            'reliableLevel': None,
            'expireTime': None,
            'deliverTime': None
        }

    # Local improve this new route
    localReq = localTSPPD(
        reqs = reqsWithoutR,
        initActions = removedActions,
        curTime = curTime)

    # Check if it is a valid bundle ===========================================
    beneficialFlag = evalBeneficial(reqsWithoutR, localReq['actions'], curTime)['beneficialFlag']
    if (not beneficialFlag):
        return {
            'dist': None,
            'actions': None,
            'potential': None,
            'reliableLevel': None,
            'expireTime': None,
            'deliverTime': None
        }

    # Calculate expire time ===================================================
    remainHours = {}
    deliverTime = {}
    if (beneficialFlag):
        travelDist = evalAccDist(reqsWithoutR, localReq['actions'], curTime)['travelDist']        
        for r in reqsWithoutR:
            deliverTime[r] = curTime + int(24 * travelDist[r] / config["SETTING_MILESPERDAY"])
            remainHours[r] = slaCheckRemainTime(
                arrivalTime = reqsWithoutR[r]['arrivalTime'], 
                curTime = curTime, 
                reqDist = reqsWithoutR[r]['travelDist'],
                travelDist = travelDist[r])['remainHours']
    expireTime = curTime + min(remainHours.values())

    dist = localReq['dist']
    actions = localReq['actions']
    potential = forecastBundlePotential(actions)['potential']
    return {
        'dist': dist,
        'actions': actions,
        'potential': potential,
        'reliableLevel': 'Ins',
        'expireTime': expireTime,
        'deliverTime': deliverTime
    }
