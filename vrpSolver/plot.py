import matplotlib.pyplot as plt
import random

from .msg import *
from .color import *

def plotGantt(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    gantt:      "List of dictionaries, in the following format\
                    [{\
                        'entityID': entityID, \
                        'timeWindow': [startTime, endTime], \
                        'desc': (optional) description of the window,\
                        'trackID': (optional, default as '0') parallel gantt chart ID, \
                        'color': (optional, default as 'Random') color, \
                        'style': (optional, default as 'solid') 'solid' \
                    }, ... , \
                    {\
                        'entityID': entityID, \
                        'timeStamps': [timeStamp1, timeStamp2, ..., timeStampN], \
                        'desc': [List of descriptions, correspond to `timeStamps`],\
                        'trackID': (optional, default as '0') parallel gantt chart ID, \
                        'color': (optional, default as 'Random') color, \
                        'style': (optional, default as 'solid') 'solid' \
                    }, ... , \
                    {\
                        'entityID': entityID, \
                        'timePin': timePin, \
                        'desc': desc, \
                        'color': (optional, default as 'Random') color \
                    }]\
                " = None,
    gridFlag:   "True if turn on the grid as background" = True,
    labelFlag:  "True if add label of entities on Gantt chart" = True,
    entities:   "The Gantt chart will be drawn in this order, if None, takes the entities in `gantt`" = None,
    forceStartEndFlag: "Force the start time and end time of gantt" = False,
    startTime:  "Start time of Gantt, default to be 0, if None, use the earliest time in `gantt`" = 0,
    endTime:    "End time of Gantt, default to be None, if None, use the latest time in `gantt`" = None,
    width:      "Width of the figure" = 12,
    height:     "Height of the figure" = 5,
    pinWidth:   "Width for pins" = 2,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Given a Gantt dictionary, plot Gantt":

    # Check for required fields ===============================================
    if (gantt == None):
        print(ERROR_MISSING_GANTT)
        return

    # Pre-calculate ===========================================================
    realStart = None
    realEnd = None
    entList = []
    for g in gantt:
        if (g['entityID'] not in entList):
            entList.append(g['entityID'])
        if ('timeWindow' in g):
            if (realStart == None or realStart > g['timeWindow'][0]):
                realStart = g['timeWindow'][0]
            if (realEnd == None or realEnd < g['timeWindow'][1]):
                realEnd = g['timeWindow'][1]
        elif ('timeStamps' in g):
            if (realStart == None or realStart > g['timeStamps'][0]):
                realStart = g['timeStamps'][0]
            if (realEnd == None or realEnd < g['timeStamps'][-1]):
                realEnd = g['timeStamps'][-1]
    if (entities != None):
        for g in entities:
            if (g not in entList):
                print(ERROR_INCOR_GANTT_MISSENT)
                return
    else:
        entities = [i for i in entList]

    # Check overwritten fields ================================================
    if (startTime != None):
        if (startTime > realStart and not forceStartEndFlag):
            print(WARNING_OVERWRITE_STARTTIME)
            startTime = realStart
    else:
        startTime = realStart

    if (endTime != None):
        if (endTime < realEnd and not forceStartEndFlag):
            print(WARNING_OVERWRITE_ENDTIME)
            endTime = realEnd
    else:
        endTime = realEnd

    # Check for multiple tracks ===============================================
    numTrack = 0
    trackPos = {}
    enableMultipleTrackFlag = False
    for g in gantt:
        if (g['entityID'] in entities):
            if ('trackID' in g):
                enableMultipleTrackFlag = True
                break
    if (enableMultipleTrackFlag):
        for g in gantt:
            if (g['entityID'] in entities):
                if (g['entityID'] not in trackPos):
                    trackPos[g['entityID']] = {}
                if (g['trackID'] not in trackPos[g['entityID']]):
                    trackPos[g['entityID']][g['trackID']] = []
    else:
        for g in gantt:
            if (g['entityID'] in entities and g['entityID'] not in trackPos):
                trackPos[g['entityID']] = {g['entityID']: []}


    # If no based matplotlib figure, define fig size ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        fig.set_figheight(height)
        fig.set_figwidth(width)
        ax.set_xlim(startTime, endTime + (endTime - startTime) * 0.05)
        ax.set_ylim(0, len(entities) + 0.2)

    # Set axis ================================================================
    ax.set_yticks([i + 0.5 for i in range(len(entities))])
    entities.reverse()
    ax.set_yticklabels(entities)
    entities.reverse()
    ax.set_xlabel("Time")

    # Grids ===================================================================
    if (gridFlag):
        ax.grid(b = True, linestyle=':')

    # Loop through `gantt` and draw gantt =====================================
    for g in gantt:
        if (g['entityID'] in entities):
            

            bottom = len(entities) - entities.index(g['entityID']) - 0.9
            top = 0
            if (labelFlag == True):
                top = len(entities) - entities.index(g['entityID']) - 0.5
            else:
                top = len(entities) - entities.index(g['entityID']) - 0.25
            if ('timeWindow' in g):
                s = g['timeWindow'][0]
                e = g['timeWindow'][1]
                x = [s, s, e, e, s]
                y = [bottom, top, top, bottom, bottom]
                ax.plot(x, y, color = 'black', linewidth = 1)
                if (g['color'] != 'random'):
                    ax.fill(x, y, color = g['color'])
                else:
                    rndColor = colorRandom()
                    ax.fill(x, y, color = rndColor)
                if (labelFlag == True):
                    ax.annotate(g['desc'], (s, top + 0.1))
                if (g['style'] != 'solid'):
                    ax.fill(x, y, hatch = g['style'], fill=False)
            elif ('timeStamps' in g):
                for i in range(len(g['timeStamps']) - 1):
                    s = g['timeStamps'][i]
                    e = g['timeStamps'][i + 1]
                    x = [s, s, e, e, s]
                    y = [bottom, top, top, bottom, bottom]
                    ax.plot(x, y, color = 'black', linewidth = 1)
                    ax.annotate(g['desc'][i], (s, top + 0.1))
                    if (g['color'] != 'random'):
                        ax.fill(x, y, color = g['color'])
                    else:
                        rndColor = colorRandom()
                        ax.fill(x, y, color = rndColor)
                    if (g['style'] != 'solid'):
                        ax.fill(x, y, hatch = g['style'], fill=False)
                if (labelFlag == True):
                    ax.annotate(g['desc'][-1], (g['timeStamps'][-1], top + 0.1))
            elif ('timePin' in g):
                s = g['timePin'] - pinWidth
                m = g['timePin']
                e = g['timePin'] + pinWidth
                x = [m, s, e, m]
                y1 = [top, top + 0.1, top + 0.1, top]
                y2 = [bottom, bottom - 0.1, bottom - 0.1, bottom]
                if (g['color'] != 'random'):
                    ax.fill(x, y1, color = g['color'])
                    ax.fill(x, y2, color = g['color'])
                    ax.plot([m, m], [bottom, top], linewidth = 2, color = g['color'])
                else:
                    rndColor = colorRandom()
                    ax.fill(x, y1, color = rndColor)
                    ax.fill(x, y2, color = rndColor)
                    ax.plot([m, m], [bottom, top], linewidth = 2, color = rndColor)
                if (labelFlag == True):
                    ax.annotate(g['desc'], (m, top + 0.15))

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)

    return fig, ax

def plotNodes(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    color:      "1) String 'Random', or\
                 2) String, color" = 'Random',
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Draw nodes":

    # Check for required fields ===============================================
    if (nodes == None):
        print(ERROR_MISSING_NODES)
        return

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for i in nodes:
            allX.append(nodes[i]['loc'][0])
            allY.append(nodes[i]['loc'][1])
        xMin = min(allX) - 0.25
        xMax = max(allX) + 0.25
        yMin = min(allY) - 0.25
        yMax = max(allY) + 0.25
        xSpan = None
        ySpan = None
        if (xMax - xMin > yMax - yMin):
            xSpan = 20
            ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
        else:
            xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
            ySpan = 20
        fig.set_figheight(ySpan)
        fig.set_figwidth(xSpan)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Draw nodes ==============================================================
    x = []
    y = []
    for node in nodes:
        x.append(nodes[node]['loc'][0])
        y.append(nodes[node]['loc'][1])
        ax.annotate(node, (nodes[node]['loc'][0], nodes[node]['loc'][1]))
    ax.plot(x, y, 'ro')

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)

    return fig, ax

def plotArcs(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    arcs:       "List of 2-tuples, arcs in format of [(nodeID1, nodeID2), ...]" = None,
    arrowFlag:  "Boolean, whether or not add arrows to route" = True,
    arrowHeadwidth: "Arrow head width" = 0.1,
    arrowHeadlength: "Arrow head length" = 0.5,
    color:      "1) String 'Random', or\
                 2) String, color" = 'Random',
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Draw arcs":

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for i in arcs:
            allX.append(i[0][0])
            allY.append(i[0][1])
            allX.append(i[1][0])
            allY.append(i[1][1])
        xMin = min(allX) - 0.25
        xMax = max(allX) + 0.25
        yMin = min(allY) - 0.25
        yMax = max(allY) + 0.25
        xSpan = None
        ySpan = None
        if (xMax - xMin > yMax - yMin):
            xSpan = 20
            ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
        else:
            xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
            ySpan = 20
        fig.set_figheight(ySpan)
        fig.set_figwidth(xSpan)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Draw arcs ===============================================================
    for arc in arcs:
        lx = [nodeLoc[arc[0]][0], nodeLoc[arc[1]][0]]
        ly = [nodeLoc[arc[0]][1], nodeLoc[arc[1]][1]]
        if (color == 'Random'):
            rndColor = colorRandom()
            plt.plot(lx, ly, color = rndColor)
        else:
            plt.plot(lx, ly, color = color)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)

    return fig, ax

def plotRoutes(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    seq:        "List of nodeIDs, visiting sequences" = None,
    arrowFlag:  "Boolean, whether or not add arrows to route" = True,
    arrowHeadwidth: "Arrow head width" = 0.1,
    arrowHeadlength: "Arrow head length" = 0.5,
    color:      "1) String 'Random', or\
                 2) String, color" = 'Random',
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Draw a route, e.g., TSP solution":

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for i in seq:
            allX.append(nodes[i]['loc'][0])
            allY.append(nodes[i]['loc'][1])
        xMin = min(allX) - 0.25
        xMax = max(allX) + 0.25
        yMin = min(allY) - 0.25
        yMax = max(allY) + 0.25
        xSpan = None
        ySpan = None
        if (xMax - xMin > yMax - yMin):
            xSpan = 20
            ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
        else:
            xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
            ySpan = 20
        fig.set_figheight(ySpan)
        fig.set_figwidth(xSpan)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Draw Seq ================================================================
    if (color == 'Random'):
        color = colorRandom()
    if (seq != None and len(seq) > 0):
        lx = []
        ly = []
        for s in range(len(seq) - 1):
            lx.append(nodes[seq[s]]['loc'][0])
            ly.append(nodes[seq[s]]['loc'][1])
        lx.append(nodes[seq[len(seq) - 1]]['loc'][0])
        ly.append(nodes[seq[len(seq) - 1]]['loc'][1])
        if (color == 'Random'):
            ax.plot(lx, ly, color=rndColor)
        else:
            ax.plot(lx, ly, color=color)

    # Draw arrows =============================================================
    if (seq != None and len(seq) > 0):
        for s in range(len(seq) - 1):
            x1 = nodes[seq[s]]['loc'][0]
            y1 = nodes[seq[s]]['loc'][1]
            x2 = nodes[seq[s + 1]]['loc'][0]
            y2 = nodes[seq[s + 1]]['loc'][1]
            dx = x2 - x1
            dy = y2 - y1
            # if (color == 'Random'):
            #     ax.arrow(x=x1, y=y1, dx=dx, dy=dy, linewidth=1, color=rndColor)
            #     ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=1, head_width=arrowHeadwidth, head_length=arrowHeadlength, color=rndColor)
            ax.arrow(x=x1, y=y1, dx=dx, dy=dy, linewidth=1, color=color)
            if (arrowFlag):
                ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=1, head_width=arrowHeadwidth, head_length=arrowHeadlength, color=color)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)

    return fig, ax

def plotReqs(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    reqs:       "Dictionary, includes all pickup and delivery requests, in the format of \
                {\
                    'reqID': {\
                        'pickup': pickupZipCode,\
                        'delivery': deliveryZipCode,\
                        'size': itemSize, \
                    }, \
                    'reqID': ... \
                }" = None,
    color:      "1) String 'Random', or \
                 2) String of color" = 'Random',
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Plot delivery requests":

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for r in reqs:
            allX.append(nodes[reqs[r]['pickup']]['loc'][0])
            allX.append(nodes[reqs[r]['delivery']]['loc'][0])
            allY.append(nodes[reqs[r]['pickup']]['loc'][1])
            allY.append(nodes[reqs[r]['delivery']]['loc'][1])
        xMin = min(allX) - 0.25
        xMax = max(allX) + 0.25
        yMin = min(allY) - 0.25
        yMax = max(allY) + 0.25
        xSpan = None
        ySpan = None
        if (xMax - xMin > yMax - yMin):
            xSpan = 20
            ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
        else:
            xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
            ySpan = 20
        fig.set_figheight(ySpan)
        fig.set_figwidth(xSpan)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Each request will be plot as arrow ======================================
    for r in reqs:
        x1 = nodes[reqs[r]['pickup']]['loc'][0]
        y1 = nodes[reqs[r]['pickup']]['loc'][1]
        x2 = nodes[reqs[r]['delivery']]['loc'][0]
        y2 = nodes[reqs[r]['delivery']]['loc'][1]
        dx = x2 - x1
        dy = y2 - y1

        # Define color
        if (color == 'Random'):
            rndColor = colorRandom()
            ax.arrow(x=x1, y=y1, dx=dx, dy=dy, linewidth=1, color=rndColor)
            ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=1, head_width=0.1, head_length=0.2, color=rndColor)
        else:
            ax.arrow(x=x1, y=y1, dx=dx, dy=dy, linewidth=1, color=color)
            ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=1, head_width=0.1, head_length=0.2, color=color)    
        ax.annotate(r, (x1 + dx / 2, y1 + dy / 2))
    
    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)

    return fig, ax

def plotActions(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    actions:    "List of 3-tuples, in format of [('pickup', reqID, nodeID), ...]" = [],
    color:      "1) String 'Random', or \
                 2) String of color" = 'Random',
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Plot delivery requests":

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for r in actions:
            allX.append(nodes[r[2]]['loc'][0])
            allY.append(nodes[r[2]]['loc'][1])
        xMin = min(allX) - 0.25
        xMax = max(allX) + 0.25
        yMin = min(allY) - 0.25
        yMax = max(allY) + 0.25
        xSpan = None
        ySpan = None
        if (xMax - xMin > yMax - yMin):
            xSpan = 20
            ySpan = 20 * ((yMax - yMin) / (xMax - xMin))
        else:
            xSpan = 20 * ((xMax - xMin) / (yMax - yMin))
            ySpan = 20
        fig.set_figheight(ySpan)
        fig.set_figwidth(xSpan)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Draw actions ============================================================
    route = []
    for act in actions:
        route.append(act[2])
    fig, ax = plotRoutes(
        fig = fig,
        ax = ax,
        nodes = nodes,
        seq = route,
        color = color,
        arrowFlag = True,
        saveFigPath = saveFigPath)

    return fig, ax
