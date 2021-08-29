import matplotlib.pyplot as plt
import random

from .msg import *
from .color import *

# [Constructing]
def plotTimeWindows(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    tws:        "List of time windows" = None, 
    color:      "1) String 'Random', or\
                 2) String, color" = 'Random',
    height:     "Height of the figure" = 2,
    width:      "Width of the figure" = 10,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None    
    ) -> "Given a list of time windows, plot it":
    
    # If there is no based matplotlib figure, define it =======================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        earliest = None
        latest = None
        for tw in tws:
            if (tw[0] != None):
                if (earliest == None or earliest > tw[0]):
                    earliest = tw[0]
            if (tw[0] != None):
                if (latest == None or latest < tw[1]):
                    latest = tw[1]
        if (earliest == None):
            pass

    return fig, ax    

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
                        'desc': (optional) [List of descriptions, correspond to `timeStamps`],\
                        'trackID': (optional, default as '0') parallel gantt chart ID, \
                        'color': (optional, default as 'Random') color, \
                        'style': (optional, default as 'solid') 'solid' \
                    }, ... , ]\
                " = None,
    gridFlag:   "True if turn on the grid as background" = True,
    labelFlag:  "True if add label of entities on Gantt chart" = True,
    linewidth:  "The width of Gantt block borders" = 1,
    entities:   "The Gantt chart will be drawn in this order, if None, takes the entities in `gantt`" = None,
    startTime:  "Start time of Gantt, default to be 0, if None, use the earliest time in `gantt`" = 0,
    endTime:    "End time of Gantt, default to be None, if None, use the latest time in `gantt`" = None,
    showTail:   "Show the latest time of all gantt blocks" = True,
    width:      "Width of the figure" = 12,
    height:     "Height of the figure" = 5,
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
    if (startTime != None and startTime > realStart):
        startTime = realStart
    else:
        startTime = realStart

    if (endTime != None and endTime < realEnd):
        endTime = realEnd
    else:
        endTime = realEnd

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
                ax.plot(x, y, color = 'black', linewidth = linewidth)
                if (g['color'] != 'random'):
                    ax.fill(x, y, color = g['color'], linewidth = linewidth)
                else:
                    rndColor = colorRandom()
                    ax.fill(x, y, color = rndColor, linewidth = linewidth)
                if (labelFlag == True):
                    ax.annotate(g['desc'], (s, top + 0.1))
                if (g['style'] != 'solid'):
                    ax.fill(x, y, hatch = g['style'], fill=False, linewidth = linewidth)
            elif ('timeStamps' in g):
                for i in range(len(g['timeStamps']) - 1):
                    s = g['timeStamps'][i]
                    e = g['timeStamps'][i + 1]
                    x = [s, s, e, e, s]
                    y = [bottom, top, top, bottom, bottom]
                    ax.plot(x, y, color = 'black', linewidth = linewidth)
                    ax.annotate(g['desc'][i], (s, top + 0.1))
                    if (g['color'] != 'random'):
                        ax.fill(x, y, color = g['color'], linewidth = linewidth)
                    else:
                        rndColor = colorRandom()
                        ax.fill(x, y, color = rndColor, linewidth = linewidth)
                    if (g['style'] != 'solid'):
                        ax.fill(x, y, hatch = g['style'], fill=False, linewidth = linewidth)
                if (labelFlag == True):
                    ax.annotate(g['desc'][-1], (g['timeStamps'][-1], top + 0.1))

    # Show time span ==========================================================
    if (showTail):
        xTicks = list(ax.get_xticks())
        xTicks.append(realEnd)
        ax.set_xticks(xTicks)

    # Fix height if fig, ax are not provided ==================================
    if (fig == None or ax == None):
        fig.set_figheight(5 * len(entities))

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
    xMin:       "min of x-axis" = None,
    xMax:       "max of x-axis" = None,
    yMin:       "min of y-axis" = None,
    yMax:       "max of y-axis" = None,
    edgeWidth:  "Width on the edge" = 0.5,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Draw nodes":

    # Check for required fields ===============================================
    if (nodes == None):
        print(ERROR_MISSING_NODES)
        return

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for i in nodes:
            allX.append(nodes[i]['loc'][0])
            allY.append(nodes[i]['loc'][1])
        if (xMin == None):
            xMin = min(allX) - edgeWidth
        if (xMax == None):
            xMax = max(allX) + edgeWidth
        if (yMin == None):
            yMin = min(allY) - edgeWidth
        if (yMax == None):
            yMax = max(allY) + edgeWidth
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
    linewidth:  "Width of arcs" = 1,
    arrowFlag:  "Boolean, whether or not add arrows to route" = True,
    arrowHeadwidth: "Arrow head width" = 0.1,
    arrowHeadlength: "Arrow head length" = 0.2,
    color:      "1) String 'Random', or\
                 2) String, color" = 'Random',
    xMin:       "min of x-axis" = None,
    xMax:       "max of x-axis" = None,
    yMin:       "min of y-axis" = None,
    yMax:       "max of y-axis" = None,
    edgeWidth:  "Width on the edge" = 0.5,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Draw arcs":

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for i in arcs:
            allX.append(nodes[i[0]]['loc'][0])
            allY.append(nodes[i[0]]['loc'][1])
            allX.append(nodes[i[1]]['loc'][0])
            allY.append(nodes[i[1]]['loc'][1])
        if (xMin == None):
            xMin = min(allX) - edgeWidth
        if (xMax == None):
            xMax = max(allX) + edgeWidth
        if (yMin == None):
            yMin = min(allY) - edgeWidth
        if (yMax == None):
            yMax = max(allY) + edgeWidth
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
        x1 = nodes[arc[0]]['loc'][0]
        y1 = nodes[arc[0]]['loc'][1]
        x2 = nodes[arc[1]]['loc'][0]
        y2 = nodes[arc[1]]['loc'][1]
        dx = x2 - x1
        dy = y2 - y1
        if (color == 'Random'):
            rndColor = colorRandom()
            plt.plot([x1, x2], [y1, y2], color = rndColor)
            if (arrowFlag):
                ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=linewidth, head_width=arrowHeadwidth, head_length=arrowHeadlength, color=rndColor)
        else:
            plt.plot([x1, x2], [y1, y2], color = color)
            if (arrowFlag):
                ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=linewidth, head_width=arrowHeadwidth, head_length=arrowHeadlength, color=color)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)

    return fig, ax

def plotSeq(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y)}, \
                        nodeID2: {'loc': (x, y)}, \
                        ... \
                    }" = None, 
    seq:        "List of nodeIDs" = None,
    linewidth:  "Width of arcs" = 1,
    arrowFlag:  "Boolean, whether or not add arrows to route" = True,
    arrowHeadwidth: "Arrow head width" = 0.1,
    arrowHeadlength: "Arrow head length" = 0.2,
    color:      "1) String 'Random', or\
                 2) String, color" = 'Random',
    xMin:       "min of x-axis" = None,
    xMax:       "max of x-axis" = None,
    yMin:       "min of y-axis" = None,
    yMax:       "max of y-axis" = None,
    edgeWidth:  "Width on the edge" = 0.5,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None

    ) -> "Draw a route, e.g., TSP solution":    
    # Create arcs =============================================================
    arcs = []
    for i in range(len(seq) - 1):
        arcs.append([seq[i], seq[i + 1]])

    # Color ===================================================================
    if (color == 'Random'):
        color = colorRandom()

    # Call plotArcs ===========================================================
    fig, ax = plotArcs(
        fig = fig,
        ax = ax,
        nodes = nodes,
        arcs = arcs,
        linewidth = linewidth,
        arrowFlag = arrowFlag,
        arrowHeadwidth = arrowHeadwidth,
        arrowHeadlength = arrowHeadlength,
        color = color,
        xMin = xMin,
        xMax = xMax,
        yMin = yMin,
        yMax = yMax,
        edgeWidth = edgeWidth,
        saveFigPath = saveFigPath)

    return fig, ax

