import matplotlib.pyplot as plt
import random

from .msg import *
from .color import *
from .geometry import *

def plotGrid(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    gridColRow: "Number of columns, and number of rows" = (None, None),
    barriers:   "List of blocking grids" = [],
    gridSize:   "Size of the grid" = 1,
    gridBackColor: "Background color of grids" = None,
    gridEdgeColor: "Edge color of grids" = 'black',
    gridOpacity: "Opacity of grids background colors" = 1,
    barrierBackColor: "Background color of barriers" = 'gray',
    barrierEdgeColor: "Edge color of barriers" = 'black',
    barrierOpacity: "Opacity of barriers" = 0.5,
    barrierBackStyle: "Background style of barriers" = '///',
    ) -> "Plot a grid with barriers":

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        xMin = (0 - 0.5) * gridSize
        yMin = (0 - 0.5) * gridSize
        xMax = (gridColRow[0] + 0.5) * gridSize
        yMax = (gridColRow[1] + 0.5) * gridSize
        fig.set_figwidth(xMax)
        fig.set_figheight(yMax)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)
        plt.axis('off')

    for col in range(gridColRow[0]):
        for row in range(gridColRow[1]):
            poly = [
                [col * gridSize, row * gridSize], 
                [col * gridSize + gridSize, row * gridSize],
                [col * gridSize + gridSize, row * gridSize + gridSize],
                [col * gridSize, row * gridSize + gridSize],
                [col * gridSize, row * gridSize]
            ]
            if ((col, row) not in barriers):
                plotPoly(
                    fig = fig,
                    ax = ax,
                    poly = poly,
                    edgeColor = gridEdgeColor,
                    fillColor = gridBackColor,
                    opacity = gridOpacity)
            else:
                plotPoly(
                    fig = fig,
                    ax = ax,
                    poly = poly,
                    edgeColor = barrierEdgeColor,
                    fillColor = barrierBackColor,
                    opacity = barrierOpacity,
                    fillStyle = barrierBackStyle)

    return fig, ax

def plotGridPath(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    gridColRow: "Number of columns, and number of rows" = (None, None),
    gridSize:   "Size of the grid" = 1,
    path:       "The sequences of visiting grids, a list of coordinates" = None,
    pathColor:  "The color of path" = 'Random',
    pathWidth:  "The width of path" = 3,
    markerSize: "Size of starting/ending points" = 15
    ) -> "Plot the path on the grid":

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        xMin = (0 - 0.5) * gridSize
        yMin = (0 - 0.5) * gridSize
        xMax = (gridColRow[0] + 0.5) * gridSize
        yMax = (gridColRow[1] + 0.5) * gridSize
        fig.set_figwidth(xMax)
        fig.set_figheight(yMax)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)
        plt.axis('off')

    # Path color ==============================================================
    if (pathColor == 'Random'):
        pathColor = colorRandom()

    # Plot the origin/destination =============================================
    ax.plot(path[0][0] * gridSize + gridSize / 2, path[0][1] * gridSize + gridSize / 2, marker = 'o', markersize= markerSize, color = pathColor)
    ax.plot(path[-1][0] * gridSize + gridSize / 2, path[-1][1] * gridSize + gridSize / 2, marker = 's', markersize= markerSize, color = pathColor)
    for i in range(len(path) - 1):
        x = [path[i][0] * gridSize + gridSize / 2, path[i + 1][0] * gridSize + gridSize / 2]
        y = [path[i][1] * gridSize + gridSize / 2, path[i + 1][1] * gridSize + gridSize / 2]
        ax.plot(x, y, color = pathColor, linewidth = pathWidth)
    return fig, ax

def plotPoly(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    poly:       "Polygon to be plot" = None,
    xyReverseFlag: "Reverse x, y. Usually use for (lat, lon)" = False,
    linewidth:  "Width of arcs" = 1,
    edgeColor:  "1) String 'Random', or\
                 2) String, color" = 'Random',
    fillColor:  "1) (default) String 'Random', or\
                 2) None, no fill, or\
                 3) String, color" = None,
    fillStyle:  "Background style, None if no style, '///' for shadow" = None,
    opacity:    "Opacity of filled area" = 0.5,
    figSize:    "Size of the figure, in (width, height)" = (5, 5), 
    xMin:       "min of x-axis" = None,
    xMax:       "max of x-axis" = None,
    yMin:       "min of y-axis" = None,
    yMax:       "max of y-axis" = None,
    edgeWidth:  "Width on the edge" = 0.5,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Draw a route, e.g., TSP solution":    

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for pt in poly:
            if (not xyReverseFlag):
                allX.append(pt[0])
                allY.append(pt[1])
            else:
                allX.append(pt[1])
                allY.append(pt[0])
        if (xMin == None):
            xMin = min(allX) - edgeWidth
        if (xMax == None):
            xMax = max(allX) + edgeWidth
        if (yMin == None):
            yMin = min(allY) - edgeWidth
        if (yMax == None):
            yMax = max(allY) + edgeWidth
        if (figSize == None):
            if (xMax - xMin > yMax - yMin):
                width = 5
                height = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                width = 5 * ((xMax - xMin) / (yMax - yMin))
                height = 5
        else:
            (width, height) = figSize
        fig.set_figwidth(width)
        fig.set_figheight(height)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Get the x, y list =======================================================
    x = []
    y = []
    for pt in poly:
        if (not xyReverseFlag):
            x.append(pt[0])
            y.append(pt[1])
        else:
            x.append(pt[1])
            y.append(pt[0])
    if (not xyReverseFlag):
        x.append(poly[0][0])
        y.append(poly[0][1])
    else:
        x.append(poly[0][1])
        y.append(poly[0][0])        

    # Plot ====================================================================
    if (edgeColor == 'Random'):
        edgeColor = colorRandom()
    if (fillColor == None):
        ax.plot(x, y, color = edgeColor, linewidth = linewidth)
    else:
        if (fillColor == 'Random'):
            fillColor = colorRandom()
        ax.fill(x, y, facecolor=fillColor, edgecolor=edgeColor, hatch=fillStyle, linewidth=linewidth, alpha=opacity)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)

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
                        'color': (optional, default as 'random') color, \
                        'style': (optional, default as 'solid') 'solid' \
                    }, ... , \
                    {\
                        'entityID': entityID, \
                        'timeStamps': [timeStamp1, timeStamp2, ..., timeStampN], \
                        'desc': (optional) [List of descriptions, correspond to `timeStamps`],\
                        'trackID': (optional, default as '0') parallel gantt chart ID, \
                        'color': (optional, default as 'random') color, \
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
    figSize:    "Size of the figure, in (width, height)" = (12, 5),
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
        fig.set_figheight(figSize[1])
        fig.set_figwidth(figSize[0])
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
    xyReverseFlag: "Reverse x, y. Usually use for (lat, lon)" = False,
    color:      "Decide the color of nodes if the 'color' tag is not in `nodes` \
                 1) String 'Random', or\
                 2) String, color" = 'Random',
    figSize:    "Size of the figure, in (width, height)" = (5, 5), 
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
            if (not xyReverseFlag):
                allX.append(nodes[i]['loc'][0])
                allY.append(nodes[i]['loc'][1])
            else:
                allX.append(nodes[i]['loc'][1])
                allY.append(nodes[i]['loc'][0])
        if (xMin == None):
            xMin = min(allX) - edgeWidth
        if (xMax == None):
            xMax = max(allX) + edgeWidth
        if (yMin == None):
            yMin = min(allY) - edgeWidth
        if (yMax == None):
            yMax = max(allY) + edgeWidth
        if (figSize == None):
            if (xMax - xMin > yMax - yMin):
                width = 5
                height = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                width = 5 * ((xMax - xMin) / (yMax - yMin))
                height = 5
        else:
            (width, height) = figSize
        fig.set_figwidth(width)
        fig.set_figheight(height)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Draw nodes ==============================================================
    for n in nodes:
        # Define color --------------------------------------------------------
        nodeColor = None
        if ('color' in nodes[n]):
            nodeColor = nodes[n]['color']
        elif (color == 'Random'):
            nodeColor = colorRandom()
        else:
            nodeColor = color

        # plot nodes ----------------------------------------------------------
        x = None
        y = None
        if (not xyReverseFlag):
            x = nodes[n]['loc'][0]
            y = nodes[n]['loc'][1]
        else:
            x = nodes[n]['loc'][1]
            y = nodes[n]['loc'][0]

        ax.plot(x, y, marker = 'o', color = nodeColor)
        if ('label' not in nodes[n]):
            ax.annotate(n, (x, y))
        else:
            ax.annotate(nodes[n]['label'], (x, y))

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
    xyReverseFlag: "Reverse x, y. Usually use for (lat, lon)" = False,
    arcs:       "List of 2-tuples, arcs in format of [(nodeID1, nodeID2), ...]" = None,
    linewidth:  "Width of arcs" = 1,
    arrowFlag:  "Boolean, whether or not add arrows to route" = True,
    arrowHeadwidth: "Arrow head width" = 0.1,
    arrowHeadlength: "Arrow head length" = 0.2,
    color:      "1) String 'Random', or\
                 2) String, color" = 'Random',
    figSize:    "Size of the figure, in (width, height)" = (5, 5), 
    xMin:       "min of x-axis" = None,
    xMax:       "max of x-axis" = None,
    yMin:       "min of y-axis" = None,
    yMax:       "max of y-axis" = None,
    edgeWidth:  "Width on the edge" = 0.5,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None
    ) -> "Draw arcs":

    # Check for required fields ===============================================
    if (nodes == None):
        print(ERROR_MISSING_NODES)
        return

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for i in arcs:
            if (not xyReverseFlag):
                allX.append(nodes[i[0]]['loc'][0])
                allY.append(nodes[i[0]]['loc'][1])
                allX.append(nodes[i[1]]['loc'][0])
                allY.append(nodes[i[1]]['loc'][1])
            else:
                allX.append(nodes[i[0]]['loc'][1])
                allY.append(nodes[i[0]]['loc'][0])
                allX.append(nodes[i[1]]['loc'][1])
                allY.append(nodes[i[1]]['loc'][0])
        if (xMin == None):
            xMin = min(allX) - edgeWidth
        if (xMax == None):
            xMax = max(allX) + edgeWidth
        if (yMin == None):
            yMin = min(allY) - edgeWidth
        if (yMax == None):
            yMax = max(allY) + edgeWidth
        if (figSize == None):
            if (xMax - xMin > yMax - yMin):
                width = 5
                height = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                width = 5 * ((xMax - xMin) / (yMax - yMin))
                height = 5
        else:
            (width, height) = figSize
        fig.set_figwidth(width)
        fig.set_figheight(height)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Draw arcs ===============================================================
    for arc in arcs:
        if (not xyReverseFlag):
            x1 = nodes[arc[0]]['loc'][0]
            y1 = nodes[arc[0]]['loc'][1]
            x2 = nodes[arc[1]]['loc'][0]
            y2 = nodes[arc[1]]['loc'][1]
        else:
            x1 = nodes[arc[0]]['loc'][1]
            y1 = nodes[arc[0]]['loc'][0]
            x2 = nodes[arc[1]]['loc'][1]
            y2 = nodes[arc[1]]['loc'][0]
        dx = x2 - x1
        dy = y2 - y1
        if (color == 'Random'):
            rndColor = colorRandom()
            ax.plot([x1, x2], [y1, y2], color = rndColor)
            if (arrowFlag):
                ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=linewidth, head_width=arrowHeadwidth, head_length=arrowHeadlength, color=rndColor)
        else:
            ax.plot([x1, x2], [y1, y2], color = color)
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
    xyReverseFlag: "Reverse x, y. Usually use for (lat, lon)" = False,
    seq:        "List of nodeIDs" = None,
    linewidth:  "Width of arcs" = 1,
    arrowFlag:  "Boolean, whether or not add arrows to route" = True,
    arrowHeadwidth: "Arrow head width" = 0.1,
    arrowHeadlength: "Arrow head length" = 0.2,
    color:      "1) String 'Random', or\
                 2) String, color" = 'Random',
    figSize:    "Size of the figure, in (width, height)" = (5, 5), 
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
        xyReverseFlag = xyReverseFlag,
        arcs = arcs,
        linewidth = linewidth,
        arrowFlag = arrowFlag,
        arrowHeadwidth = arrowHeadwidth,
        arrowHeadlength = arrowHeadlength,
        color = color,
        figSize = figSize,
        xMin = xMin,
        xMax = xMax,
        yMin = yMin,
        yMax = yMax,
        edgeWidth = edgeWidth,
        saveFigPath = saveFigPath)

    return fig, ax

