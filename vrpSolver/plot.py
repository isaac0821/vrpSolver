import matplotlib.pyplot as plt
import random

from .color import *
from .geometry import *
from .msg import *
from .weather import *

def plotGrid(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    gridColRow: "Number of columns, and number of rows" = (None, None),
    barriers:   "List of blocking grids" = [],
    labeling:   "Additional labeling of grids to support different color or annotation, in format of \
                    {\
                        (coordX, coordY): { \
                            'color': color, \
                            'annotation': annotation \
                        }\
                    }" = None,
    gridSize:   "Size of the grid" = 1,
    gridBackColor: "Background color of grids" = None,
    gridEdgeColor: "Edge color of grids" = 'black',
    gridOpacity: "Opacity of grids background colors" = 1,
    barrierBackColor: "Background color of barriers" = 'gray',
    barrierEdgeColor: "Edge color of barriers" = 'black',
    barrierOpacity: "Opacity of barriers" = 0.5,
    barrierBackStyle: "Background style of barriers" = '///',
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
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
                plotPolygon(
                    fig = fig,
                    ax = ax,
                    poly = poly,
                    edgeColor = gridEdgeColor,
                    fillColor = gridBackColor,
                    opacity = gridOpacity)
            else:
                plotPolygon(
                    fig = fig,
                    ax = ax,
                    poly = poly,
                    edgeColor = barrierEdgeColor,
                    fillColor = barrierBackColor,
                    opacity = barrierOpacity,
                    fillStyle = barrierBackStyle)

    # Add labels ==============================================================
    if (labeling != None):
        for (col, row) in labeling:
            if ('color' in labeling[(col, row)]):
                plotPolygon(
                    fig = fig,
                    ax = ax,
                    poly = poly,
                    edgeColor = gridEdgeColor,
                    fillColor = labeling[(col, row)]['color'],
                    opacity = gridOpacity)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotGridPath(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    gridColRow: "Number of columns, and number of rows" = (None, None),
    barriers:   "List of blocking grids, needed if plotGridFlag is True" = [],
    gridSize:   "Size of the grid" = 1,
    plotGridFlag: "True if plot the grid as background, false otherwise" = True,
    path:       "The sequences of visiting grids, a list of coordinates" = None,
    pathColor:  "The color of path" = 'Random',
    pathWidth:  "The width of path" = 3,
    markerSize: "Size of starting/ending points" = 15,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
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

    # Plot grid background ====================================================
    if (plotGridFlag):
        fig, ax = plotGrid(
            fig = fig,
            ax = ax,
            gridColRow = gridColRow,
            barriers = barriers,
            gridSize = gridSize)

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

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotWinds(
    fig:        "fig" = None,
    ax:         "ax" = None,
    width:      "Width of the figure" = None,
    height:     "Height of the figure" = 3,
    winds:      "List of dictionary, wind speed/wind direct in time sequence, time starts from 0\
                [{\
                    'startTime': startTime,\
                    'endTime': endTime,\
                    'windSpd': windSpd, # In [m/s]\
                    'windDeg': winDeg, # In [Degree]\
                }, ...]" = None,
    startOfDay: "Start time of the day in [h]" = 0,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
    ) -> "Given a list of wind, plot the wind direction and wind speed":

    # Draw frame ==============================================================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        if (height != None):
            fig.set_figheight(height)
        if (width != None):
            fig.set_figwidth(width)
        else:
            fig.set_figwidth(len(winds) * 0.5 + 0.25)
    ax.set_xlabel("Time of the day [h]")
    ax.set_ylabel("Wind speed [m/s]")
    
    # Draw wind speeds ========================================================
    ticks = [winds[0]['startTime'] / 3600 + startOfDay]
    for w in range(len(winds)):
        ticks.append(winds[w]['endTime'] / 3600 + startOfDay)
        ax.plot(
            [winds[w]['startTime'] / 3600 + startOfDay, winds[w]['endTime'] / 3600 + startOfDay], 
            [winds[w]['windSpd'], winds[w]['windSpd']],
            color = 'black')
    for w in range(1, len(winds)):
        ax.plot(
            [winds[w - 1]['endTime'] / 3600 + startOfDay, winds[w]['startTime'] / 3600 + startOfDay], 
            [winds[w - 1]['windSpd'], winds[w]['windSpd']],
            color = 'black')

    # Draw wind direction =====================================================
    for w in range(len(winds)):
        arrowMidX = (winds[w]['startTime'] + winds[w]['endTime']) / 7200 + startOfDay
        arrowMidY = winds[w]['windSpd'] + 0.7
        dx = math.sin(math.radians(winds[w]['windDeg'])) * 0.25
        dy = math.cos(math.radians(winds[w]['windDeg'])) * 0.25
        ax.arrow(
            x = arrowMidX - dx / 2, 
            y = arrowMidY - dy / 2, 
            dx = dx, 
            dy = dy, 
            linewidth=2, head_width=0.15, head_length=0.32)

    # Set ticks ===============================================================
    ax.set_xticks(ticks)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotCloudsInTime(
    fig:        "fig" = None,
    ax:         "ax" = None,
    figSize:    "Size of the figure, in (width, height)" = (5, 5), 
    xMin:       "min of x-axis" = None,
    xMax:       "max of x-axis" = None,
    yMin:       "min of y-axis" = None,
    yMax:       "max of y-axis" = None,
    edgeWidth:  "Width on the edge" = 0.05,
    clouds:     "A list of clouds" = None,
    nodes:      "nodes with locations coordinates" = None,
    timeStamp:  "Time stamps of the frame" = None,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
    ) -> "Given a time stamp, plot the locations of clouds and customers":

    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []

        # locs of nodes -------------------------------------------------------
        for n in nodes:
            allX.append(nodes[n]['loc'][1])
            allY.append(nodes[n]['loc'][0])

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

    # Plot nodes ==============================================================
    fig, ax = plotNodes(
        fig = fig,
        ax = ax,
        nodes = nodes,
        color = 'green',
        xyReverseFlag = True)

    # Plot clouds =============================================================
    for c in clouds:
        currentCloudPosition = getCloudCurrentPosition(c, timeStamp)
        if (currentCloudPosition != None):
            fig, ax = plotPolygon(
                fig = fig,
                ax = ax,
                edgeColor = 'black',
                fillColor = 'gray',
                fillStyle = '///',
                opacity = 0.8,
                poly = currentCloudPosition,
                xyReverseFlag = True)
    plt.close(fig)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotRoadNetwork(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    roadNetwork: "A road network dictionary" = None,
    linewidth:  "Width of arcs" = 1,
    color:      "1) String 'Random', or\
                 2) String, color" = 'black',
    figSize:    "Size of the figure, in (width, height)" = (5, 5), 
    xMin:       "min of x-axis" = None,
    xMax:       "max of x-axis" = None,
    yMin:       "min of y-axis" = None,
    yMax:       "max of y-axis" = None,
    edgeWidth:  "Width on the edge" = 0.005,
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
    ) -> "Draw a set of polylines, usually for plotting road network": 

    # FIXME: In future, we might want to distinguish roads by max speed or show the names of roads
    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for road in roadNetwork:
            for pt in roadNetwork[road]['line']:
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
    for road in roadNetwork:
        x = []
        y = []
        for pt in roadNetwork[road]['line']:
            x.append(pt[1])
            y.append(pt[0])
        edgeColor = color
        if (color == 'Random'):
                edgeColor = colorRandom()
        ax.plot(x, y, color = edgeColor, linewidth = linewidth)

    plt.close(fig)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotPolygon(
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
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
    ) -> "Draw a polygon, e.g., cloud":

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
    plt.close(fig)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotNodes(
    fig:        "Based matplotlib figure object" = None,
    ax:         "Based matplotlib ax object" = None,
    nodes:      "Dictionary, returns the coordinate of given nodeID, \
                    {\
                        nodeID1: {'loc': (x, y), 'marker': 'r', 'color': 'red', 'size': 3}, \
                        nodeID2: {'loc': (x, y), 'marker': 'r', 'color': 'red', 'size': 3}, \
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
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
    ) -> "Draw nodes":

    # Check for required fields ===============================================
    if (nodes == None):
        msgError(ERROR_MISSING_NODES)
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

        # Define marker and marker size ---------------------------------------
        nodeMarker = 'o'
        if ('marker' in nodes[n]):
            nodeMarker = nodes[n]['marker']
        nodeMarkersize = None
        if ('markersize' in nodes[n]):
            nodeMarkersize = nodes[n]['markersize']

        # plot nodes ----------------------------------------------------------
        x = None
        y = None
        if (not xyReverseFlag):
            x = nodes[n]['loc'][0]
            y = nodes[n]['loc'][1]
        else:
            x = nodes[n]['loc'][1]
            y = nodes[n]['loc'][0]
        if (nodeMarker == None):
            ax.plot(x, y, color = nodeColor, marker = nodeMarker)
        else:
            ax.plot(x, y, color = nodeColor, marker = nodeMarker, markersize = nodeMarkersize)
        if ('label' not in nodes[n]):
            ax.annotate(n, (x, y))
        else:
            ax.annotate(nodes[n]['label'], (x, y))

    plt.close(fig)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

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
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
    ) -> "Draw arcs":

    # Check for required fields ===============================================
    if (nodes == None):
        msgError(ERROR_MISSING_NODES)
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
    plt.close(fig)

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

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
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
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
    entities:   "1) None, takes the entities in `gantt` \
                 2) List of strings, the Gantt chart will be drawn in this order, \
                 3) List of lists (strings), the Gantt chart will be drawn in groups" = None,
    startTime:  "Start time of Gantt, default to be 0, if None, use the earliest time in `gantt`" = 0,
    endTime:    "End time of Gantt, default to be None, if None, use the latest time in `gantt`" = None,
    showTail:   "Show the latest time of all gantt blocks" = True,
    figSize:    "Size of the figure, in (width, height)" = (12, 5),
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
    ) -> "Given a Gantt dictionary, plot Gantt":

    # Check for required fields ===============================================
    if (gantt == None):
        msgError(ERROR_MISSING_GANTT)
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

    # Check overwritten fields ================================================
    if (startTime != None):
        startTime = startTime
    else:
        startTime = realStart

    if (endTime != None):
        endTime = endTime
    else:
        endTime = realEnd

    # Arrange entities ========================================================              
    if (entities != None):
        # Check inputs
        groupFlag = False
        for e in entities:
            if (type(e) == list):
                groupFlag = True
                break
        # If the type of entities is List of Lists, the entities are grouped
        if (groupFlag == True):
            groupEntities = []
            for el in entities:
                if (type(el) != list):
                    msgError(ERROR_INCOR_GANTT_ENTITYGROUP)
                    return
                for e in el:
                    groupEntities.append(e)
                groupEntities.append(None)
            groupEntities = groupEntities[:-1] # Remove the last None
            entities = [i for i in groupEntities]
        for e in entities:
            if (e not in entList and e != None):
                msgError(ERROR_INCOR_GANTT_MISSENT)
                return
    elif (entities == None):
        entities = [i for i in entList]

    # If no based matplotlib figure, define fig size ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        fig.set_figheight(figSize[1])
        fig.set_figwidth(figSize[0])
        ax.set_xlim(startTime, endTime + (endTime - startTime) * 0.05)
        ax.set_ylim(0, len(entities) + 0.2)

    # Set axis ================================================================
    entities.reverse()
    yticks = []
    pos = 0.5
    for i in range(len(entities)):
        yticks.append(pos)
        if (entities[i] != None):
            pos += 1
        else:
            pos += 0.5
    ax.set_yticks(yticks)   
    ax.set_yticklabels(entities)
    entities.reverse()
    ax.set_xlabel("Time")

    # Grids ===================================================================
    if (gridFlag):
        ax.grid(b = True, linestyle=':')

    # Loop through `gantt` and draw gantt =====================================
    for g in gantt:
        if (g['entityID'] in entities):
            bottom = yticks[len(yticks) - 1 - entities.index(g['entityID'])] + 0.4 - 0.8
            top = 0
            if (labelFlag == True):
                top = yticks[len(yticks) - 1 - entities.index(g['entityID'])] + 0.20
            else:
                top = yticks[len(yticks) - 1 - entities.index(g['entityID'])] + 0.4
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
        fig.set_figheight(5 * max(pos))

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)
    return fig, ax

def plotStep(
    fig:        "Based matplotlib figure object" = None, 
    ax:         "Based matplotlib ax object" = None,
    step:       "List of dictionaries, in the following format\
                    [{\
                        'resID': resource ID, \
                        'timeStamp': list of time stamps, \
                        'useLevel': number of resource that are been used after correspond time stamp, \
                        'color': (optional, default as 'random') color, \
                        'style': (optional, default as 'solid') 'solid' \
                    }, ... , ]\
                " = None,
    stepInt:    "Interval of stat" = 1,
    showPercentageFlag: "True if show the percentage of time for each level" = True,
    gridFlag:   "True if turn on the grid as background" = True,
    labelFlag:  "True if add label of entities on Gantt chart" = True,
    linewidth:  "The width of step block borders" = 1,
    entities:   "1) None, takes the entities in `step` \
                 2) List of Strings, the step chart will be drawn in this order" = None,
    startTime:  "Start time of step, default to be 0, if None, use the earliest time in `step`" = 0,
    endTime:    "End time of step, default to be None, if None, use the latest time in `step`" = None,
    showTail:   "Show the latest time of all step blocks" = True,
    figSize:    "Size of the figure, in (width, height)" = (12, 5),
    saveFigPath:"1) None, if not exporting image, or \
                 2) String, the path for exporting image" = None,
    showFig:    "True if shows the figure in environment such as Jupyter Notebook, \
                 recommended to turn off if generate a batch of images" = True
    ) -> "Given a step dictionary, plot step chart":

    # Check for required fields ===============================================
    if (step == None):
        msgError(ERROR_MISSING_STEP)
        return

    # Pre-calculate ===========================================================
    realStart = None
    realEnd = None
    resList = []
    for st in step:
        if (st['resID'] not in resList):
            resList.append(st['resID'])

        for t in st['timeStamp']:
            if (realStart == None or realStart > st['timeStamp'][0]):
                realStart = st['timeStamp'][0]
            if (realEnd == None or realEnd < st['timeStamp'][-1]):
                realEnd = st['timeStamp'][-1]
    if (entities != None):
        for e in entities:
            if (e not in resList):
                msgError(ERROR_INCOR_GANTT_MISSENT)
                return
    elif (entities == None):
        entities = [i for i in resList]

    # Check overwritten fields ================================================
    if (startTime != None):
        startTime = startTime
    else:
        startTime = realStart

    if (endTime != None):
        endTime = endTime
    else:
        endTime = realEnd

    # If no based matplotlib figure, define fig size ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        fig.set_figwidth(figSize[0])
        ax.set_xlim(startTime, endTime + (endTime - startTime) * 0.05)

    # Set axis ================================================================
    rangeByEntity = {}
    for st in step:
        rangeByEntity[st['resID']] = [min(st['useLevel']), max(st['useLevel'])]

    entities.reverse()
    # yticks = []

    btm = {}
    yticks = [0]
    yticklabels = [0]
    accBtm = 0
    for i in range(len(entities)):
        btm[entities[i]] = accBtm
        if (accBtm + rangeByEntity[entities[i]][0] != yticks[-1]):
            yticks.append(accBtm + rangeByEntity[entities[i]][0])
            yticklabels.append(rangeByEntity[entities[i]][0])
        for j in range(math.ceil(rangeByEntity[entities[i]][0]), math.floor(rangeByEntity[entities[i]][1])):
            yticks.append(accBtm + j)
            yticklabels.append(j)

        accBtm += rangeByEntity[entities[i]][1] - rangeByEntity[entities[i]][0]
        if (rangeByEntity[entities[i]][1] != yticklabels[-1]):
            yticks.append(accBtm)
            yticklabels.append(rangeByEntity[entities[i]][1])
        accBtm += 1 # Space between step charts

    ax.set_yticks(yticks)
    # entities.reverse()
    ax.set_yticklabels(yticklabels)
    entities.reverse()
    ax.set_xlabel("Time")

    if (fig == None and ax == None):
        ax.set_ylim(0, accBtm)

    # Grids ===================================================================
    if (gridFlag):
        ax.grid(b = True, linestyle=':')

    # Loop through `step` and draw step chart =================================
    for st in step:
        if (st['resID'] in resList):
            bottom = btm[st['resID']]
            polyX = [realStart]
            polyY = [bottom]
            polyX = [st['timeStamp'][0]]
            polyY = [bottom]
            for i in range(len(st['timeStamp']) - 1):
                polyX.append(st['timeStamp'][i])
                polyY.append(bottom + st['useLevel'][i] - rangeByEntity[st['resID']][0])
                polyX.append(st['timeStamp'][i + 1])
                polyY.append(bottom + st['useLevel'][i] - rangeByEntity[st['resID']][0])
            polyX.append(st['timeStamp'][-1])
            polyY.append(bottom + st['useLevel'][-1])
            if (st['timeStamp'][-1] >= realEnd):
                polyX.append(realEnd)
                polyY.append(bottom)

            ax.plot(polyX, polyY, color = 'black', linewidth = linewidth)            
            if ('color' in st and st['color'] != 'random'):
                ax.fill(polyX, polyY, color = st['color'], linewidth = linewidth)
                ax.annotate(st['resID'], (0, bottom + 0.1))
            else:
                rndColor = colorRandom()
                ax.fill(polyX, polyY, color = str(rndColor), linewidth = linewidth)
                ax.annotate(st['resID'], (0, bottom + 0.1), color = 'black')
            if ('style' in st and st['style'] != 'solid'):
                ax.fill(polyX, polyY, hatch = st['style'], fill=False, linewidth = linewidth)

    # Show percentage =========================================================
    if (showPercentageFlag):
        # Get stat
        stStatScale = {}
        for st in step:
            if (st['resID'] in resList):
                stStatScale[st['resID']] = []
                for r in range(0, math.ceil(rangeByEntity[st['resID']][1] - rangeByEntity[st['resID']][0]), stepInt):
                    rangeStat = [rangeByEntity[st['resID']][0] + r * stepInt, rangeByEntity[st['resID']][0] + (r + 1) * stepInt]
                    stStatScale[st['resID']].append({
                        'rangeStat': rangeStat,
                        'totalLength': 0
                    })
                for i in range(len(st['timeStamp']) - 1):
                    lengthOfLevel = st['timeStamp'][i + 1] - st['timeStamp'][i]
                    for j in range(len(stStatScale[st['resID']])):
                        if (st['useLevel'][i] > stStatScale[st['resID']][j]['rangeStat'][0] and st['useLevel'][i] <= stStatScale[st['resID']][j]['rangeStat'][1]):
                            stStatScale[st['resID']][j]['totalLength'] += lengthOfLevel
        # Plot Stat
        for st in step:
            if (st['resID'] in resList):
                for stat in stStatScale[st['resID']]:
                    note = "[" + str(stat['rangeStat'][0]) + ", " + str(stat['rangeStat'][1]) + "]: " + str(round(stat['totalLength'] / realEnd * 100, 2)) + "%"
                    ax.annotate(note, (ax.get_xlim()[1] + 5, 0.1 + btm[st['resID']] + stat['rangeStat'][0]))

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
    if (not showFig):
        plt.close(fig)
    return fig, ax
