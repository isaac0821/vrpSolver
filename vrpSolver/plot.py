import matplotlib.pyplot as plt

from .common import *
from .color import *
from .msg import *
from .province import *
from .error import *

# History =====================================================================
# 20230518 - `plotNodes()` now will plot the neighborhood of nodes
# 20230624 - Rename functions `plotArcs()`, `plotLocSeq()`, `plotNodeSeq()`
# =============================================================================

def plotNodes(
    nodes: dict, 
    nodeColor: str = 'Random',
    nodeMarker: str = 'o',
    nodeMarkersize: float = 1,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    boundingBox: tuple[int|float|None, int|float|None, int|float|None, int|float|None] = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    """Draw nodes

    Parameters
    ----------
    nodes: dictionary, required
        The coordinates and other attributions of the nodes to be plotted, in the following format::
            >>> nodes = {
            ...     nodeID1: {
            ...         'loc': (x, y),
            ...         'neighbor': poly, # Optional, indicate if need to display the neighborhood
            ...         'marker': 'r',    # Optional, default as 'o'
            ...         'markersize': 2,  # Optional, default as None
            ...         'color': 'red',   # Optional, default as 'Random'
            ...         'size': 3,        # Optional, default as 3
            ...         'fontsize': 3,    # Optional, default as 3
            ...     }, # ...
            ... }
    nodeColor: str, optional, default 'Random'
        Alternative option. If 'color' is provided in `node`, this will be ignored.
    neighborColor: str, optional, default 'gray'
        If nodes have 'neighbor' label, will plot the neighbor area in this color
    neighborOpacity: float, optional, default 0.5
        The opacity of neighborhood.
    neighborEntranceColor: str, optional, default 'black', 
        The color of neighbor entrance
    xyReverseFlag: bool, optional, default False
        True if need to reverse the x, y coordinates, e.g., plot for (lat, lon)
    fig: matplotlib object, optional, defaut None
        `fig` and `ax` indicates the matplotlib object to plot on, if not provided, plot in a new figure
    ax: matplotlib object, optional, default None
        See `fig`
    figSize: 2-tuple, optional, default as (None, 5)
        Size of the figure in (width, height). If width or height is set to be None, it will be auto-adjusted.
    boundingBox: 4-tuple, optional, default as (None, None, None, None)
        (xMin, xMax, yMin, yMax), defines four boundaries of the figure
    saveFigPath: string, optional, default as None
        The path for exporting image if provided
    showFig: bool, optional, default as True
        True if show the figure in Juypter Notebook environment

    Returns
    -------
    fig, ax: matplotlib.pyplot object
    """

    # Check for required fields ===============================================
    if (nodes == None):
        raise MissingParameterError(ERROR_MISSING_NODES)

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        # Adjust bounding box
        (xMin, xMax, yMin, yMax) = boundingBox
        if (xMin == None or xMax == None or yMin == None or yMax == None):
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
                xMin = min(allX) - 0.1 * abs(max(allX) - min(allX))
            if (xMax == None):
                xMax = max(allX) + 0.1 * abs(max(allX) - min(allX))
            if (yMin == None):
                yMin = min(allY) - 0.1 * abs(max(allY) - min(allY))
            if (yMax == None):
                yMax = max(allY) + 0.1 * abs(max(allY) - min(allY))
        # Adjust width and height
        width = 0
        height = 0
        if (figSize == None or (figSize[0] == None and figSize[1] == None)):
            if (xMax - xMin > yMax - yMin):
                width = 5
                height = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                width = 5 * ((xMax - xMin) / (yMax - yMin))
                height = 5
        elif (figSize != None and figSize[0] != None and figSize[1] == None):
            width = figSize[0]
            height = figSize[0] * ((yMax - yMin) / (xMax - xMin))
        elif (figSize != None and figSize[0] == None and figSize[1] != None):
            width = figSize[1] * ((xMax - xMin) / (yMax - yMin))
            height = figSize[1]
        else:
            (width, height) = figSize

        if (isinstance(fig, plt.Figure)):
            fig.set_figwidth(width)
            fig.set_figheight(height)
            ax.set_xlim(xMin, xMax)
            ax.set_ylim(yMin, yMax)

    # Draw nodes ==============================================================
    for n in nodes:
        # Define color --------------------------------------------------------
        color = None
        if ('color' in nodes[n]):
            color = nodes[n]['color']
        elif (nodeColor == 'Random'):
            color = colorRandom()
        else:
            color = nodeColor

        # Define marker and marker size ---------------------------------------
        if ('marker' in nodes[n]):
            nodeMarker = nodes[n]['marker']
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
        if (nodeMarkersize == None):
            ax.plot(x, y, color = color, marker = nodeMarker)
        else:
            ax.plot(x, y, color = color, marker = nodeMarker, markersize = nodeMarkersize)
        if ('label' not in nodes[n]):
            lbl = n
        else:
            lbl = nodes[n]['label']

        ha = 'left'
        if ('ha' in nodes[n]):
            ha = nodes[n]['ha']
        va = 'top'
        if ('va' in nodes[n]):
            va = nodes[n]['va']
        ax.annotate(lbl, (x, y), ha=ha, va=va)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None and isinstance(fig, plt.Figure)):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotArcs(
    arcs: dict,
    arcColor: str = 'Random',
    arcWidth: float = 1.0,
    arrowFlag: bool = True,
    arrowHeadWidth: float = 2.0,
    arrowHeadLength: float = 3.0,
    startColor: str = 'black',
    endColor: str = 'black',
    bothEndSize: int|float = 2.0,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    boundingBox: tuple[int|float|None, int|float|None, int|float|None, int|float|None] = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):
    
    """Draw arcs

    Parameters
    ----------

    arcs: dict, required
        A set of arcs, each arc is defined by two points
    arcColor: string, optional, default 'Random'
        Color of arcs
    arcWidth: float, optional, default 1
        Width of arcs
    arrowFlag: bool, optional, default True
        True if plot arrows
    arrowHeadWidth: float, optional, default 0.1
        Width of arrow head
    arrowHeadLength: float, optional, default 0.2
        Length of arrow head
    fig: matplotlib object, optional, default None
        `fig` and `ax` indicates the matplotlib object to plot on, if not provided, plot in a new figure
    ax: matplotlib object, optional, default None
        See `fig`
    xyReverseFlag: bool, optional, default False
        True if need to reverse the x, y coordinates, e.g., plot for (lat, lon)
    figSize: 2-tuple, optional, default as (None, 5)
        Size of the figure in (width, height). If width or height is set to be None, it will be auto-adjusted.
    boundingBox: 4-tuple, optional, default as (None, None, None, None)
        (xMin, xMax, yMin, yMax), defines four boundaries of the figure
    saveFigPath: string, optional, default as None
        The path for exporting image if provided
    showFig: bool, optional, default as True
        True if show the figure in Juypter Notebook environment

    Returns
    -------
    fig, ax: matplotlib.pyplot object
    """

    # Check for required fields ===============================================
    if (arcs == None):
        raise MissingParameterError("ERROR: Missing required field `arcs`.")

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        # Adjust bounding box
        (xMin, xMax, yMin, yMax) = boundingBox
        if (xMin == None or xMax == None or yMin == None or yMax == None):
            allX = []
            allY = []
            for i in arcs:
                if (not xyReverseFlag):
                    allX.append(arcs[i]['arc'][0][0])
                    allX.append(arcs[i]['arc'][1][0])
                    allY.append(arcs[i]['arc'][0][1])
                    allY.append(arcs[i]['arc'][1][1])
                else:
                    allX.append(arcs[i]['arc'][0][1])
                    allX.append(arcs[i]['arc'][1][1])
                    allY.append(arcs[i]['arc'][0][0])
                    allY.append(arcs[i]['arc'][1][0])
            
            if (xMin == None):
                xMin = min(allX) - 0.1 * abs(max(allX) - min(allX))
            if (xMax == None):
                xMax = max(allX) + 0.1 * abs(max(allX) - min(allX))
            if (yMin == None):
                yMin = min(allY) - 0.1 * abs(max(allY) - min(allY))
            if (yMax == None):
                yMax = max(allY) + 0.1 * abs(max(allY) - min(allY))
        # Adjust width and height
        width = 0
        height = 0
        if (figSize == None or (figSize[0] == None and figSize[1] == None)):
            if (xMax - xMin > yMax - yMin):
                width = 5
                height = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                width = 5 * ((xMax - xMin) / (yMax - yMin))
                height = 5
        elif (figSize != None and figSize[0] != None and figSize[1] == None):
            width = figSize[0]
            height = figSize[0] * ((yMax - yMin) / (xMax - xMin))
        elif (figSize != None and figSize[0] == None and figSize[1] != None):
            width = figSize[1] * ((xMax - xMin) / (yMax - yMin))
            height = figSize[1]
        else:
            (width, height) = figSize

        if (isinstance(fig, plt.Figure)):
            fig.set_figwidth(width)
            fig.set_figheight(height)
            ax.set_xlim(xMin, xMax)
            ax.set_ylim(yMin, yMax)

    # Draw arcs ===============================================================
    for i in arcs:
        x1, y1, x2, y2 = 0, 0, 0, 0
        if (not xyReverseFlag):
            x1 = arcs[i]['arc'][0][0]
            x2 = arcs[i]['arc'][1][0]
            y1 = arcs[i]['arc'][0][1]
            y2 = arcs[i]['arc'][1][1]
        else:
            x1 = arcs[i]['arc'][0][1]
            x2 = arcs[i]['arc'][1][1]
            y1 = arcs[i]['arc'][0][0]
            y2 = arcs[i]['arc'][1][0]
        dx = x2 - x1
        dy = y2 - y1
        if (arcColor == 'Random'):
            rndColor = colorRandom()
            ax.plot([x1, x2], [y1, y2], color = rndColor)
            if (arrowFlag):
                ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=arcWidth, head_width=arrowHeadWidth, head_length=arrowHeadLength, color=rndColor)
        else:
            ax.plot([x1, x2], [y1, y2], color = arcColor)
            if (arrowFlag):
                ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=arcWidth, head_width=arrowHeadWidth, head_length=arrowHeadLength, color=arcColor)

        ax.plot(x1, y1, color = startColor, marker = 'o', markersize = bothEndSize)
        ax.plot(x2, y2, color = endColor, marker = 'o', markersize = bothEndSize)
        if ('label' not in arcs[i]):
            lbl = i
        else:
            lbl = arcs[i]['label']
        ha = 'left'
        if ('ha' in arcs[i]):
            ha = arcs[i]['ha']
        va = 'top'
        if ('va' in arcs[i]):
            va = arcs[i]['va']
        ax.annotate(lbl, (x1 + dx / 2, y1 + dy / 2), ha=ha, va=va)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None and isinstance(fig, plt.Figure)):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotLocSeq(
    locSeq: list[pt],
    lineColor: str = 'Random',
    lineWidth: float = 1.0,
    arrowFlag: bool = True,
    arrowHeadWidth: float = 0.1,
    arrowHeadLength: float = 0.2,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    boundingBox: tuple[int|float|None, int|float|None, int|float|None, int|float|None] = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):
    
    """Given a list of coordinates, plot a open polyline by sequences.

    Parameters
    ----------

    locSeq: list[pt], required
        A list of coordinates to form a sequence
    lineColor: string, optional, default 'Random'
        Color of lines
    lineWidth: float, optional, default 1
        Width of lines
    arrowFlag: bool, optional, default True
        True if plot arrows
    arrowHeadWidth: float, optional, default 0.1
        Width of arrow head
    arrowHeadLength: float, optional, default 0.2
        Length of arrow head
    xyReverseFlag: bool, optional, default False
        True if need to reverse the x, y coordinates, e.g., plot for (lat, lon)
    fig: matplotlib object, optional, default None
        `fig` and `ax` indicates the matplotlib object to plot on, if not provided, plot in a new figure
    ax: matplotlib object, optional, default None
        See `fig`
    figSize: 2-tuple, optional, default as (None, 5)
        Size of the figure in (width, height). If width or height is set to be None, it will be auto-adjusted.
    boundingBox: 4-tuple, optional, default as (None, None, None, None)
        (xMin, xMax, yMin, yMax), defines four boundaries of the figure
    saveFigPath: string, optional, default as None
        The path for exporting image if provided
    showFig: bool, optional, default as True
        True if show the figure in Juypter Notebook environment

    Returns
    -------
    fig, ax: matplotlib.pyplot object
    """

    # Check for required fields ===============================================
    if (locSeq == None):
        raise MissingParameterError("ERROR: Missing required field `locSeq`.")

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for i in locSeq:
            if (not xyReverseFlag):
                allX.append(i[0])
                allY.append(i[1])
            else:
                allX.append(i[1])
                allY.append(i[0])
        (xMin, xMax, yMin, yMax) = boundingBox
        if (xMin == None):
            xMin = min(allX) - 0.1 * abs(max(allX) - min(allX))
        if (xMax == None):
            xMax = max(allX) + 0.1 * abs(max(allX) - min(allX))
        if (yMin == None):
            yMin = min(allY) - 0.1 * abs(max(allY) - min(allY))
        if (yMax == None):
            yMax = max(allY) + 0.1 * abs(max(allY) - min(allY))
        width = 0
        height = 0
        if (figSize == None or (figSize[0] == None and figSize[1] == None)):
            if (xMax - xMin > yMax - yMin):
                width = 5
                height = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                width = 5 * ((xMax - xMin) / (yMax - yMin))
                height = 5
        elif (figSize != None and figSize[0] != None and figSize[1] == None):
            width = figSize[0]
            height = figSize[0] * ((yMax - yMin) / (xMax - xMin))
        elif (figSize != None and figSize[0] == None and figSize[1] != None):
            width = figSize[1] * ((xMax - xMin) / (yMax - yMin))
            height = figSize[1]
        else:
            (width, height) = figSize

        if (isinstance(fig, plt.Figure)):
            fig.set_figwidth(width)
            fig.set_figheight(height)
            ax.set_xlim(xMin, xMax)
            ax.set_ylim(yMin, yMax)

    # Draw lines ===============================================================
    if (lineColor == 'Random'):
        lineColor = colorRandom() 
    for i in range(len(locSeq) - 1):
        line = [locSeq[i], locSeq[i + 1]]
        if (not xyReverseFlag):
            x1 = line[0][0]
            y1 = line[0][1]
            x2 = line[1][0]
            y2 = line[1][1]
        else:
            x1 = line[0][1]
            y1 = line[0][0]
            x2 = line[1][1]
            y2 = line[1][0]
        dx = x2 - x1
        dy = y2 - y1
        ax.plot([x1, x2], [y1, y2], color = lineColor, linewidth=lineWidth)
        if (arrowFlag):
            ax.arrow(x=x1, y=y1, dx=dx / 2, dy=dy / 2, linewidth=lineWidth, head_width=arrowHeadWidth, head_length=arrowHeadLength, color=lineColor)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None and isinstance(fig, plt.Figure)):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotNodeSeq(
    nodes: dict, 
    nodeSeq: list[int|str],
    lineColor: str = 'Random',
    lineWidth: float = 1,
    arrowFlag: bool = True,
    arrowHeadWidth: float = 0.1,
    arrowHeadLength: float = 0.2,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    boundingBox: tuple[int|float|None, int|float|None, int|float|None, int|float|None] = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    """Given a `nodes` dictionary and a sequence of node IDs, plot a route that visits each node by IDs.

    Parameters
    ----------

    nodes: dictionary, required
        The coordinates and other attributions of the nodes to be plotted, in the following format::
            >>> nodes = {
            ...     nodeID1: {
            ...         'loc': (x, y),
            ...         'marker': 'r',    # Optional, default as 'o'
            ...         'markersize': 2,  # Optional, default as None
            ...         'color': 'red',   # Optional, default as 'Random'
            ...         'size': 3,        # Optional, default as 3
            ...         'fontsize': 3,    # Optional, default as 3
            ...         'neighbor': poly, # Optional, indicate if need to display the neighborhood
            ...     }, # ...
            ... }
    nodeSeq: list[int|str], required
        A list of nodeIDs which will form a visiting sequence
    arcColor: string, optional, default 'Random'
        Color of arcs
    arcWidth: float, optional, default 1
        Width of arcs
    arrowFlag: bool, optional, default True
        True if plot arrows
    arrowHeadWidth: float, optional, default 0.1
        Width of arrow head
    arrowHeadLength: float, optional, default 0.2
        Length of arrow head
    xyReverseFlag: bool, optional, default False
        True if need to reverse the x, y coordinates, e.g., plot for (lat, lon)
    fig: matplotlib object, optional, default None
        `fig` and `ax` indicates the matplotlib object to plot on, if not provided, plot in a new figure
    ax: matplotlib object, optional, default None
        See `fig`
    figSize: 2-tuple, optional, default as (None, 5)
        Size of the figure in (width, height). If width or height is set to be None, it will be auto-adjusted.
    boundingBox: 4-tuple, optional, default as (None, None, None, None)
        (xMin, xMax, yMin, yMax), defines four boundaries of the figure
    saveFigPath: string, optional, default as None
        The path for exporting image if provided
    showFig: bool, optional, default as True
        True if show the figure in Juypter Notebook environment

    Returns
    -------
    fig, ax: matplotlib.pyplot object
    """

    # Create arcs =============================================================
    if (nodeSeq == None):
        raise MissingParameterError("ERROR: Missing required field `nodeSeq`.")
    arcs = []
    for i in range(len(nodeSeq) - 1):
        arcs.append([nodeSeq[i], nodeSeq[i + 1]])

    # Color ===================================================================
    if (lineColor == 'Random'):
        lineColor = colorRandom()

    # Call plotArcs ===========================================================

    locSeq = []
    for n in nodeSeq:
        if (n not in nodes):
            raise UnsupportedInputError("ERROR: Cannot find %s in nodes" % n)
        locSeq.append(nodes[n]['loc'])

    fig, ax = plotLocSeq(
        fig = fig,
        ax = ax,
        locSeq = locSeq,
        lineColor = lineColor,
        lineWidth = lineWidth,
        arrowFlag = arrowFlag,
        arrowHeadWidth = arrowHeadWidth,
        arrowHeadLength = arrowHeadLength,
        xyReverseFlag = xyReverseFlag,
        figSize = figSize,
        boundingBox = boundingBox,
        showAxis = showAxis,
        saveFigPath = saveFigPath,
        showFig = showFig)

    return fig, ax

def plotPolygon(
    poly: poly, 
    edgeWidth: float = 0.5,
    edgeColor: str|None = 'Random',
    fillColor: str|None = None,
    fillStyle: str = "///",
    opacity: float = 0.5,
    xyReverseFlag: bool = False,    
    fig = None,
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    boundingBox: tuple[int|float|None, int|float|None, int|float|None, int|float|None] = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    """Draw arcs

    Parameters
    ----------

    poly: poly, required, default None
        A polygon to be plotted
    edgeWidth: float, optional, default 0.5
        Width of the edge
    edgeColor: string, optional, default 'Random'
        Color of the edge
    fillColor: string, optional, default None
        Color filled in the polygon
    fillStyle: string, optional, default "///"
        Style filled in the polygon
    opacity: float, optional, default 0.5
        Opacity of the polygon
    xyReverseFlag: bool, optional, default False
        True if need to reverse the x, y coordinates, e.g., plot for (lat, lon)
    fig: matplotlib object, optional, defaut None
        `fig` and `ax` indicates the matplotlib object to plot on, if not provided, plot in a new figure
    ax: matplotlib object, optional, default None
        See `fig`
    figSize: 2-tuple, optional, default as (None, 5)
        Size of the figure in (width, height). If width or height is set to be None, it will be auto-adjusted.
    boundingBox: 4-tuple, optional, default as (None, None, None, None)
        (xMin, xMax, yMin, yMax), defines four boundaries of the figure
    saveFigPath: string, optional, default as None
        The path for exporting image if provided
    showFig: bool, optional, default as True
        True if show the figure in Juypter Notebook environment

    Returns
    -------
    fig, ax: matplotlib.pyplot object
    """

    # Check for required fields ===============================================
    if (poly == None):
        raise MissingParameterError("ERROR: Missing required field `poly`.")

    # If no based matplotlib figure provided, define boundary =================
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
        (xMin, xMax, yMin, yMax) = boundingBox
        if (xMin == None):
            xMin = min(allX) - 0.1 * abs(max(allX) - min(allX))
        if (xMax == None):
            xMax = max(allX) + 0.1 * abs(max(allX) - min(allX))
        if (yMin == None):
            yMin = min(allY) - 0.1 * abs(max(allY) - min(allY))
        if (yMax == None):
            yMax = max(allY) + 0.1 * abs(max(allY) - min(allY))
        width = 0
        height = 0
        if (figSize == None or (figSize[0] == None and figSize[1] == None)):
            if (xMax - xMin > yMax - yMin):
                width = 5
                height = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                width = 5 * ((xMax - xMin) / (yMax - yMin))
                height = 5
        elif (figSize != None and figSize[0] != None and figSize[1] == None):
            width = figSize[0]
            height = figSize[0] * ((yMax - yMin) / (xMax - xMin))
        elif (figSize != None and figSize[0] == None and figSize[1] != None):
            width = figSize[1] * ((xMax - xMin) / (yMax - yMin))
            height = figSize[1]
        else:
            (width, height) = figSize

        if (isinstance(fig, plt.Figure)):
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

    if (abs(x[-1] - x[0]) > CONST_EPSILON):
        x.append(x[0])
        y.append(y[0])

    # Plot ====================================================================
    if (edgeColor == 'Random'):
        edgeColor = colorRandom()
    if (fillColor == None):
        ax.plot(x, y, color = edgeColor, linewidth = edgeWidth)
    else:
        if (fillColor == 'Random'):
            ax.fill(x, y, facecolor=colorRandom(), edgecolor=edgeColor, hatch=fillStyle, linewidth=edgeWidth, alpha=opacity)
        else:
            ax.fill(x, y, facecolor=fillColor, edgecolor=edgeColor, hatch=fillStyle, linewidth=edgeWidth, alpha=opacity)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None and isinstance(fig, plt.Figure)):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotArcSegs(
    arcSegs: list[arcSeg],
    arcSegColor: str = 'Random',
    arcSegWidth: float = 1.0,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    boundingBox: tuple[int|float|None, int|float|None, int|float|None, int|float|None] = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ): 
    return fig, ax

def plotArcPoly():
    return fig, ax

def plotProvinceMap(
    country: str = 'U.S.',
    province: list[str]|str = [],
    edgeWidth: float = 0.5,
    edgeColor: str = 'Random',
    fillColor: str|None = None,
    fillStyle: str = "///",
    opacity: float = 0.5,   
    fig = None,
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    """Draw arcs

    Parameters
    ----------


    country: string, required, default 'U.S.'
        Country of the province
    province: string | list[string], required, default ['New York']
        A province or a list of provinces to be plotted.
    edgeWidth: float, optional, default 0.5
        Width of the edge
    edgeColor: string, optional, default 'Random'
        Color of the edge
    fillColor: string, optional, default None
        Color filled in the polygon
    fillStyle: string, optional, default "///"
        Style filled in the polygon
    opacity: float, optional, default 0.5
        Opacity of the polygon
    xyReverseFlag: bool, optional, default False
        True if need to reverse the x, y coordinates, e.g., plot for (lat, lon)
    fig: matplotlib object, optional, defaut None
        `fig` and `ax` indicates the matplotlib object to plot on, if not provided, plot in a new figure
    ax: matplotlib object, optional, default None
        See `fig`
    figSize: 2-tuple, optional, default as (None, 5)
        Size of the figure in (width, height). If width or height is set to be None, it will be auto-adjusted.    
    saveFigPath: string, optional, default as None
        The path for exporting image if provided
    showFig: bool, optional, default as True
        True if show the figure in Juypter Notebook environment

    Returns
    -------
    fig, ax: matplotlib.pyplot object
    """

    if (fig == None or ax == None):
        fig, ax = plt.subplots()
    
    if (type(province) == str):
        province = [province]
    for prv in province:
        prvPoly = []
        if (country == 'U.S.'):
            if prv in usState:
                prvPoly = usState[prv]
            elif prv in usStateAbbr:
                prvPoly = usState[usStateAbbr[prv]]
        else:
            raise VrpSolverNotAvailableError("Error: %s is not included yet, please stay tune." % country) 
        fig, ax = plotPolygon(
            fig = fig,
            ax = ax,
            poly = prvPoly,
            edgeWidth = edgeWidth,
            edgeColor = edgeColor,
            fillColor = fillColor,
            fillStyle = fillStyle,
            opacity = opacity,
            xyReverseFlag = True,
            figSize=figSize,
            saveFigPath=saveFigPath,
            showAxis = showAxis,
            showFig=showFig)

    return fig, ax

def plotRoadNetwork(
    roads: dict,
    roadWidth: dict[str,float]|float = {
            'motorway': 2,
            'motorway_link': 2,
            'truck': 1.5,
            'truck_link': 1.5,
            'primary': 1.2,
            'primary_link': 1.2,
            'secondary': 1,
            'secondary_link': 1,
            'tertiary': 1,
            'tertiary_link': 1,
            'residential': 0.8,
            'others': 0.8
        },
    roadColors: dict[str,str]|str = {
            'motorway': 'red',
            'motorway_link': 'red',
            'truck': 'orange',
            'truck_link': 'orange',
            'primary': 'orange',
            'primary_link': 'orange',
            'secondary': 'orange',
            'secondary_link': 'orange',
            'tertiary': 'orange',
            'tertiary_link': 'orange',
            'residential': 'green',
            'others': 'gray'
        },
    roadShowFlags: list[str, bool]|str|bool = 'All',
    bldColors: dict[str,str]|str = {
            'building': 'yellow',
            'commercial': 'yellow',
            'residential': 'green',
            'house': 'green',
            'static_caravan': 'green',
            'industrial': 'orange',
            'manufacture': 'orange'
        },
    bldShowFlags: list[str, bool]|str|bool = False,
    fig = None,
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    boundingBox: tuple[int|float|None, int|float|None, int|float|None, int|float|None] = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ): 

    """Plot road network (and buildings) using given OSM transformed dictionary

    Parameters
    ----------

    roads: dict, required
        The road network dictionary, including the geometry shape. In the following format::
            >>> roads = {
            ...     'road': {roadID: {'class': class, 'shape': shape}}
            ...     'building': {buildingID: {'type': type, 'shape': shape}}
            ... }
    roadWidth: dict[str, float]|float, optional
        The width of roads. If a dictionary is given, will use the width in the dictionary, otherwise will be default as 1. If a float is given will use the float value.
    roadColor: dict[str, str]|str, optional
        The color of roads, works in the same way as `roadWidth`
    roadShowFlags: dict[str, bool]|str|bool, optional, 'All'
        Whether or not show some types of roads, works in the same way as `roadWidth`
    bldColor: dict[str, str]|str, optional
        The color of buildings, works in the same way as `roadWidth`
    bldShowFlags: dict[str, bool]|str|bool, optional
        Whether or not show some types of buildings, works in the same way as `roadWidth`
    fig: matplotlib object, optional, defaut None
        `fig` and `ax` indicates the matplotlib object to plot on, if not provided, plot in a new figure
    ax: matplotlib object, optional, default None
        See `fig`
    figSize: 2-tuple, optional, default as (None, 5)
        Size of the figure in (width, height). If width or height is set to be None, it will be auto-adjusted.
    boundingBox: 4-tuple, optional, default as (None, None, None, None)
        (xMin, xMax, yMin, yMax), defines four boundaries of the figure
    saveFigPath: string, optional, default as None
        The path for exporting image if provided
    showFig: bool, optional, default as True
        True if show the figure in Juypter Notebook environment
    """

    # FIXME: In future, we might want to distinguish roads by max speed or show the names of roads
    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for pt in roads['boundary']:
            allX.append(pt[1])
            allY.append(pt[0])
        (xMin, xMax, yMin, yMax) = boundingBox
        if (xMin == None):
            xMin = min(allX) - 0.1 * abs(max(allX) - min(allX))
        if (xMax == None):
            xMax = max(allX) + 0.1 * abs(max(allX) - min(allX))
        if (yMin == None):
            yMin = min(allY) - 0.1 * abs(max(allY) - min(allY))
        if (yMax == None):
            yMax = max(allY) + 0.1 * abs(max(allY) - min(allY))
        width = 0
        height = 0
        if (figSize == None or (figSize[0] == None and figSize[1] == None)):
            if (xMax - xMin > yMax - yMin):
                width = 5
                height = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                width = 5 * ((xMax - xMin) / (yMax - yMin))
                height = 5
        elif (figSize != None and figSize[0] != None and figSize[1] == None):
            width = figSize[0]
            height = figSize[0] * ((yMax - yMin) / (xMax - xMin))
        elif (figSize != None and figSize[0] == None and figSize[1] != None):
            width = figSize[1] * ((xMax - xMin) / (yMax - yMin))
            height = figSize[1]
        else:
            (width, height) = figSize

        if (isinstance(fig, plt.Figure)):
            fig.set_figwidth(width)
            fig.set_figheight(height)
            ax.set_xlim(xMin, xMax)
            ax.set_ylim(yMin, yMax)

    # Plot roads ==============================================================
    for road in roads['road']:
        if (roadShowFlags == 'All' or roadShowFlags == True 
            or (type(roadShowFlags) == list 
                and roads['road'][road]['class'] in roadShowFlags)):
            x = []
            y = []
            for pt in roads['road'][road]['shape']:
                x.append(pt[1])
                y.append(pt[0])
            color = None
            if (roadColors == 'Random'):
                color = colorRandom()
            elif (type(roadColors) == str):
                color = roadColors
            elif (type(roadColors) == dict):
                if (roads['road'][road]['class'] in roadColors):
                    color = roadColors[roads['road'][road]['class']]
                else:
                    color = 'gray'
            rw = 1
            if (type(roadWidth) == dict):
                if (roads['road'][road]['class'] in roadWidth):
                    rw = roadWidth[roads['road'][road]['class']]
            else:
                rw = roadWidth
            ax.plot(x, y, color = color, linewidth = rw)

    # Plot buildings ==========================================================
    for building in roads['building']:
        if (bldShowFlags == 'All' or bldShowFlags == True
            or (type(bldShowFlags) == list
                and roads['building'][building]['type'] in bldShowFlags)):
            x = []
            y = []
            color = None
            for pt in roads['building'][building]['shape']:
                x.append(pt[1])
                y.append(pt[0])
            if (bldColors == 'Random'):
                color = colorRandom()
            elif (type(bldColors) == str):
                color = bldColors
            elif (type(bldColors) == dict):
                if (roads['building'][building]['type'] in bldColors):
                    color = bldColors[roads['building'][building]['type']]
                else:
                    color = 'gray'
            ax.fill(x, y, facecolor=color)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotGantt(
    gantt: list[dict],
    group: list[dict] | None = None, 
    phase: list[dict] | None = None,
    xlabel: str = "Time",
    ylabel: str = "",
    ganttHeight: float = 0.8,
    ganttLinewidth: float = 1.0,
    ganttBlockHeight: float = 1.0,
    groupSeparatorHeight: float = 0.5,
    groupSeparatorStyle: str = "-",
    groupSeparatorLinewidth: float = 0.5,
    groupSeparatorColor: str = 'gray',
    phaseSeparatorStyle: str = "--",
    phaseSeparatorLinewidth: float = 0.5,
    phaseSeparatorColor: str = 'gray',
    startTime: float = 0,
    endTime: float|None = None,
    fig = None, 
    ax = None,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (12, 5), 
    showGridFlag: bool = False,
    saveFigPath: str|None = None,
    showTailFlag: bool = True,
    showFig: bool = True
    ) -> "Given a Gantt dictionary, plot Gantt":

    """Given a Gantt dictionary, plot a Gantt chart

    Parameters
    ----------
    gantt: list of dictionaries, required
        A list of dictionaries, each represents a gantt block, in the following format\
            >>> gantt = [{
            ...     'entityID': entityID, 
            ...     'timeWindow': [startTime, endTime], 
            ...     'desc': (optional) description of the window,
            ...     'color': (optional, default as 'random') color, 
            ...     'style': (optional, default as 'solid') 'solid' 
            ... }, ... ]
    group: list of dictionaries, optional, default []
        Groups of entityIDs, in the following format
            >>> group = [{
            ...     'title': groupTitle,
            ...     'entities': [entityID1, entityID2, ...],
            ...     'showFlag': True,
            ...     'backgroundColor': 'gray',
            ...     'backgroundOpacity': 0.5
            ... }, ... ]
    phase: list of dictionaries, optional, default []
        A list of phases, in the following format
            >>> phase = [{
            ...     'title': phaseTitle,
            ...     'timeWindow': [startTime, endTime],
            ...     'backgroundColor': 'gray',
            ...     'backgroundOpacity': 0.5
            ... }, ... ]
    xlabel: string, optional, default "Time"
        The label of x-axis
    ylabel: string, optional, default ""
        The label of y-axis
    fig: matplotlib object, optional, default None
        `fig` and `ax` indicates the matplotlib object to plot on, if not provided, plot in a new figure
    ax: matplotlib object, optional, default None
        See `fig`
    figSize: 2-tuple, optional, default as (None, 5)
        Size of the figure in (width, height). If width or height is set to be None, it will be auto-adjusted.    
    showGridFlag: bool, optional, default as False
        True if turn the grids on.
    saveFigPath: string, optional, default as None
        The path for exporting image if provided
    showFig: bool, optional, default as True
        True if show the figure in Juypter Notebook environment
    """

    # Check for required fields ===============================================
    if (gantt == None):
        raise MissingParameterError(ERROR_MISSING_GANTT)

    # Pre-calculate ===========================================================
    realStart = None
    realEnd = None
    for g in gantt:
        if ('entityID' not in g or 'timeWindow' not in g):
            raise MissingParameterError("ERROR: Missing 'entityID' and 'timeWindow' in gantt.")

        if (realStart == None or realStart > g['timeWindow'][0]):
            realStart = g['timeWindow'][0]
        if (realEnd == None or realEnd < g['timeWindow'][1]):
            realEnd = g['timeWindow'][1]
    if (startTime == None):
        startTime = realStart
    if (endTime == None):
        endTime = realEnd

    # Arrange entities ========================================================
    if (group == None or group == []):
        group = [{
            'title': "",
            'entities': [],
            'showFlag': True,
            'backgroundColor': None,
            'backgroundOpacity': 0
        }]
        for g in gantt:
            if (g['entityID'] not in group[0]['entities']):
                group[0]['entities'].append(g['entityID'])

    entities = []
    for g in group:
        if (g['showFlag']):
            entities.extend(g['entities'])
            entities.append(None)
    if (entities[-1] == None):
        entities = entities[:-1]
    entities.reverse()

    totalHeight = len([i for i in entities if i != None]) * ganttBlockHeight + len([i for i in entities if i == None]) * groupSeparatorHeight
    yAcc = 0.5 * ganttBlockHeight
    yticks = [0.5 * ganttBlockHeight]
    for i in range(1, len(entities)):        
        if (entities[i - 1] != None and entities[i] != None):
            yAcc += ganttBlockHeight
            yticks.append(yAcc)
        elif (entities[i - 1] == None or entities[i] == None):
            yAcc += 0.5 * groupSeparatorHeight + 0.5 * ganttBlockHeight
            yticks.append(yAcc)

    # If no based matplotlib figure, define fig size ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        width = 0
        height = 0
        if (figSize == None or (figSize[0] == None and figSize[1] == None)):
            width = 5 * ((endTime - startTime) / totalHeight)
            height = 5
        elif (figSize != None and figSize[0] != None and figSize[1] == None):
            width = figSize[0]
            height = figSize[0] * (totalHeight / (endTime - startTime))
        elif (figSize != None and figSize[0] == None and figSize[1] != None):
            width = figSize[1] * ((endTime - startTime) / totalHeight)
            height = figSize[1]
        else:
            (width, height) = figSize
        if (isinstance(fig, plt.Figure)):
            fig.set_figwidth(width)
            fig.set_figheight(height)
            ax.set_xlim(startTime, endTime + (endTime - startTime) * 0.05)
            ax.set_ylim(0, totalHeight)
    ax.set_yticks(yticks)   
    ax.set_yticklabels(entities)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    # Plot groups =============================================================
    groupColor = []
    groupOpacity = []
    groupTitle = []
    for i in range(len(group)):
        if (group[i]['showFlag']):
            if (group[i]['backgroundColor'] != 'Random'):
                groupColor.append(group[i]['backgroundColor'])
            else:
                groupColor.append(colorRandom())
            groupOpacity.append(group[i]['backgroundOpacity'])
            groupTitle.append(group[i]['title'])
    groupStyle = []
    yS = 0
    yE = 0
    for i in range(len(entities)):
        if (entities[i] == None):
            yE = yticks[i]
            groupStyle.append([yS, yE])
            yS = yticks[i]
    groupStyle.append([yS, yticks[-1] + 0.5 * ganttBlockHeight])
    for g in range(len(groupStyle)):
        s = groupStyle[g][0]
        e = groupStyle[g][1]
        ax.fill([startTime, startTime, endTime + (endTime - startTime) * 0.05, endTime + (endTime - startTime) * 0.05, startTime], 
            [s, e, e, s, s], color=groupColor[g], alpha=groupOpacity[g])
        ax.annotate(groupTitle[g], (endTime + (endTime - startTime) * 0.05, e), ha='right', va='top')

    # Plot phases =============================================================
    phaseSep = []
    phaseAnchor = []
    phaseTitle = []
    if (phase != None and phase != []):
        if (phaseSeparatorStyle != None):
            for i in range(1, len(phase)):
                sep = phase[i - 1]['timeWindow'][1] + (phase[i]['timeWindow'][0] - phase[i - 1]['timeWindow'][1]) / 2
                ax.plot([sep, sep], [0, totalHeight], color = phaseSeparatorColor, linewidth = phaseSeparatorLinewidth, linestyle = phaseSeparatorStyle)
            ax.plot([phase[-1]['timeWindow'][1], phase[-1]['timeWindow'][1]], [0, totalHeight], color = phaseSeparatorColor, linewidth = phaseSeparatorLinewidth, linestyle = phaseSeparatorStyle)
        for i in range(len(phase)):
            if (phase[i]['backgroundColor'] != None):
                ax.fill([phase[i]['timeWindow'][0], phase[i]['timeWindow'][0], phase[i]['timeWindow'][1], phase[i]['timeWindow'][1], phase[i]['timeWindow'][0]], 
                    [0, totalHeight, totalHeight, 0, 0], 
                    color=phase[i]['backgroundColor'] if phase[i]['backgroundColor'] != 'Random' else colorRandom(), 
                    alpha=phase[i]['backgroundOpacity'])
                ax.annotate(phase[i]['title'], (phase[i]['timeWindow'][0] + (phase[i]['timeWindow'][1] - phase[i]['timeWindow'][0]) / 2, totalHeight), ha='center', va='bottom')

    # Plot each gantt block ===================================================
    for i in range(len(entities)):
        if (entities[i] != None):
            center = yticks[i]
            top = None
            bottom = yticks[i] - ganttHeight / 2
            ent = [g for g in gantt if g['entityID'] == entities[i]]
            for g in ent:
                s = g['timeWindow'][0]
                e = g['timeWindow'][1]
                x = [s, s, e, e, s]
                if ('desc' not in g or ('descPosition' in g and g['descPosition'] == 'Inside')):
                    top = yticks[i] + ganttHeight / 2
                else:
                    top = yticks[i] + ganttHeight / 4
                y = [bottom, top, top, bottom, bottom]
                ax.plot(x, y, color = 'black', linewidth = ganttLinewidth)
                if ('color' in g or g['color'] != 'random'):
                    ax.fill(x, y, color = g['color'], linewidth = ganttLinewidth)
                else:
                    rndColor = colorRandom()
                    ax.fill(x, y, color = rndColor, linewidth = ganttLinewidth)
                if ('desc' in g):
                    if ('descPosition' not in g or g['descPosition'] == 'Top'):
                        ax.annotate(g['desc'], (s + ganttHeight / 8, top + ganttHeight / 8))
                    elif ('descPosition' in g and g['descPosition'] == 'Inside'):
                        ax.annotate(g['desc'], (s + ganttHeight / 8, bottom + ganttHeight / 8))
                if (g['style'] != 'solid'):
                    ax.fill(x, y, hatch = g['style'], fill=False, linewidth = ganttLinewidth)
        elif (groupSeparatorStyle != None):
            center = yticks[i]
            ax.plot((startTime, endTime + (endTime - startTime) * 0.05), (center, center), linestyle = groupSeparatorStyle, linewidth = groupSeparatorLinewidth, color = groupSeparatorColor)

    # Show time span ==========================================================
    if (showTailFlag):
        xTicks = list(ax.get_xticks())
        xTicks.append(realEnd)
        ax.set_xticks(xTicks)  

    # Grids ===================================================================
    if (showGridFlag):
        ax.grid(b = True, linestyle=':')

    # Fix height if fig, ax are not provided ==================================
    if (fig == None or ax == None):
        fig.set_figheight(5 * max(pos))

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax
