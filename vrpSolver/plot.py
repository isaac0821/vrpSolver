import matplotlib.pyplot as plt

from .common import *
from .color import *
from .msg import *
from .province import *
from .error import *
from .geometry import *

# History =====================================================================
# 20230518 - `plotNodes()` now will plot the neighborhood of nodes
# 20230624 - Rename functions `plotArcs()`, `plotLocSeq()`, `plotNodeSeq()`
# =============================================================================

def plotNodes(
    nodes: dict, 
    locFieldName = 'loc',
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
                    allX.append(nodes[i][locFieldName][0])
                    allY.append(nodes[i][locFieldName][1])
                else:
                    allX.append(nodes[i][locFieldName][1])
                    allY.append(nodes[i][locFieldName][0])        
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
            x = nodes[n][locFieldName][0]
            y = nodes[n][locFieldName][1]
        else:
            x = nodes[n][locFieldName][1]
            y = nodes[n][locFieldName][0]
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
    arcFieldName = 'arc',
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
                    allX.append(arcs[i][arcFieldName][0][0])
                    allX.append(arcs[i][arcFieldName][1][0])
                    allY.append(arcs[i][arcFieldName][0][1])
                    allY.append(arcs[i][arcFieldName][1][1])
                else:
                    allX.append(arcs[i][arcFieldName][0][1])
                    allX.append(arcs[i][arcFieldName][1][1])
                    allY.append(arcs[i][arcFieldName][0][0])
                    allY.append(arcs[i][arcFieldName][1][0])
            
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
            x1 = arcs[i][arcFieldName][0][0]
            x2 = arcs[i][arcFieldName][1][0]
            y1 = arcs[i][arcFieldName][0][1]
            y2 = arcs[i][arcFieldName][1][1]
        else:
            x1 = arcs[i][arcFieldName][0][1]
            x2 = arcs[i][arcFieldName][1][1]
            y1 = arcs[i][arcFieldName][0][0]
            y2 = arcs[i][arcFieldName][1][0]
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
    locFieldName: str = 'loc',
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
        locSeq.append(nodes[n][locFieldName])

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

def plotPoly(
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

    """Draw a polygon

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

def plotCircle(
    center: pt, 
    radius: float,
    lod: int = 30,
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

    polyCircle = [[
        center[0] + radius * math.sin(2 * d * math.pi / lod),
        center[1] + radius * math.cos(2 * d * math.pi / lod),
    ] for d in range(lod + 1)]

    fig, ax = plotPoly(
        poly = polyCircle,
        edgeWidth = edgeWidth,
        edgeColor = edgeColor,
        fillColor = fillColor,
        fillStyle = fillStyle,
        opacity = opacity,
        xyReverseFlag = xyReverseFlag,
        fig = fig,
        ax = ax,
        figSize = figSize,
        boundingBox = boundingBox,
        showAxis = showAxis,
        saveFigPath = saveFigPath,
        showFig = showFig,
    )
    return fig, ax
    
def plotPolys(
    polys: dict,
    polyFieldName: str = 'poly',
    showAnchorFlag: bool = True,
    anchorFieldName: str = 'anchor',
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

    # Sanity check ============================================================
    if (type(polys) != dict):
        raise UnsupportedInputError("ERROR: `polys` is a dictionary, to plot an individual polygon, please use plotPoly() instead.")

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for p in polys:
            for i in polys[p][polyFieldName]:
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

    # First, plot polygons ====================================================
    for p in polys:
        fig, ax = plotPoly(
            fig = fig,
            ax = ax,
            poly = polys[p][polyFieldName],
            edgeWidth = edgeWidth,
            edgeColor = edgeColor,
            fillColor = fillColor,
            fillStyle = fillStyle,
            opacity = opacity,
            xyReverseFlag = xyReverseFlag,
            figSize = figSize)

    # Next, plot anchors ======================================================
    if (showAnchorFlag):
        for p in polys:
            if ('label' not in polys[p]):
                lbl = p
            else:
                lbl = polys[p]['label']
            ha = 'center'
            va = 'center'
            ct = polys[p]['anchor']
            if (xyReverseFlag):
                ct = [ct[1], ct[0]]
            ax.annotate(lbl, ct, ha=ha, va=va)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None and isinstance(fig, plt.Figure)):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

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
    
    prvPoly = {}
    if (type(province) == str):
        province = [province]
    for prv in province:
        if (country == 'U.S.'):
            if prv in usState:
                ct = ptPolyCenter(usState[prv])
                prvPoly[prv] = {
                    'anchor': ct,
                    'poly': usState[prv]
                }                
            elif prv in usStateAbbr:
                ct = ptPolyCenter(usState[usStateAbbr[prv]])
                prvPoly[prv] = {
                    'anchor': ct,
                    'poly': usState[usStateAbbr[prv]]
                }
        else:
            raise VrpSolverNotAvailableError("Error: %s is not included yet, please stay tune." % country) 
    fig, ax = plotPolys(
        fig = fig,
        ax = ax,
        polys = prvPoly,
        showAnchorFlag = True,
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

def aniRouting(
    timeRange: tuple[int, int],
    nodes: dict|None=None,
    locFieldName: str = 'loc',
    vehicles: dict|None = None,
    polys: dict|None = None,
    polyFieldName = 'poly',
    nodeColor: str = 'black',
    nodeMarker: str = 'o',
    nodeMarkersize: float = 2,
    nodeNeighborColor: str|None = 'gray',
    nodeNeighborOpacity: float = 0.5,
    vehColor: str = 'blue',
    vehPathColor: str|None = 'gray',
    vehTraceColor: str|None = 'orange',
    vehTraceTime: float|None = None,
    vehShowSpeedFlag: bool = True,
    vehShowNoteFlag: bool = True,
    polyEdgeColor: str = 'black',
    polyEdgeWidth: float = 1,
    polyFillColor: str|None = 'gray',
    polyFillStyle: str|None = '///',
    polyOpacity: float = 0.5,
    speed: int = 1,
    fps: int = 1,
    repeatFlag: bool = True,
    xyReverseFlag: bool = False,
    figSize: list[int|float|None] | tuple[int|float|None, int|float|None] = (None, 5), 
    boundingBox: tuple[int|float|None, int|float|None, int|float|None, int|float|None] = (None, None, None, None),
    aniSavePath: str|None = None,
    aniSaveDPI: int = 300
    ):

    """Given nodes, vehicles, static polygons, and dynamic polygons, create animation

    Parameters
    ----------

    nodes: dictionary, optional, default None
        A dictionary of nodes.
            >>> nodes[nID] = {
            ...     'loc': [x, y],
            ...     'neighbor': poly, # A polygon indicating neighborhood
            ...     'direction': direction, # Moving direction
            ...     'speed': speed, # Moving speed
            ...     'marker': 'r',    # Optional, default as 'o'
            ...     'markersize': 2,  # Optional, default as None
            ...     'color': 'red',   # Optional, default as 'Random'
            ...     'size': 3,        # Optional, default as 3
            ...     'fontsize': 3,    # Optional, default as 3
            ... }
    vehicles: dictionary, optional, default None
        A dictionary of vehicle. 
            >>> vehicles[vID] = {
            ...     'vehicle': vehicleName,
            ...     'seq': [], # A sequence of visiting locations
            ...     'timeStamp': [], # A sequence of time stamps when visiting each location
            ...     'note': [], # Note to show during each leg of sequence
            ...     'color': vehColor,
            ...     'pathColor': pathColor,
            ...     'traceColor': traceColor,
            ...     'traceTime': traceTime
            ... }
    polys: dictionary, optional, default None
        A dictionary indicating polygons that are dynamic/static in the animation
            >>> polys[pID] = {
            ...     'anchor': [x, y], # The anchor of the polygon
            ...     'poly': [pt1, pt2], # A sequence of extreme points, coordinates are relative to 'anchor'
            ...     'direction': direction, # Moving direction
            ...     'speed': speed, # Moving speed,
            # ...     'clockwise': None, # True if clockwise, False if counter-clockwise, None if not rotating
            # ...     'rSpeed': radSpeed, # Rad speed,
            ...     'timeRange': [ts, te], # Time range of the movement
            ...     'edgeColor': color, # Edge color
            ...     'fillColor': color, # Fill color
            ...     'fillStyle': "///", # Fill style
            ...     'opacity': 0.5, # opacity
            ... }

    """

    # Check for required fields ===============================================
    if (nodes == None and vehicles == None and polys == None):
        raise MissingParameterError("ERROR: Need to provide entities for animation.")

    # If no based matplotlib figure provided, define boundary =================
    fig, ax = plt.subplots()
    allX = []
    allY = []
    if (nodes != None):
        for i in nodes:
            if (not xyReverseFlag):
                allX.append(nodes[i][locFieldName][0])
                allY.append(nodes[i][locFieldName][1])
            else:
                allX.append(nodes[i][locFieldName][1])
                allY.append(nodes[i][locFieldName][0])
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

    # Styling
    nodeStyle = {}
    if (nodes != None):
        for nID in nodes:
            nodeStyle[nID] = {}
            if (nodeColor == 'Random'):
                nodeStyle[nID]['nodeColor'] = vrpSolver.colorRandom()
            elif (nodeColor != None):
                nodeStyle[nID]['nodeColor'] = nodeColor
            elif ('nodeColor' in nodes[nID]):
                nodeStyle[nID]['nodeColor'] = nodes[nID]['nodeColor']
            
            if (nodeMarker != None):
                nodeStyle[nID]['nodeMarker'] = nodeMarker
            elif ('nodeMarker' in nodes[nID]):
                nodeStyle[nID]['nodeMarker'] = nodes[nID]['nodeMarker']

            if (nodeMarkersize != None):
                nodeStyle[nID]['nodeMarkersize'] = nodeMarkersize
            elif ('nodeMarkersize' in nodes[nID]):
                nodeStyle[nID]['nodeMarkersize'] = nodes[nID]['nodeMarkersize']

            if (nodeNeighborColor == 'Random'):
                nodeStyle[nID]['nodeNeighborColor'] = vrpSolver.colorRandom()
            elif (nodeNeighborColor != None):
                nodeStyle[nID]['nodeNeighborColor'] = nodeNeighborColor
            elif ('nodeNeighborColor' in nodes[nID]):
                nodeStyle[nID]['nodeNeighborColor'] = nodes[nID]['nodeNeighborColor']

    polyStyle = {}
    if (polys != None):
        for pID in polys:
            polyStyle[pID] = {}
            if (polyEdgeColor == 'Random'):
                polyStyle[pID]['edgeColor'] = vrpSolver.colorRandom()
            elif (polyEdgeColor != None):
                polyStyle[pID]['edgeColor'] = polyEdgeColor
            elif ('edgeColor' in polys[pID]):
                polyStyle[pID]['edgeColor'] = polys[pID]['edgeColor']

            if (polyEdgeWidth != None):
                polyStyle[pID]['edgeWidth'] = polyEdgeWidth
            elif ('edgeWidth' in polys[pID]):
                polyStyle[pID]['edgeWidth'] = polys[pID]['edgeWidth']

            if (polyFillColor == 'Random'):
                polyStyle[pID]['fillColor'] = vrpSolver.colorRandom()
            elif (polyFillColor != None):
                polyStyle[pID]['fillColor'] = polyFillColor
            elif ('fillColor' in polys[pID]):
                polyStyle[pID]['fillColor'] = polys[pID]['fillColor']

            if (polyFillStyle != None):
                polyStyle[pID]['fillStyle'] = polyFillStyle
            elif ('fillStyle' in polys[pID]):
                polyStyle[pID]['fillStyle'] = polys[pID]['fillStyle']

            if (polyOpacity != None):
                polyStyle[pID]['opacity'] = polyOpacity
            elif ('opacity' in polys[pID]):
                polyStyle[pID]['opacity'] = polys[pID]['opacity']

    vehicleStyle = {}
    if (vehicles != None):
        for vID in vehicles:
            vehicleStyle[vID] = {}
            if (vehColor == 'Random'):
                vehicleStyle[vID]['vehColor'] = vrpSolver.colorRandom()
            elif (vehColor != None):
                vehicleStyle[vID]['vehColor'] = vehColor
            elif ('color' in vehicles[vID]):
                vehicleStyle[vID]['vehColor'] = vehicles[vID]['color']

            if (vehPathColor == 'Random'):
                vehicleStyle[vID]['pathColor'] = vrpSolver.colorRandom()
            elif (vehPathColor != None):
                vehicleStyle[vID]['pathColor'] = vehPathColor
            elif ('pathColor' in vehicles[vID]):
                vehicleStyle[vID]['pathColor'] = vehicles[vID]['pathColor']

            if (vehTraceColor == 'Random'):
                vehicleStyle[vID]['traceColor'] = vrpSolver.colorRandom()
            elif (vehTraceColor != None):
                vehicleStyle[vID]['traceColor'] = vehTraceColor
            elif ('traceColor' in vehicles[vID]):
                vehicleStyle[vID]['traceColor'] = vehicles[vID]['traceColor']

            if (vehTraceTime != None):
                vehicleStyle[vID]['traceTime'] = vehTraceTime
            elif ('traceTime' in vehicles[vID]):
                vehicleStyle[vID]['traceTime'] = vehicles[vID]['traceTime']
            else:
                vehicleStyle[vID]['traceTime'] = None

    def animate(t):
        ax.clear()
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

        # Clock
        clock = timeRange[0] + t * speed / (fps * 1.0)
        ax.set_title("Clock: %s[s]" % round(clock, 2))

        # Plot static/dynamic polygons
        if (polys != None):
            # Plot each polygon -----------------------------------------------
            for pID in polys:
                pX = []
                pY = []
                for p in polys[pID][polyFieldName]:
                    pt = None
                    if ('direction' in polys[pID] and 'speed' in polys[pID]):
                        if (clock < polys[pID]['timeRange'][0]):
                            pt = p
                        elif (clock < polys[pID]['timeRange'][1]):
                            pt = vrpSolver.ptInDistXY(p, polys[pID]['direction'], polys[pID]['speed'] * clock)
                        else:
                            pt = vrpSolver.ptInDistXY(p, polys[pID]['direction'], polys[pID]['speed'] * (polys[pID]['timeRange'][1] - polys[pID]['timeRange'][0]))
                    else:
                        pt = p
                    if (not xyReverseFlag):
                        pX.append(pt[0])
                        pY.append(pt[1])
                    else:
                        pX.append(pt[1])
                        pY.append(pt[0])

                if ('fillColor' not in polyStyle[pID]):
                    ax.plot(pX, pY, color=polyStyle[pID]['edgeColor'], linewidth=polyStyle[pID]['edgeWidth'])
                elif ('fillStyle' not in polyStyle[pID]):
                    ax.fill(pX, pY, facecolor=polyStyle[pID]['fillColor'], edgecolor=polyStyle[pID]['edgeColor'], 
                        linewidth=polyStyle[pID]['edgeWidth'], alpha=polyStyle[pID]['opacity'])
                else:
                    ax.fill(pX, pY, facecolor=polyStyle[pID]['fillColor'], edgecolor=polyStyle[pID]['edgeColor'], 
                        hatch=polyStyle[pID]['fillStyle'], linewidth=polyStyle[pID]['edgeWidth'], alpha=polyStyle[pID]['opacity'])

        # Plot nodes
        if (nodes != None):
            # Plot neighborhood of each node ----------------------------------
            for nID in nodes:
                plotNodeFlag = False
                if ('timeWindow' in nodes[nID]):
                    if (nodes[nID]['timeWindow'][0] <= clock <= nodes[nID]['timeWindow'][1]):
                        plotNodeFlag = True
                else:
                    plotNodeFlag = True
                if (plotNodeFlag and 'neighbor' in nodes[nID] and nodeNeighborColor != None):
                    curNeighbor = [i for i in nodes[nID]['neighbor']]
                    neiX = []
                    neiY = []
                    for i in range(len(curNeighbor)):
                        neiPt = None
                        if ('direction' in nodes[nID] and 'speed' in nodes[nID]):
                            if (clock < nodes[nID]['timeRange'][0]):
                                neiPt = nodes[nID]['neighbor'][i]
                            elif (clock < nodes[nID]['timeRange'][1]):
                                neiPt = vrpSolver.ptInDistXY(nodes[nID]['neighbor'][i], nodes[nID]['direction'], nodes[nID]['speed'] * clock)
                            else:
                                neiPt = vrpSolver.ptInDistXY(nodes[nID]['neighbor'][i], nodes[nID]['direction'], nodes[nID]['speed'] * (nodes[nID]['timeRange'][1] - nodes[nID]['timeRange'][0]))
                        else:
                            neiPt = nodes[nID]['neighbor'][i]
                        if (not xyReverseFlag):
                            neiX.append(neiPt[0])
                            neiY.append(neiPt[1])
                        else:
                            neiX.append(neiPt[1])
                            neiY.append(neiPt[0])
                    ax.plot(neiX, neiY, color = 'black', linewidth = 1)
                    if ('nodeNeighborColor' in nodeStyle[nID]):
                        ax.fill(neiX, neiY, facecolor = nodeStyle[nID]['nodeNeighborColor'], edgecolor = 'black', hatch = '///', linewidth = 1, alpha = nodeNeighborOpacity)
                
            # Plot the location of nodes --------------------------------------
            for nID in nodes:
                plotNodeFlag = False
                if ('timeWindow' in nodes[nID]):
                    if (nodes[nID]['timeWindow'][0] <= clock <= nodes[nID]['timeWindow'][1]):
                        plotNodeFlag = True
                else:
                    plotNodeFlag = True
                if (plotNodeFlag):
                    curLoc = None
                    if ('direction' in nodes[nID] and 'speed' in nodes[nID]):
                        if (clock < nodes[nID]['timeRange'][0]):
                            curLoc = nodes[nID][locFieldName]
                        elif (clock < nodes[nID]['timeRange'][1]):
                            curLoc = vrpSolver.ptInDistXY(nodes[nID][locFieldName], nodes[nID]['direction'], nodes[nID]['speed'] * clock)
                        else:
                            curLoc = vrpSolver.ptInDistXY(nodes[nID][locFieldName], nodes[nID]['direction'], nodes[nID]['speed'] * (nodes[nID]['timeRange'][1] - nodes[nID]['timeRange'][0]))
                    else:
                        curLoc = nodes[nID][locFieldName]
                    x = None
                    y = None
                    if (not xyReverseFlag):
                        x = curLoc[0]
                        y = curLoc[1]
                    else:
                        x = curLoc[1]
                        y = curLoc[0]
                    ax.plot(x, y, color = nodeStyle[nID]['nodeColor'], marker = nodeStyle[nID]['nodeMarker'], markersize = nodeStyle[nID]['nodeMarkersize'])
                    if ('label' not in nodes[nID]):
                        lbl = nID
                    else:
                        lbl = nodes[nID]['label']
                    ax.annotate(lbl, (x, y))

        # Plot vehicle
        if (vehicles != None):
            for vID in vehicles:
                # Plot path ---------------------------------------------------
                if ('pathColor' in vehicleStyle[vID]):
                    pathX = []
                    pathY = []
                    if (not xyReverseFlag):
                        pathX = [i[0] for i in vehicles[vID]['seq']]
                        pathY = [i[1] for i in vehicles[vID]['seq']]
                    else:
                        pathX = [i[1] for i in vehicles[vID]['seq']]
                        pathY = [i[0] for i in vehicles[vID]['seq']]
                    ax.plot(pathX, pathY, color=vehicleStyle[vID]['pathColor'], linewidth = 1)

                # Plot trace --------------------------------------------------
                if ('traceColor' in vehicleStyle[vID]):
                    ts = timeRange[0]
                    te = clock
                    if (vehicleStyle[vID]['traceTime'] != None):
                        ts = max(te - vehicleStyle[vID]['traceTime'], timeRange[0])

                    if (ts < te):
                        trace = vrpSolver.traceInTimedSeq(
                            seq = vehicles[vID]['seq'],
                            timeStamp = vehicles[vID]['timeStamp'],
                            ts = ts,
                            te = te)

                        if (len(trace) > 0 and 'traceColor' in vehicleStyle[vID]):
                            traceX = []
                            traceY = []
                            if (not xyReverseFlag):
                                traceX = [i[0] for i in trace]
                                traceY = [i[1] for i in trace]
                            else:
                                traceX = [i[1] for i in trace]
                                traceY = [i[0] for i in trace]
                            ax.plot(traceX, traceY, color=vehicleStyle[vID]['traceColor'], linewidth = 1)

                # Plot vehicle ------------------------------------------------
                curLoc = vrpSolver.locInTimedSeq(
                    seq = vehicles[vID]['seq'],
                    timeStamp = vehicles[vID]['timeStamp'],
                    t = clock)
                ax.plot(curLoc[0], curLoc[1], color = vehicleStyle[vID]['vehColor'], marker = 'o', markersize = 4)
                if ('vehicle' not in vehicles[vID]):
                    lbl = vID
                else:
                    lbl = vehicles[vID]['vehicle']

                if (vehShowSpeedFlag):
                    curSpd = vrpSolver.speedInTimedSeq(seq = vehicles[vID]['seq'], timeStamp = vehicles[vID]['timeStamp'], t = clock)
                    lbl += " %s[m/s]" % (round(curSpd, 2))
                
                if (vehShowNoteFlag and 'note' in vehicles[vID]):
                    note = ""
                    if (clock <= vehicles[vID]['timeStamp'][0]):
                        note = vehicles[vID]['note'][0]
                    elif (clock >= vehicles[vID]['timeStamp'][-1]):
                        note = vehicles[vID]['note'][-1]
                    else:
                        for i in range(len(vehicles[vID]['timeStamp']) - 1):
                            if (vehicles[vID]['timeStamp'][i] <= clock < vehicles[vID]['timeStamp'][i + 1]):
                                note = vehicles[vID]['note'][i]
                    ax.annotate(lbl + "\n" + note, (curLoc[0], curLoc[1]))
                else:
                    ax.annotate(lbl, (curLoc[0], curLoc[1]))

    ani = FuncAnimation(
        fig, animate, 
        frames=int((timeRange[1] / (speed * 1.0) - timeRange[0] / (speed * 1.0)) * fps), 
        interval= int(1000 / fps), 
        repeat = repeatFlag)

    if (aniSavePath):
        ani.save("%s.gif" % aniSavePath, dpi=aniSaveDPI,
         writer=PillowWriter(fps=fps))

    return ani
