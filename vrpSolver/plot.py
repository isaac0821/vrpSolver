import matplotlib.pyplot as plt

from matplotlib import rcParams
# rcParams['font.family'] = 'SimSun'
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

from .common import *
from .color import *
from .msg import *
from .province import *
from .geometry import *

# History =====================================================================
# 20230518 - `plotNodes()` now will plot the neighborhood of nodes
# 20230624 - Rename functions `plotArcs()`, `plotLocSeq()`, `plotNodeSeq()`
# 20231022 - Refine `aniRouting()`
# =============================================================================

def plotLocs(
    locs: list[pt],
    locColor: str = 'Random',
    locMarker: str = 'o',
    locMarkersize: float = 1,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    """Plot locations on a figure

    Parameters
    ----------
    locs: list of pt, Required
        A list of locations to be plotted
    locColor: str, Optional, Default as 'Random'
        The color of locations to be plotted, 'Random' if the color is randomized
    


    """

    # Check for required fields ===============================================
    if (locs == None):
        raise MissingParameterError("ERROR: Missing required field `locs`.")

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        # Adjust bounding box
        (xMin, xMax, yMin, yMax) = boundingBox
        if (xMin == None or xMax == None or yMin == None or yMax == None):
            allX = []
            allY = []
            for i in locs:
                if (not xyReverseFlag):
                    allX.append(i[0])
                    allY.append(i[1])
                else:
                    allX.append(i[1])
                    allY.append(i[0])        
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

    # Draw locs ==============================================================
    for i in locs:
        # Define color --------------------------------------------------------
        color = None
        if (locColor == 'Random'):
            color = colorRandom()
        else:
            color = locColor

        # plot nodes ----------------------------------------------------------
        x = None
        y = None
        if (not xyReverseFlag):
            x = i[0]
            y = i[1]
        else:
            x = i[1]
            y = i[0]
        if (locMarkersize == None):
            ax.plot(x, y, color = color, marker = locMarker)
        else:
            ax.plot(x, y, color = color, marker = locMarker, markerSize = locMarkersize)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None and isinstance(fig, plt.Figure)):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax
    
def plotNodes(
    nodes: dict, 
    locFieldName = 'loc',
    nodeColor: str = 'Random',
    nodeMarker: str = 'o',
    nodeMarkersize: float = 1,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    """
    Draw nodes

    Parameters
    ----------
    nodes: dictionary, required
        A `nodes` dictionary. See :ref:`nodes` for reference.
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
    fig, ax
        matplotlib.pyplot object
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
        if ('markerSize' in nodes[n]):
            nodeMarkersize = nodes[n]['markerSize']

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
            ax.plot(x, y, color = color, marker = nodeMarker, markerSize = nodeMarkersize)
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
    arcStartLocFieldName = 'startLoc',
    arcEndLocFieldName = 'endLoc',
    arcColor: str = 'Random',
    arcWidth: float = 1.0,
    arcLabel: str = None,
    arcStyle: str = 'solid',
    arcDashes: tuple = (5, 2),
    arrowFlag: bool = True,
    arrowPosition: float = 0.5,
    arrowHeadWidth: float = 2.0,
    arrowHeadLength: float = 3.0,
    startColor: str = 'black',
    endColor: str = 'black',
    bothEndSize: int|float = 2.0,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):
    
    """
    Draw arcs

    Parameters
    ----------

    arcs: dict, required
        An `arcs` dictionary, each arc is defined by two points. See :ref:`arcs` for reference.
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
    fig, ax
        matplotlib.pyplot object
    """

    # Check for required fields ===============================================
    if (arcs == None):
        raise MissingParameterError("ERROR: Missing required field `arcs`.")
    for i in arcs:
        if (arcFieldName not in arcs[i] and arcStartLocFieldName not in arcs[i] and arcEndLocFieldName not in arcs[i]):
            raise MissingParameterError("ERROR: Cannot find arc field in given `arcs` - %s" % i)

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        # Adjust bounding box
        (xMin, xMax, yMin, yMax) = boundingBox
        if (xMin == None or xMax == None or yMin == None or yMax == None):
            allX = []
            allY = []
            if (arcFieldName in arcs):
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
            else:
                for i in arcs:
                    if (not xyReverseFlag):
                        allX.append(arcs[i][arcStartLocFieldName][0])
                        allY.append(arcs[i][arcStartLocFieldName][1])
                        allX.append(arcs[i][arcEndLocFieldName][0])
                        allY.append(arcs[i][arcEndLocFieldName][1])
                    else:
                        allX.append(arcs[i][arcStartLocFieldName][1])
                        allY.append(arcs[i][arcStartLocFieldName][0])
                        allX.append(arcs[i][arcEndLocFieldName][1])
                        allY.append(arcs[i][arcEndLocFieldName][0])
                        
                                    
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
        if (arcFieldName in arcs[i]):        
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
        elif (arcStartLocFieldName in arcs[i] and arcEndLocFieldName in arcs[i]):
            if (not xyReverseFlag):
                x1 = arcs[i][arcStartLocFieldName][0]
                y1 = arcs[i][arcStartLocFieldName][1]
                x2 = arcs[i][arcEndLocFieldName][0]
                y2 = arcs[i][arcEndLocFieldName][1]
            else:
                x1 = arcs[i][arcStartLocFieldName][1]
                y1 = arcs[i][arcStartLocFieldName][0]
                x2 = arcs[i][arcEndLocFieldName][1]
                y2 = arcs[i][arcEndLocFieldName][0]
        else:
            raise MissingParameterError("ERROR: Cannot find corresponded `arcFieldName` or (`arcStartLocFieldName` and `arcEndLocFieldName`) in given `arcs` structure.")
        dx = x2 - x1
        dy = y2 - y1
        if (arcColor == 'Random'):
            rndColor = colorRandom()
            ax.plot([x1, x2], [y1, y2], color = rndColor, linewidth=arcWidth, linestyle = arcStyle, dashes = arcDashes)
            if (arrowFlag):
                deg = headingXY([x1, y1], [x2, y2])
                ptC = [x1 + (x2 - x1) * arrowPosition, y1 + (y2 - y1) * arrowPosition]
                ptM = ptInDistXY(ptC, direction = deg + 180, dist = arrowHeadLength / 2)
                ptH = ptInDistXY(ptC, direction = deg, dist = arrowHeadLength / 2)
                pt1 = ptInDistXY(ptM, direction = deg + 90, dist = arrowHeadWidth / 2)
                pt2 = ptInDistXY(ptM, direction = deg - 90, dist = arrowHeadWidth / 2)
                ax.fill([ptH[0], pt1[0], pt2[0]], [ptH[1], pt1[1], pt2[1]], facecolor=rndColor, edgecolor=rndColor, linewidth=0)
                # ax.arrow(x=x1, y=y1, dx=dx * arrowPosition, dy=dy * arrowPosition, linewidth=arcWidth, head_width=arrowHeadWidth, head_length=arrowHeadLength, color=rndColor)
        else:
            ax.plot([x1, x2], [y1, y2], color = arcColor, linewidth=arcWidth, linestyle = arcStyle, dashes = arcDashes)
            if (arrowFlag):
                deg = headingXY([x1, y1], [x2, y2])
                ptC = [x1 + (x2 - x1) * arrowPosition, y1 + (y2 - y1) * arrowPosition]
                ptM = ptInDistXY(ptC, direction = deg + 180, dist = arrowHeadLength / 2)
                ptH = ptInDistXY(ptC, direction = deg, dist = arrowHeadLength / 2)
                pt1 = ptInDistXY(ptM, direction = deg + 90, dist = arrowHeadWidth / 2)
                pt2 = ptInDistXY(ptM, direction = deg - 90, dist = arrowHeadWidth / 2)
                ax.fill([ptH[0], pt1[0], pt2[0]], [ptH[1], pt1[1], pt2[1]], facecolor=arcColor, edgecolor=arcColor, linewidth=0)

        ax.plot(x1, y1, color = startColor, marker = 'o', markerSize = bothEndSize)
        ax.plot(x2, y2, color = endColor, marker = 'o', markerSize = bothEndSize)
        if (arcLabel == None and 'label' not in arcs[i]):
            lbl = i
        elif (arcLabel == None):
            lbl = arcs[i]['label']
        else:
            lbl = arcLabel
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
    lineStyle: str = 'solid',
    lineDashes: tuple = (5, 2),
    nodeColor: str = 'black',
    nodeMarkersize: float = 1,
    arrowFlag: bool = True,
    arrowPosition: float = 0.5,
    arrowHeadWidth: float = 0.1,
    arrowHeadLength: float = 0.2,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
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

    arcs = {}
    for i in range(len(locSeq) - 1):
        arcs[i] = {'arc': [locSeq[i], locSeq[i + 1]]}

    # Color ===================================================================
    if (lineColor == 'Random'):
        lineColor = colorRandom()

    fig, ax = plotArcs(
        fig = fig,
        ax = ax,
        arcs = arcs,
        arcFieldName = 'arc',
        arcColor = lineColor,
        arcWidth = lineWidth,
        arcStyle = lineStyle,
        arcDashes = lineDashes if lineStyle == 'dashed' else (None, None),
        arcLabel = "",
        arrowFlag = arrowFlag,
        arrowPosition = arrowPosition,
        arrowHeadWidth = arrowHeadWidth,
        arrowHeadLength = arrowHeadLength,
        startColor = nodeColor,
        endColor = nodeColor ,
        bothEndSize = nodeMarkersize,
        xyReverseFlag = xyReverseFlag,
        figSize = figSize,
        boundingBox = boundingBox,
        showAxis = showAxis,
        saveFigPath = saveFigPath,
        showFig = showFig)
    return fig, ax

def plotNodeSeq(
    nodes: dict, 
    nodeSeq: list[int|str],
    locFieldName: str = 'loc',
    lineColor: str = 'Random',
    lineWidth: float = 1,
    lineStyle: str = 'solid',
    lineDashes: tuple = (5, 2),
    nodeColor: str = 'black',
    nodeMarkersize: float = 1,
    arrowFlag: bool = True,
    arrowPosition: float = 0.5,
    arrowHeadWidth: float = 0.1,
    arrowHeadLength: float = 0.2,
    xyReverseFlag: bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    """Given a `nodes` dictionary and a sequence of node IDs, plot a route that visits each node by IDs.

    Parameters
    ----------
    nodes: dict, required
        A `node` dictionary. See :ref:`nodes` for reference.
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
    # Call plotArcs ===========================================================
    for n in nodeSeq:
        if (n not in nodes):
            raise UnsupportedInputError("ERROR: Cannot find %s in nodes" % n)

    arcs = {}
    for i in range(len(nodeSeq) - 1):
        arcs[i] = {'arc': [nodes[nodeSeq[i]][locFieldName], nodes[nodeSeq[i + 1]][locFieldName]]}

    # Color ===================================================================
    if (lineColor == 'Random'):
        lineColor = colorRandom()

    fig, ax = plotArcs(
        fig = fig,
        ax = ax,
        arcs = arcs,
        arcFieldName = 'arc',
        arcColor = lineColor,
        arcWidth = lineWidth,
        arcStyle = lineStyle,
        arcDashes = lineDashes if lineStyle == 'dashed' else (None, None),
        arcLabel = "",
        arrowFlag = arrowFlag,
        arrowPosition = arrowPosition,
        arrowHeadWidth = arrowHeadWidth,
        arrowHeadLength = arrowHeadLength,
        startColor = nodeColor,
        endColor = nodeColor ,
        bothEndSize = nodeMarkersize,
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
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
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
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
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

def plotPolygons(
    polygons: dict,
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
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    # Sanity check ============================================================
    if (type(polygons) != dict):
        raise UnsupportedInputError("ERROR: `polygons` is a dictionary, to plot an individual polygon, please use plotPoly() instead.")

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        for p in polygons:
            for i in polygons[p][polyFieldName]:
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
    for p in polygons:
        fig, ax = plotPoly(
            fig = fig,
            ax = ax,
            poly = polygons[p][polyFieldName],
            edgeWidth = edgeWidth,
            edgeColor = edgeColor,
            fillColor = fillColor,
            fillStyle = fillStyle,
            opacity = opacity,
            xyReverseFlag = xyReverseFlag,
            figSize = figSize)

    # Next, plot anchors ======================================================
    if (showAnchorFlag):
        for p in polygons:
            if ('label' not in polygons[p]):
                lbl = p
            else:
                lbl = polygons[p]['label']
            ha = 'center'
            va = 'center'
            ct = polygons[p]['anchor']
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
    fig, ax = plotPolygons(
        fig = fig,
        ax = ax,
        polygons = prvPoly,
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

def plotRoads(
    roads: dict,
    roadBoundary: poly|polys = None,
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
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ): 
    # FIXME: In future, we might want to distinguish roads by max speed or show the names of roads
    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        if (roadBoundary != None):
            try:    
                for poly in roadBoundary:
                    for pt in poly:
                        allX.append(float(pt[1]))
                        allY.append(float(pt[0]))
            except:
                for pt in roadBoundary:
                    allX.append(float(pt[1]))
                    allY.append(float(pt[0]))
        else:
            for road in roads:
                for pt in roads[road]['shape']:
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
    for road in roads:
        if (roadShowFlags == 'All' or roadShowFlags == True 
            or (type(roadShowFlags) == list 
                and roads[road]['class'] in roadShowFlags)):
            x = []
            y = []
            for pt in roads[road]['shape']:
                x.append(pt[1])
                y.append(pt[0])
            color = None
            if (roadColors == 'Random'):
                color = colorRandom()
            elif (type(roadColors) == str):
                color = roadColors
            elif (type(roadColors) == dict):
                if (roads[road]['class'] in roadColors):
                    color = roadColors[roads[road]['class']]
                else:
                    color = 'gray'
            rw = 1
            if (type(roadWidth) == dict):
                if (roads[road]['class'] in roadWidth):
                    rw = roadWidth[roads[road]['class']]
            else:
                rw = roadWidth
            ax.plot(x, y, color = color, linewidth = rw)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotBuildings(
    buildings: dict,
    buildingBoundary: poly|polys = None,
    buildingColors: dict[str,str]|str = {
            'building': 'yellow',
            'commercial': 'yellow',
            'residential': 'green',
            'house': 'green',
            'static_caravan': 'green',
            'industrial': 'orange',
            'manufacture': 'orange'
        },
    buildingShowFlags: list[str, bool]|str|bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ): 
    # FIXME: In future, we might want to distinguish roads by max speed or show the names of roads
    # If no based matplotlib figure, define boundary ==========================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        allX = []
        allY = []
        if (buildingBoundary != None):
            try:    
                for poly in buildingBoundary:
                    for pt in poly:
                        allX.append(float(pt[1]))
                        allY.append(float(pt[0]))
            except:
                for pt in buildingBoundary:
                    allX.append(float(pt[1]))
                    allY.append(float(pt[0]))

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

    # Plot buildings ==========================================================
    for building in buildings:
        if (buildingShowFlags == 'All' or buildingShowFlags == True
            or (type(buildingShowFlags) == list
                and buildings[building]['type'] in buildingShowFlags)):
            x = []
            y = []
            color = None
            for pt in buildings[building]['shape']:
                x.append(pt[1])
                y.append(pt[0])
            if (buildingColors == 'Random'):
                color = colorRandom()
            elif (type(buildingColors) == str):
                color = buildingColors
            elif (type(buildingColors) == dict):
                if (buildings[building]['type'] in buildingColors):
                    color = buildingColors[buildings[building]['type']]
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
    timedSeqFieldName: str = 'timedSeq',
    timeWindowFieldName: str = 'timeWindow',
    nodeColor: str = 'black',
    nodeMarker: str = 'o',
    nodeMarkersize: float = 2,
    vehicles: dict|None = None,
    vehTimedSeqFieldName: str = 'timedSeq',
    vehLabelFieldName: str = 'label',
    vehColor: str = 'blue',
    vehMarker: str = '^',
    vehMarkersize: float = 5,
    vehPathColor: str|None = 'gray',
    vehPathWidth: float|int|None = 3,
    vehTraceColor: str|None = 'orange',
    vehTraceWidth: float|int|None = 3,
    vehTraceShadowTime: float|None = None,
    vehSpdShowLabelFlag: bool = True,
    vehSpdShowArrowFlag: bool = True,
    vehSpdArrowLength: float = 5,
    vehShowNoteFlag: bool = True,
    polygons: dict|None = None,
    polyFieldName = 'poly',    
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
            ...     'markerSize': 2,  # Optional, default as None
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
            ...     'traceShadowTime': traceShadowTime
            ... }
    polygons: dictionary, optional, default None
        A dictionary indicating polygons that are dynamic/static in the animation
            >>> polygons[pID] = {
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
    if (nodes == None and vehicles == None and polygons == None):
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

    fig.set_figwidth(width) 
    fig.set_figheight(height)

    # Styling =================================================================
    nodeStyle = {}
    if (nodes != None):
        for nID in nodes:
            nodeStyle[nID] = {}
            if (nodeColor == None or nodeColor == 'Random'):
                nodeStyle[nID]['nodeColor'] = colorRandom()
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

    polyStyle = {}
    if (polygons != None):
        for pID in polygons:
            polyStyle[pID] = {}
            if (polyEdgeColor == None or polyEdgeColor == 'Random'):
                polyStyle[pID]['edgeColor'] = colorRandom()
            elif (polyEdgeColor != None):
                polyStyle[pID]['edgeColor'] = polyEdgeColor
            elif ('edgeColor' in polygons[pID]):
                polyStyle[pID]['edgeColor'] = polygons[pID]['edgeColor']

            if (polyEdgeWidth != None):
                polyStyle[pID]['edgeWidth'] = polyEdgeWidth
            elif ('edgeWidth' in polygons[pID]):
                polyStyle[pID]['edgeWidth'] = polygons[pID]['edgeWidth']

            if (polyFillColor == 'Random'):
                polyStyle[pID]['fillColor'] = colorRandom()
            elif (polyFillColor != None):
                polyStyle[pID]['fillColor'] = polyFillColor
            elif ('fillColor' in polygons[pID]):
                polyStyle[pID]['fillColor'] = polygons[pID]['fillColor']

            if (polyFillStyle != None):
                polyStyle[pID]['fillStyle'] = polyFillStyle
            elif ('fillStyle' in polygons[pID]):
                polyStyle[pID]['fillStyle'] = polygons[pID]['fillStyle']

            if (polyOpacity != None):
                polyStyle[pID]['opacity'] = polyOpacity
            elif ('opacity' in polygons[pID]):
                polyStyle[pID]['opacity'] = polygons[pID]['opacity']

    vehicleStyle = {}
    if (vehicles != None):
        for vID in vehicles:
            vehicleStyle[vID] = {}

            if (vehColor == 'Random'):
                vehicleStyle[vID]['vehColor'] = colorRandom()
            elif (vehColor != None):
                vehicleStyle[vID]['vehColor'] = vehColor
            elif ('color' in vehicles[vID]):
                vehicleStyle[vID]['vehColor'] = vehicles[vID]['color']

            if (vehMarker != None):
                vehicleStyle[vID]['vehMarker'] = vehMarker
            elif ('marker' in vehicles[vID]):
                vehicleStyle[vID]['vehMarker'] = vehicles[vID]['marker']

            if (vehMarkersize != None):
                vehicleStyle[vID]['vehMarkersize'] = vehMarkersize
            elif ('markerSize' in vehicles[vID]):
                vehicleStyle[vID]['vehMarkersize'] = vehicles[vID]['markerSize']

            if (vehPathColor == 'Random'):
                vehicleStyle[vID]['pathColor'] = colorRandom()
            elif (vehPathColor != None):
                vehicleStyle[vID]['pathColor'] = vehPathColor
            elif ('pathColor' in vehicles[vID]):
                vehicleStyle[vID]['pathColor'] = vehicles[vID]['pathColor']

            if (vehPathWidth != None):
                vehicleStyle[vID]['pathWidth'] = vehPathWidth
            elif ('pathWidth' in polygons[pID]):
                vehicleStyle[vID]['pathWidth'] = vehicles[vID]['pathWidth']

            if (vehTraceColor == 'Random'):
                vehicleStyle[vID]['traceColor'] = colorRandom()
            elif (vehTraceColor != None):
                vehicleStyle[vID]['traceColor'] = vehTraceColor
            elif ('traceColor' in vehicles[vID]):
                vehicleStyle[vID]['traceColor'] = vehicles[vID]['traceColor']

            if (vehTraceWidth != None):
                vehicleStyle[vID]['traceWidth'] = vehTraceWidth
            elif ('traceWidth' in polygons[pID]):
                vehicleStyle[vID]['traceWidth'] = vehicles[vID]['traceWidth']

            if (vehTraceShadowTime != None):
                vehicleStyle[vID]['traceShadowTime'] = vehTraceShadowTime
            elif ('traceShadowTime' in vehicles[vID]):
                vehicleStyle[vID]['traceShadowTime'] = vehicles[vID]['traceShadowTime']
            else:
                vehicleStyle[vID]['traceShadowTime'] = None

    def animate(t):
        ax.clear()
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

        # Clock
        clock = timeRange[0] + t * speed / (fps * 1.0)
        ax.set_title("Clock: %s[s]" % round(clock, 2))

        # Plot static/dynamic polygons
        if (polygons != None):
            # Plot each polygon -----------------------------------------------
            for pID in polygons:
                # 判定此时刻是否需要绘制poly
                plotPolyFlag = False
                if (timeWindowFieldName not in polygons[pID] or polygons[pID][timeWindowFieldName][0] <= clock <= polygons[pID][timeWindowFieldName][1]):
                    plotPolyFlag = True

                # 每个Poly的坐标轮廓
                pX = []
                pY = []
                if (plotPolyFlag):
                    for p in polygons[pID][polyFieldName]:
                        pt = None
                        if ('direction' in polygons[pID] and 'speed' in polygons[pID]):
                            if (clock < polygons[pID][timeWindowFieldName][0]):
                                pt = p
                            elif (clock < polygons[pID][timeWindowFieldName][1]):
                                pt = ptInDistXY(p, polygons[pID]['direction'], polygons[pID]['speed'] * clock)
                            else:
                                pt = ptInDistXY(p, polygons[pID]['direction'], polygons[pID]['speed'] * (polygons[pID]['timeRange'][1] - polygons[pID]['timeRange'][0]))
                        else:
                            pt = p

                        if (not xyReverseFlag):
                            pX.append(pt[0])
                            pY.append(pt[1])
                        else:
                            pX.append(pt[1])
                            pY.append(pt[0])

                # Plot polygons with styling
                if (plotPolyFlag):
                    if ('fillColor' not in polyStyle[pID]):
                        ax.plot(pX, pY, 
                            color=polyStyle[pID]['edgeColor'], 
                            linewidth=polyStyle[pID]['edgeWidth'])
                    elif ('fillStyle' not in polyStyle[pID]):
                        ax.fill(pX, pY, 
                            facecolor=polyStyle[pID]['fillColor'], 
                            edgecolor=polyStyle[pID]['edgeColor'], 
                            linewidth=polyStyle[pID]['edgeWidth'], 
                            alpha=polyStyle[pID]['opacity'])
                    else:
                        ax.fill(pX, pY, 
                            facecolor=polyStyle[pID]['fillColor'], 
                            edgecolor=polyStyle[pID]['edgeColor'], 
                            hatch=polyStyle[pID]['fillStyle'], 
                            linewidth=polyStyle[pID]['edgeWidth'], 
                            alpha=polyStyle[pID]['opacity'])

        # Plot nodes
        if (nodes != None):
            # Plot the location of nodes --------------------------------------
            for nID in nodes:
                plotNodeFlag = False
                if (timeWindowFieldName not in nodes[nID] or nodes[nID][timeWindowFieldName][0] <= clock <= nodes[nID][timeWindowFieldName][1]):
                    plotNodeFlag = True

                x = None
                y = None       
                curLoc = None         
                if (plotNodeFlag):
                    if (timedSeqFieldName not in nodes[nID]):
                        curLoc = nodes[nID][locFieldName]
                    else:
                        curSnap = snapInTimedSeq(
                            timedSeq = nodes[nID][timedSeqFieldName],
                            t = clock)
                        curLoc = curSnap['loc']

                if (not xyReverseFlag):
                    x = curLoc[0]
                    y = curLoc[1]
                else:
                    x = curLoc[1]
                    y = curLoc[0]
                # Styling of each node
                if (plotNodeFlag):
                    ax.plot(x, y, 
                        color = nodeStyle[nID]['nodeColor'], 
                        marker = nodeStyle[nID]['nodeMarker'], 
                        markerSize = nodeStyle[nID]['nodeMarkersize'])
                    if ('label' not in nodes[nID]):
                        lbl = nID
                    else:
                        lbl = nodes[nID]['label']
                    ax.annotate(lbl, (x, y))

        # Plot vehicle
        if (vehicles != None):
            # Plot each vehicle -----------------------------------------------
            for vID in vehicles:
                # Plot path
                if ('pathColor' in vehicleStyle[vID]):
                    pathX = []
                    pathY = []
                    if (not xyReverseFlag):
                        pathX = [vehicles[vID][vehTimedSeqFieldName][i][0][0] for i in range(len(vehicles[vID][vehTimedSeqFieldName]))]
                        pathY = [vehicles[vID][vehTimedSeqFieldName][i][0][1] for i in range(len(vehicles[vID][vehTimedSeqFieldName]))]
                    else:
                        pathX = [vehicles[vID][vehTimedSeqFieldName][i][0][1] for i in range(len(vehicles[vID][vehTimedSeqFieldName]))]
                        pathY = [vehicles[vID][vehTimedSeqFieldName][i][0][0] for i in range(len(vehicles[vID][vehTimedSeqFieldName]))]
                    ax.plot(pathX, pathY, color=vehicleStyle[vID]['pathColor'], linewidth = 1)

                # Plot trace
                if ('traceColor' in vehicleStyle[vID]):
                    ts = timeRange[0]
                    te = clock
                    if (vehicleStyle[vID]['traceShadowTime'] != None):
                        ts = max(te - vehicleStyle[vID]['traceShadowTime'], timeRange[0])

                    if (ts < te):
                        trace = traceInTimedSeq(
                            timedSeq = vehicles[vID][vehTimedSeqFieldName],
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

                # Plot vehicle
                curSnap = snapInTimedSeq(
                    timedSeq = vehicles[vID][vehTimedSeqFieldName],
                    t = clock)
                curLoc = curSnap['loc']
                ax.plot(curLoc[0], curLoc[1], 
                    color = vehicleStyle[vID]['vehColor'], 
                    marker = vehicleStyle[vID]['vehMarker'], 
                    markerSize = vehicleStyle[vID]['vehMarkersize'])

                if (vehLabelFieldName not in vehicles[vID]):
                    lbl = str(vID)
                else:
                    lbl = vehicles[vID][vehLabelFieldName]

                if (vehSpdShowLabelFlag):
                    curSpd = curSnap['speed']
                    lbl += " %s[m/s]" % (round(curSpd, 2))
                
                if (vehShowNoteFlag and 'note' in vehicles[vID]):
                    note = ""
                    if (clock <= vehicles[vID]['timedSeq'][0][1]):
                        note = vehicles[vID]['note'][0]
                    elif (clock >= vehicles[vID]['timedSeq'][-1][1]):
                        note = vehicles[vID]['note'][-1]
                    else:
                        for i in range(len(vehicles[vID]['timedSeq']) - 1):
                            if (vehicles[vID]['timedSeq'][i][1] <= clock < vehicles[vID]['timedSeq'][i + 1][1]):
                                note = vehicles[vID]['note'][i]
                                break
                    ax.annotate(lbl + "\n" + note, (curLoc[0], curLoc[1]))
                else:
                    ax.annotate(lbl, (curLoc[0], curLoc[1]))

    ani = FuncAnimation(
        fig, animate, 
        frames=int((timeRange[1] / (speed * 1.0) - timeRange[0] / (speed * 1.0)) * fps), 
        interval= int(1000 / fps), 
        repeat = repeatFlag)

    if (aniSavePath):
        ani.save("%s.gif" % aniSavePath, dpi=aniSaveDPI, writer=PillowWriter(fps=fps))

    return ani