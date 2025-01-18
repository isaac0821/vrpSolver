import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from .common import *
from .color import *
from .msg import *
from .province import *
from .geometry import *
from .travel import *

def findFigSize(boundingBox, width = None, height = None, latLonFlag = False):
    """
    Given a bounding box, a width(or height), returns the height(or width) of the figure

    Parameters
    ----------

    boundingBox: 4-tuple, required
        The bounding box of the figure
    width: float|None, optional, default as None
        The desired width of the figure
    height: float|None, optional, default as None
        The desired height of the figure

    Returns
    -------
    float, float
        The (width, height) proportional to bounding box

    """
    (xMin, xMax, yMin, yMax) = boundingBox

    w = None
    h = None

    if (not latLonFlag):
        if (width == None and height == None):
            if (xMax - xMin > yMax - yMin):
                w = 5
                h = 5 * ((yMax - yMin) / (xMax - xMin))
            else:
                w = 5 * ((xMax - xMin) / (yMax - yMin))
                h = 5
        elif (width != None and height == None):
            w = width
            h = width * ((yMax - yMin) / (xMax - xMin))
        elif (width == None and height != None):
            w = height * ((xMax - xMin) / (yMax - yMin))
            h = height
        else:
            w = width
            h = height
    else:
        heightDelta = distLatLon(
            (xMin, yMin + (yMax - yMin) / 2), 
            (xMax, yMin + (yMax - yMin) / 2))
        widthDelta = distLatLon(
            (xMin + (xMax - xMin) / 2, yMin),
            (xMin + (xMax - xMin) / 2, yMax))

        if (width == None and height == None):
            if (widthDelta > heightDelta):
                w = 5
                h = 5 * heightDelta / widthDelta
            else:
                w = 5 * widthDelta / heightDelta
                h = 5
        elif (width != None and height == None):
            w = width
            h = width * heightDelta / widthDelta
        elif (width == None and height != None):
            w = height * widthDelta / heightDelta
            h = height
        else:
            w = width
            h = height

    return w, h

def plotLocs(
    locs: list[pt],
    locColor: str = 'Random',
    locMarker: str = 'o',
    locMarkerSize: float = 1,
    latLonFlag: bool = False,
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
    locs: list of pt, required
        A list of locations to be plotted
    locColor: str, optional, default as 'Random'
        The color of locations to be plotted, 'Random' if the color is randomized.
    locMarker: str, optional, default as 'o'
        The shape of the marker
    locMarkerSize: str, optional, default as 1
        The size of the marker
    latLonFlag: bool, optional, default as False
        Reverse the x, y, (x, y) => (y, x). Used in plotting (lat, lon) coordinates.
    fig: matplotlib object, optional, default as None
        If fig and ax are provided, 


    """

    # Check for required fields ===============================================
    if (locs == None):
        raise MissingParameterError("ERROR: Missing required field `locs`.")

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig, ax = plt.subplots()
        boundingBox = findBoundingBox(
            boundingBox = boundingBox, 
            pts = locs,
            latLonFlag = latLonFlag)
        (xMin, xMax, yMin, yMax) = boundingBox
        (width, height) = findFigSize(boundingBox, figSize[0], figSize[1], latLonFlag)
        fig.set_figwidth(width)
        fig.set_figheight(height)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

    # Draw locs ==============================================================
    for i in locs:
        # Define color --------------------------------------------------------
        color = None
        if (locColor == 'Random'):
            color = rndColor()
        else:
            color = locColor

        # plot nodes ----------------------------------------------------------
        x = None
        y = None
        if (not latLonFlag):
            x = i[0]
            y = i[1]
        else:
            x = i[1]
            y = i[0]
        if (locMarkerSize == None):
            ax.plot(x, y, color = color, marker = locMarker)
        else:
            ax.plot(x, y, color = color, marker = locMarker, markersize = locMarkerSize)

    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None and isinstance(fig, plt.Figure)):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

    return fig, ax

def plotLocs3D(
    locs3D: list[pt3D],
    locColor: str = 'Random',
    locMarker: str = 'o',
    locMarkerSize: float = 1,
    latLonFlag: bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox3D = (None, None, None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    # Check for required fields ===============================================
    if (locs3D == None):
        raise MissingParameterError("ERROR: Missing required field `locs3D`.")

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        boundingBox3D = findBoundingBox3D(
            boundingBox3D = boundingBox3D, 
            pts3D = locs3D)
        (xMin, xMax, yMin, yMax, zMin, zMax) = boundingBox3D
        # (width, height) = findFigSize(boundingBox3D, figSize[0], figSize[1], latLonFlag)
        # fig.set_figwidth(width)
        # fig.set_figheight(height)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)
        ax.set_zlim(zMin, zMax)

    # Draw locs ==============================================================
    for i in locs3D:
        # Define color --------------------------------------------------------
        color = None
        if (locColor == 'Random'):
            color = rndColor()
        else:
            color = locColor

        # plot nodes ----------------------------------------------------------        
        x = None
        y = None
        z = i[2]
        if (not latLonFlag):
            x = i[0]
            y = i[1]
        else:
            x = i[1]
            y = i[0]

        if (locMarkerSize == None):
            ax.scatter3D(x, y, z, c = color, marker = locMarker)
        else:
            ax.scatter3D(x, y, z, c = color, marker = locMarker, s = locMarkerSize)

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
    nodeMarkerSize: float = 1,
    neighborFieldName = 'neighbor',
    neighborColor: str|None = 'gray',
    neighborOpacity: float = 0.5,
    neighborFillStyle: str = '///',
    latLonFlag: bool = False,
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
    locFieldName: str, optional, default as 'loc'
        The key value in `nodes` indicating the location of each node.
    nodeColor: str, optional, default 'Random'
        Alternative option. If 'color' is provided in `nodes`, this field will be ignored.
    nodeMarker: str, optional, default 'o'
        Alternative option for node marker. If 'marker' is provided in `nodes`, this field will be ignored.
    nodeMarker: str, optional, default 'o'
        Alternative option for node marker size. If 'markerSize' is provided in `nodes`, this field will be ignored.
    neighborColor: str, optional, default 'gray'
        If nodes have 'neighbor' label, will plot the neighbor area in this color
    neighborOpacity: float, optional, default 0.5
        The opacity of neighborhood.
    neighborFillStyle: str, optional, default '///'
        The fill style of the neighborhood.
    latLonFlag: bool, optional, default False
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
        boundingBox = findBoundingBox(
            boundingBox = boundingBox, 
            nodes = nodes, 
            locFieldName = locFieldName, 
            latLonFlag = latLonFlag)
        (xMin, xMax, yMin, yMax) = boundingBox
        (width, height) = findFigSize(boundingBox, figSize[0], figSize[1], latLonFlag)
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
            color = rndColor()
        else:
            color = nodeColor

        # Define marker and marker size ---------------------------------------
        if ('marker' in nodes[n]):
            nodeMarker = nodes[n]['marker']
        if ('markerSize' in nodes[n]):
            nodeMarkerSize = nodes[n]['markerSize']

        # plot nodes ----------------------------------------------------------
        x = None
        y = None
        if (not latLonFlag):
            x = nodes[n][locFieldName][0]
            y = nodes[n][locFieldName][1]
        else:
            x = nodes[n][locFieldName][1]
            y = nodes[n][locFieldName][0]
        if (nodeMarkerSize == None):
            ax.plot(x, y, color = color, marker = nodeMarker)
        else:
            ax.plot(x, y, color = color, marker = nodeMarker, markersize = nodeMarkerSize)
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

    # Draw node neighborhoods =================================================
    for n in nodes:
        # If node has neighbor, plot the neighbor color first -----------------
        if ('neiShape' in nodes[n] and nodes[n]['neiShape'] == 'Poly'):
            fig, ax = plotPoly(
                fig = fig,
                ax = ax,
                poly = nodes[n][neighborFieldName],
                edgeWidth = 1,
                edgeColor = 'black',
                fillColor = neighborColor,
                opacity = neighborOpacity,
                latLonFlag = latLonFlag,
                showAxis = showAxis,
                fillStyle = neighborFillStyle)
        if ('neiShape' in nodes[n] and nodes[n]['neiShape'] == 'Isochrone'):
            for poly in nodes[n][neighborFieldName]:
                fig, ax = plotPoly(
                    fig = fig,
                    ax = ax,
                    poly = poly,
                    edgeWidth = 1,
                    edgeColor = 'black',
                    fillColor = neighborColor,
                    opacity = neighborOpacity / len(nodes[n][neighborFieldName]),
                    latLonFlag = latLonFlag,
                    showAxis = showAxis,
                    fillStyle = neighborFillStyle)            
        if ('neiShape' in nodes[n] and 'radius' in nodes[n] and nodes[n]['neiShape'] == 'Circle'):
            neiPoly = circleByCenterXY(
                center = nodes[n][locFieldName],
                radius = nodes[n]['radius'])
            fig, ax = plotPoly(
                fig = fig,
                ax = ax,
                poly = neiPoly,
                edgeWidth = 1,
                edgeColor = 'black',
                fillColor = neighborColor,
                opacity = neighborOpacity,
                latLonFlag = latLonFlag,
                showAxis = showAxis,
                fillStyle = neighborFillStyle)

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
    neighborOpacity: float = 0.5,
    neighborEntWidth: float = 1.2,
    neighborEntColor: str = 'black',
    neighborBtwWidth: float = 1,
    neighborBtwColor: str = 'gray',
    neighborFillStyle: str = '///',
    latLonFlag: bool = False,
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
    latLonFlag: bool, optional, default False
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
        boundingBox = findBoundingBox(
            boundingBox = boundingBox, 
            arcs = arcs, 
            arcFieldName = arcFieldName,
            arcStartLocFieldName = arcStartLocFieldName,
            arcEndLocFieldName = arcEndLocFieldName, 
            latLonFlag = latLonFlag)
        (xMin, xMax, yMin, yMax) = boundingBox
        (width, height) = findFigSize(boundingBox, figSize[0], figSize[1], latLonFlag)
        if (isinstance(fig, plt.Figure)):
            fig.set_figwidth(width)
            fig.set_figheight(height)
            ax.set_xlim(xMin, xMax)
            ax.set_ylim(yMin, yMax)

    # Draw arcs ===============================================================
    for i in arcs:
        x1, y1, x2, y2 = 0, 0, 0, 0
        if (arcFieldName in arcs[i]):        
            if (not latLonFlag):
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
            if (not latLonFlag):
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
            rndRGB = rndColor()
            ax.plot([x1, x2], [y1, y2], color = rndRGB, linewidth=arcWidth, linestyle = arcStyle, dashes = arcDashes)
            if (arrowFlag):
                deg = headingXY([x1, y1], [x2, y2])
                ptC = [x1 + (x2 - x1) * arrowPosition, y1 + (y2 - y1) * arrowPosition]
                ptM = ptInDistXY(ptC, direction = deg + 180, dist = arrowHeadLength / 2)
                ptH = ptInDistXY(ptC, direction = deg, dist = arrowHeadLength / 2)
                pt1 = ptInDistXY(ptM, direction = deg + 90, dist = arrowHeadWidth / 2)
                pt2 = ptInDistXY(ptM, direction = deg - 90, dist = arrowHeadWidth / 2)
                ax.fill([ptH[0], pt1[0], pt2[0]], [ptH[1], pt1[1], pt2[1]], facecolor=rndRGB, edgecolor=rndRGB, linewidth=0)
                # ax.arrow(x=x1, y=y1, dx=dx * arrowPosition, dy=dy * arrowPosition, linewidth=arcWidth, head_width=arrowHeadWidth, head_length=arrowHeadLength, color=rndRGB)
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

        ax.plot(x1, y1, color = startColor, marker = 'o', markersize = bothEndSize)
        ax.plot(x2, y2, color = endColor, marker = 'o', markersize = bothEndSize)
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

    # Draw arcs neighborhoods =================================================
    for i in arcs:
        if ('neiA' in arcs[i]):
            fig, ax = plotPoly(
                fig = fig,
                ax = ax,
                poly = arcs[i]['neiA'],
                edgeWidth = neighborEntWidth,
                edgeColor = 'black',
                fillColor = neighborEntColor,
                opacity = neighborOpacity,
                latLonFlag = latLonFlag,
                showAxis = showAxis,
                fillStyle = neighborFillStyle)
        if ('neiB' in arcs[i]):
            fig, ax = plotPoly(
                fig = fig,
                ax = ax,
                poly = arcs[i]['neiB'],
                edgeWidth = neighborEntWidth,
                edgeColor = 'black',
                fillColor = neighborEntColor,
                opacity = neighborOpacity,
                latLonFlag = latLonFlag,
                showAxis = showAxis,
                fillStyle = neighborFillStyle)
        if ('neiBtw' in arcs[i]):
            fig, ax = plotPoly(
                fig = fig,
                ax = ax,
                poly = arcs[i]['neiBtw'],
                edgeWidth = neighborBtwWidth,
                edgeColor = 'black',
                fillColor = neighborBtwColor,
                opacity = neighborOpacity,
                latLonFlag = latLonFlag,
                showAxis = showAxis,
                fillStyle = neighborFillStyle)

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
    nodeMarkerSize: float = 1,
    arrowFlag: bool = True,
    arrowPosition: float = 0.5,
    arrowHeadWidth: float = 0.1,
    arrowHeadLength: float = 0.2,
    latLonFlag: bool = False,
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
    latLonFlag: bool, optional, default False
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
        if (not is2PtsSame(locSeq[i], locSeq[i + 1])):
            arcs[i] = {'arc': [locSeq[i], locSeq[i + 1]]}

    # Color ===================================================================
    if (lineColor == 'Random'):
        lineColor = rndColor()

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
        bothEndSize = nodeMarkerSize,
        latLonFlag = latLonFlag,
        figSize = figSize,
        boundingBox = boundingBox,
        showAxis = showAxis,
        saveFigPath = saveFigPath,
        showFig = showFig)
    return fig, ax

def plotLocSeq3D(
    locSeq3D: list[pt3D],
    lineColor: str = 'Random',
    lineWidth: float = 1.0,
    lineStyle: str = 'solid',
    lineDashes: tuple = (5, 2),
    nodeColor: str = 'black',
    nodeMarkerSize: float = 1,
    latLonFlag: bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox3D = (None, None, None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

    # Check for required fields ===============================================
    if (locSeq3D == None):
        raise MissingParameterError("ERROR: Missing required field `locSeq3D`.")

    # If no based matplotlib figure provided, define boundary =================
    if (fig == None or ax == None):
        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        boundingBox3D = findBoundingBox3D(
            boundingBox3D = boundingBox3D, 
            pts3D = locSeq3D)
        (xMin, xMax, yMin, yMax, zMin, zMax) = boundingBox3D
        # (width, height) = findFigSize(boundingBox3D, figSize[0], figSize[1], latLonFlag)
        # fig.set_figwidth(width)
        # fig.set_figheight(height)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)
        ax.set_zlim(zMin, zMax)

    x = []
    y = []
    z = []
    for loc in locSeq3D:
        if (not latLonFlag):
            x.append(loc[0])
            y.append(loc[1])
        else:
            x.append(loc[1])
            y.append(loc[0])
        z.append(loc[2])

    color = None
    if (lineColor == 'Random'):
        color = rndColor()
    else:
        color = lineColor

    if (lineStyle == 'dashed'):
        ax.plot(x, y, z, c = color, linewidth=lineWidth, linestyle = 'dashed', dashes = lineDashes)
    else:
        ax.plot(x, y, z, c = color, linewidth=lineWidth, linestyle = lineStyle)
            
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
    lineStyle: str = 'solid',
    lineDashes: tuple = (5, 2),
    nodeColor: str = 'black',
    nodeMarkerSize: float = 1,
    arrowFlag: bool = True,
    arrowPosition: float = 0.5,
    arrowHeadWidth: float = 0.1,
    arrowHeadLength: float = 0.2,
    latLonFlag: bool = False,
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
    latLonFlag: bool, optional, default False
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
        lineColor = rndColor()

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
        bothEndSize = nodeMarkerSize,
        latLonFlag = latLonFlag,
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
    latLonFlag: bool = False,
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
    latLonFlag: bool, optional, default False
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
        boundingBox = findBoundingBox(
            boundingBox = boundingBox, 
            poly = poly,
            latLonFlag = latLonFlag)

        (xMin, xMax, yMin, yMax) = boundingBox
        width = None
        height = None
        if (figSize == None or figSize[0] == None or figSize[1] == None):
            (width, height) = findFigSize(boundingBox, figSize[0], figSize[1], latLonFlag)
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
        if (not latLonFlag):
            x.append(pt[0])
            y.append(pt[1])
        else:
            x.append(pt[1])
            y.append(pt[0])
    if (not latLonFlag):
        x.append(poly[0][0])
        y.append(poly[0][1])
    else:
        x.append(poly[0][1])
        y.append(poly[0][0])

    if (abs(x[-1] - x[0]) > ERRTOL['distPt2Pt']):
        x.append(x[0])
        y.append(y[0])

    # Plot ====================================================================
    if (edgeColor == 'Random'):
        edgeColor = rndColor()
    if (fillColor == None):
        ax.plot(x, y, color = edgeColor, linewidth = edgeWidth)
    else:
        if (fillColor == 'Random'):
            ax.fill(x, y, facecolor=rndColor(), edgecolor=edgeColor, hatch=fillStyle, linewidth=edgeWidth, alpha=opacity)
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

def plotPolyhedron3D(
    polyhedron: poly, 
    edgeWidth: float = 0.5,
    edgeColor: str|None = 'Random',
    fillColor: str|None = None,
    fillStyle: str = "///",
    opacity: float = 0.5,
    latLonFlag: bool = False,
    fig = None,
    ax = None,
    figSize = (None, 5), 
    boundingBox = (None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True
    ):

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
    latLonFlag: bool = False,
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
        latLonFlag = latLonFlag,
        fig = fig,
        ax = ax,
        figSize = figSize,
        boundingBox = boundingBox,
        showAxis = showAxis,
        saveFigPath = saveFigPath,
        showFig = showFig,
    )
    return fig, ax

def plotCone3D(
    cone: dict,
    lod: int = 30,
    tanAlpha: float|None = None,    
    fig = None,
    ax = None,
    boundingBox3D = (None, None, None, None, None, None),
    showAxis: bool = True,
    saveFigPath: str|None = None,
    showFig: bool = True):

    if (fig == None or ax == None):
        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        boundingBox3D = findBoundingBox3D(
            boundingBox3D = boundingBox3D, 
            cone = cone)
        (xMin, xMax, yMin, yMax, zMin, zMax) = boundingBox3D
        # (width, height) = findFigSize(boundingBox3D, figSize[0], figSize[1], latLonFlag)
        # fig.set_figwidth(width)
        # fig.set_figheight(height)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)
        ax.set_zlim(zMin, zMax)

    x0 = cone['center'][0]
    y0 = cone['center'][1]
    z0 = cone['center'][2] if len(cone['center']) >= 3 else 0
    zMax = cone['maxHeight']
    if ('tanAlpha' in cone):
        tanAlpha = cone['tanAlpha']
    rMax = zMax * tanAlpha

    r = np.linspace(0, rMax, lod)
    theta = np.linspace(0, 2 * np.pi, lod)
    r, theta = np.meshgrid(r, theta)

    X = r * np.sin(theta) + x0
    Y = r * np.cos(theta) + y0

    def f(x, y):
        return np.sqrt((x - x0)**2 + (y - y0)**2) / tanAlpha + z0

    Z = f(X, Y)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap = cm.jet, alpha = 0.1,
                edgecolor='black');
    
    # Axis on and off =========================================================
    if (not showAxis):
        plt.axis('off')

    # Save figure =============================================================
    if (saveFigPath != None and isinstance(fig, plt.Figure)):
        fig.savefig(saveFigPath)
    if (not showFig):
        plt.close(fig)

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
    latLonFlag: bool = False,
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
        boundingBox = findBoundingBox(
            boundingBox = boundingBox, 
            polygons = polygons,
            anchorFieldName = anchorFieldName,
            polyFieldName = polyFieldName,
            latLonFlag = latLonFlag)
        (xMin, xMax, yMin, yMax) = boundingBox
        (width, height) = findFigSize(boundingBox, figSize[0], figSize[1], latLonFlag)
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
            latLonFlag = latLonFlag,
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
            if (latLonFlag):
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

def plotMapChina(
    regionList: list,
    readFromLocalFlag: bool = True,
    localDirectory: str = "D:/Zoo/gull/geoVeRoPy/data/map",
    saveToLocalFlag: bool = True,
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

    multiPoly = []

    # 一会儿要绘制的边界线们
    for region in regionList:
        multiPoly.extend(findPolysChina(
            province = region[0],
            city = region[1],
            district = region[2],
            readFromLocalFlag = readFromLocalFlag,
            localDirectory = localDirectory,
            saveToLocalFlag = saveToLocalFlag))

    # 计算范围大小
    x = []
    y = []
    for polys in multiPoly:
        for poly in polys:
            for pt in poly:
                x.append(pt[0])
                y.append(pt[1])
    boundingBox = (min(x), max(x), min(y), max(y))
    figSize = findFigSize(boundingBox, figSize[0], figSize[1], True)

    # 读取本地的
    for polys in multiPoly:
        for poly in polys:
            fig, ax = plotPoly(
                fig = fig,
                ax = ax,
                poly = poly,
                edgeWidth = edgeWidth,
                edgeColor = edgeColor,
                fillColor = fillColor,
                fillStyle = fillStyle,
                opacity = opacity,
                latLonFlag = True,
                boundingBox = boundingBox,
                figSize = figSize,
                saveFigPath = saveFigPath,
                showAxis = showAxis,
                showFig = showFig)

    return fig, ax

def plotMapUS(
    stateList: list[str]|str = [],
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
    if (type(stateList) == str):
        stateList = [stateList]
    for prv in stateList:
        if prv in prov:
            ct = ptPolyCenter(prov[prv]['shape'])
            prvPoly[prv] = {
                'anchor': ct,
                'label': prov[prv]['abbr'],
                'poly': prov[prv]['shape'],
            }
        else:
            for k in prov:
                if (prv == prov[k]['abbr']):
                    ct = ptPolyCenter(prov[k]['shape'])
                    prvPoly[prv] = {
                        'anchor': ct,
                        'label': prov[k]['abbr'],
                        'poly': prov[k]['shape'],
                    }

    fig, ax = plotPolygons(
        fig = fig,
        ax = ax,
        polygons = prvPoly,
        showAnchorFlag = False,
        edgeWidth = edgeWidth,
        edgeColor = edgeColor,
        fillColor = fillColor,
        fillStyle = fillStyle,
        opacity = opacity,
        latLonFlag = True,
        figSize = figSize,
        saveFigPath = saveFigPath,
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
                        allX.append(float(pt[0]))
                        allY.append(float(pt[1]))
            except:
                for pt in roadBoundary:
                    allX.append(float(pt[0]))
                    allY.append(float(pt[1]))
        else:
            for road in roads:
                for pt in roads[road]['shape']:
                    allX.append(pt[0])
                    allY.append(pt[1])

        (xMin, xMax, yMin, yMax) = boundingBox
        if (xMin == None):
            xMin = min(allX) - 0.05 * abs(max(allX) - min(allX))
        if (xMax == None):
            xMax = max(allX) + 0.05 * abs(max(allX) - min(allX))
        if (yMin == None):
            yMin = min(allY) - 0.05 * abs(max(allY) - min(allY))
        if (yMax == None):
            yMax = max(allY) + 0.05 * abs(max(allY) - min(allY))
        boundingBox = (xMin, xMax, yMin, yMax)
        figSize = findFigSize(boundingBox, figSize[0], figSize[1], True)
        (width, height) = figSize

        fig.set_figwidth(width)
        fig.set_figheight(height)
        ax.set_xlim(yMin, yMax)
        ax.set_ylim(xMin, xMax)

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
                color = rndColor()
            elif (type(roadColors) == str):
                color = roadColors
            elif (type(roadColors) == dict):
                if (roads[road]['class'] in roadColors):
                    color = roadColors[roads[road]['class']]
                else:
                    color = roadColors['others']
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
                color = rndColor()
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
                groupColor.append(rndColor())
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
                    color=phase[i]['backgroundColor'] if phase[i]['backgroundColor'] != 'Random' else rndColor(), 
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
                    rndRGB = rndColor()
                    ax.fill(x, y, color = rndRGB, linewidth = ganttLinewidth)
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
