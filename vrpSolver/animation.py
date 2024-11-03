import matplotlib.pyplot as plt

from matplotlib import rcParams
# rcParams['font.family'] = 'SimSun'
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from .common import *
from .color import *
from .msg import *
from .province import *
from .geometry import *

def aniRouting(
    timeRange: tuple[int, int],
    # Nodes -------------------------------------------------------------------
    nodes: dict|None = None,
    locFieldName: str = 'loc',
    nodeTimedSeqFieldName: str = 'timedSeq',
    nodeTimeWindowFieldName: str = 'timeWindow',
    nodeColor: str = 'black',
    nodeMarker: str = 'o',
    nodeMarkerSize: float = 2,
    # Vehicles ----------------------------------------------------------------
    vehicles: dict|None = None,
    vehTimedSeqFieldName: str = 'timedSeq',
    vehLabelFieldName: str = 'label',
    vehColor: str = 'blue',
    vehMarker: str = '^',
    vehMarkerSize: float = 5,
    vehPathColor: str|None = 'gray',
    vehPathWidth: float|int|None = 3,
    vehTraceColor: str|None = 'orange',
    vehTraceWidth: float|int|None = 3,
    vehTraceShadowTime: float|None = None,
    vehSpdShowLabelFlag: bool = True,
    vehSpdShowArrowFlag: bool = True,
    vehSpdArrowLength: float = 5,
    vehShowNoteFlag: bool = True,
    # Polygons ----------------------------------------------------------------
    polygons: dict|None = None,
    polyAnchorFieldName: str = None,
    polyTimedSeqFieldName: str = None,
    polyTimeWindowFieldName: str = None,
    polyFieldName = 'poly',    
    polyEdgeColor: str = 'black',
    polyEdgeWidth: float = 1,
    polyFillColor: str|None = 'gray',
    polyFillStyle: str|None = '///',
    polyOpacity: float = 0.5,
    # Settings ----------------------------------------------------------------
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

            if (nodeMarkerSize != None):
                nodeStyle[nID]['nodeMarkerSize'] = nodeMarkerSize
            elif ('nodeMarkerSize' in nodes[nID]):
                nodeStyle[nID]['nodeMarkerSize'] = nodes[nID]['nodeMarkerSize']

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

            if (vehMarkerSize != None):
                vehicleStyle[vID]['vehMarkerSize'] = vehMarkerSize
            elif ('markerSize' in vehicles[vID]):
                vehicleStyle[vID]['vehMarkerSize'] = vehicles[vID]['markerSize']

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
                if (polyTimeWindowFieldName not in polygons[pID] or polygons[pID][polyTimeWindowFieldName][0] <= clock <= polygons[pID][polyTimeWindowFieldName][1]):
                    plotPolyFlag = True

                # 每个Poly的坐标轮廓
                pX = []
                pY = []
                if (plotPolyFlag):
                    for p in polygons[pID][polyFieldName]:
                        pt = None
                        if (polyAnchorFieldName in polygons[pID] and polyTimedSeqFieldName in polygons[pID]):
                            # 如果还没开始动，在原点不动
                            if (clock < polygons[pID][polyTimedSeqFieldName][0][1]):
                                pt = p
                            # 如果到达终点了，在终点不动
                            elif (clock > polygons[pID][polyTimedSeqFieldName][-1][1]):
                                dx = (polygons[pID][polyTimedSeqFieldName][-1][0][0] - polygons[pID][polyAnchorFieldName][0])
                                dy = (polygons[pID][polyTimedSeqFieldName][-1][0][1] - polygons[pID][polyAnchorFieldName][1])
                                pt = (p[0] + dx, p[1] + dy)
                            else: 
                                curSnap = snapInTimedSeq(
                                    timedSeq = polygons[pID][polyTimedSeqFieldName],
                                    t = clock)
                                dx = (curSnap['loc'][0] - polygons[pID][polyAnchorFieldName][0])
                                dy = (curSnap['loc'][1] - polygons[pID][polyAnchorFieldName][1])
                                pt = (p[0] + dx, p[1] + dy)
                        elif ('direction' in polygons[pID] and 'speed' in polygons[pID]):
                            if (clock < polygons[pID][polyTimeWindowFieldName][0]):
                                pt = p
                            elif (clock < polygons[pID][polyTimeWindowFieldName][1]):
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
                if (nodeTimeWindowFieldName not in nodes[nID] or nodes[nID][nodeTimeWindowFieldName][0] <= clock <= nodes[nID][nodeTimeWindowFieldName][1]):
                    plotNodeFlag = True

                x = None
                y = None       
                curLoc = None         
                if (plotNodeFlag):
                    if (nodeTimedSeqFieldName not in nodes[nID]):
                        curLoc = nodes[nID][locFieldName]
                    else:
                        curSnap = snapInTimedSeq(
                            timedSeq = nodes[nID][nodeTimedSeqFieldName],
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
                        markersize = nodeStyle[nID]['nodeMarkerSize'])
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
                    markersize = vehicleStyle[vID]['vehMarkerSize'])

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