import math
import networkx as nx
import shapely
from shapely.geometry import mapping

import gurobipy as grb

from .common import *
from .geometry import *
from .travel import *
from .polyTour import *
from .obj2Obj import *

def solveCETSP(
    startLoc: pt,
    endLoc: pt,
    nodes: dict, # Index from 1
    neighbor: str = "CircleXY",
    timeLimit: int | None = None,
    **kwargs
    ) -> dict | None:

    # WARNING: This is a freeze version, do no edit, origin version is in geoVeRo

    # Sanity check ============================================================
    if (nodes == None or type(nodes) != dict):
        raise MissingParameterError(ERROR_MISSING_NODES)

    if (neighbor in ['CircleXY', 'Poly', 'CircleLatLon']):
        if (neighbor == 'CircleXY'):
            if ('radius' not in kwargs and 'radiusFieldName' not in kwargs):
                raise MissingParameterError("ERROR: Must provide an uniform radius as `radius` or the field name of radius as `radiusFieldName`.")
        elif (neighbor == 'Poly'):
            if ('polyFieldName' not in kwargs):
                warnings.warn("WARNING: `polyFieldName` is not provided, set to be default as `poly`.")
                kwargs['polyFieldName'] = 'poly'
        elif (neighbor == 'CircleLatLon'):
            if ('radiusMeter' not in kwargs):
                raise MissingParameterError("ERROR: Missing `radiusMeter`.")
    else:
        raise UnsupportedInputError("ERROR: Neighborhood type is not supported")

    # Model initialization ====================================================
    CETSP = grb.Model("CETSP")
    CETSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        CETSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    convergence = []
    repSeqHis = {}
    degenCuts = []

    cutCount = {
        'subtour': 0,
        'gbc': 0,
        'degen': 0,
    }

    startTime = datetime.datetime.now()

    # Define sets
    startID = 0
    endID = len(nodes) + 1
    nodeAll = [i for i in range(0, len(nodes) + 2)]

    # Create neighborhoods ====================================================
    startLocMercator = None
    endLocMercator = None
    if (neighbor == "CircleXY"):
        for i in nodes:
            nodes[i]['neighbor'] = circleByCenterXY(
                center = nodes[i]['loc'], 
                radius = kwargs['radius'] if 'radius' in kwargs else nodes[i][kwargs['radiusFieldName']],
                lod = 60)
    elif (neighbor == "CircleLatLon"):
        for i in nodes:
            radius = kwargs['radiusMeter'] if 'radiusMeter' in kwargs else nodes[i][kwargs['radiusFieldName']]
            polyLatLon = circleByCenterLatLon(
                center = nodes[i]['loc'],
                radius = radius)
            polyXYMercator = polyLatLon2XYMercator(polyLatLon)
            nodes[i]['polyXYMercator'] = [pt for pt in polyXYMercator]
        startLocMercator = ptLatLon2XYMercator(startLoc)
        endLocMercator = ptLatLon2XYMercator(endLoc)

    # Parameters ==============================================================
    # anchor starts from depotLoc, in between are a list of circles, ends with depotLoc
    allX = []
    allY = []
    if (neighbor == "CircleXY"):
        allX = [startLoc[0]]
        allY = [startLoc[1]]
        for i in nodes:
            r = None
            if ('radius' in kwargs):
                r = kwargs['radius']
            elif ('radiusFieldName' in kwargs):
                r = nodes[i][kwargs['radiusFieldName']]
            allX.append(nodes[i]['loc'][0] - r)
            allX.append(nodes[i]['loc'][0] + r)
            allY.append(nodes[i]['loc'][1] - r)
            allY.append(nodes[i]['loc'][1] + r)
        allX.append(endLoc[0])
        allY.append(endLoc[1])
    elif (neighbor == "Poly"):
        allX = [startLoc[0]]
        allY = [startLoc[1]]
        for i in nodes:
            for p in nodes[i][kwargs['polyFieldName']]:
                allX.append(p[0])
                allY.append(p[1])
        allX.append(endLoc[0])
        allY.append(endLoc[1])
    elif (neighbor == "CircleLatLon"):
        startLocMercator = ptLatLon2XYMercator(startLoc)
        endLocMercator = ptLatLon2XYMercator(endLoc)
        allX = [startLocMercator[0]]
        allY = [startLocMercator[1]]
        for i in nodes:
            for p in nodes[i]['polyXYMercator']:
                allX.append(p[0])
                allY.append(p[1])
        allX.append(endLocMercator[0])
        allY.append(endLocMercator[1])

    lbX = min(allX) - 1
    lbY = min(allY) - 1
    ubX = max(allX) + 1
    ubY = max(allY) + 1

    zBar = {}
    if (neighbor == "CircleXY"):
        tau, _ = matrixDist(
            nodes = nodes, 
            edges = 'Euclidean', 
            locFieldName = 'loc')
        tauStart, _, _, _ = vectorDist(
            loc = startLoc,
            nodes = nodes,
            edges = 'Euclidean',
            locFieldName = 'loc')
        tauEnd, _, _, _ = vectorDist(
            loc = endLoc,
            nodes = nodes,
            edges = 'Euclidean',
            locFieldName = 'loc')
        for i in nodes:
            for j in nodes:
                if (i != j):
                    sumR = (2 * kwargs['radius']) if 'radius' in kwargs else nodes[i][kwargs['radiusFieldName']] + nodes[j][kwargs['radiusFieldName']]
                    zBar[i, j] = max(tau[i, j] - sumR, 0)
        for i in nodes:
            startR = kwargs['radius'] if 'radius' in kwargs else nodes[i][kwargs['radiusFieldName']]
            zBar[startID, i] = max(tauStart[i] - startR, 0)
            zBar[i, startID] = max(tauStart[i] - startR, 0)
            
            endR = kwargs['radius'] if 'radius' in kwargs else nodes[i][kwargs['radiusFieldName']]
            zBar[i, endID] = max(tauEnd[i] - endR, 0)
            zBar[endID, i] = max(tauEnd[i] - endR, 0)
        zBar[endID, startID] = 0
        zBar[startID, endID] = 0

    elif (neighbor == "Poly"):
        for i in nodes:
            for j in nodes:
                if (i != j):
                    zBar[i, j] = distPoly2Poly(nodes[i][kwargs['polyFieldName']], nodes[j][kwargs['polyFieldName']])
        for i in nodes:
            zBar[startID, i] = distPt2Poly(startLoc, nodes[i][kwargs['polyFieldName']])
            zBar[i, startID] = distPt2Poly(startLoc, nodes[i][kwargs['polyFieldName']])
            zBar[endID, i] = distPt2Poly(endLoc, nodes[i][kwargs['polyFieldName']])
            zBar[i, endID] = distPt2Poly(endLoc, nodes[i][kwargs['polyFieldName']])
        zBar[endID, startID] = 0
        zBar[startID, endID] = 0

    elif (neighbor == "CircleLatLon"):
        for i in nodes:
            for j in nodes:
                if (i != j):
                    zBar[i, j] = distPoly2Poly(nodes[i]['polyXYMercator'], nodes[j]['polyXYMercator'])
        for i in nodes:
            zBar[startID, i] = distPt2Poly(startLocMercator, nodes[i]['polyXYMercator'])
            zBar[i, startID] = distPt2Poly(startLocMercator, nodes[i]['polyXYMercator'])
            zBar[endID, i] = distPt2Poly(endLocMercator, nodes[i]['polyXYMercator'])
            zBar[i, endID] = distPt2Poly(endLocMercator, nodes[i]['polyXYMercator'])
        zBar[endID, startID] = 0
        zBar[startID, endID] = 0

    # Decision variables ======================================================
    # e[i,j] == 1 if disk i is visited instantly prior to j
    e = {}
    for i in nodeAll:
        for j in nodeAll:
            if (i != j):
                e[i, j] = CETSP.addVar(
                    vtype = grb.GRB.BINARY, 
                    obj = zBar[i, j],
                    name = 'e_%s_%s' % (i, j))

    theta = CETSP.addVar(
        vtype = grb.GRB.CONTINUOUS, 
        obj = 1,
        name = 'theta')

    # Objective function ======================================================
    CETSP.modelSense = grb.GRB.MINIMIZE
    CETSP.Params.lazyConstraints = 1
    CETSP.update()

    # TSP constraints =========================================================
    for i in nodeAll:
        CETSP.addConstr(grb.quicksum(e[i, j] for j in nodeAll if i != j) == 1)
    for i in nodeAll:
        CETSP.addConstr(grb.quicksum(e[j, i] for j in nodeAll if i != j) == 1)
    CETSP.addConstr(e[startID, endID] == 0)
    CETSP.addConstr(e[endID, startID] == 1)

    CETSP._e = e
    def GBDCutInfo(coeff, dvSeq, note) -> str:
        cutInfo = "Add %s: - theta >= " % note
        for i in range(len(dvSeq)):
            cutInfo += "%s * e[%s, %s] + " % (
                round((coeff[i] - zBar[dvSeq[i]]), 3), dvSeq[i][0], dvSeq[i][1]) 
        cutInfo = cutInfo[:-3] 
        return cutInfo

    def callbackCuts(model, where):
        if (where == grb.GRB.Callback.MIPSOL):
            TSPFeasibleFlag = True

            # Get visiting sequence -------------------------------------------
            e_sol = model.cbGetSolution(model._e)
            G = nx.Graph()
            for (i, j) in e.keys():
                if (e_sol[i, j] > 0.9):
                    G.add_edge(i, j)

            # Subtour elimination ---------------------------------------------
            components = [list(c) for c in nx.connected_components(G)]
            for component in components:
                if (len(component) < len(nodes) + 2):
                    TSPFeasibleFlag = False
                    model.cbLazy(grb.quicksum(e[i, j] for i in component for j in component if i != j) <= len(component) - 1)
                    cutCount['subtour'] += 1

            # Benders cut -----------------------------------------------------      
            if (TSPFeasibleFlag):
                # 还原访问顺序
                repArc = []
                for (i, j) in e.keys():
                    if (e_sol[i, j] > 0.9):
                        repArc.append([i, j])
                repSeq = []
                currentNode = startID
                repSeq.append(currentNode)
                while (len(repArc) > 0):
                    for i in range(len(repArc)):
                        if (repArc[i][0] == currentNode):
                            currentNode = repArc[i][1]
                            repSeq.append(currentNode)
                            repArc.pop(i)
                            break
                repSeq = repSeq[:-1]

                writeLog("\n")
                writeLog(list2String(repSeq))
                writeLog(hyphenStr())

                # Check if current solution is degenerated/duplicated
                BendersFlag = True
                if (tuple(repSeq) in repSeqHis):
                    writeLog("Duplicated repSeq: " + list2String(repSeq)) 
                    BendersFlag = False

                if (BendersFlag):
                    p2p = None
                    if (neighbor == "CircleXY"):
                        circles = []
                        for i in repSeq:
                            if (i != startID and i != endID):
                                circles.append({
                                    'center': nodes[i]['loc'],
                                    'radius': kwargs['radius'] if 'radius' in kwargs else nodes[i][kwargs['radiusFieldName']]
                                })
                        p2p = circle2CirclePath(startPt = startLoc, endPt = endLoc, circles = circles)
                    elif (neighbor == "Poly"):
                        polys = []
                        for i in repSeq:
                            if (i != startID and i != endID):
                                polys.append(nodes[i][kwargs['polyFieldName']])
                        p2p = poly2PolyPath(startPt = startLoc, endPt = endLoc, polys = polys)
                    elif (neighbor == "CircleLatLon"):
                        polys = []
                        for i in repSeq:
                            if (i != startID and i != endID):
                                polys.append(nodes[i]['polyXYMercator'])
                        p2p = poly2PolyPath(startPt = startLocMercator, endPt = endLocMercator, polys = polys)
                    
                    minDist = float('inf')
                    for seq in repSeqHis:
                        if (repSeqHis[seq] < minDist):
                            minDist = repSeqHis[seq]

                    # Cut preparing
                    residual = p2p['dist']
                    for i in range(len(repSeq) - 1):
                        residual -= zBar[repSeq[i], repSeq[i + 1]]
                    writeLog("Residual: %s" % residual)

                    polyFieldName = None
                    if (neighbor == "CircleXY"):
                        polyFieldName = 'neighbor'
                    elif (neighbor == "Poly"):
                        polyFieldName = kwargs['polyFieldName']
                    elif (neighbor == "CircleLatLon"):
                        polyFieldName = 'polyXYMercator'
                    mileage = polyPath2Mileage(repSeq, p2p['path'], nodes, polyFieldName = polyFieldName)
                    repSeqHis[tuple([i for i in repSeq])] = p2p['dist']

                    # Add default cut ---------------------------------------------
                    coeff = []
                    dvSeq = []
                    for i in range(len(repSeq) - 1):
                        coeff.append(distEuclideanXY(p2p['path'][i], p2p['path'][i + 1])['dist'])
                        dvSeq.append((repSeq[i], repSeq[i + 1]))
                    writeLog(GBDCutInfo(coeff, dvSeq, "default cut"))
                    model.cbLazy(theta >= grb.quicksum(
                        e[dvSeq[i]] * (coeff[i] - zBar[dvSeq[i]]) for i in range(len(dvSeq))))
                    cutCount['gbc'] += 1

                    objBound = model.cbGet(grb.GRB.Callback.MIPSOL_OBJBND)
                    objIncum = model.cbGet(grb.GRB.Callback.MIPSOL_OBJBST)
                    timePassed = round((datetime.datetime.now() - startTime).total_seconds(), 2)
                    writeLog("Time Pass: " + str(timePassed) + "[s]"
                        + "\nCut: subtour - %s, gbc - %s" % (cutCount['subtour'], cutCount['gbc'])
                        + "\nSo far Best Dist: " + str(minDist)
                        + "\nCurrSol Dist: " + str(p2p['dist'])
                        + "\nCurrObj: " + str(objIncum)
                        + "\nCurrBound: " + str(objBound))
                    convergence.append((p2p['dist'], objIncum, objBound, timePassed))             

    # TSP with no callback ====================================================
    CETSP.update()
    CETSP.optimize(callbackCuts)

    # Reconstruct solution ====================================================
    ofv = None
    seq = []
    activeX = []
    solType = None
    gap = None
    lb = None    
    ub = None
    runtime = None

    ofv = CETSP.getObjective().getValue()
    for i, j in e: 
        if (e[i, j].x > 0.5):
            activeX.append([i, j])
    currentNode = startID
    seq.append(currentNode)
    while (len(activeX) > 0):
        for i in range(len(activeX)):
            if (activeX[i][0] == currentNode):
                currentNode = activeX[i][1]
                seq.append(currentNode)
                activeX.pop(i)
                break
    seq = seq[:-1]

    p2p = None
    if (neighbor == "CircleXY"):
        circles = []
        for i in seq:
            if (i != startID and i != endID):
                circles.append({
                    'center': nodes[i]['loc'],
                    'radius': kwargs['radius'] if 'radius' in kwargs else nodes[i][kwargs['radiusFieldName']]
                })
        p2p = circle2CirclePath(startPt = startLoc, endPt = endLoc, circles = circles)
    elif (neighbor == "Poly"):
        polys = []
        for i in seq:
            if (i != startID and i != endID):
                polys.append(nodes[i][kwargs['polyFieldName']])
        p2p = poly2PolyPath(startPt = startLoc, endPt = endLoc, polys = polys)
    elif (neighbor == "CircleLatLon"):
        polys = []
        for i in seq:
            if (i != startID and i != endID):
                polys.append(nodes[i]['polyXYMercator'])
        p2pLatLon = poly2PolyPath(startPt = startLocMercator, endPt = endLocMercator, polys = polys)
        p2p = {
            'dist': p2pLatLon['dist'],
            'path': []
        }
        for i in range(len(p2pLatLon['path'])):
            p2p['path'].append(ptXY2LatLonMercator(p2pLatLon['path'][i]))

    if (CETSP.status == grb.GRB.status.OPTIMAL):
        solType = 'IP_Optimal'
        gap = 0
        lb = ofv
        ub = ofv
        runtime = CETSP.Runtime
    elif (CETSP.status == grb.GRB.status.TIME_LIMIT):
        solType = 'IP_TimeLimit'
        gap = CETSP.MIPGap
        lb = CETSP.ObjBoundC
        ub = CETSP.ObjVal
        runtime = CETSP.Runtime
    else:
        return None

    return {
        'ofv': ofv,
        'dist': p2p['dist'],
        'seq': seq,
        'path': p2p['path'],
        'gap': gap,
        'solType': solType,
        'lowerBound': lb,
        'upperBound': ub,
        'runtime': runtime,
        'convergence': convergence,
        'cutCount': cutCount
    }
