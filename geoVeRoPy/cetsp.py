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
    nodes: dict, 
    neighbor: str = "CircleXY",
    algo: str = "GeneticAlgorithm",
    timeLimit: int | None = None,
    **kwargs):

    # WARNING: This is a freeze version, do no edit, improved version is in geoVeRoPyPri
    # Sanity check ============================================================
    if (nodes == None or type(nodes) != dict):
        raise MissingParameterError(ERROR_MISSING_NODES)

    if (neighbor in ['CircleXY', 'PolyXY', 'CircleLatLon', 'PolyLatLon']):
        if (neighbor == 'CircleXY'):
            if ('radius' not in kwargs and 'radiusFieldName' not in kwargs):
                raise MissingParameterError("ERROR: Must provide an uniform radius as `radius` or the field name of radius as `radiusFieldName`.")
        elif (neighbor == 'PolyXY'):
            if ('polyFieldName' not in kwargs):
                warnings.warn("WARNING: `polyFieldName` is not provided, set to be default as `poly`.")
                kwargs['polyFieldName'] = 'poly'
        elif (neighbor == 'CircleLatLon'):
            if ('radiusMeter' not in kwargs):
                raise MissingParameterError("ERROR: Missing `radiusMeter`.")
        elif (neighbor == 'PolyLatLon'):
            if ('polyFieldName' not in kwargs):
                warnings.warn("WARNING: `polyFieldName` is not provided, set to be default as `poly`.")
                kwargs['polyFieldName'] = 'poly'
    else:
        raise UnsupportedInputError("ERROR: Neighborhood type is not supported")


    return

def _solveCETSPGeralizedBendersDecomp(
    startLoc: pt,
    endLoc: pt,
    nodes: dict, # Index from 1
    neighbor: str = "CircleXY",
    timeLimit: int | None = None,
    **kwargs
    ) -> dict | None:

    # Model initialization ====================================================
    CETSP = grb.Model("CETSP")
    CETSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        CETSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    convergence = []
    repSeqHis = {}
    degenCuts = []

    mapType = 'XY'
    if (neighbor == 'CircleXY' or neighbor == 'PolyXY'):
        mapType = 'XY'
    elif (neighbor == 'CircleLatLon' or neighbor == 'PolyLatLon'):
        mapType = 'LatLon'

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
    elif (neighbor == "PolyLatLon"):
        for i in nodes:
            nodes[i]['polyXYMercator'] = []
            for pt in nodes[i][kwargs['polyFieldName']]:
                ptMercator = polyLatLon2XYMercator(pt)
                nodes[i]['polyXYMercator'].append(ptMercator)
        startLocMercator = ptLatLon2XYMercator(startLoc)
        endLocMercator = ptLatLon2XYMercator(endLoc)
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
    elif (neighbor == "PolyXY"):
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
    elif (neighbor == "PolyXY"):
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
                    elif (neighbor == "PolyXY"):
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
                    elif (neighbor == "PolyLatLon"):
                        # FIXME: Working
                        polys = []
                        for i in repSeq:
                            if (i != startID and i != endID):
                                polys.append(ptLatLon2XYMercator(nodes[i]))
                    
                    minDist = float('inf')
                    for seq in repSeqHis:
                        if (repSeqHis[seq] < minDist):
                            minDist = repSeqHis[seq]

                    # Cut preparing
                    residual = p2p['dist']
                    for i in range(len(repSeq) - 1):
                        residual -= zBar[repSeq[i], repSeq[i + 1]]
                    writeLog("Residual: %s" % residual)

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
    elif (neighbor == "PolyXY"):
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

def _solveCETSPGeneticAlgorithm(
    startLoc: pt,
    endLoc: pt,
    nodes: dict, # Index from 1
    neighbor: str = "Circle",
    timeLimit: int | None = None,
    popSize: int = None,
    neighRatio: dict = {},
    stop: dict = {},
    **kwargs
    ) -> dict | None:

    class chromosomeCETSP:
        def __init__(self, startLoc, endLoc, nodes, seq):
            # NOTE: seq以depotID开始和结束
            # NOTE: 每个seq都需要补全为一条合法的cetsp路径
            # Complete logic:
            # STEP 1: 记录得到所有的转折点
            # STEP 2: 记录所有未经过的点，计算距离，若没有未经过点，结束
            # STEP 3: 将最近的未经过点插入，转入STEP 2
            
            # 记录nodes的信息
            self.startLoc = startLoc
            self.endLoc = endLoc
            self.nodes = nodes
            self.initSeq = [i for i in seq]

            # 原始输入的seq
            self.seq = geoVeRoPy.Ring()
            for i in seq:
                n = geoVeRoPy.RingNode(i)
                self.seq.append(n)
            self.seq.rehead(0)

            # 转折点列表
            self.turning = []
            # 穿越点列表
            self.trespass = []
            # 未访问点及距离，暂存用，最终需要为空
            self.dist2NotInclude = {}

            # 补全
            self.seq2Path()

        def seq2Path(self):
            # 需要先得到一组turn point
            circles = []
            seqTra = [n.key for n in self.seq.traverse()]
            seqTra.append(0)
            for i in range(1, len(seqTra) - 1):
                circles.append({
                    'center': self.nodes[seqTra[i]]['loc'],
                    'radius': self.nodes[seqTra[i]]['radius']
                })
            c2c = geoVeRoPy.circle2CirclePath(
                startPt = self.startLoc,
                endPt = self.endLoc,
                circles = circles,
                algo = 'SOCP')
            degen = geoVeRoPy.seqRemoveDegen(seq = c2c['path'])

            # 找turn point/trespass point
            # NOTE: 这里的trespass point不完整，足够避免重复计算了
            self.turning = []
            self.trespass = []
            for i in range(len(degen['aggNodeList'])):
                if (degen['removedFlag'][i] == False):
                    self.turning.extend([seqTra[k] for k in degen['aggNodeList'][i]])
                else:
                    self.trespass.extend([seqTra[k] for k in degen['aggNodeList'][i]])
            self.path = degen['newSeq']
            self.dist = c2c['dist']

            # 少进行一次circle2CirclePath()
            initPathFlag = True

            completeFlag = False
            # 快速判断当前的seq是否已经覆盖，如果已经覆盖，则要求seq长度为self.nodes + 1
            if (len(seqTra) == len(self.nodes) + 1):
                completeFlag = True

            # 现在开始补齐        
            while (not completeFlag):
                completeFlag = True
                # print(self.turning, self.trespass, self.dist2NotInclude)
                
                # 先按照turnpoint构造一个路径
                if (initPathFlag):
                    # 如果是第一次循环，使用之前计算的c2c和degen，不重复计算
                    initPathFlag = False
                else:
                    circles = []
                    for i in range(1, len(self.turning) - 1):
                        circles.append({
                            'center': self.nodes[self.turning[i]]['loc'],
                            'radius': self.nodes[self.turning[i]]['radius']
                        })
                    # 得到seq对应路径
                    c2c = geoVeRoPy.circle2CirclePath(
                        startPt = self.startLoc,
                        endPt = self.endLoc,
                        circles = circles,
                        algo = 'SOCP')
                    degen = geoVeRoPy.seqRemoveDegen(seq = c2c['path'])

                    self.trespass = []
                    self.path = degen['newSeq']
                    self.dist = c2c['dist']

                # 判断剩余的点是否为trespass点
                self.dist2NotInclude = {}
                for i in self.nodes:
                    if (i not in self.turning and i not in self.trespass):
                        res = geoVeRoPy.distPt2Seq(
                            pt = self.nodes[i]['loc'], 
                            seq = degen['newSeq'],
                            closedFlag = True,
                            detailFlag = True)
                        if (res['dist'] <= self.nodes[i]['radius']):
                            self.trespass.append(i)
                        else:
                            self.dist2NotInclude[i] = res
                            completeFlag = False

                # 将一个最近的点加入turn point
                if (not completeFlag):
                    closestIdx = None
                    closestDist = float('inf')
                    closestInsertID = None
                    for i in self.dist2NotInclude:
                        if (self.dist2NotInclude[i]['dist'] < closestDist):
                            closestIdx = i
                            closestInsertID = self.dist2NotInclude[i]['nearestIdx'][1]
                    self.turning.insert(closestInsertID, closestIdx)

            # 更新seq为新的turning point
            self.seq = geoVeRoPy.Ring()
            for i in range(len(self.turning) - 1):
                n = geoVeRoPy.RingNode(self.turning[i])
                self.seq.append(n)
            self.seq.rehead(0)

    def swap(chromo, idxI):
        seq = [i.key for i in chromo.seq.traverse()]
        if (idxI < len(seq) - 1):
            seq[idxI], seq[idxI + 1] = seq[idxI + 1], seq[idxI]
        else:
            seq[idxI], seq[0] = seq[0], seq[idxI]
        return chromosomeCETSP(startLoc, endLoc, nodes, seq)

    def exchange(chromo, idxI, idxJ):
        seq = [i.key for i in chromo.seq.traverse()]
        seq[idxI], seq[idxJ] =  seq[idxJ], seq[idxI]
        return chromosomeCETSP(startLoc, endLoc, nodes, seq)

    def rotate(chromo, idxI, idxJ):
        seq = [i.key for i in chromo.seq.traverse()]
        if (idxI > idxJ):
            idxI, idxJ = idxJ, idxI
        newSeq = [seq[i] for i in range(idxI)]
        newSeq.extend([seq[idxJ - i] for i in range(idxJ - idxI + 1)])
        newSeq.extend([seq[i] for i in range(idxJ + 1, len(seq))])
        return chromosomeCETSP(startLoc, endLoc, nodes, newSeq)
    
    def crossover(chromo1, chromo2, idx1I, idx1J, idx2I, idx2J):
        # 原始序列
        seq1 = [i.key for i in chromo1.seq.traverse()]
        seq2 = [i.key for i in chromo2.seq.traverse()]

        # 把idxI和idxJ排个序换一下，保证idxI在前面
        if (idx1I > idx1J):
            idx1I, idx1J = idx1J, idx1I
        if (idx2I > idx2J):
            idx2I, idx2J = idx2J, idx2I
        # Before:
        # = = = idx1I.prev idx1I - - - idx1J idx1J.next = = =
        # = = = idx2I.prev idx2I - - - idx2J idx2J.next = = =
        # After:
        # = = = idx1I.prev idx2I - - - idx2J idx1J.next = = =
        # = = = idx2I.prev idx1I - - - idx1J idx2J.next = = =
        
        # 构造新序列
        newSeq1 = [seq2[i] for i in range(idx2I, idx2J)]
        for i in range(idx1J, len(seq1)):
            if (seq1[i] not in newSeq1):
                newSeq1.append(seq1[i])
        for i in range(idx1I):
            if (seq1[i] not in newSeq1):
                newSeq1.append(seq1[i])
        if (0 not in newSeq1):
            newSeq1.append(0)

        newSeq2 = [seq1[i] for i in range(idx1I, idx1J)]
        for i in range(idx2J, len(seq2)):
            if (seq2[i] not in newSeq2):
                newSeq2.append(seq2[i])
        for i in range(idx2I):
            if (seq2[i] not in newSeq2):
                newSeq2.append(seq2[i])
        if (0 not in newSeq2):
            newSeq2.append(0)

        newChromo1 = chromosomeCETSP(startLoc, endLoc, nodes, newSeq1)
        newChromo2 = chromosomeCETSP(startLoc, endLoc, nodes, newSeq2)
        return newChromo1, newChromo2

    # Initialize ==============================================================
    dashboard = {
        'bestOfv': float('inf'),
        'bestChromo': None
    }
    startTime = datetime.datetime.now()
    convergence = []

    # Initialize population by randomization ==================================
    popObj = []
    for k in range(popSize):
        seq = [i for i in nodes]
        seq.append(0)
        random.shuffle(seq)
        popObj.append(chromosomeCETSP(startLoc, endLoc, nodes, seq))

    for chromo in popObj:
        if (chromo.dist < dashboard['bestOfv']):
            dashboard['bestOfv'] = chromo.dist
            dashboard['bestSeq'] = [i.key for i in chromo.seq.traverse()]
            dashboard['bestChromo'] = chromo

    contFlag = True

    iterTotal = 0
    iterNoImp = 0

    while (contFlag):
        # Crossover and create offspring
        while (len(popObj) <= (int)((1 + neighRatio['crossover']) * popSize)):
            # Randomly select two genes, the better gene has better chance to have offspring            
            rnd1 = None
            rnd2 = None
            while (rnd1 == None or rnd2 == None or rnd1 == rnd2):
                coeff = []
                for i in range(len(popObj)):
                    coeff.append(1 / (popObj[i].dist - 0.8 * dashboard['bestOfv']))
                rnd1 = rndPick(coeff)
                rnd2 = rndPick(coeff)
            # Randomly select a window from the first chromo and the second chromo
            [idx1I, idx1J] = random.sample([i for i in range(popObj[rnd1].seq.count)], 2)
            [idx2I, idx2J] = random.sample([i for i in range(popObj[rnd2].seq.count)], 2)
            newSeq1, newSeq2 = crossover(popObj[rnd1], popObj[rnd2], idx1I, idx1J, idx2I, idx2J)

            popObj.append(newSeq1)
            popObj.append(newSeq2)

        # Mutation
        # NOTE: will always keep the worse outcome
        # swap
        numSwap = (int)(neighRatio['swap'] * popSize)
        for k in range(numSwap):
            rnd = random.randint(0, len(popObj) - 1)            
            idx = random.randint(0, popObj[rnd].seq.count - 1)
            popObj[rnd] = swap(popObj[rnd], idx)

        # exchange
        numExchange = (int)(neighRatio['exchange'] * popSize)
        for k in range(numExchange):
            rnd = random.randint(0, len(popObj) - 1)
            [idxI, idxJ] = random.sample([i for i in range(popObj[rnd].seq.count)], 2)
            while (abs(idxJ - idxI) <= 2
                or idxI == 0 and idxJ == popObj[rnd].seq.count - 1
                or idxI == popObj[rnd].seq.count - 1 and idxJ == 0):
                [idxI, idxJ] = random.sample([i for i in range(popObj[rnd].seq.count)], 2)
            popObj[rnd] = exchange(popObj[rnd], idxI, idxJ)

        # rotate
        numRotate = (int)(neighRatio['rotate'] * popSize)
        for k in range(numRotate):
            rnd = random.randint(0, len(popObj) - 1)
            [idxI, idxJ] = random.sample([i for i in range(popObj[rnd].seq.count)], 2)
            while (abs(idxJ - idxI) <= 2
                or idxI == 0 and idxJ == popObj[rnd].seq.count - 1
                or idxI == popObj[rnd].seq.count - 1 and idxJ == 0):
                [idxI, idxJ] = random.sample([i for i in range(popObj[rnd].seq.count)], 2)
            popObj[rnd] = rotate(popObj[rnd], idxI, idxJ)

        # Tournament
        while (len(popObj) > popSize):
            # Randomly select two genes
            rnd1 = None
            rnd2 = None
            while (rnd1 == None or rnd2 == None or rnd1 == rnd2):
                rnd1 = random.randint(0, len(popObj) - 1)
                rnd2 = random.randint(0, len(popObj) - 1)
            # kill the loser
            if (popObj[rnd1].dist > popObj[rnd2].dist):
                del popObj[rnd1]
            else:
                del popObj[rnd2]

        # Update dashboard
        newOfvFound = False
        for chromo in popObj:
            if (chromo.dist < dashboard['bestOfv']):
                newOfvFound = True
                dashboard['bestOfv'] = chromo.dist
                dashboard['bestSeq'] = [i.key for i in chromo.seq.traverse()]
                dashboard['bestChromo'] = chromo
        print(hyphenStr())
        print("Iter: ", iterTotal, "\tofv: ", dashboard['bestOfv'])
        if (newOfvFound):
            iterNoImp = 0
        else:
            iterNoImp += 1
        iterTotal += 1

        # Check stopping criteria
        if ('numNoImproveIter' in stop):
            if (iterNoImp > stop['numNoImproveIter']):
                contFlag = False
                break
        if ('numIter' in stop):
            if (iterTotal > stop['numIter']):
                contFlag = False
                break
        if ('runtime' in stop):
            if ((datetime.datetime.now() - startTime).total_seconds() > stop['runtime']):
                contFlag = False
                break

    return {
        'ofv': dashboard['bestOfv'],
        'seq': dashboard['bestSeq'],
        'path': dashboard['bestChromo'].path,
        'runtime': runtime,
        'convergence': convergence,
    }
