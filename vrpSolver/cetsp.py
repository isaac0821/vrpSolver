import math
import networkx as nx
import shapely
from shapely.geometry import mapping

# MILP Solver
AVAIL_SOLVER = 'Gurobi'
env = None
try:
    import gurobipy as grb
except:
    raise ("ERROR: Gurobi is not found.")

from .common import *
from .geometry import *
from .travel import *

def mileage2BothEndCut(repSeq, mileage: list[dict]):
    # An example ==============================================================
    # o: depot
    # x: 转折点
    # -: 穿透点
    #  o       -       -         x         -       -       x       o
    # (0) --- (7) --- (4) -- (1, 3, 8) -- (5) --- (2) --- (6) --- (9)
    #  0
    #          1
    #                  2      
    #          4
    #                  5                           5
    #                            6
    #                                      7      (7)
    #                                     (8)     (8)       8
    #                                              9
    #                                      10              
    #                                                              11
    # All-first-cut:
    # 0 => 1 => 2 => 6 => 7 => 7 => 8 => 11
    # 1 e[0, 7] + 1 e[7, 4] + 4 e[4, 1] + 0 e[1, 3] + 0 e[3, 8] + 1 e[8, 5] + 0 e[5, 2] + 1 e[2, 6] + 3 e[6, 9]
    # All-last-cut:
    # 0 => 4 => 5 => 6 => 8 => 8 => 8 => 11
    # 4 e[0, 7] + 1 e[7, 4] + 1 e[4, 1] + 0 e[1, 3] + 0 e[3, 8] + 2 e[8, 5] + 0 e[5, 2] + 0 e[2, 6] + 3 e[6, 9]

    # 记录是否穿透了某个区域
    degenFlag = False

    aggNodeList = []
    mileageByStep = []
    for i in range(len(mileage)):
        aggNodeList.append(mileage[i]['polyID'])
        mileageByStep.append(mileage[i]['mileage'])
        if (type(mileage[i]['mileage']) == list and len(mileage[i]['mileage']) > 0):
            degenFlag = True
    totalMileage = mileage[-1]['mileage']

    if (degenFlag == False):
        return []

    # Tailoring mileage =======================================================
    # NOTE: 向前和向后对齐mileage
    # NOTE: 如果裁剪时发现区间不相交，则返回不可行
    # Forwarding - 裁剪mileage区间的后半
    for i in range(len(mileage) - 1):
        # 穿透点对应的是mileage的区间,转折点则不需要
        if (mileage[i]['type'] == 'Intersect'):
            for j in range(i + 1, len(mileage)):
                # 只有穿透点才可能继续往后搜索
                if (mileage[j]['type'] == 'Intersect'):
                    # 如果下一个mileage区间的结束点在这个区间内,这个区间内真正可达的范围将被裁剪
                    if (mileageByStep[j][-1] < mileageByStep[i][-1]):
                        mileageByStep[i][-1] = mileageByStep[j][-1]
                        if (mileageByStep[j][-1] < mileageByStep[i][0]):
                            return []
                # 如果接下来是个转折点就只需要最后比对一次,然后结束
                else:
                    if (mileageByStep[j] < mileageByStep[i][-1]):
                        mileageByStep[i][-1] = mileageByStep[j]
                        if (mileageByStep[j] < mileageByStep[i][0]):
                            return []
                    break
    # Backwarding
    for i in range(1, len(mileage)):
        # 穿透点对应的是mileage的区间,转折点不需要
        if (mileage[i]['type'] == 'Intersect'):
            for j in range(1, i + 1):
                # 只有穿透点才可能继续往前搜索
                if (mileage[i - j]['type'] == 'Intersect'):
                    if (mileageByStep[i - j][0] > mileageByStep[i][0]):
                        mileageByStep[i][0] = mileageByStep[i - j][0]
                        if (mileageByStep[i - j][0] > mileageByStep[i][-1]):
                            return []
                else:
                    if (mileageByStep[i - j] > mileageByStep[i][0]):
                        mileageByStep[i][0] = mileageByStep[i - j]
                        if (mileageByStep[i - j] > mileageByStep[i][-1]):
                            return []
                    break

    # Decision variable sequence ==============================================
    dvSeq = []
    for i in range(1, len(repSeq)):
        dvSeq.append((repSeq[i - 1], repSeq[i]))
    
    # Two cuts ================================================================
    allFirstCut = []
    allFirstDegen = False
    acc = 0
    for i in range(1, len(aggNodeList)):
        # 如果不是重合点
        if (type(aggNodeList[i]) != list):
            # 如果不是转折点，all-first-cut取第一个点的mileage
            if (type(mileageByStep[i]) == list):
                allFirstCut.append(mileageByStep[i][0] - acc)
                acc = mileageByStep[i][0]
            # 如果是转折点，all-first-cut只能取该点处的mileage
            else:
                allFirstCut.append(mileageByStep[i] - acc)
                acc = mileageByStep[i]
        # 如果是重合点
        else:
            # 只留下第一个点的距离和之前点的距离
            allFirstCut.append(mileageByStep[i] - acc)
            acc = mileageByStep[i]
            # 其他的距离都设为0
            for k in range(len(aggNodeList[i]) - 1):
                allFirstCut.append(0)
                allFirstDegen = True

    allLastCut = []
    allLastDegen = False
    acc = 0
    for i in range(1, len(aggNodeList)):
        # 如果不是重合点
        if (type(aggNodeList[i]) != list):
            # 如果不是转折点，all-last-cut取第二个点的mileage
            if (type(mileageByStep[i]) == list):
                allLastCut.append(mileageByStep[i][1] - acc)
                acc = mileageByStep[i][1]
            # 如果是转折点，all-last-cut只能取该点处的mileage
            else:
                allLastCut.append(mileageByStep[i] - acc)
                acc = mileageByStep[i]
        # 如果是重合点
        else:
            # 只留下第一个点的距离和之前点的距离
            allLastCut.append(mileageByStep[i] - acc)
            acc = mileageByStep[i]
            # 其他的距离都设为0
            for k in range(len(aggNodeList[i]) - 1):
                allLastCut.append(0)
                allLastDegen = True

    return [
        {'coeff': allFirstCut, 'dvSeq': dvSeq, 'degenFlag': allFirstDegen},
        {'coeff': allLastCut, 'dvSeq': dvSeq, 'degenFlag': allLastDegen},
    ]

def solveCETSP(
    startLoc: pt,
    endLoc: pt,
    nodes: dict,
    radius: float | str, 
    cutSetting: dict,
    timeLimit: int | None = None,
    ) -> dict | None:

    # Model initialization ====================================================
    CETSP = grb.Model("CETSP")
    CETSP.setParam('OutputFlag', 0)
    if (timeLimit != None):
        CETSP.setParam(grb.GRB.Param.TimeLimit, timeLimit)
    convergence = []
    repSeqHis = []
    GBDCuts = {}

    startTime = datetime.datetime.now()

    # Define sets
    startID = 0
    endID = len(nodes) + 1
    nodeAll = [i for i in range(0, len(nodes) + 2)]

    # Parameters ==============================================================
    # anchor starts from depotLoc, in between are a list of circles, ends with depotLoc
    allX = [startLoc[0]]
    allY = [startLoc[1]]
    for i in nodes:
        r = radius if type(radius) != str else nodes[i][radius]
        allX.append(nodes[i]['loc'][0] - r)
        allX.append(nodes[i]['loc'][0] + r)
        allY.append(nodes[i]['loc'][1] - r)
        allY.append(nodes[i]['loc'][1] + r)
    allX.append(endLoc[0])
    allY.append(endLoc[1])

    lbX = min(allX) - 1
    lbY = min(allY) - 1
    ubX = max(allX) + 1
    ubY = max(allY) + 1

    # Define visiting matrix for the centers of disk
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

    zBar = {}
    for i in nodes:
        for j in nodes:
            if (i != j):
                sumR = 2 * radius if type(radius) != str else nodes[i][radius] + nodes[j][radius]
                zBar[i, j] = max(tau[i, j] - sumR, 0)
    for i in nodes:
        startR = radius if type(radius) != str else nodes[i][radius]
        zBar[startID, i] = max(tauStart[i] - startR, 0)
        zBar[i, startID] = max(tauStart[i] - startR, 0)
        
        endR = radius if type(radius) != str else nodes[i][radius]
        zBar[i, endID] = max(tauEnd[i] - endR, 0)
        zBar[endID, i] = max(tauEnd[i] - endR, 0)
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
            feasibleFlag = True

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
                    feasibleFlag = False
                    model.cbLazy(grb.quicksum(e[i, j] for i in component for j in component if i != j) <= len(component) - 1)
    
            # Benders cut -----------------------------------------------------
            if (feasibleFlag):
                repArc = []
                for (i, j) in e.keys():
                    if (e_sol[i, j] > 0.9):
                        repArc.append([i, j])

                # 还原访问顺序
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
                repSeqHis.append(repSeq)

                writeLog("\n")
                writeLog(list2String(repSeq))
                writeLog(hyphenStr())

                circles = []
                neiPolys = {}
                for i in repSeq:
                    if (i != startID and i != endID):
                        circles.append({
                            'center': nodes[i]['loc'],
                            'radius': radius if type(radius) != str else nodes[i][radius]
                        })
                        neiPolys[i] = {
                            'neighbor': circleByCenterXY(
                                center = nodes[i]['loc'], 
                                radius = radius if type(radius) != str else nodes[i][radius])
                        }

                p2p = circle2CirclePath(
                    startPt = startLoc, 
                    endPt = endLoc, 
                    circles = circles)
                mileage = polyPath2Mileage(repSeq, p2p['path'], neiPolys)

                # Add default cut ---------------------------------------------
                if ('defaultCut' in cutSetting and cutSetting['defaultCut'] == True):
                    coeff = []
                    dvSeq = []
                    for i in range(len(repSeq) - 1):
                        coeff.append(distEuclideanXY(p2p['path'][i], p2p['path'][i + 1])['dist'])
                        dvSeq.append((repSeq[i], repSeq[i + 1]))
                    writeLog(GBDCutInfo(coeff, dvSeq, "default cut"))
                    model.cbLazy(theta >= grb.quicksum(
                        e[dvSeq[i]] * (coeff[i] - zBar[dvSeq[i]]) for i in range(len(dvSeq))))
                    GBDCuts[len(GBDCuts)] = {
                        'type': 'Default',
                        'seq': repSeq,
                        'dist': p2p['dist'],
                        'zBar': sum([zBar[dvSeq[i]] for i in range(len(dvSeq))]),
                        'theta': sum([(coeff[i] - zBar[dvSeq[i]]) for i in range(len(dvSeq))]),
                        'coeff': [i for i in coeff],
                        'dvSeq': [i for i in dvSeq]
                    }

                # Add multi-cuts ----------------------------------------------
                # NOTE: All-first-cut/all-last-cut/reverse-all-first-cut/reverse-all-last-cut
                if ('multiCut' in cutSetting and cutSetting['multiCut'] == True):
                    multiCut = mileage2BothEndCut(repSeq, mileage)
                    for cut in multiCut:
                        if ('multiCutDegenOnlyFlag' in cutSetting and cutSetting['multiCutDegenOnlyFlag'] == True):
                            if (cut['degenFlag'] == True):
                                writeLog(GBDCutInfo(cut['coeff'], cut['dvSeq'], "multi-degen cut"))
                                model.cbLazy(theta >= grb.quicksum(
                                    e[cut['dvSeq'][i]] * (cut['coeff'][i] - zBar[cut['dvSeq'][i]]) for i in range(len(cut['dvSeq']))))
                                GBDCuts[len(GBDCuts)] = {
                                    'type': 'Multi-cut-Degen',
                                    'seq': repSeq,
                                    'dist': p2p['dist'],
                                    'zBar': sum([zBar[cut['dvSeq'][i]] for i in range(len(cut['dvSeq']))]),
                                    'theta': sum([(cut['coeff'][i] - zBar[cut['dvSeq'][i]]) for i in range(len(cut['dvSeq']))]),
                                    'coeff': [i for i in cut['coeff']],
                                    'dvSeq': [i for i in cut['dvSeq']]
                                }
                        else:
                            writeLog(GBDCutInfo(cut['coeff'], cut['dvSeq'], "multi cut"))
                            model.cbLazy(theta >= grb.quicksum(
                                e[cut['dvSeq'][i]] * (cut['coeff'][i] - zBar[cut['dvSeq'][i]]) for i in range(len(cut['dvSeq']))))
                            GBDCuts[len(GBDCuts)] = {
                                'type': 'Multi-cut',
                                'seq': repSeq,
                                'dist': p2p['dist'],
                                'zBar': sum([zBar[cut['dvSeq'][i]] for i in range(len(cut['dvSeq']))]),
                                'theta': sum([(cut['coeff'][i] - zBar[cut['dvSeq'][i]]) for i in range(len(cut['dvSeq']))]),
                                'coeff': [i for i in cut['coeff']],
                                'dvSeq': [i for i in cut['dvSeq']]
                            }

                objBound = model.cbGet(grb.GRB.Callback.MIPSOL_OBJBND)
                objIncum = model.cbGet(grb.GRB.Callback.MIPSOL_OBJBST)
                timePassed = round((datetime.datetime.now() - startTime).total_seconds(), 2)
                writeLog("Time Pass: " + str(timePassed) + "[s]"
                    + "\nCurrSol Dist: " + str(p2p['dist'])
                    + "\nCurrObj: " + str(model.cbGet(grb.GRB.Callback.MIPSOL_OBJ))
                    + "\nCurrObjBest: " + str(objIncum)
                    + "\nCurrBound: " + str(objBound))
                convergence.append((p2p['dist'], objIncum, objBound, timePassed))

    # TSP with no callback ====================================================
    CETSP.update()
    # CETSP.write("Origin.lp")
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
    circles = []
    for i in seq:
        if (i != startID and i != endID):
            circles.append({
                'center': nodes[i]['loc'],
                'radius': radius if type(radius) != str else nodes[i][radius]
            })
    p2p = circle2CirclePath(startPt = startLoc, endPt = endLoc, circles = circles)

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
        'numIter': len(repSeqHis),
        'numGBD': len(GBDCuts)
    }
