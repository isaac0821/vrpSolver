import gurobipy as grb
import networkx as nx

from .geometry import *
from .common import *
from .msg import *
from .ds import *

# Loc seq related =============================================================
def seqSortPts(seq: list[pt], pts: list[pt], allowNotIncludeFlag = True) -> list:
    # First, calculate accumulated dist since start for each turning point of seq
    acc = [0]
    sofar = 0
    for i in range(len(seq) - 1):
        sofar += distEuclideanXY(seq[i], seq[i + 1])['dist']
        acc.append(sofar)
    ptHeap = []
    notIncludedPts = []

    # 看看每个节点落在哪个区间
    for pt in pts:
        isOnSegFlag = False
        for i in range(len(seq) - 1):
            # 若和多个线段相交，只考虑第一次相交的情况，否则turnpoint上的始终会被记录两次
            if (isPtOnSeg(pt, [seq[i], seq[i + 1]])):
                isOnSegFlag = True
                addDist = distEuclideanXY(seq[i], pt)['dist']
                heapq.heappush(ptHeap, (acc[i] + addDist, pt))
                break
        if (not isOnSegFlag):
            if (allowNotIncludeFlag):
                notIncludedPts.append(pt)
            else:
                raise UnsupportedInputError("ERROR: %s does not belong to given loc sequence" % str(pt))

    sortedPts = []
    while (len(ptHeap) > 0):
        sortedPts.append(heapq.heappop(ptHeap)[1])

    return {
        'sortedPts': sortedPts,
        'notIncludedPts': notIncludedPts
    }

def seqSortNodes(seq: list[pt], nodes: dict, locFieldName: str = 'loc', allowNotIncludeFlag: bool = True) -> list: 
    acc = [0]
    sofar = 0
    for i in range(len(seq) - 1):
        sofar += distEuclideanXY(seq[i], seq[i + 1])['dist']
        acc.append(sofar)
    ptHeap = []
    notIncludedNodeIDs = []

    for nodeID in nodes:
        isOnSegFlag = False
        for i in range(len(seq) - 1):
            if (isPtOnSeg(nodes[nodeID][locFieldName], [seq[i], seq[i + 1]])):
                isOnSegFlag = True
                addDist = distEuclideanXY(seq[i], nodes[nodeID][locFieldName])['dist']
                print(acc[i] + addDist, nodeID)
                heapq.heappush(ptHeap, (acc[i] + addDist, nodeID))
                break
        if (not isOnSegFlag):
            if (allowNotIncludeFlag):
                notIncludedNodeIDs.append(nodeID)
            else:
                raise UnsupportedInputError("ERROR: %s does not belong to given loc sequence" % str(nodeID))
    
    sortedNodeIDs = []
    while (len(ptHeap) > 0):
        sortedNodeIDs.append(heapq.heappop(ptHeap)[1])

    return {
        'sortedNodeIDs': sortedNodeIDs,
        'notIncludedNodeIDs': notIncludedNodeIDs
    }

# Mileage of points on seq ====================================================
def mileagePt(seq: list[pt], pt: pt, error = CONST_EPSILON) -> None | float | list[float]:
    """
    Given a sequence of locs, and a point, returns the mileage of point on the seq, None if not on the sequence

    Parameters
    ----------
    seq: list[pt], Required
        A sequence of locations
    pt: pt, Required
        A point

    Return
    ------
    None | float | list[float]
        None if pt is not on the sequence, float if only one mileage, list[float] if the pt appears multiple times

    """
    seq = seqRemoveDegen(seq)['newSeq']

    m = []

    # 找到所有pt所处的线段
    acc = 0
    onSeg = []
    closest2Snap = float('inf')
    bestSnap = None
    closestSeg = None
    for i in range(len(seq) - 1):
        seg = [seq[i], seq[i + 1]]
        leng = distEuclideanXY(seq[i], seq[i + 1])['dist']
        if (isPtOnSeg(pt, seg, interiorOnly=False, error=error)):
            onSeg.append((acc, seg))
        # 如果
        elif (len(onSeg) == 0):
            snapPt = ptFoot2Line(pt, seg)
            dist2Snap = distEuclideanXY(pt, snapPt)['dist']
            if (dist2Snap < closest2Snap):
                closest2Snap = dist2Snap
                closestSeg = (acc, seg)
                bestSnap = snapPt
        acc += leng 
    # 如果pt不在任何一段线段上,那么把pt给snap到最近的线段上
    if (len(onSeg) == 0):
        onSeg.append(closestSeg)
        pt = bestSnap

    # 逐个计算mileage
    for i in range(len(onSeg)):
        dist2Ori = distEuclideanXY(pt, onSeg[i][1][0])['dist']
        newMileage = dist2Ori + onSeg[i][0]
        repeatedFlag = False
        for exist in m:
            if (abs(exist - newMileage) <= CONST_EPSILON):
                repeatedFlag = True
        if (not repeatedFlag):
            m.append(newMileage)

    if (len(m) == 0):
        return None
    elif (len(m) == 1):
        return m[0]
    else:
        return m

# Path touring through polygons ===============================================
def seqRemoveDegen(seq: list[pt], error:float=CONST_EPSILON):
    """
    Given a sequence of points, returns a subset of points that only includes turning points of the sequence. If there are multiple points overlapped at the same location, keeps one of those points.

    Parameters
    ----------
    seq: list[pt], required
        The coordinates of a sequence of points.
    error: float, optional, default as CONST_EPSILON
        Error tolerance.

    Returns
    -------
    dict
        A new dictionary, in the format of 
        
        >>> {
        ...     'newSeq': newSeq, 
        ...     'aggNodeList': aggNodeList, 
        ...     'removedFlag': removedFlag, 
        ...     'locatedSeg': locatedSeg
        ... }

        - The 'newSeq' returns a new sequence which only has turn points of the origin sequence
        - The 'aggNodeList' is a list of lists, if a point overlaps with its previous/next point, the index of both points will be aggregated into the same list.
        - The 'removedFlag' indicates whether a point is removed as a non-turning point, true if the point is not included in the 'newSeq'
        - The 'locatedSeq' returns a list of line segments, for each point removed, it returns the line segment it belongs to 

    Examples
    --------
    For the following inputs
        >>> seq = [[-1, 0], [0, 0], [0, 0], [1, 0], [2, 0], [2, 1], [1, 1], [1, 0], [1, -1]]
        >>> res = seqRemoveDegen(seq)
    The result is as follows
        >>> res = {
        ...     'newSeq': [[-1, 0], [2, 0], [2, 1], [1, 1], [1, -1]],
        ...     'aggNodeList': [[0], [1, 2], [3], [4], [5], [6], [7], [8]],
        ...     'removedFlag': [False, True, True, False, False, False, True, False],
        ...     'locatedSeg': [None,
        ...         [[-1, 0], [2, 0]],
        ...         [[-1, 0], [2, 0]],
        ...         None,
        ...         None,
        ...         None,
        ...         [[1, 1], [1, -1]],
        ...         None]
        ... }
    The result shows that, the new sequence is [[-1, 0], [2, 0], [2, 1], [1, 1], [1, -1]], in the new sequence, seq[1], seq[2], seq[6] are not included since they are not turn points.
    seq[1] and seq[2] are aggregated due to overlaps. Although seq[3] and seq[6] are overlapped, they are not aggregated because they are not neighboring. 
    For the removed points, 'locatedSeq' finds the segment they located.

    """

    # Step 1: 先按是否重合对点进行聚合  
    curLocList = [seq[0]]
    curAgg = [0]
    aggNodeList = []
    # 挤香肠算法
    for i in range(1, len(seq)):
        # 如果当前点和任意一个挤出来的点足够近，则计入
        samePtFlag = False

        for pt in curLocList:
            if (is2PtsSame(pt, seq[i], error)):
                curAgg.append(i)
                curLocList.append(seq[i])
                samePtFlag = True
                break

        # 若不重复，了结
        if (not samePtFlag):
            aggNodeList.append([k for k in curAgg])
            curAgg = [i]
            curLocList = [seq[i]]

    aggNodeList.append([k for k in curAgg])

    # Step 2: 对聚合后的点，判断是否为转折点
    # NOTE: removeFlag的长度和aggNodeList一致
    removedFlag = [False]
    for i in range(1, len(aggNodeList) - 1):
        preLoc = seq[aggNodeList[i - 1] if type(aggNodeList[i - 1]) != list else aggNodeList[i - 1][0]]
        curLoc = seq[aggNodeList[i] if type(aggNodeList[i]) != list else aggNodeList[i][0]]
        sucLoc = seq[aggNodeList[i + 1] if type(aggNodeList[i + 1]) != list else aggNodeList[i + 1][0]]

        if (is2PtsSame(preLoc, sucLoc)):
            removedFlag.append(False)
        else:
            dev = distPt2Seg(curLoc, [preLoc, sucLoc])
            if (dev <= error):
                removedFlag.append(True)
            else:
                removedFlag.append(False)
    removedFlag.append(False)

    # 得到去掉共线和重合点后的折线
    newSeq = []
    for i in range(len(aggNodeList)):
        if (removedFlag[i] == False):
            if (type(aggNodeList[i]) == list):
                newSeq.append(seq[aggNodeList[i][0]])
            else:
                newSeq.append(seq[aggNodeList[i]])

    # 对于被移除的共线点，找到其所在的线段
    locatedSeg = []
    # 把没有移除的点的序号记一下
    seqPre = []
    seqSuc = []
    # 查找移除点之前和之后一个removeFlag为False的对应aggNode，得到对应线段
    for i in range(len(removedFlag)):
        # Log the prev that is not removed
        if (removedFlag[i] == False):
            seqPre.append(i)
        else:
            seqPre.append(seqPre[-1])
        # Log the next that is not removed
        if (removedFlag[len(removedFlag) - 1 - i] == False):
            seqSuc.insert(0, len(removedFlag) - 1 - i)
        else:
            seqSuc.insert(0, seqSuc[0])

    for i in range(len(removedFlag)):
        if (removedFlag[i] == False):
            locatedSeg.append(None)
            para = []
        else:
            startLoc = None
            if (type(aggNodeList[seqPre[i]]) == list):
                startLoc = seq[aggNodeList[seqPre[i]][0]]
            else:
                startLoc = seq[aggNodeList[seqPre[i]]]
            endLoc = None
            if (type(aggNodeList[seqSuc[i]]) == list):
                endLoc = seq[aggNodeList[seqSuc[i]][0]]
            else:
                endLoc = seq[aggNodeList[seqSuc[i]]]
            locatedSeg.append([startLoc, endLoc])

    return {
        'newSeq': newSeq,
        'aggNodeList': aggNodeList,
        'removedFlag': removedFlag,
        'locatedSeg': locatedSeg
    }


def ptSetSeq2Poly(seq, polygons:dict, polyFieldName = 'polygon'):
    """Given a sequence and a dictionary of polygons, finds the intersection points between seq and polygons

    Parameters
    ----------
    seq: list[pt], required
        A list of points as a sequence
    polygons: dict, required
        A dictionary, each key is the ID of polygon, the field of polygon is in `polyFieldName`
    polyFieldName: string, optional, default 'polygon'
        Default field name for polygon

    Return
    ------
    list of dict
        A list of dictionaries, each includes the information of a segment/point, a list of polygon IDs that the segment belongs to

    """

    if (polygons == None):
        raise MissingParameterError("ERROR: Missing required field 'polygons'.")

    # First, for each leg in the seq, find the individual polygons intersect with the leg
    actions = []
    accMileage = 0

    for i in range(len(seq) - 1):
        # seg[0]的mileage更小
        seg = [seq[i], seq[i + 1]]

        # NOTE: 准备用segment tree
        # NOTE: For now use the naive way - checking the bounding box
        for pID in polygons:
            # 根据Seg和poly的相交情况做判断
            segIntPoly = intSeg2Poly(seg = seg, poly = polygons[pID][polyFieldName])
            # 如果相交得到多个部分，则分别进行处理

            ints = []
            if (type(segIntPoly) == list):
                for intPart in segIntPoly:
                    ints.append(intPart)
            else:
                ints = [segIntPoly]

            for intPart in ints:
                if (intPart['status'] == 'NoCross'):
                    # No intersection pass
                    pass

                elif (intPart['status'] == 'Cross' and intPart['intersectType'] == 'Point'):
                    intPt = intPart['intersect']

                    # 距离seg[0]有多远，加入mileage
                    m = distEuclideanXY(seg[0], intPt)['dist']
                    actions.append({
                        'loc': intPt,
                        'action': 'touch',
                        'polyID': pID,
                        'mileage': accMileage + m
                    })

                elif (intPart['status'] == 'Cross' and intPart['intersectType'] == 'Segment'):

                    # 线段和polygon的两个交点，分别判断是不是在polygon的内部
                    intPt1 = intPart['intersect'][0]
                    intPt2 = intPart['intersect'][1]
                    intPt1InnerFlag = False
                    intPt2InnerFlag = False
                    if (isPtInPoly(intPt1, polygons[pID][polyFieldName], interiorOnly=True)):
                        intPt1InnerFlag = True
                    if (isPtInPoly(intPt2, polygons[pID][polyFieldName], interiorOnly=True)):
                        intPt2InnerFlag = True

                    # 两个端点的mileage
                    m1 = distEuclideanXY(seg[0], intPt1)['dist']
                    m2 = distEuclideanXY(seg[0], intPt2)['dist']

                    if (m1 > m2):
                        m1, m2 = m2, m1
                        intPt1, intPt2 = intPt2, intPt1
                        intPt1InnerFlag, intPt2InnerFlag = intPt2InnerFlag, intPt1InnerFlag

                    if (abs(m1 - m2) <= CONST_EPSILON):
                        # 交点距离太近可能是误判
                        actions.append({
                            'loc': intPt1,
                            'action': 'touch',
                            'polyID': pID,
                            'mileage': accMileage + m1
                        })
                    else:
                        # 如果第一个点在poly内，第二个在poly外
                        # m1 => inside, m2 => leave
                        if (intPt1InnerFlag and not intPt2InnerFlag):
                            actions.append({
                                'loc': intPt2,
                                'action': 'leave',
                                'polyID': pID,
                                'mileage': accMileage + m2,
                            })

                        # 如果第一个点在poly外，第二个点在poly内
                        # m1 => enter, m2 => inside
                        elif (not intPt1InnerFlag and intPt2InnerFlag):
                            actions.append({
                                'loc': intPt1,
                                'action': 'enter',
                                'polyID': pID,
                                'mileage': accMileage + m1,
                            })

                        # 如果两点都在内部，整段折线都在内部
                        elif (intPt1InnerFlag and intPt2InnerFlag):
                            pass

                        # 如果俩点都在外面，只有可能是两个转折点都在边缘上，一进一出
                        elif (not intPt1InnerFlag and not intPt2InnerFlag):
                            actions.append({
                                'loc': intPt1,
                                'action': 'enter',
                                'polyID': pID,
                                'mileage': accMileage + m1,
                            })
                            actions.append({
                                'loc': intPt2,
                                'action': 'leave',
                                'polyID': pID,
                                'mileage': accMileage + m2,
                            })
    
        accMileage += distEuclideanXY(seg[0], seg[1])['dist']

    actions = sorted(actions, key = lambda d: d['mileage'])

    return actions

def splitOverlapSegs(seg1: line, seg2: line, belong1: list = [1], belong2: list = [2]):
    """Given two segments that are overlapped, split into segments that are not overlapped

    Parameters
    ----------
    seg1: line, required
        The first segment
    seg2: line, required
        The second segment
    belong1: list, optional, default [1]
        A sequence of indexes, which seg1 represents
    belong2: list, optional, default [2]
        A sequence of indexes, which seg2 represents

    Return
    ------
    list of dict
        A list of dictionary with split segments and belonging info.

    """

    s2s = intSeg2Seg(seg1, seg2)
    if (s2s['status'] == 'Collinear' and s2s['interiorFlag'] == True and s2s['intersectType'] == 'Segment'):
        # print("Part ", i, segSet[i], "Part ", j, segSet[j], "S2S")
        # 取出seg1和seg2的端点
        A = seg1[0]
        B = seg1[1]
        X = seg2[0]
        Y = seg2[1]

        AXB = isPtOnSeg(pt = X, seg = [A, B], interiorOnly = True)
        AYB = isPtOnSeg(pt = Y, seg = [A, B], interiorOnly = True)
        XAY = isPtOnSeg(pt = A, seg = [X, Y], interiorOnly = True)
        XBY = isPtOnSeg(pt = B, seg = [X, Y], interiorOnly = True)

        sameAX = is2PtsSame(A, X)
        sameAY = is2PtsSame(A, Y)
        sameBX = is2PtsSame(B, X)
        sameBY = is2PtsSame(B, Y)

        distAX = distEuclideanXY(A, X)['dist']
        distBX = distEuclideanXY(B, X)['dist']
        distAY = distEuclideanXY(A, Y)['dist']
        distBY = distEuclideanXY(B, Y)['dist']

        jointBelong = [k for k in belong1]
        jointBelong.extend([k for k in belong2 if k not in belong1])
        
        newSegSet = []
        # NOTE: This is so stupid, but efficient
        # Case 1: A - X - Y - B
        if (AXB and AYB and distAX < distAY):
            # print("# Case 1: A - X - Y - B")
            newSegSet.append({
                'shape': [A, X],
                'belong': [k for k in belong1],
            })
            newSegSet.append({
                'shape': [X, Y],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [Y, B],
                'belong': [k for k in belong1],
            })
        # Case 2: A - Y - X - B
        elif (AXB and AYB and distAX > distAY):
            # print("# Case 2: A - Y - X - B")
            newSegSet.append({
                'shape': [A, Y],
                'belong': [k for k in belong1],
            })
            newSegSet.append({
                'shape': [Y, X],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [X, B],
                'belong': [k for k in belong1],
            })
        # Case 3: A - X - B - Y
        elif (AXB and XBY):
            # print("# Case 3: A - X - B - Y")
            newSegSet.append({
                'shape': [A, X],
                'belong': [k for k in belong1],
            })
            newSegSet.append({
                'shape': [X, B],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [B, Y],
                'belong': [k for k in belong2],
            })
        # Case 4: A - Y - B - X
        elif (AYB and XBY):
            # print("# Case 4: A - Y - B - X")
            newSegSet.append({
                'shape': [A, Y],
                'belong': [k for k in belong1],
            })
            newSegSet.append({
                'shape': [Y, B],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [B, X],
                'belong': [k for k in belong2],
            })
        # Case 5: X - A - Y - B
        elif (XAY and AYB):
            # print("# Case 5: X - A - Y - B")
            newSegSet.append({
                'shape': [X, A],
                'belong': [k for k in belong2],
            })
            newSegSet.append({
                'shape': [A, Y],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [Y, B],
                'belong': [k for k in belong1],
            })
        # Case 6: Y - A - X - B
        elif (XAY and AXB):
            # print("# Case 6: Y - A - X - B")
            newSegSet.append({
                'shape': [Y, A],
                'belong': [k for k in belong2],
            })
            newSegSet.append({
                'shape': [A, X],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [X, B],
                'belong': [k for k in belong1],
            })
        # Case 7: X - A - B - Y
        elif (XAY and XBY and distAX < distBX):
            # print("# Case 7: X - A - B - Y")
            newSegSet.append({
                'shape': [X, A],
                'belong': [k for k in belong2],
            })
            newSegSet.append({
                'shape': [A, B],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [B, Y],
                'belong': [k for k in belong2],
            })
        # Case 8: Y - A - B - X
        elif (XAY and XBY and distAX > distBX):
            # print("# Case 8: Y - A - B - X")
            newSegSet.append({
                'shape': [Y, A],
                'belong': [k for k in belong2],
            })
            newSegSet.append({
                'shape': [A, B],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [B, X],
                'belong': [k for k in belong2],
            })
        # Case 9: A - X - (BY)
        elif (AXB and sameBY):
            # print("# Case 9: A - X - (BY)")
            newSegSet.append({
                'shape': [A, X],
                'belong': [k for k in belong1],
            })
            newSegSet.append({
                'shape': [X, B],
                'belong': [k for k in jointBelong],
            })
        # Case 10: A - Y - (BX)
        elif (AYB and sameBX):
            # print("# Case 10: A - Y - (BX)")
            newSegSet.append({
                'shape': [A, Y],
                'belong': [k for k in belong1],
            })
            newSegSet.append({
                'shape': [Y, B],
                'belong': [k for k in jointBelong],
            })
        # Case 11: X - A - (BY)
        elif (XAY and sameBY):
            # print("# Case 11: X - A - (BY)")
            newSegSet.append({
                'shape': [X, A],
                'belong': [k for k in belong2],
            })
            newSegSet.append({
                'shape': [A, B],
                'belong': [k for k in jointBelong],
            })
        # Case 12: Y - A - (XB)
        elif (XAY and sameBX):
            # print("# Case 12: Y - A - (XB)")
            newSegSet.append({
                'shape': [Y, A],
                'belong': [k for k in belong2],
            })
            newSegSet.append({
                'shape': [A, B],
                'belong': [k for k in jointBelong],
            })
        # Case 13: (AX) - Y - B
        elif (AYB and sameAX):
            # print("# Case 13: (AX) - Y - B")
            newSegSet.append({
                'shape': [A, Y],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [Y, B],
                'belong': [k for k in belong1],
            })
        # Case 14: (AX) - B - Y
        elif (XBY and sameAX):
            # print("# Case 14: (AX) - B - Y")
            newSegSet.append({
                'shape': [A, B],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [B, Y],
                'belong': [k for k in belong2],
            })
        # Case 15: (AY) - X - B
        elif (AXB and sameAY):
            # print("# Case 15: (AY) - X - B")
            newSegSet.append({
                'shape': [A, X],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [X, B],
                'belong': [k for k in belong1],
            })
        # Case 16: (AY) - B - X
        elif (XBY and sameAY):
            # print("# Case 16: (AY) - B - X")
            newSegSet.append({
                'shape': [A, B],
                'belong': [k for k in jointBelong],
            })
            newSegSet.append({
                'shape': [B, X],
                'belong': [k for k in belong2],
            })
        # Case 17: (AX) - (BY) or (AY) - (BX)
        elif (is2SegsSame([A, B], [X, Y])):
            # print("# Case 17: (AX) - (BY) or (AY) - (BX)")
            newSegSet.append({
                'shape': [A, B],
                'belong': [k for k in jointBelong],
            })
        return newSegSet
    return None

def segSetSeq2Poly(seq: list, polygons: dict, polyFieldName: str = 'polygon', seqDegenFlag: bool = False):
    """Given a sequence and a dictionary of polygons, finds the intersection between seq and polygons

    Parameters
    ----------
    seq: list[pt], required
        A list of points as a sequence
    polygons: dict, required
        A dictionary, each key is the ID of polygon, the field of polygon is in `polyFieldName`
    polyFieldName: string, optional, default 'polygon'
        Default field name for polygon
    seqDegenFlag: boolean, optional, default False
        True if the sequence has already been processed to remove degenerated points.

    Return
    ------
    list of dict
        A list of dictionaries, each includes the information of a segment/point, a list of polygon IDs that the segment belongs to

    """

    # NOTE: 简化线段，去掉穿越点，如果已经处理过了，就不用重复计算了
    if (not seqDegenFlag):
        seq = seqRemoveDegen(seq, error = 0.03)['newSeq']

    # Step 0: turnPts =========================================================
    # NOTE: 所有的转折点必须在一个多边形的边缘上，所以必须给每个turn point找到一个polygon
    # NOTE: 只有出现了转折点本身在一个多边形的边缘上才需要记录tangle的形式，否则都会作为seg的一部分
    turnPts = {}
    for i in range(len(seq)):
        pt = seq[i]
        # 记录到每个polygon的距离
        tansPolys = []
        inerPolys = []
        for p in polygons:
            d2Edge = distPt2Seq(pt = pt, seq = polygons[p][polyFieldName], closedFlag = True)
            if (d2Edge <= 0.03):
                tansPolys.append(p)
                inerPolys.append(p)
            if (isPtInPoly(pt=pt, poly=polygons[p][polyFieldName])):
                inerPolys.append(p)
        turnPts[i] = {
            'loc': seq[i],
            'tansPolys': tansPolys,
            'inerPolys': inerPolys
        }

    # Step 1: Initialize ======================================================
    segIntPoly = {}
    accMileage = 0
    for i in range(len(seq) - 1):
        segIntPoly[i] = {
            'seg': [seq[i], seq[i + 1]],
            'startMileage': accMileage,
            'intPolys': [],
            'stTangPt': None,
        }
        if (len(turnPts[i]['tansPolys']) > 0):
            segIntPoly[i]['stTangPt'] = {
                'shape': seq[i],
                'type': 'Point',
                'belong': turnPts[i]['inerPolys'],
                'mileage': accMileage
            }
        accMileage += distEuclideanXY(seq[i], seq[i + 1])['dist']
        segIntPoly[i]['endMileage'] = accMileage

    # Step 2: For each segment, gets polygon intersected ======================
    for i in segIntPoly:
        segIntPoly[i]['intMileageList'] = []
        st = segIntPoly[i]['startMileage']
        ed = segIntPoly[i]['endMileage']

        segIntPoly[i]['intMileageList'].append((st, -1, segIntPoly[i]['seg'][0], 'start'))

        for p in polygons:
            segInt = intSeg2Poly([seq[i], seq[i + 1]], polygons[p][polyFieldName])
            # 如果交出多个来的话，分别计算
            if (type(segInt) == list):
                for s in segInt:
                    # 如果交出一个线段，计算线段的起始结束mileage
                    if (s['status'] == 'Cross' and s['intersectType'] == 'Segment'):
                        int1 = s['intersect'][0]
                        int2 = s['intersect'][1]
                        d1 = distEuclideanXY(segIntPoly[i]['seg'][0], int1)['dist']
                        d2 = distEuclideanXY(segIntPoly[i]['seg'][0], int2)['dist']
                        if (d1 < d2):
                            segIntPoly[i]['intMileageList'].append((d1 + st, p, int1, 'enter'))
                            segIntPoly[i]['intMileageList'].append((d2 + st, p, int2, 'leave'))
                        else:
                            segIntPoly[i]['intMileageList'].append((d2 + st, p, int2, 'enter'))
                            segIntPoly[i]['intMileageList'].append((d1 + st, p, int1, 'leave'))
                            
                    # 如果交出一个点，计算点的mileage, 这个点不能是两边端点
                    elif (s['status'] == 'Cross' and s['intersectType'] == 'Point'):
                        intP = s['intersect']
                        if (not is2PtsSame(intP, segIntPoly[i]['seg'][0]) and not is2PtsSame(intP, segIntPoly[i]['seg'][1])):
                            dP = distEuclideanXY(segIntPoly[i]['seg'][0], intP)['dist']
                            segIntPoly[i]['intMileageList'].append((dP + st, p, intP, 'tangle'))

            # 如果交出一个线段，计算线段的起始结束mileage
            elif (segInt['status'] == 'Cross' and segInt['intersectType'] == 'Segment'):
                int1 = segInt['intersect'][0]
                int2 = segInt['intersect'][1]
                d1 = distEuclideanXY(segIntPoly[i]['seg'][0], int1)['dist']
                d2 = distEuclideanXY(segIntPoly[i]['seg'][0], int2)['dist']
                if (d1 < d2):
                    segIntPoly[i]['intMileageList'].append((d1 + st, p, int1, 'enter'))
                    segIntPoly[i]['intMileageList'].append((d2 + st, p, int2, 'leave'))
                else:
                    segIntPoly[i]['intMileageList'].append((d2 + st, p, int2, 'enter'))
                    segIntPoly[i]['intMileageList'].append((d1 + st, p, int1, 'leave'))
            # 如果交出一个点，计算点的mileage
            elif (segInt['status'] == 'Cross' and segInt['intersectType'] == 'Point'):
                intP = segInt['intersect']
                dP = distEuclideanXY(segIntPoly[i]['seg'][0], intP)['dist']
                segIntPoly[i]['intMileageList'].append((dP + st, p, intP, 'tangle'))

        # 对intMileageList进行排序
        segIntPoly[i]['intMileageList'].sort()
        segIntPoly[i]['intMileageList'].append((ed, -1, segIntPoly[i]['seg'][1], 'end'))

    # Step 3: Restore =========================================================
    segSet = []
    for i in segIntPoly:
        polyInside = []
        segIntPoly[i]['segSet'] = []
        curMileage = segIntPoly[i]['startMileage']
        curPt = segIntPoly[i]['seg'][0]
        if (segIntPoly[i]['stTangPt'] != None):
            segSet.append(segIntPoly[i]['stTangPt'])

        for k in range(len(segIntPoly[i]['intMileageList'])):
            # 下一个点
            newMileage = segIntPoly[i]['intMileageList'][k][0]
            newPt = segIntPoly[i]['intMileageList'][k][2]

            # 如果mileage不增长，则不会单独交一段出来，除非是有tangle
            if (abs(newMileage - curMileage) > CONST_EPSILON):
                segSet.append({
                    'shape': [curPt, newPt],
                    'type': 'Segment',
                    'belong': [k for k in polyInside],
                    'mileage': [curMileage, newMileage]
                    })
                curMileage = newMileage
                curPt = newPt

            if (segIntPoly[i]['intMileageList'][k][3] == 'enter'):
                polyInside.append(segIntPoly[i]['intMileageList'][k][1])

            elif (segIntPoly[i]['intMileageList'][k][3] == 'leave'):
                polyInside.remove(segIntPoly[i]['intMileageList'][k][1])

    return segSet

def polyPath2Mileage(repSeq: list, path: list[pt], nodes: dict, polyFieldName: str = 'neighbor'):
    # RECORD: 20240724 这个函数真是麻烦得要死

    # NOTE: 根据p2pPath，得到足够多的子问题信息以生成cut
    # NOTE: repSeq的第一项和最后一项应该是startLoc和endLoc对应的ID

    # Step 1: 先把转折点找出来
    # NOTE: error取值不能太小，因为用的是30边形拟合的圆 + poly2Poly，导致误差其实还蛮大的
    # NOTE: 这个有问题，但是不好调...error太大了不行，太小了也不行
    degenPath = seqRemoveDegen(path, error = 0.015)

    # Step 2: 按照转折点，找到路径与每个poly的合法相交部分
    # 需要返回的字典列表
    mileage = []
    turnPtMileageAcc = 0

    # 针对聚合后的aggNodeList上的点逐个分情况讨论，
    # NOTE: 这里
    for i in range(len(degenPath['removedFlag'])):
        # 接下来分情况讨论：
        # NOTE: 单独/重合 => 在该坐标上有一个解还是多个解
        # NOTE: 转折/穿透 => path访问该poly的时候是相切还是相交
        # Case 1: 单独转折点 - len(aggNode) == 1 and removeFlag == False
        # Case 2: 重合转折点 - len(aggNode) >  1 and removeFlag == False
        # Case 3: 单独穿透点 - len(aggNode) == 1 and removeFlag == True
        # Case 4: 重合穿透点 - len(aggNode) >  1 and removeFlag == True
        aggNode = degenPath['aggNodeList'][i]

        # 转折点的情形
        # 此时涉及的线段包括两条，[lastTurnPt, curTurnPt]和[curTurnPt, nextTurnPt]
        if (degenPath['removedFlag'][i] == False):
            # 找到上个转折点
            lastTurnPt = None
            if (i >= 1):
                # 向前找到前面最后一个False的值的位置，对应的index在aggNodeList里
                # NOTE: 第一个值肯定为False，所以一定能找到               
                for lastTurnIdx in range(i):
                    if (degenPath['removedFlag'][i - 1 - lastTurnIdx] == False):
                        lastTurnPt = path[degenPath['aggNodeList'][i - 1 - lastTurnIdx][0]]
                        break

            # 当前的转折点和进一步累加的里程
            curTurnPt = path[aggNode[0]]
            distAdd2Acc = 0
            if (lastTurnPt != None):
                distAdd2Acc = distEuclideanXY(lastTurnPt, curTurnPt)['dist']

            # 找到下个转折点
            nextTurnPt = None
            if (i <= len(degenPath['removedFlag']) - 2):
                # 找到后续的removedFlag里第一个为False的值，对应的index在aggNodeList里
                # NOTE: 最后一个值肯定为False，所以一定能找到
                for nextTurnIdx in range(i, len(degenPath['removedFlag'])):
                    if (degenPath['removedFlag'][nextTurnIdx] == False):
                        nextTurnPt = path[degenPath['aggNodeList'][nextTurnIdx][0]]
                        break
            
            # Case 1: 单独转折点
            # NOTE: 生成一个Touch的Single点，该点的mileage为lastTurnPt到该点的距离
            # NOTE: 最简单情形
            if (len(aggNode) == 1):
                mileage.append({
                    'polyID': repSeq[aggNode[0]],
                    'type': 'Touch',
                    'intersect': curTurnPt,
                    'mileage': turnPtMileageAcc + distAdd2Acc
                })
            
            # Case 2: 重合转折点
            # NOTE: 这种情况下，要区分每个重合在此处的转折点与neighborhood是相交还是相切
            # NOTE: 对于相交的，返回mileage的范围，对于相切的，则视作转折点
            # NOTE: 注意，至少一个是转折点
            # NOTE: 找到转折点对应的nodeID，在那之前和之后的分别属于上一段和下一段
            elif (len(aggNode) > 1):
                # 基本流程：
                # Step 1: 找到aggNode里相切点的连续子集，也就是第一个相切点到最后一个相切点之间
                # Step 2: 相切连续子集的前面与上一段线段相交
                # Step 3: 相切连续子集的后面与后一段线段相交
                
                # 先区分出是相切还是相交
                intType = []
                for k in range(len(aggNode)):
                    # 如果是出发点或返回点，则必须是相切点
                    if (aggNode[k] == repSeq[0] or aggNode[k] == repSeq[-1]):
                        intType.append("Point")
                    # 判断转折点是在poly的内部还是边缘
                    # NOTE: 如果转折点是在
                    else:                        
                        insideInteriorFlag = isPtInPoly(
                            pt = curTurnPt,
                            poly = nodes[repSeq[aggNode[k]]][polyFieldName],
                            interiorOnly = True)
                        if (insideInteriorFlag):
                            insideInteriorFlag = not (distPt2Seq(
                                pt = curTurnPt,
                                seq = nodes[repSeq[aggNode[k]]][polyFieldName],
                                closedFlag = True) <= 0.01)
                        if (insideInteriorFlag):
                            intType.append("Seg")
                        else:
                            intType.append("Point")
                # print("intType", intType)

                # 先区分出相切前和相切后
                # NOTE: 先从前往后，再从后往前，分别标记出Pre和Post
                intPhase = ["Tang" for i in range(len(aggNode))]
                for k in range(len(intType)):
                    if (intType[k] == "Seg"):
                        intPhase[k] = "Prev"
                    else:
                        break
                for k in range(len(intType)):
                    if (intType[len(intType) - 1 - k] == "Seg"):
                        intPhase[len(intType) - 1 - k] = "Post"
                    else:
                        break
                # print("intPhase", intPhase)

                # 按照分区确定mileage
                prevInt = [] # 与上一段相交
                tangInt = [] # 相切
                postInt = [] # 与下一段相交
                for k in range(len(aggNode)):
                    if (intPhase[k] == "Prev"):
                        prevInt.append(aggNode[k])
                    elif (intPhase[k] == "Tang"):
                        tangInt.append(aggNode[k])
                    elif (intPhase[k] == "Post"):
                        postInt.append(aggNode[k])
                # print(prevInt, tangInt, postInt)

                # 与上一段相交的部分
                for k in range(len(prevInt)):
                    neiIntSeg = intSeg2Poly(
                        [lastTurnPt, curTurnPt], 
                        nodes[repSeq[prevInt[k]]][polyFieldName])
                    # NOTE: 这段必须处理成线段，如果是距离很短，那就是数值问题
                    if (neiIntSeg['intersectType'] == 'Segment'):
                        [loc1, loc2] = neiIntSeg['intersect']
                        dist1 = distEuclideanXY(loc1, lastTurnPt)['dist']
                        dist2 = distEuclideanXY(loc2, lastTurnPt)['dist']
                        if (dist1 < dist2):
                            mileage.append({
                                'polyID': repSeq[prevInt[k]],
                                'type': 'Intersect',
                                'intersect': [loc1, loc2],
                                'mileage': [turnPtMileageAcc + dist1, turnPtMileageAcc + distAdd2Acc]
                            })
                        else:
                            mileage.append({
                                'polyID': repSeq[prevInt[k]],
                                'type': 'Intersect',
                                'intersect': [loc2, loc1],
                                'mileage': [turnPtMileageAcc + dist2, turnPtMileageAcc + distAdd2Acc]
                            })
                    elif (neiIntSeg['intersectType'] == 'Point'):
                        loc = neiIntSeg['intersect']
                        mileage.append({
                            'polyID': repSeq[postInt[k]],
                            'type': 'Touch',
                            'intersect': loc,
                            'mileage': turnPtMileageAcc + distAdd2Acc
                            })
                    else:
                        print("ERROR: 说好的相交了呢...")
                        print("Intersection", neiIntSeg)
                        print(repSeq)
                        print(path)
                        print(degenPath)
                        print(prevInt, tangInt, postInt, intType)

                # 相切的部分
                if (len(tangInt) == 1):
                    mileage.append({
                        'polyID': repSeq[tangInt[0]],
                        'type': 'Touch',
                        'intersect': curTurnPt,
                        'mileage': turnPtMileageAcc + distAdd2Acc
                    })
                else:
                    mileage.append({
                        'polyID': [repSeq[tangInt[k]] for k in range(len(tangInt))],
                        'type': 'Touch',
                        'intersect': curTurnPt,
                        'mileage': turnPtMileageAcc + distAdd2Acc
                    })

                # 与下一段相交的部分
                for k in range(len(postInt)):
                    neiIntSeg = intSeg2Poly(
                        [curTurnPt, nextTurnPt], 
                        nodes[repSeq[postInt[k]]][polyFieldName])
                    # NOTE: 这段必须处理成线段，如果是距离很短，那就是数值问题
                    if (neiIntSeg['intersectType'] == 'Segment'):
                        [loc1, loc2] = neiIntSeg['intersect']
                        dist1 = distEuclideanXY(loc1, curTurnPt)['dist']
                        dist2 = distEuclideanXY(loc2, curTurnPt)['dist']
                        if (dist1 < dist2):
                            mileage.append({
                                'polyID': repSeq[postInt[k]],
                                'type': 'Intersect',
                                'intersect': [loc1, loc2],
                                'mileage': [turnPtMileageAcc + distAdd2Acc, 
                                    turnPtMileageAcc + distAdd2Acc + dist2]
                                })
                        else:
                            mileage.append({
                                'polyID': repSeq[postInt[k]],
                                'type': 'Intersect',
                                'intersect': [loc2, loc1],
                                'mileage': [turnPtMileageAcc + distAdd2Acc, 
                                    turnPtMileageAcc + distAdd2Acc + dist1]
                                })
                    elif (neiIntSeg['intersectType'] == 'Point'):
                        loc = neiIntSeg['intersect']
                        mileage.append({
                            'polyID': repSeq[postInt[k]],
                            'type': 'Touch',
                            'intersect': loc,
                            'mileage': turnPtMileageAcc + distAdd2Acc
                            })
                    else:
                        print("ERROR: 说好的相交了呢...")
                        print("Intersection", neiIntSeg)
                        print(repSeq)
                        print(path)
                        print(degenPath)
                        print(prevInt, tangInt, postInt, intType)

            # 更新一下累积到lastTurnPt的累积mileage
            turnPtMileageAcc += distAdd2Acc
        
        # 穿透点的情形
        elif (degenPath['removedFlag'][i] == True):
            # NOTE: 穿透点一般而言将找到一条相交的线段，也有可能极小情况下实际上是相切的，那种情况得处理成转折点

            lastTurnPt = degenPath['locatedSeg'][i][0]

            # Case 3: 单独穿透点
            if (len(aggNode) == 1):
                # 穿透点所在的线段
                inclNode = aggNode[0]
                neiIntSeg = intSeg2Poly(
                    degenPath['locatedSeg'][i], 
                    nodes[repSeq[inclNode]][polyFieldName])
                
                # Case 3.1: 最正常的情况，path穿过poly，相交为一个线段
                if (neiIntSeg['intersectType'] == 'Segment'):
                    intSeg = neiIntSeg['intersect']
                    loc1 = intSeg[0]
                    loc2 = intSeg[1]
                    dist1 = distEuclideanXY(loc1, lastTurnPt)['dist']
                    dist2 = distEuclideanXY(loc2, lastTurnPt)['dist']
                    if (dist1 < dist2):
                        mileage.append({
                            'polyID': repSeq[inclNode],
                            'type': 'Intersect',
                            'intersect': [loc1, loc2],
                            'mileage': [turnPtMileageAcc + dist1, turnPtMileageAcc + dist2]
                        })
                    else:
                        mileage.append({
                            'polyID': repSeq[inclNode],
                            'type': 'Intersect',
                            'intersect': [loc2, loc1],
                            'mileage': [turnPtMileageAcc + dist2, turnPtMileageAcc + dist1]
                        })

                # Case 3.2: 特殊情况下，如果穿透点+单独点为neighbor的切点，此时把穿透点处理成转折点
                elif (neiIntSeg['intersectType'] == 'Point'):
                    tangLoc = neiIntSeg['intersect']
                    dist = distEuclideanXY(tangLoc, lastTurnPt)['dist']
                    mileage.append({
                        'polyID': repSeq[inclNode],
                        'type': 'Touch',
                        'intersect': tangLoc,
                        'mileage': turnPtMileageAcc + dist,
                        'info': 'Tangent'
                    })
                    lastTurnPt = tangLoc
                    turnPtMileageAcc += dist
                
                # Case 3.No: 正常情况下这个分支不应该存在，但是实际上因为精度的问题就是会出现
                # NOTE: 处理成相切点，相切处为线段上离poly最近点
                elif (neiIntSeg['intersectType'] == None):
                    tangLoc = nearestPtLine2Poly(
                        degenPath['locatedSeg'][i], 
                        nodes[repSeq[inclNode]][polyFieldName])['ptOnLine']
                    dist = distEuclideanXY(tangLoc, lastTurnPt)['dist']
                    mileage.append({
                        'polyID': repSeq[inclNode],
                        'type': 'Touch',
                        'intersect': tangLoc,
                        'mileage': turnPtMileageAcc + dist,
                        'info': 'TangentError'
                    })
                    lastTurnPt = tangLoc
                    turnPtMileageAcc += dist
                    warnings.warn("WARNING: Numerical issue when calculating mileage.")
            
            # Case 4: 重合穿透点
            # FIXME: 这部分代码要好好走查一下
            # NOTE: 这个情况很复杂，如果存在至少一个相切的情况，那么该点实际上是转折点，且是多重转折点
            # NOTE: 需要挨个确认是否是相切点，如果是相切点，按相切点处理（聚合在一起），如果不是相切点，依次计算mileage
            else:
                tangFlag = False
                tangLoc = None
                neiIntSet = []

                for k in degenPath['aggNodeList'][i]:
                    # 穿透点所在的线段
                    neiInt = intSeg2Poly(
                        degenPath['locatedSeg'][i], 
                        nodes[repSeq[k]][polyFieldName])
                    neiIntSet.append((repSeq[k], neiInt))
                    if (neiInt['intersectType'] == 'Point'):
                        tangFlag = True
                        tangLoc = neiInt['intersect']
                    elif (neiInt['intersectType'] == None):
                        tangFlag = True
                        tangLoc = nearestPtLine2Poly(
                            degenPath['locatedSeg'][i], 
                            nodes[repSeq[k]][polyFieldName])['ptOnLine']

                if (tangFlag == False):
                    for intSeg in neiIntSet:
                        loc1 = intSeg[1]['intersect'][0]
                        loc2 = intSeg[1]['intersect'][1]
                        dist1 = distEuclideanXY(loc1, lastTurnPt)['dist']
                        dist2 = distEuclideanXY(loc2, lastTurnPt)['dist']
                        if (dist1 < dist2):
                            mileage.append({
                                'polyID': intSeg[0],
                                'type': 'Intersect',
                                'intersect': [loc1, loc2],
                                'mileage': [turnPtMileageAcc + dist1, turnPtMileageAcc + dist2]
                            })
                        else:
                            mileage.append({
                                'polyID': intSeg[0],
                                'type': 'Intersect',
                                'intersect': [loc2, loc1],
                                'mileage': [turnPtMileageAcc + dist2, turnPtMileageAcc + dist1]
                            })

                # Case 4.1: 特殊情况下，如果穿透点+重合点为neighbor的切点，此时把穿透点处理成转折点
                # NOTE: 这个目前很罕见，但是应该也可以生成对应的算例
                else:
                    dist = distEuclideanXY(tangLoc, lastTurnPt)['dist']
                    mileage.append({
                        'polyID': [repSeq[k] for k in degenPath['aggNodeList'][i]],
                        'type': 'Touch',
                        'intersect': tangLoc,
                        'mileage': turnPtMileageAcc + dist,
                        'info': 'Tangent'
                    })
                    lastTurnPt = tangLoc
                    turnPtMileageAcc += dist

    return mileage

def serviceTimeCETSP(
    seq: list[pt], 
    polygons: dict, 
    maxSpeed: float, 
    demand: float,
    polyFieldName: str = 'polygon',
    mode: str = "AccuOverlap"):

    seq = seqRemoveDegen(
        seq = seq, 
        error = 0.03)['newSeq']

    # 每个Seg所属的polygon的集合: segSet[seg]['belong']
    segSet = segSetSeq2Poly(
        seq = seq, 
        polygons = polygons, 
        polyFieldName = 'neighbor',
        seqDegenFlag = True)

    # 补充计算下每段的长度
    for i in range(len(segSet)):
        if (segSet[i]['type'] == 'Segment'):
            segSet[i]['length'] = distEuclideanXY(segSet[i]['shape'][0], segSet[i]['shape'][1])['dist']
        else:
            segSet[i]['length'] = 0

    # 配置demand service time
    for p in polygons:
        if ('demandST' not in polygons[p]):
            polygons[p]['demandST'] = demand
        polygons[p]['serviceTime'] = 0

    # 每个polygon内的seg的集合
    for p in polygons:
        polygons[p]['segSet'] = []
    for i in range(len(segSet)):
        for p in segSet[i]['belong']:
            polygons[p]['segSet'].append(i)

    ST = grb.Model("ServiceTime")
    ST.setParam('LogToConsole', 0)
    ST.setParam('OutputFlag', 0)
    ST.modelSense = grb.GRB.MINIMIZE

    # Decision variables ======================================================
    # csp - Binary, 1 if segment s is to serve polygon p
    csp = {}
    hasGRB = []
    for s in range(len(segSet)):
        for p in polygons:
            if (p in segSet[s]['belong']):
                csp[s, p] = ST.addVar(vtype = grb.GRB.BINARY)
                hasGRB.append((s, p))
            else:
                csp[s, p] = 0

    # Travel time on each segment s
    ts = {}
    for s in range(len(segSet)):
        ts[s] = ST.addVar(vtype = grb.GRB.CONTINUOUS, obj = 1)

    # Aux variables to ensure continuous visit
    hL = {}
    hR = {}
    for s in range(len(segSet)):
        for p in polygons:
            hL[s, p] = ST.addVar(vtype = grb.GRB.BINARY)
            hR[s, p] = ST.addVar(vtype = grb.GRB.BINARY)

    # Constraints =============================================================
    # 每个seg的时间至少是ts[s]
    for s in range(len(segSet)):
        ST.addConstr(ts[s] >= segSet[s]['length'] / maxSpeed)
    # 每个poly的服务时长需要满足
    for p in polygons:
        ST.addConstr(grb.quicksum(csp[s, p] * ts[s] for s in polygons[p]['segSet']) >= polygons[p]['demandST'])

    # Accu v.s. Cont ==========================================================
    if (mode == "ContOverlap" or mode == "ContNonOverlap"):
        for p in polygons:
            for s in range(len(segSet)):
                ST.addConstr(hL[s, p] + hR[s, p] <= 1)
            for s in range(len(segSet)):
                for k in range(len(segSet)):
                    if (k < s):
                        ST.addConstr(hL[s, p] >= csp[k, p] - len(segSet) * csp[s, p])
                    elif (k > s):
                        ST.addConstr(hR[s, p] >= csp[k, p] - len(segSet) * csp[s, p])
    # Overlap v.s. nonOverlap =================================================
    if (mode == "AccuNonOverlap" or mode == "ContNonOverlap"):
        for s in range(len(segSet)):
            ST.addConstr(grb.quicksum(csp[s, p] for p in segSet[s]['belong']) <= 1)

    ST.update()
    ST.optimize()

    # Rebuild solution ========================================================
    ofv = None
    if (ST.status == grb.GRB.status.OPTIMAL):
        ofv = ST.getObjective().getValue()

        # 计算每个点上的时间
        timedSeq = []
        for s in range(len(segSet)):
            segSet[s]['duration'] = ts[s].x

        # 每段的Note
        noteSeq = []
        for s in range(len(segSet)):
            objs = []
            for (k, p) in csp:
                if (k == s and (k, p) in hasGRB and csp[k, p].x > 0.9):
                    objs.append(p)
            noteSeq.append(
                "Obj: " + list2String(objs) + "\n" + 
                "Dur: " + str(round(ts[s].x, 2)) + "\n" + 
                "Len: " + str(round(segSet[s]['length'], 2)))

        # 出入每个边界点的时间戳
        accTravelTime = 0
        for i in range(len(segSet)):
            if (type(segSet[i]['mileage']) == list):
                timedSeq.append((segSet[i]['shape'][0], accTravelTime))
                accTravelTime += segSet[i]['duration']
            else:
                timedSeq.append((segSet[i]['shape'], accTravelTime))
                accTravelTime += segSet[i]['duration']
        if (type(segSet[len(segSet)-1]['mileage']) == list):
            timedSeq.append((segSet[len(segSet)-1]['shape'][1], accTravelTime))
        else:
            timedSeq.append((segSet[len(segSet)-1]['shape'], accTravelTime))

        # 记录每个polygon累计服务时间
        for (s, p) in csp:
            if ((s, p) in hasGRB and csp[s, p].x > 0.9):
                polygons[p]['serviceTime'] += segSet[s]['duration']

    return {
        'ofv': ofv,
        'timedSeq': timedSeq,
        'note': noteSeq
    }
