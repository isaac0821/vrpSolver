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

def segSetSeq2Poly(seq: list, polygons: dict, polyFieldName: str = 'polygon', tangPts: dict = None, polyInt: dict = None):
    # 简化线段,去掉穿越点
    seq = seqRemoveDegen(seq, error = 0.03)['newSeq']

    # NOTE: 仅是加速用
    if (polyInt == None):
        polyInt = {}
        for i in polygons:
            for j in polygons:
                if (i < j and isPolyIntPoly(
                    poly1 = polygons[i][polyFieldName], 
                    poly2 = polygons[j][polyFieldName])):
                    polyInt[i, j] = True
                    polyInt[j, i] = True

    # Step 0: tangPts =========================================================
    # NOTE: 所有的转折点必须在一个多边形的边缘上
    if (tangPts == None):


    # Step 1: For each polygon, get intersections =============================
    for p in polygons:
        # 用来记录可以分配时间的经停点/经过线段
        polygons[p]['intersect'] = []
        # 由于计算精度的问题,加上简化线段时损失的精度,可能会导致path和polygon没有交点,属于提前需要把seq预设的represent point记录上来
        if (tangPts != None and p in tangPts):
            polygons[p]['tmpPoint'] = [{
                'status': 'Cross',
                'intersectType': 'Point',
                'intersect': tangPts[p],
                'interiorFlag':  True # 其实可能是在boundary上,但是不重要
            }]
        else:
            polygons[p]['tmpPoint'] = []
        for i in range(len(seq) - 1):
            segInt = intSeg2Poly([seq[i], seq[i + 1]], polygons[p][polyFieldName])
            # 对于Nonconvex的poly，这部分不应该出现，但是实测下来Shapely还是不够可靠
            if (type(segInt) == list):
                for s in segInt:
                    if (s['intersectType'] == 'Point'):
                        polygons[p]['tmpPoint'].append(s)
                    else:
                        length = distEuclideanXY(s['intersect'][0], s['intersect'][1])['dist']
                        if (length > 0.01):
                            polygons[p]['intersect'].append(s)
                        else:
                            k = {
                                'status': 'Cross',
                                'intersectType': 'Point',
                                'intersect': s['intersect'][0],
                                'interiorFlag': False
                            }
                            polygons[p]['tmpPoint'].append(k)
            # 对于相交的情形，先存起来，再判断相交情况，特别是点，对于所有的点，判断是不是已经存在
            elif (segInt['status'] == 'Cross'):
                if (segInt['intersectType'] == 'Point'):
                    polygons[p]['tmpPoint'].append(segInt)
                else:
                    length = distEuclideanXY(segInt['intersect'][0], segInt['intersect'][1])['dist']
                    if (length > 0.01):
                        polygons[p]['intersect'].append(segInt)
                    else:
                        k = {
                            'status': 'Cross',
                            'intersectType': 'Point',
                            'intersect': segInt['intersect'][0],
                            'interiorFlag': False
                        }
                        polygons[p]['tmpPoint'].append(k)

    for i in polygons:
        # 对于所有的相交项，判断重合性，主要判断Point是不是重合了已经
        for pt in polygons[i]['tmpPoint']:
            existFlag = False
            for obj in polygons[i]['intersect']:
                if (obj['intersectType'] == 'Point'):
                    if (is2PtsSame(pt1 = pt['intersect'], pt2 = obj['intersect'], error = 0.03)):
                        existFlag = True
                        break
                else:
                    if (isPtOnSeg(pt = pt['intersect'], seg = obj['intersect'], error = 0.03)):
                        existFlag = True
                        break
            if (not existFlag):
                polygons[i]['intersect'].append(pt)
        polygons[i].pop('tmpPoint')

    # Step 2: Initialize ======================================================
    # 把所有polygon内的seg/pt原封不动的放进来
    segSet = {}
    segIDAcc = 0
    for i in polygons:
        for s in polygons[i]['intersect']:
            segSet[segIDAcc] = {
                'shape': s['intersect'],
                'type': s['intersectType'],
                'belong': [i],
                'split': False
            }
            segIDAcc += 1

    # Step 3: Split the segments ==============================================
    canSplitFlag = True
    while (canSplitFlag):
        canSplitFlag = False

        newSegSet = []
        for i in segSet:
            canSplitBetweenIJFlag = False  
            for j in segSet:
                # 两个seg/pt看看能不能相交
                if (i < j and segSet[i]['split'] == False and segSet[j]['split'] == False):
                    # 分别看看i和j所属的poly有没有可能相交,如果有可能,才尝试,不可能相交就没必要了
                    trySplitFlag = False
                    for iBelong in segSet[i]['belong']:
                        for jBelong in segSet[j]['belong']:
                            if ((iBelong, jBelong) in polyInt and polyInt[iBelong, jBelong] == True):
                                trySplitFlag = True

                    if (trySplitFlag):
                        # 如果i和j都是segment
                        if (segSet[i]['type'] == 'Segment' and segSet[j]['type'] == 'Segment'):
                            s2s = intSeg2Seg(segSet[i]['shape'], segSet[j]['shape'])
                            if (s2s['status'] == 'Collinear' and s2s['interiorFlag'] == True and s2s['intersectType'] == 'Segment'):
                                segOverlapSeg = splitOverlapSegs(
                                    seg1 = segSet[i]['shape'],
                                    seg2 = segSet[j]['shape'],
                                    belong1 = [k for k in segSet[i]['belong']],
                                    belong2 = [k for k in segSet[j]['belong']])
                                if (segOverlapSeg != None):
                                    for k in range(len(segOverlapSeg)):
                                        segOverlapSeg[k]['type'] = 'Segment'
                                        segOverlapSeg[k]['split'] = False                                
                                newSegSet.extend(segOverlapSeg)
                                
                                segSet[i]['split'] = True
                                segSet[j]['split'] = True
                                canSplitBetweenIJFlag = True
                                break

                        # 如果i是segment, j是point
                        elif (segSet[i]['type'] == 'Segment' and segSet[j]['type'] == 'Point'):
                            if (isPtOnSeg(pt = segSet[j]['shape'], seg = segSet[i]['shape'], interiorOnly = True)):
                                # print("Part ", i, segSet[i], "Part ", j, segSet[j], "S2P")
                                newSegSet.append({
                                    'shape': [segSet[i]['shape'][0], segSet[j]['shape']],
                                    'type': 'Segment',
                                    'belong': [k for k in segSet[i]['belong']],
                                    'split': False
                                })
                                newSegSet.append({
                                    'shape': segSet[j]['shape'],
                                    'type': 'Point',
                                    'belong': [k for k in segSet[j]['belong']],
                                    'split': False
                                })
                                newSegSet.append({
                                    'shape': [segSet[j]['shape'], segSet[i]['shape'][1]],
                                    'type': 'Segment',
                                    'belong': [k for k in segSet[i]['belong']],
                                    'split': False
                                })                                                                
                                segSet[i]['split'] = True
                                segSet[j]['split'] = True
                                canSplitBetweenIJFlag = True
                                break

                        # 如果i是point, j是segment
                        elif (segSet[i]['type'] == 'Point' and segSet[j]['type'] == 'Segment'):
                            if (isPtOnSeg(pt = segSet[i]['shape'], seg = segSet[j]['shape'], interiorOnly = True)):
                                # print("Part ", i, segSet[i], "Part ", j, segSet[j], "P2S")
                                newSegSet.append({
                                    'shape': [segSet[j]['shape'][0], segSet[i]['shape']],
                                    'type': 'Segment',
                                    'belong': [k for k in segSet[j]['belong']],
                                    'split': False
                                })
                                newSegSet.append({
                                    'shape': segSet[i]['shape'],
                                    'type': 'Point',
                                    'belong': [k for k in segSet[i]['belong']],
                                    'split': False
                                })
                                newSegSet.append({
                                    'shape': [segSet[i]['shape'], segSet[j]['shape'][1]],
                                    'type': 'Segment',
                                    'belong': [k for k in segSet[j]['belong']],
                                    'split': False
                                })                                
                                segSet[i]['split'] = True
                                segSet[j]['split'] = True
                                canSplitBetweenIJFlag = True
                                break

                        # 如果i和j都是point
                        elif (segSet[i]['type'] == 'Point' and segSet[j]['type'] == 'Point'):
                            # 若两点重合, 则合并两者的belong, 否则不需要进行处理
                            if (is2PtsSame(pt1 = segSet[i]['shape'], pt2 = segSet[j]['shape'])):
                                # print("Part ", i, segSet[i], "Part ", j, segSet[j], "P2P")
                                jointBelong = [k for k in segSet[i]['belong']]
                                jointBelong.extend([k for k in segSet[j]['belong'] if k not in segSet[i]['belong']])
                                newSegSet.append({
                                    'shape': segSet[i]['shape'],
                                    'type': 'Point',
                                    'belong': jointBelong,
                                    'split': False
                                })                                
                                segSet[i]['split'] = True
                                segSet[j]['split'] = True
                                canSplitBetweenIJFlag = True
                                break
            if (canSplitBetweenIJFlag):                
                break

        if (len(newSegSet) > 0):
            # print("newSegSet", newSegSet)
            # print(hyphenStr())
            canSplitFlag = True

        for k in newSegSet:
            segSet[segIDAcc] = k
            segIDAcc += 1
        
        segSet = {k: v for k, v in segSet.items() if (v['split'] == False)}
    
    # NOTE: 到这里是对的
    # return segSet

    # Step 4: Sort segSet by mileage ==========================================
    # NOTE: 每个segSet内的元素应该是两两不相交
    heapSegSet = []
    sortedSegSet = {}
    # NOTE: 先得到每个segSet内元素的mileage
    for i in segSet:
        if (segSet[i]['type'] == 'Point'):
            mileage = mileagePt(seq, segSet[i]['shape'], error = 0.03)
            # 找不到mileage是真没办法...
            # NOTE: 最差的办法就是做投影然后求mileage
            if (mileage == None):
                raise RuntimeError("ERROR: Precision error")
            # 这个好麻烦...
            # NOTE: 如果仅仅是[0, end]的话还好说, 属于不应该出现在这里
            elif (type(mileage) == list):
                raise RuntimeError("ERROR: Depot is on the edge of one of polygons")
            else:
                pass
            # print(i,  mileage)
            segSet[i]['mileage'] = mileage
            heapq.heappush(heapSegSet, (mileage, segSet[i]))
        else:
            mileage1 = mileagePt(seq, segSet[i]['shape'][0], error = 0.03)
            mileage2 = mileagePt(seq, segSet[i]['shape'][1], error = 0.03)
            # NOTE: 真的好讨厌啊...
            if (mileage1 == None):
                raise RuntimeError("ERROR: Precision error")
            elif (type(mileage1) == list):
                # 这边专门针对[0, end]这种特殊情形处理
                if (type(mileage2) != list and mileage2 < distEuclideanXY(seq[0], seq[1])['dist']):
                    mileage1 = min(mileage1)
                elif (type(mileage2) != list and mileage2 > distEuclideanXY(seq[0], seq[1])['dist']):
                    mileage1 = max(mileage1)
            else:
                pass
            if (mileage2 == None):
                raise RuntimeError("ERROR: Precision error")
            elif (type(mileage2) == list):
                # 这边专门针对[0, end]这种特殊情形处理
                if (type(mileage1) != list and mileage1 < distEuclideanXY(seq[0], seq[1])['dist']):
                    mileage2 = min(mileage2)
                elif (type(mileage1) != list and mileage1 > distEuclideanXY(seq[0], seq[1])['dist']):
                    mileage2 = max(mileage2)
            else:
                pass
            # print(i, mileage1, mileage2)
            segSet[i]['mileage'] = [min(mileage1, mileage2), max(mileage1, mileage2)]
            heapq.heappush(heapSegSet, ((mileage1 + mileage2) / 2, segSet[i]))

    curMileage = 0
    while(len(heapSegSet) > 0):
        nextSeg = heapq.heappop(heapSegSet)[1]
        nextSegStart = nextSeg['mileage'] if type(nextSeg['mileage']) != list else nextSeg['mileage'][0]
        nextSegEnd = nextSeg['mileage'] if type(nextSeg['mileage']) != list else nextSeg['mileage'][1]

        if (abs(curMileage - nextSegStart) <= CONST_EPSILON):
            sortedSegSet[len(sortedSegSet)] = nextSeg
        else:
            sortedSegSet[len(sortedSegSet)] = {
                'shape': [ptInSeqMileage(seq, curMileage), ptInSeqMileage(seq, nextSegStart)],
                'type': 'Segment',
                'belong': [],
                'mileage': [curMileage, nextSegStart],
                'split': False
            }
            sortedSegSet[len(sortedSegSet)] = nextSeg
        curMileage = nextSegEnd

    return sortedSegSet
