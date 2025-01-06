import gurobipy as grb
import networkx as nx

from .geometry import *
from .common import *
from .msg import *
from .ds import *

# Distance between timed-objects ==============================================
# [Constructing]
def distTimedSeg2TimedSeg(timedSeg1: list[tuple[pt, float]], timedSeg2: list[tuple[pt, float]], startTime: float = 0):
    # 先区分出timedSeg哪段在前，哪段在后
    if (timedSeg1[0][1] < timedSeg2[0][1]):
        (sx1, sy1) = timedSeg1[0][0]
        (ex1, ey1) = timedSeg1[1][0]
        st1 = timedSeg1[0][1]
        et1 = timedSeg1[1][1]
        (sx2, sy2) = timedSeg2[0][0]
        (ex2, ey2) = timedSeg2[1][0]
        st2 = timedSeg2[0][1]
        et2 = timedSeg2[1][1]
    else:
        (sx1, sy1) = timedSeg2[0][0]
        (ex1, ey1) = timedSeg2[1][0]
        st1 = timedSeg2[0][1]
        et1 = timedSeg2[1][1]
        (sx2, sy2) = timedSeg1[0][0]
        (ex2, ey2) = timedSeg1[1][0]
        st2 = timedSeg1[0][1]
        et2 = timedSeg1[1][1]

    # 分三类情况讨论求解
    minDist = None
    if (st1 <= et1 <= st2 <= et2):
        # t \in [st1, et1] 
        # 1 move, 2 stay
        pt2 = (sx2, sy2)
        seg1 = [(sx1, sy1), (ex1, ey1)]
        minDist1 = distPt2Seg(pt2, seg1)

        # t \in [et1, st2]
        # 1 stay, 2 stay
        minDist2 = distEuclideanXY((ex1, ey1), (ex2, ey2))

        # t \in [st2, et2]
        # 1 stay, 2 move
        pt1 = (ex1, ey1)
        seg2 = [(sx2, sy2), (ex2, ey2)]
        minDist3 = distPt2Seg(pt1, seg2)

        return min([minDist1, minDist2, minDist3])

    elif (st1 <= st2 <= et1 <= et2):
        # t \in [st1, st2]
        # 1 move, 2 stay
        seg1AtS2 = snapInTimedSeq([((sx1, sy1), st1), ((ex1, ey1), et1)], st2)['loc']
        seg1 = [(sx1, ex1), seg1AtS2]
        pt2 = (sx2, sy2)
        minDist1 = distPt2Seg(pt2, seg1)

        # t \in [st2, et1]
        seg2AtE1 = snapInTimedSeq([((sx2, sy2), st2), ((ex2, ey2), et2)], et1)['loc']
        seg1 = [seg1AtS2, (ex1, ey1)]
        seg2 = None
        # 1 move, 2 move

        # t \in [et1, et2]
        # 1 stay, 2 move

        return min([minDist1, minDist2, minDist3])

    elif (st1 <= st2 <= et2 <= et1):
        # t \in [st1, st2]
        # 1 move, 2 stay

        # t \in [st2, et2]
        # 1 move, 2 move

        # t \in [et2, et1]
        # 1 move, 2 stay

        return min([minDist1, minDist2, minDist3])

def timedSeg2Vec(timedSeg):
    dt = timedSeg[1][1] - timedSeg[0][1]
    dx = timedSeg[1][0][0] - timedSeg[0][0][0]
    dy = timedSeg[1][0][1] - timedSeg[0][0][1]
    l = math.sqrt(dx ** 2 + dy ** 2)
    lx = dx * dt / l 
    ly = dy * dt / l
    return timedSeg[0][0], (lx, ly)

def distVec2Vec(pt1, vec1, pt2, vec2, earliest:None|float = None, latest:None|float = None):
    # NOTE: 俩点同时都在动
    # NOTE: 这段代码目前先用gurobi来做，之后要换成O(1)代入公式
    x1, y1 = pt1
    x2, y2 = pt2
    vx1, vy1 = vec1
    vx2, vy2 = vec2

    model = grb.Model("SOCP")
    model.setParam("OutputFlag", 0)

    # Decision variables ======================================================
    d = model.addVar(vtype=grb.GRB.CONTINUOUS, obj=1)
    t = model.addVar(vtype=grb.GRB.CONTINUOUS)
    dx = model.addVar(vtype=grb.GRB.CONTINUOUS, lb = -float('inf'))
    dy = model.addVar(vtype=grb.GRB.CONTINUOUS, lb = -float('inf'))

    # Constraints =============================================================
    model.addConstr(dx == (x1 - x2) + t * (vx1 - vx2))
    model.addConstr(dy == (y1 - y2) + t * (vy1 - vy2))
    model.addConstr(d ** 2 >= dx ** 2 + dy ** 2)
    if (earliest != None):
        model.addConstr(t >= earliest)
    if (latest != None):
        model.addConstr(t <= latest)

    # Optimize ================================================================
    model.modelSense = grb.GRB.MINIMIZE
    model.optimize()
    minDist = model.getObjective().getValue()

    return minDist

def travelVec2Vec(pt1, vec1, pt2, vec2, speed, earliest:None|float = None, latest:None|float = None):
    # NOTE: 俩点同时都在动
    # NOTE: 这段代码目前先用gurobi来做，之后要换成O(1)代入公式
    x1, y1 = pt1
    x2, y2 = pt2
    vx1, vy1 = vec1
    vx2, vy2 = vec2

    model = grb.Model("SOCP")
    model.setParam("OutputFlag", 0)

    # Decision variables ======================================================
    d = model.addVar(vtype=grb.GRB.CONTINUOUS, obj=1)

    t1 = model.addVar(vtype=grb.GRB.CONTINUOUS)
    t2 = model.addVar(vtype=grb.GRB.CONTINUOUS)

    dx = model.addVar(vtype=grb.GRB.CONTINUOUS, lb = -float('inf'))
    dy = model.addVar(vtype=grb.GRB.CONTINUOUS, lb = -float('inf'))

    # Constraints =============================================================
    model.addConstr(dx == x1 - x2 + t1 * vx1 - t2 * vx2)
    model.addConstr(dy == y1 - y2 + t1 * vy1 - t2 * vy2)
    model.addConstr(d ** 2 >= dx ** 2 + dy ** 2)
    model.addConstr(t2 == t1 + d * (1 / speed))

    if (earliest != None):
        model.addConstr(t1 >= earliest)
    if (latest != None):
        model.addConstr(t2 <= latest)

    # Optimize ================================================================
    model.modelSense = grb.GRB.MINIMIZE
    model.optimize()

    if (model.status == grb.GRB.status.OPTIMAL):
        minDist = model.getObjective().getValue()
        timedSeg = [
            ((x1 + t1.x * vx1, y1 + t1.x * vy1), t1.x), 
            ((x2 + t2.x * vx2, y2 + t2.x * vy2), t2.x)
        ]
        return {
            'minDist': minDist,
            'minTime': (t2.x - t1.x),
            'timedSeg': timedSeg
        }
    else:
        return None

def earliest2Vec(oriPt, speed, movPt, movVec):
    # NOTE: 这俩同时都在动
    oriPtX, oriPtY = oriPt
    movPtX, movPtY = movPt
    movVX, movVY = movVec

    A = (movVX ** 2 + movVY ** 2 - speed ** 2)
    B = 2 * (movVX * (movPtX - oriPtX) + movVY * (movPtY - oriPtY))
    C = (movPtX - oriPtX) ** 2 + (movPtY - oriPtY) ** 2

    # 求根，如有两个解，均输出
    delta = B ** 2 - 4 * A * C

    if (delta < 0):
        return None
    else:
        t1 = (- B + math.sqrt(delta)) / (2 * A)
        t2 = (- B - math.sqrt(delta)) / (2 * A)

        if (max([t1, t2]) < 0):
            return None
        elif (min([t1, t2]) < 0 and max([t1, t2]) >= 0):
            meetTime = max([t1, t2])
            meetPtX = movPtX + meetTime * movVX
            meetPtY = movPtY + meetTime * movVY
            meetPt = (meetPtX, meetPtY)
            return (meetPt, meetTime)
        else: 
            meetTime = min([t1, t2])
            meetPtX = movPtX + meetTime * movVX
            meetPtY = movPtY + meetTime * movVY
            meetPt = (meetPtX, meetPtY)
            return (meetPt, meetTime)
