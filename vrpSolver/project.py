import math
import relation
import geometry

def projSeg2LineXY(
	seg:        "Line segment to be projected to line" = None,
	line:       "Line to project to" = None,
	vec:        "Vector of projecting" = None
	) -> "Given a line segment, project it to a line using a given vector":

	line4SegEnd1 = [seg[0], [seg[0][0] + vec[0], seg[0][1] + vec[1]]]
	line4SegEnd2 = [seg[1], [seg[1][0] + vec[0], seg[1][1] + vec[1]]]

	projPt1 = intLine2Line(line4SegEnd1, line)
	projPt2 = intLine2Line(line4SegEnd2, line)

	return {
		'shadowSegOnLine': [projPt1, projPt2]
	}

def projSeg2SegXY(
	seg1:       "Line segment to be projected to seg2" = None,
	seg2:       "Line segment that is being projected to" = None,
	vec:        "Vector of projecting" = None
	) -> "Given a line segment, project it to another line segment, which might give us \
		1) shadowOnSeg2, the part on the line segment being projected \
		2) projFromSeg1, the part of seg1 that has been projected":

	line4Seg1End1 = [seg1[0], [seg1[0][0] + vec[0], seg1[0][1] + vec[1]]]
	line4Seg1End2 = [seg1[1], [seg1[1][0] + vec[0], seg1[1][1] + vec[1]]]
	line4Seg2End1 = [seg2[0], [seg2[0][0] + vec[0], seg2[0][1] + vec[1]]]
	line4Seg2End2 = [seg2[1], [seg2[1][0] + vec[0], seg2[1][1] + vec[1]]]

	projSeg1End1OnLine2 = intLine2Line(line4Seg1End1, seg2)
	projSeg1End2OnLine2 = intLine2Line(line4Seg1End2, seg2)
	projSeg2End1OnSeg1 = intLine2Line(line4Seg2End1, seg1)
	projSeg2End2OnSeg1 = intLine2Line(line4Seg2End2, seg1)

	seg1End1ProjOnSeg2 = isPtOnSeg(projSeg1End1OnLine2, seg2)
	seg1End2ProjOnSeg2 = isPtOnSeg(projSeg1End2OnLine2, seg2)

	# Case 1: both projPt on seg2
	if (seg1End1ProjOnSeg2 and seg1End2ProjOnSeg2):
		return {
			'shadowOnSeg2': [projSeg1End1OnLine2, projSeg1End2OnLine2],
			'projFromSeg1': seg1
		}

	# Case 2: projPt of end1 of seg1 is on seg2, projPt of end2 of seg1 is not on seg2
	elif (seg1End1ProjOnSeg2 and not seg1End2ProjOnSeg2):
		# Case 2.1: projPt of end1 of seg2 is on seg1
		if (isPtOnSeg(projSeg2End1OnSeg1)):
			return {
				'shadowOnSeg2': [seg2[0], projSeg1End1OnLine2],
				'projFromSeg1': [seg1[0], projSeg2End1OnSeg1]
			}
		else:
			return {
				'shadowOnSeg2': [projSeg1End1OnLine2, seg2[1]],
				'projFromSeg1': [projSeg2End1OnSeg1, seg1[1]]
			}

	# Case 3: projPt of end1 of seg1 is not on seg2, projPt of end1 of seg1 is on seg2

	# Case 4: both projPt not on seg2, seg2 not projecting on seg1

	# Case 5: both projPt not on seg2, seg2 entirely projecting on seg1

