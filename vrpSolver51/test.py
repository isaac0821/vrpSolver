
# [Constructing]
def threeOpt(
    seq:        "A given sequence of vehicle route, assuming this route is feasible", 
    tau:        "Traveling cost matrix", 
    i:          "Index in the sequence, 0 <= i <= len(seq) - 2", 
    j:          "Index in the sequence, 0 <= j <= len(seq) - 2, i.Next() != j",
    k:          "Index in the sequence, 0 <= k <= len(seq) - 2, j.Next() != k",
    oldC:       "Total cost before swapping, if we need to calculate the new cost, input this value, otherwise newC will be None" = None,
    asymFlag:   "True if asymmetric" = False
    ) -> "Reconnect arc between i, i+1 and j, j+1, if it is applicable":

    # Before: ... --> nI -> nINext --> ...abcd... --> nJ     -> nJNext --> ...efgh... --> nK     -> nKNext --> ...
    # After1: ... --> nI -> nINext --> ...abcd... --> nJ     -> nK     --> ...hgfe... --> nJNext -> nKNext --> ...
    # After2: ... --> nI -> nJ     --> ...dcba... --> nINext -> nJNext --> ...efgh... --> nK     -> nKNext --> ...
    # After3: ... --> nI -> nJ     --> ...dcba... --> nINext -> nK     --> ...hgfe... --> nJNext -> nKNext --> ...
    # After4: ... --> nI -> nJNext --> ...efgh... --> nK     -> nINext --> ...abcd... --> nJ     -> nKNext --> ...
    # After5: ... --> nI -> nJNext --> ...efgh... --> nK     -> nJ     --> ...dcba... --> nINext -> nKNext --> ...
    # After6: ... --> nI -> nK     --> ...hgfe... --> nJNext -> nINext --> ...abcd... --> nJ     -> nKNext --> ...
    # After7: ... --> nI -> nK     --> ...hgfe... --> nJNext -> nJ     --> ...dcba... --> nINext -> nKNext --> ...
    N = len(seq)
    nIPrev = seq[iterSeq(N, i, 'prev')]
    nI = seq[i]
    nINext = seq[iterSeq(N, i, 'next')]
    nJPrev = seq[iterSeq(N, j, 'prev')]
    nJ = seq[j]
    nJNext = seq[iterSeq(N, j, 'next')]
    nKPrev = seq[iterSeq(N, k, 'prev')]
    nK = seq[k]
    nKNext = seq[iterSeq(N, k, 'next')]

    revSeq = [seq[len(seq) - i - 1] for i in range(len(seq))]
    costABCD = calSeqCostMatrix(tau, seq, iterSeq(N, i, 'next'), j)
    costEFGH = calSeqCostMatrix(tau, seq, iterSeq(N, j, 'next'), k)
    costDCBA = calSeqCostMatrix(tau, revSeq, N - j, N - iterSeq(N, i, 'next'))
    costHGFE = calSeqCostMatrix(tau, revSeq, N - k, N - iterSeq(N, j, 'next'))

    # Check feasibility
    availDCBA = True
    availHGFE = True
    if (asymFlag):
    	for m in range(N - j, N - iterSeq(N, i, 'next')):
    		if ((revSeq[m], revSeq[m + 1]) not in tau):
    			availDCBA = False
    			break
    	for m in range(N - k, N - iterSeq(N, j, 'next')):
    		if ((revSeq[m], revSeq[m + 1] not in tau)):
    			availHGFE = False
    			break

    # Seq 1
    newSeq1 = []
    newSeq1.extend([seq[m] for m in range(i + 1)])
    newSeq1.extend([seq[j - m] for m in range(j - i)])
    newSeq1.extend([seq[m] for m in range(j + 1, len(seq))])

    newSeq1 = []
    newSeq1.extend([seq[m] for m in range(j + 1)])


    print(seq)
    print(i, j, k)
    print(nI, nJ, nK)
    print(newSeq1)
    print(newSeq2)
    print(newSeq3)
    print(newSeq4)
    print(newSeq5)
    print(newSeq6)
    print(newSeq7)

    return