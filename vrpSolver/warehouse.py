
def createWarehouseLayout(
    shelfCol: "Integer, Number of columns in the warehouse layout" = 5,
    shelfRow: "Integer, Number of rows in warehouse layout" = 10,
    shelfWidth: "Integer, Width of shelf" = 1,
    shelfDeep:  "Integer, Number of deeps of the shelf" = 1,
    aisleWidth: "Integer, Width of the aisles" = 1,
    layoutType: "1) String, 'Traditional' or\
                 2) String, 'CrossAisle'" = 'CrossAisle',
    RSPointType: "1) String, 'BottomMiddle' or\
                 2) String, 'LeftMiddle' or\
                 3) String, 'Corner'" = 'BottomMiddle'
    ) -> "Create the travel network and the grid":

    # Initialize ==============================================================
    dicLayout = {}

    # Metadata ================================================================
    dicLayout['meta'] = {
        'shelfCol': shelfCol,
        'shelfRow': shelfRow,
        'shelfWidth': shelfWidth,
        'shelfDeep': shelfDeep,
        'aisleWidth': aisleWidth,
        'layoutType': layoutType,
        'RSPointType': RSPointType
    }

    # Create grids ============================================================
    dicLayout['grid'] = {}
    numCol = aisleWidth + (shelfWidth * shelfDeep + aisleWidth) * shelfCol
    numRow = aisleWidth * 2 + shelfWidth * shelfRow + aisleWidth * (1 if layoutType == 'CrossAisle' else 0)
    dicLayout['grid']['colRow'] = (numCol, numRow)
    dicLayout['grid']['barriers'] = []
    if (layoutType == 'Traditional'):
        for i in range(1, shelfRow + 1):
            for j in range(shelfCol):
                for d in range(shelfDeep):
                    dicLayout['grid']['barriers'].append((aisleWidth + (shelfDeep * shelfWidth + aisleWidth) * j + shelfWidth * d, i))
    elif (layoutType == 'CrossAisle'):
        for i in range(1, int((shelfRow + 1) / 2) + 1):
            for j in range(shelfCol):
                for d in range(shelfDeep):
                    dicLayout['grid']['barriers'].append((aisleWidth + (shelfDeep * shelfWidth + aisleWidth) * j + shelfWidth * d, i))
        for i in range(int((shelfRow + 1) / 2) + 2, shelfRow + 2):
            for j in range(shelfCol):
                for d in range(shelfDeep):
                    dicLayout['grid']['barriers'].append((aisleWidth + (shelfDeep * shelfWidth + aisleWidth) * j + shelfWidth * d, i))

    # Create RS location ======================================================
    RSX = None
    RSY = None
    if (RSPointType == 'BottomMiddle'):
        RSX = (aisleWidth + (aisleWidth + shelfDeep * shelfWidth) * shelfCol) / 2
        RSY = aisleWidth / 2
        dicLayout['RSGridID'] = (int((aisleWidth + (aisleWidth + shelfDeep * shelfWidth) * shelfCol) / 2), 0)
    elif (RSPointType == 'LeftMiddle'):
        RSX = aisleWidth / 2
        RSY = (aisleWidth * 2 + shelfWidth * shelfRow + (aisleWidth if layoutType == 'CrossAisle' else 0)) / 2
        dicLayout['RSGridID'] = (0, int((aisleWidth * 2 + shelfWidth * shelfRow + (aisleWidth if layoutType == 'CrossAisle' else 0)) / 2))
    elif (RSPointType == 'Corner'):
        RSX = aisleWidth / 2
        RSY = aisleWidth / 2
        dicLayout['RSGridID'] = (0, 0)
    dicLayout['RS'] = [RSX, RSY]

    # Create shelves ==========================================================
    dicLayout['shelf'] = {}
    shelfID = 1
    for col in range(shelfCol):
        for row in range(shelfRow):
            for deep in range(shelfDeep):
                topLeftX = aisleWidth + col * (aisleWidth + shelfDeep * shelfWidth) + deep * shelfWidth
                topTopY = None
                if (layoutType == 'Traditional'):
                    topTopY = aisleWidth + shelfWidth + row * shelfWidth
                elif (layoutType == 'CrossAisle'):
                    if (row < shelfRow / 2):
                        topTopY = aisleWidth + shelfWidth + row * shelfWidth
                    else:
                        topTopY = 2 * aisleWidth + shelfWidth + row * shelfWidth
                poly = [[topLeftX, topTopY], 
                        [topLeftX, topTopY - shelfWidth], 
                        [topLeftX + shelfWidth, topTopY - shelfWidth], 
                        [topLeftX + shelfWidth, topTopY]]
                dicLayout['shelf'][shelfID] = {
                    'poly': poly,
                    'rowID': row,
                    'colID': col,
                    'deepID': deep
                }
                shelfID += 1

    # For the boundary ========================================================
    dicLayout['boundary'] = {}
    dicLayout['boundary']['ploy'] = [
        [0, 0], 
        [aisleWidth + (aisleWidth + shelfDeep * shelfWidth) * shelfCol, 0],
        [aisleWidth + (aisleWidth + shelfDeep * shelfWidth) * shelfCol, aisleWidth * 2 + shelfWidth * shelfRow + (aisleWidth if layoutType == 'CrossAisle' else 0)],
        [0, aisleWidth * 2 + shelfWidth * shelfRow + (aisleWidth if layoutType == 'CrossAisle' else 0)]        
    ]   

    return dicLayout