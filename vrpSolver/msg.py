
# NOTICE ======================================================================
# The file name `msg.py` does NOT stands for Monosodium glutamate. And PLEASE,
# PLEASE, PLEASE do NOT put any MSG into your dishes.

# Warning message - Overwrite inputted field
WARNING_OVERWRITE_STARTTIME         = "WARNING: Field `startTime` has been overwritten"
WARNING_OVERWRITE_ENDTIME           = "WARNING: Field `endTime` has been overwritten"

# Error message - Missing required field
ERROR_MISSING_N                     = "ERROR: Missing required field `N`, `N` is an integer indicating number of vertices"
ERROR_MISSING_NODES                 =("ERROR: Missing required field `nodes`. The format of `nodes` is \n"
                                    + "    { \n "
                                    + "        nodeID1: {'loc': (x, y)}, \n"
                                    + "        nodeID2: {'loc': (x, y)}, \n"
                                    + "        ... \n"
                                    + "    }")
ERROR_MISSING_PLAINNODES            = "ERROR: Missing required field `plainNodes`, which can be generated using `rndPlainNodes()`"
ERROR_MISSING_DISTRARGS             = "ERROR: Missing required field `distrArgs`"
ERROR_MISSING_DISTRARGS_UNISQ       = "ERROR: `xRange`, `yRange` are required in `distrArgs`"
ERROR_MISSING_DISTRARGS_UNICC       = "ERROR: `radius`, `centerLoc` are required in `distrArgs`"
ERROR_MISSING_DISTRARGS_UNIPOLY     = "ERROR: `poly` or `polys` is required in `distrArgs`"
ERROR_MISSING_DISTRARGS_UNIROADNETWORK = "ERROR: `roadNetwork` is required in `distrArgs`"
ERROR_MISSING_DISTRARGS_RING        = "ERROR: `numNodes`, `radius`, `degOffset` are required in `distrArgs`"
ERROR_MISSING_DISTRARGS_CLUSTER     = "ERROR: (`numCluster`, `xRange`, `yRange`) or `centroidLocs` are required in `distrArgs`"
ERROR_MISSING_TIMESPAN              = "ERROR: Missing required field `timeSpan`"
ERROR_MISSING_TWKNOWNARGS           = "ERROR: Missing required field `twKnownArgs`"
ERROR_MISSING_TWKNOWNARGS_TW        = "ERROR: Missing `timeWindows` in `twKnownArgs"
ERROR_MISSING_TWPERIODICARGS        = "ERROR: Missing required field `twPeriodicArgs`"
ERROR_MISSING_TWPERIODICARGS_CT     = "ERROR: Missing `cycleTime` in `twPeriodicArgs`"
ERROR_MISSING_TWPERIODICARGS_ST     = "ERROR: Missing `startTimeInCycle` in `twPeriodicArgs`"
ERROR_MISSING_TWPERIODICARGS_ET     = "ERROR: Missing `endTimeInCycle` in `twPeriodicArgs`"
ERROR_MISSING_TWRANDOMARGS          = "ERROR: Missing required field `twRandomArgs`"
ERROR_MISSING_TWRANDOMARGS_FG       = "ERROR: Missing `indFlag` in `twRandomArgs`"
ERROR_MISSING_TWRANDOMARGS_BTW      = "ERROR: Missing `avgBtwDur` in `twRandomArgs`"
ERROR_MISSING_TWRANDOMARGS_TW       = "ERROR: Missing `avgTWDur` in `twRandomArgs`"
ERROR_MISSING_GANTT                 =("ERROR: Missing required field `gantt`. The format of `gantt` is \n"
                                    + "    [{\n"
                                    + "        'entityID': entityID, \n"
                                    + "        'timeWindow': [startTime, endTime], \n"
                                    + "        'desc': description of the window,\n"
                                    + "        'color': color, \n"
                                    + "        'style': 'solid' \n"
                                    + "    }, ... , \n"
                                    + "    {\n"
                                    + "        'entityID': entityID, \n"
                                    + "        'timeStamps': [timeStamp1, timeStamp2, ..., timeStampN], \n"
                                    + "        'desc': [List of descriptions, correspond to `timeStamps`],\n"
                                    + "        'color': color, \n"
                                    + "        'style': 'solid' \n"
                                    + "    }]\n")
ERROR_ZERO_VECTOR                   = "ERROR: Segment or Ray should be defined using two different points"

# Error message - Incorrect input option
ERROR_OPTS_TWTYPE                   = "ERROR: Incorrect `twType`, valid options includes %s" % ("")
ERROR_OPTS_DISTR                    = "ERROR: Incorrect `distr`, valid options includes %s" % ("")

ERROR_OPTS_SHORTESTPATH_ALGO        = "ERROR: There exists negative weight arcs, cannot choose this algorithm"

# Error message - Incorrect input contents
ERROR_INCOR_NODEIDS                 = "ERROR: Incorrect `nodeIDs`, options are 1) 'All', 2) A list of node IDs"
ERROR_INCOR_TAU                     = "ERROR: Incorrect `edges`, options are 'Euclidean', 'LatLon', 'Grid' (additional info needed to create grid), and travel matrix."
ERROR_INCOR_GANTT_MISSENT           = "ERROR: Missing entity in `gantt`"
ERROR_INCOR_DISTARG                 =("ERROR: Incorrect `distr` for customer locations, options are \n"
                                    + "1) 'uniformSquare', \n"
                                    + "2) 'uniformCircleXY', \n"
                                    + "3) 'uniformCircleLatLon', \n"
                                    + "4) 'uniformPoly', \n"
                                    + "5) 'uniformPolys', \n"
                                    + "6) 'uniformRoadNetwork', \n"
                                    + "7) 'clusterXY', \n"
                                    + "8) 'clusterLatLon'")