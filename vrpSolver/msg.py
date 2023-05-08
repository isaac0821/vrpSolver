
# NOTICE ======================================================================
# The file name `msg.py` does NOT stands for Monosodium glutamate. And PLEASE,
# PLEASE, PLEASE do NOT put any MSG into your dishes.
# =============================================================================

# Warning message - Overwrite inputted field
WARNING_OVERWRITE_STARTTIME             = "WARNING: Field `startTime` has been overwritten"
WARNING_OVERWRITE_ENDTIME               = "WARNING: Field `endTime` has been overwritten"

# Error message - Missing required field
ERROR_MISSING_N                         = "ERROR: Missing required field `N`, `N` is an integer indicating number of vertices"
ERROR_MISSING_NODES                     =("ERROR: Missing required field `nodes`. The format of `nodes` is \n"
                                        + "    { \n "
                                        + "        nodeID1: {'loc': (x, y)}, \n"
                                        + "        nodeID2: {'loc': (x, y)}, \n"
                                        + "        ... \n"
                                        + "    }")
ERROR_MISSING_NODES_DISTR               =("ERROR: Missing required field `distr`, `distri` can be chosen from the following options\n"
                                        + "1) (default) Uniformly sample from a square on the Euclidean space\n"
                                        + "    distr = {\n"
                                        + "        'method': 'uniformSquareXY', \n"
                                        + "        'xRange': (0, 100), # A 2-tuple with minimum/maximum range of x, default as (0, 100), \n"
                                        + "        'yRange': (0, 100), # A 2-tuple with minimum/maximum range of y, default as (0, 100), \n"
                                        + "    }\n"
                                        + "2) Uniformly sample from a given polygon on the Euclidean space\n"
                                        + "    distr = {\n"
                                        + "        'method': 'uniformPolyXY', \n"
                                        + "        'polyXY': poly, # polygon of the area, (no holes)\n"
                                        + "        'polyXYs': polys, # alternative option for 'polyXY', as a list of polygons \n"
                                        + "    }\n"
                                        + "3) Uniformly sample from a circle on the Euclidean space\n"
                                        + "    distr = {\n"
                                        + "        'method': 'uniformCircleXY',\n"
                                        + "        'centerXY': (0, 0), # centering location, default as (0, 0), \n"
                                        + "        'radius': 100, # radius of the circle , default as 100\n"
                                        + "    }\n"
                                        + "4) Uniformly sample from a given polygon by lat/lon\n"
                                        + "    distr = {\n"
                                        + "        'method': 'uniformPolyLatLon', \n"
                                        + "        'polyLatLon': polygon of the area, (no holes)\n"
                                        + "        'polyLatLons': alternative option for 'polyLatLon', as a list of polygons \n"
                                        + "    }\n"
                                        + "5) Uniformly sample from a given circle by lat/lon,\n"
                                        + "    distr = {\n"
                                        + "        'method': 'uniformCircleLatLon', \n"
                                        + "        'centerLatLon': required, centering location in lat/lon, \n"
                                        + "        'radiusInMeters': radius of the circle in meters \n"
                                        + "    }\n"
                                        + "6) Uniformly generate from a given polygon on a road network\n"
                                        + "    distr = {\n"
                                        + "        'method': 'roadNetworkPolyLatLon'\n"
                                        + "        'roadNetwork': list of arcs that can be sampled \n"
                                        + "        'polyLatLon': nodes should generated within the polygon, if not provided, will consider the entire network, \n"
                                        + "        'roadClass': list of classes that can be sampled \n"
                                        + "    }\n"
                                        + "7) Uniformly generate from a given circle on a road network\n"
                                        + "    distr = {\n"
                                        + "        'method': 'roadNetworkCircleLatLon', \n"
                                        + "        'roadNetwork': list of arcs that can be sampled \n"
                                        + "        'centerLatLon': [lat, lon], \n"
                                        + "        'radiusInMeters': radius in [m] \n"
                                        + "        'roadClass': list of classes that can be sampled\n"
                                        + "    }\n")
ERROR_MISSING_NODES_DISTR_POLYXY        = "ERROR: Missing required key `polyXY` or `polyXYs` in field `distr`, which indicates a polygon / a list of polygons in the Euclidean space"
ERROR_MISSING_NODES_DISTR_POLYLATLON    = "ERROR: Missing required key `polyLatLon` or `polyLatLons` in field `distr`, which indicates a polygon / a list of polygons in lat/lon format"
ERROR_MISSING_NODES_DISTR_ROADNETWORK   = "ERROR: Missing required key `roadNetwork` in field `distr`. Need to provide the road network where the nodes are generated."


ERROR_MISSING_GANTT                     =("ERROR: Missing required field `gantt`. The format of `gantt` is \n"
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
    
ERROR_ZERO_VECTOR                       = "ERROR: Segment or Ray should be defined using two different points"

# Error message - Incorrect input contents
ERROR_INCOR_NODEIDS                     = "ERROR: Incorrect `nodeIDs`, options are 1) 'All', 2) A list of node IDs"
ERROR_INCOR_TAU                         = "ERROR: Incorrect `edges`, options are 'Euclidean', 'LatLon', 'Grid' (additional info needed to create grid), and travel matrix."
ERROR_INCOR_GANTT_MISSENT               = "ERROR: Missing entity in `gantt`"
ERROR_INCOR_GANTT_ENTITYGROUP           = "ERROR: Incorrect `entities`, options are 1) None (by default), 2) List of strings of entity IDs, or 3) List of list, each list is a list of strings of entity IDs."

# Error message - Time windows related
ERROR_TW_OVERLAP                        = "ERROR: Overlapped time windows."
