
# Information message - 

# Warning message - Overwrite inputted field
WARNING_OVERWRITE_STARTTIME     = "WARNING: Field `startTime` has been overwritten"
WARNING_OVERWRITE_ENDTIME       = "WARNING: Field `endTime` has been overwritten"

# Error message - Missing required field
ERROR_MISSING_N                 = "ERROR: Missing required field `N`, `N` is an integer indicating number of vertices"
ERROR_MISSING_NODES             =("ERROR: Missing required field `nodes`. The format of `nodes` is \n"
                                + "    { \n "
                                + "        nodeID1: {'loc': (x, y)}, \n"
                                + "        nodeID2: {'loc': (x, y)}, \n"
                                + "        ... \n"
                                + "    }")
ERROR_MISSING_TIMESPAN          = "ERROR: Missing required field `timeSpan`"
ERROR_MISSING_TWKNOWNARGS       = "ERROR: Missing required field `twKnownArgs`"
ERROR_MISSING_TWKNOWNARGS_TW    = "ERROR: Missing `timeWindows` in `twKnownArgs"
ERROR_MISSING_TWPERIODICARGS    = "ERROR: Missing required field `twPeriodicArgs`"
ERROR_MISSING_TWPERIODICARGS_CT = "ERROR: Missing `cycleTime` in `twPeriodicArgs`"
ERROR_MISSING_TWPERIODICARGS_ST = "ERROR: Missing `startTimeInCycle` in `twPeriodicArgs`"
ERROR_MISSING_TWPERIODICARGS_ET = "ERROR: Missing `endTimeInCycle` in `twPeriodicArgs`"
ERROR_MISSING_TWRANDOMARGS      = "ERROR: Missing required field `twRandomArgs`"
ERROR_MISSING_TWRANDOMARGS_FG   = "ERROR: Missing `indFlag` in `twRandomArgs`"
ERROR_MISSING_TWRANDOMARGS_BTW  = "ERROR: Missing `avgBtwDur` in `twRandomArgs`"
ERROR_MISSING_TWRANDOMARGS_TW   = "ERROR: Missing `avgTWDur` in `twRandomArgs`"
ERROR_MISSING_GANTT             =("ERROR: Missing required field `gantt`. The format of `gantt` is \n"
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

# Error message - Incorrect input option
ERROR_OPTS_TWTYPE               = "ERROR: Incorrect `twType`, valid options includes %s"

# Error message - Incorrect input contents
ERROR_INCOR_GANTT_MISSENT       = "ERROR: Missing entity in `gantt`"