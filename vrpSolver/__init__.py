__version__ = "0.0.30"
__author__ = "Lan Peng"

# History =====================================================================
# v0.0.1  - 04202021 - Initial upload
# v0.0.2  - 04282021 - Add `plotGantt()` function to plot Gantt chart
# v0.0.3  - 05022021 - `plotGantt()` add descriptions
# v0.0.4  - 05042021 - Reorganize structure of package
#                    - Remove `saTSP()` prepare for rewrite
#                    - Add `rndTimeWindowsNodes()` function to generate nodes 
#                      with time windows
# v0.0.5  - 05042021 - Minor fixes
# v0.0.6  - 05052021 - Add `timePin` option to `plotGantt()` function
# v0.0.7  - 05262021 - Gantt enables force drawing time window
# v0.0.8  - 07182021 - Massive change of structure
#                    - Add `heuTSP()` for TSP heuristic
#                    - Add time window functions
# v0.0.9  - 07212021 - Fix some functions in `timeWindows.py` and `plotGantt()`
# v0.0.10 - 08252021 - Sweep algorithm for TSP
# v0.0.11 - 09082021 - Jarvis algorithm for convex hull
# v0.0.12 - 09082021 - Fix convex hull algorithm
# v0.0.13 - 09122021 - More options for `rndPlainNodes()`
# v0.0.14 - 09132021 - Reconstruct
# v0.0.15 - 09292021 - k-Nearest for TSP
# v0.0.16 - 10102021 - Minor fixes
# v0.0.17 - 10102021 - Minor fixes
# v0.0.18 - 10312021 - Minor fixes
# v0.0.19 - 11102021 - A bunch of `geometry.py` functions for finding the intersection
#                      between line segments, rays, polygons, etc.
# v0.0.20 - 11152021 - Add `ptXY2LatLonMercator()`, `ptLatLon2XYMercator()` for projection
#                    - Add `twMovingPtInsidePolyXY()`, `twMovingPtInsidePolyLatLon()` for 
#                      time windows of a moving object that is shadowed by another moving object
# v0.0.21 - 11172021 - Add `getCloseNeighbor()` to find cluster of nodes that all close to each others
# v0.0.22 - 11192021 - Add `getCentroid()` and minor fixes
# v0.0.23 - 11192021 - Minor fixes
# v0.0.24 - 11242021 - Add `rectInWidthLengthOrientationXY()` and `rectInWidthLengthOrientationLatLon()`
# v0.0.25 - 11252021 - Add `plotPoly()`
# v0.0.26 - 12012021 - Add cheapest insertion for TSP
# v0.0.27 - 12062021 - Minor fixes
# v0.0.28 - 12062021 - Minor fixes
# v0.0.29 - 12062021 - Minor fixes
# v0.0.30 - 12122021 - Add `gridPathFinding()` with A* algorithm, `plotGrid()`, `plotGridPath()`
# =============================================================================

# Constants and messages
from .const import *
from .msg import *

# Basic modules
from .common import *
from .plot import *
from .instance import *
from .timeWindows import *
from .color import *
from .vector import *

# Geometry
from .geometry import *
from .relation import *

# Graph/network algorithms
from .graph import *

# Parallel machine scheduling problem
from .PMS import *

# Warehouse
from .warehouse import *

# TSP
from .ipTSP import *
from .heuTSP import *

# VRP
# from .ipVRP import *	# Debugging