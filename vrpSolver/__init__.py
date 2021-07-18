__version__ = "0.0.8"
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

# Graph/network algorithms
from .mst import *
from .matching import *
from .shortestPath import *

# TSP
from .TSP import *

# VRP and its variants
from .consVRP import *
from .ipCVRP import *
from .cgCVRPTW import *