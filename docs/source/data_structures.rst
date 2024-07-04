Data Structures
===============

*vrpSolver* uses dictionaries to store different objects such as `nodes`, `arcs`, `roads`, etc. One of the main reasons for doing so is that dictionaries do not require additional Python packages. Additionally, dictionaries can have flexible definitions of the keys, so that the data used in *vrpSolver* can be used in other Python packages that the user is using. In addition, two easy-to-use functions are provided (in `common.py`) for exporting and importing dictionaries. They are :func:`~vrpSolver.common.saveDictionary()` and :func:`~vrpSolver.common.loadDictionary()`.


`nodes`
-------

The `nodes` dictionary has to contain the location information of a node, additionally, it may include other information as follows:

- 'loc': The location of the nodes in (x, y) or (lat, lon). In some functions, the key of 'loc' may be redirected by the `locFieldName` parameter.
- 'demand': The demand of each node.
- 'timeWindow': The available time window of each node.
- 'neighbor': Polygons that are surrounding nodes. In the format of a list of points. This is particularly useful for the "close-enough" problems.
- 'color': Used in :func:`~vrpSolver.plot.plotNodes()`, indicating the color of each node.

`arcs`
------


`roads`
-------

`tau`
-----

`path`
------
