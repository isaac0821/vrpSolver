.. _dictionaries:

Formatted Dictionaries
======================

*geoVeRoPy* uses dictionaries to store different objects such as `nodes`, `arcs`, etc. One of the main reasons for doing so is that dictionaries do not require additional Python packages. Additionally, dictionaries can have flexible definitions of the keys, so that the data used in *geoVeRoPy* can be used in other Python packages that the user is using. In addition, two easy-to-use functions are provided (in `common.py`) for exporting and importing dictionaries. They are :func:`~geoVeRoPy.common.saveDictionary()` and :func:`~geoVeRoPy.common.loadDictionary()`.

.. _nodes:

`nodes`
-------

The `nodes` dictionary has to contain the location information of a node, additionally, the dictionary may include other information such as demands, time windows, etc.

- 'loc': The location of the nodes in (x, y) or (lat, lon). In some functions, the key of 'loc' may be redirected using `locFieldName` argument.
- 'demand': The demand of each node. The key of 'demand' may be redirected using `demandFieldName` argument.
- 'timeWindow': The available time window of each node. The key of 'timeWindow' may be redirected using `timeWindowFieldName` argument.
- 'neighbor': Polygons that are surrounding nodes. In the format of a list of points. This is particularly useful for the "close-enough" problems.
- 'timedSeq': A list of 3-tuples, in the format of (x, y, time). This particularly useful for the "moving target" problems.
- 'color': Used in :func:`~geoVeRoPy.plot.plotNodes()`, indicating the color of each node.
- 'marker': Used in :func:`~geoVeRoPy.plot.plotNodes()`, indicating the marker of each node.
- 'markerSize': Used in :func:`~geoVeRoPy.plot.plotNodes()`, indicating the size of the marker of each node.
- 'label': Used in :func:`~geoVeRoPy.plot.plotNodes()`, if provided, a label will be placed at the node location.
- 'ha': Used in :func:`~geoVeRoPy.plot.plotNodes()`, horizontal alignment of the label, options are 'left', 'center', and 'right'.
- 'va': Used in :func:`~geoVeRoPy.plot.plotNodes()`, vertical alignment of the label, options are 'top', 'middle', and 'bottom'.

Who creates/modifies it?

- :func:`~geoVeRoPy.instance.rndNodes()` creates `nodes`.
- :func:`~geoVeRoPy.instance.rndNodeNeighbors()` modifies `nodes` by adding `neighbor` field.

Who uses it?

- :func:`~geoVeRoPy.instance.rndNodeNeighbors()` must use an existing `nodes` dictionary as input.
- :func:`~geoVeRoPy.geometry.matrixDist()` takes `nodes` as input, calculates traveling matrix between nodes.
- :func:`~geoVeRoPy.plot.plotNodes()` takes `nodes` as input, returns a matplotlib figure with nodes plotted.
- :func:`~geoVeRoPy.plot.plotNodeSeq()` takes `nodes` as input, additionally, this function requires a list of node IDs to plot a sequence of nodes. The function returns a matplotlib figure.
- :func:`~geoVeRoPy.tsp.solveTSP()` takes `nodes` as input, finds the TSP route.

.. _arcs:

`arcs`
------

The `arcs` dictionary has to contain the location of both ends.

- 'arc': The location of both ends of the arc. In some functions, the key of 'arc' may be redirected using `arcFieldName` argument.
- 'color': Used in :func:`~geoVeRoPy.plot.plotArcs()`, indicating the color of each arc.
- 'label': Used in :func:`~geoVeRoPy.plot.plotArcs()`, if provided, a label will be placed at middle of the arc.
- 'ha': Used in :func:`~geoVeRoPy.plot.plotArcs()`, horizontal alignment of the label, options are 'left', 'center', and 'right'.
- 'va': Used in :func:`~geoVeRoPy.plot.plotArcs()`, vertical alignment of the label, options are 'top', 'middle', and 'bottom'.

Who creates/modifies it?

- :func:`~geoVeRoPy.instance.rndArcs()` creates `arcs`

Who uses it?

- :func:`~geoVeRoPy.plot.plotArcs()` takes `arcs` as input, creates a matplotlib figure with arcs plotted.

.. _polygons:

`polygons`
----------

The `polygon` dictionary defines the information of polygons. This dictionary is rarely used, in most scenarios, a :ref:`polys` is used instead. The `polys` data type is a list of :ref:`poly` s, which is a list of points (in (x, y) format, or [x, y] format) that defines a polygon.

- 'anchor': A point within the polygon to "represent" the polygon. In some functions, the key of 'anchor' may be redirected using `anchorFieldName` argument.
- 'poly': A polygon surrounding 'anchor'.

Who creates/modifies it?

- :func:`~geoVeRoPy.instance.rndPolys()` creates `polygons` if returnAsListFlag is set to be False, otherwise, this function creates `polys`.

Who uses it?

- :func:`~geoVeRoPy.plot.plotPolygons()` takes `polygons` as input, creates a matplotlib figure with polygons plotted.

.. .. _vehicles:

.. `vehicles`
.. ----------

.. The `vehicles` dictionary defines the information of different vehicles, which can be used for animation.

.. - 'speed': The speed of the vehicle in [m/s]
.. - 