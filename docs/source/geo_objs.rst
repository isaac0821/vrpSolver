.. _geo_objs:

Geometry Objects
================

*vrpSolver* uses lists and tuples to represent simple geometry objects such as points, lines, line segments, rays, line segment sequence, polygons and etc.

.. _pt:

`pt`
----

A `pt` data type defines a point/location in Euclidean space or in lat/lon. It can be defined as a 2-tuple, such as (10, 10), or a 2-list, such as [10, 10]. In most cases, these two formats are interchangeable.

.. _line:

`line`
------

A `line` data type is a list of two :ref:`pt` s, the format is [:ref:`pt`, :ref:`pt`]. It can define lines, line segments, or rays. It will be considered differently when it is used in different functions. For example, for
	>>> s = [(0, 0), (5, 5)]
s is considered as a line in the function :func:`~vrpSolver.geometry.isPtOnLine()`
	>>> res = vrpSolver.isPtOnLine(pt = (1, 1), line = s)
s is considered as a line segment in the function :func:`~vrpSolver.geometry.isPtOnSeg()`
	>>> res = vrpSolver.isPtOnSeg(pt = (1, 1), seg = s)
s is considered as a ray in the function :func:`~vrpSolver.geometry.isPtOnRay()`
	>>> res = vrpSolver.isPtOnRay(pt = (1, 1), ray = s)

.. _poly:

`poly`
------

A `poly` data type is a list of :ref:`pt` s, the format is [:ref:`pt`, :ref:`pt`, ..., :ref:`pt`]. It can define line segment sequences and polygons. It is more often used as representing polygons that does not have holes.

.. _polys:

`polys`
-------

A `polys` data type is a list of polys, which can define a list of polygons. The `polys` only contains the coordinate information of polygons, and neglect all other information such as the identifier, colors, styles, etc. If a function requires an input as 'polys', a `polys` will need to be provided, otherwise, if the function asks for 'polygons', a :ref:`polygons` will be needed as it may need those additional information. The `polys` are more often used in this package then :ref:`polygons`.