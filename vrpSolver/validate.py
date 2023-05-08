
import warnings

from .common import *

# A function to validate if a list is a legal representation of a polygon
def validPoly(poly):
    # Check if its a list
    if (type(poly) != list):
        raise InvalidPolygonError("Incorrect format of polygon")
    # Check if its a list of 2-tuples/2-lists
    elif (len(poly) <= 2):
        raise InvalidPolygonError("Polygon need to has at least 2 extreme points")
    else:
        for i in poly:
            if (type(i) is not tuple and type(i) is not list and len(i) != 2):
                raise InvalidPolygonError("Incorrect format of polygon, each extreme point in the list should be a 2-tuple/2-list")
    # TODO: Check if self-cross, etc...
    return True

def validPolys(polys):
    # Check if its a list
    if (type(polys) != list):
        return False
    else:
        for p in polys:
            if (validPoly(p) == False):
                return False
    return True