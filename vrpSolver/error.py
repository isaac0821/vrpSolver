class UnsupportedInputError(Exception):
    pass

class ZeroVectorError(Exception):
    pass

class InvalidPolygonError(Exception):
    pass

class EmptyError(Exception):
    pass

class MissingParameterError(Exception):
    # Missing required input parameter(s)
    pass

class KeyExistError(Exception):
    pass

class KeyNotExistError(Exception):
    pass

class OutOfRangeError(Exception):
    pass

class NonOverlapIntervalTreeError(Exception):
    pass