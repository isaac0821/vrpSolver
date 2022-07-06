import random

def colorRandom(
    ) -> "Returns a random color code":
    color = "#%06x" % random.randint(0, 0xFFFFFF)
    return color
