import random

def colorRandom(
    baseColor:  "Random color based on a main color, e.g. different depth of red" = None
    ) -> "Returns a random color code":
    color = "#%06x" % random.randint(0, 0xFFFFFF)
    return color

def colorScale(
	) -> "":

	return lstColor