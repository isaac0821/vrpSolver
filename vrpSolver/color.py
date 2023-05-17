import random
import numpy as np

def colorScale(baseColor = 'R'):
    color = ""
    if (baseColor == 'R'):
        color = "#FFFF00"
    elif (baseColor == 'G'):
        color = "#FFFF00"
    elif (baseColor == 'B'):
        color = "#00FFFF"
    return color

def hex2RGB(colorHex):
    colorHex = colorHex.lstrip('#')
    lv = len(colorHex)
    return tuple(int(colorHex[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def getBWText(hexColor):
    (r, g, b) = hex2RGB(hexColor)     
    m = np.array([
        [0.29900, -0.147108,  0.614777],
        [0.58700, -0.288804, -0.514799],
        [0.11400,  0.435912, -0.099978]
    ])
     
    yuv = np.dot(np.array([[[r, g, b]]]), m)
    yuv[:, :, 1:] += 0.5
    nR = yuv[0][0][0]
    nG = yuv[0][0][1]
    nB = yuv[0][0][2]

    gl = 0.299 * nR + 0.587 * nG + 0.114 * nB
    if (gl >= 192):
        return 'black'
    else:
        return 'white'

def RGB2Hex(colorRGB):
    return '#%02x%02x%02x' % colorRGB

def colorRandom(
    ) -> "Returns a random color code":
    color = "#%06x" % random.randint(0, 0xFFFFFF)
    return color

def colorComplement(hexColor):
    (r, g, b) = hex2RGB(hexColor)
    def hilo(a, b, c):
        if c < b: b, c = c, b
        if b < a: a, b = b, a
        if c < b: b, c = c, b
        return a + c
    k = hilo(r, g, b)
    (revR, revG, revB) = tuple(k - u for u in (r, g, b))
    return RGB2Hex((revR, revG, revB))
