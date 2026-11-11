# Routine to get the colors for MOKAS
# 

import numpy as np
import matplotlib.colors as mpl_colors
from matplotlib import cm
from colorsys import hsv_to_rgb, hls_to_rgb
import bokeh.palettes as palettes

# Load ral colors
# See http://www.ralcolor.com/
def get_ral_colors():
    ral_colors = []
    with open("ral_color_selected2.txt") as f:
        rals = f.readlines()
    ral_colors = [r.split()[-3:] for r in rals]
    ral_colors = [[int(rgb) for rgb in r.split()[-3:]] for r in rals]
    ral_colors = np.array(ral_colors)
    ral_colors = np.concatenate((ral_colors, ral_colors / 2))
    ral_colors = np.concatenate((ral_colors, ral_colors))
    return ral_colors

def get_liza_colors(color='green', whiter=0.6):
    black = (0,0,0)
    if color == 'green':
        c10 = mpl_colors.hex2color("#204A20")
        c11 = mpl_colors.hex2color("#548e72")
        c01 = mpl_colors.hex2color("#377D2C")
        c00 = mpl_colors.hex2color("#1B2E34")
        #c00 = [0.5*(c+wc) for c in c00]
    elif color == 'red':
        c00 = mpl_colors.hex2color("#770811")
        c01 = mpl_colors.hex2color("#727681")
        c10 = mpl_colors.hex2color("#827A74")
        c11 = mpl_colors.hex2color("#4e443c")
    elif color == 'blue':
        c00 = mpl_colors.hex2color("#1C1F39")
        c01 = mpl_colors.hex2color("#CBA465")
        c10 = mpl_colors.hex2color("#7E6647")
        #c11 = mpl_colors.hex2color("#356975") # It's too dark
        c11 = mpl_colors.hex2color("#98BEC7")
        if whiter is not None:
            c00 = [0.5*(c+whiter) for c in c00]    
    if whiter is not None:
        c10 = [0.5*(c+whiter) for c in c10]
        c01 = [0.5*(c+whiter) for c in c01]
        c11 = [0.5*(c+whiter) for c in c11]

    clrs = [black, c00, c01, c10, c11]
    return clrs, mpl_colors.ListedColormap(clrs,'liza_'+color)


def get_cmap(N, cmap='hsv', norm=False):
    """Returns a function that maps each index 
    in 0, 1, ... N-1 to a distinct RGB color.
    http://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib
    """
    color_norm  = mpl_colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap=cmap) 
    pColor = [scalar_map.to_rgba(i)[:3] for i in range(N)]
    if norm:
        norm_factor = 1
    else:
        norm_factor = 255
    map_index_to_rgb_color = [[col[0]*norm_factor,col[1]*norm_factor,col[2]*norm_factor] for col in pColor]
    return map_index_to_rgb_color

def get_colors(num_colors, palette='pastel', norm=False, 
                visualization_library='mpl'):
    """
    for bokeh
    possible palettes: magma, inferno, viridis, cividis

    """
    black = np.array(3*[0])
    white = np.array(3*[255])
    if palette == 'hue':
        colors=[]
        for i in np.arange(0., 360., 360. / num_colors):
            hue = i/360.
            lightness = (50 + np.random.rand() * 10)/100.
            saturation = (90 + np.random.rand() * 10)/100.
            col = hls_to_rgb(hue, lightness, saturation)
            c = [int(col[0]*255),int(col[1]*255),int(col[2]*255)]
            colors.append(c)
        colors = np.random.permutation(colors)
    elif palette == 'pastel':
        colors = (np.random.randint(0, 256, (num_colors,3)) + white) / 2
    else:
        if visualization_library == 'mpl':
            colors = get_cmap(num_colors, palette, norm)
            #colors = np.vstack((black,colors))
        elif visualization_library == 'bokeh':
            _num_colors = num_colors
            if _num_colors > 256:
                num_colors = 256
            if palette == 'magma':
                colors = palettes.magma(num_colors)
            elif palette == 'inferno':
                colors = palettes.inferno(num_colors)
            elif palette == 'cividis':
                colors = palettes.cividis(num_colors)
            elif palette == 'viridis':
                colors = paletters.viridis(num_colors)
            # Check if larger than 256
            colors = (_num_colors//256) * colors + colors[:_num_colors%256]

    return colors

def getKoreanColors(i, n):
    """
    Make a palette in the korean style
    """
    n = float(i)/float(n)*3.
    R = (n<=1.)+ (2.-n)*(n>1.)*(n<=2.)
    G = n*(n<=1.)+ (n>1.)*(n<=2.)+(3.-n)*(n>2.)
    B = (n-1.)*(n>=1.)*(n<2.)+(n>=2.)
    R, G, B = [int(i*255) for i in [R,G,B]]
    return R,G,B

def getPalette(n, palette='ral', noSwitchColor='white', koreanPalette=None):
    """
    get the color palette
    n: int
        n. of points 
    """
    #print("You are using the {} palette".format(palette))
    if type(palette) is not type('str'):
        return palette

    white = np.array([255,255,255]) 
    if koreanPalette is None:
        # Prepare the Korean Palette
        koreanPalette = np.array([getKoreanColors(i, n) 
                                        for i in range(n)])

    if palette == 'korean':
        pColor = koreanPalette
    elif palette == 'randomKorean':
        pColor = np.random.permutation(koreanPalette)
    elif palette == 'random':
        pColor = np.random.randint(0, 256, koreanPalette.shape)
    elif palette == 'pastel':
        pColor = (np.random.randint(0, 256, koreanPalette.shape) + white) / 2
    elif palette == 'randomHue':
        # Use equally spaced colors in the HUE weel, and
        # then randomize
        pColor = [hsv_to_rgb(j/float(n),1, np.random.uniform(0.75,1)) for j in range(n)]
        pColor = np.random.permutation(pColor)
    elif palette == 'hue':
        # Use equally spaced colors in the HUE weel
        pColor = [hsv_to_rgb(j/float(n),1, 1) 
                             for j in range(n)]
    elif palette == 'randomRal':
        ral_colors = get_ral_colors()[:n]
        pColor = np.random.permutation(ral_colors)
    elif palette == 'ral':
        pColor = get_ral_colors()[:n]
    else:
        try:
            pColor = get_cmap(n, palette)
            if palette == 'coolwarm':
                pColor = pColor[::-1]
        except:
            print("Color not available, use pastel instead")
            pColor = (np.random.randint(0, 256, koreanPalette.shape) + white) / 2


    if noSwitchColor == 'black':
        noSwitchColorValue = 3*[0]
    elif noSwitchColor == 'white':
        noSwitchColorValue = 3*[255]
    elif noSwitchColor == 'gray':
        noSwitchColorValue = 3*[125]
    elif noSwitchColor == 'lightgray':
        noSwitchColorValue = 3*[220]
    else:
        print("No color, assuming black")
        noSwitchColorValue = 3*[0]
    return np.concatenate(([noSwitchColorValue], pColor))/255.  


if __name__ == "__main__":
    cw = get_cmap(134, 'coolwarm')
    print(cw)

