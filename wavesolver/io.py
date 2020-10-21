import os
import shutil
from wavesolver.linalg import *

# ----------------------------------------------------------------------------

# COLORS
clist = ['limegreen', 'firebrick', 'slateblue', 'brown', 'chocolate',
         'coral', 'hotpink', 'royalblue', 'mediumaquamarine', 'darkslategrey']
cxkcd = ['xkcd:electric blue', 'xkcd:barbie pink', 'xkcd:emerald',
         'xkcd:dandelion', 'xkcd:barney purple', 'xkcd:salmon', 'xkcd:blurple']
cxkcd = ['xkcd:electric blue', 'xkcd:barbie pink', 'xkcd:emerald',
         'xkcd:dandelion', 'xkcd:pale orange', 'xkcd:grey blue',
         'xkcd:blurple']

# teal light brown mauve grey blue
cgreys = ['black', 'dimgrey', 'darkgrey', 'lightgrey', 'gainsboro',
          'whitesmoke']
cdarks = ['seagreen', 'navy', 'indigo', 'maroon', 'darkorange', 'olive']
clights = ['mediumseagreen', 'mediumblue', 'darkviolet', 'firebrick',
           'orange', 'yellow']

# ----------------------------------------------------------------------------


class io():
    """Input/Ouput Configuration"""
    def __init__(self, path='',
                 movie=False,
                 save=True,
                 sep='/',
                 id=0,
                 name='',
                 format='.png',
                 plotCoord='th',
                 plotx=None,
                 ploty=None,
                 Lmax=21,
                 plotRefEigs=False):
        self.path = path
        self.sep = sep
        self.movie = movie
        self.save = save
        self.id = id
        self.name = name
        self.plotCoord = plotCoord
        self.format = format
        self.plotx = plotx
        self.ploty = ploty
        self.Lmax = Lmax
        self.plotRefEigs = plotRefEigs

# ----------------------------------------------------------------------------


def getopts(argv):
    """ Get command line arguments (made redundant by argparse package) """
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            if argv[1]:
                val = argv[1]
            opts[argv[0]] = val  # Add key and value to the dictionary.
        # Reduce the argument list by copying it starting from index 1.
        argv = argv[1:]
    return opts
# ----------------------------------------------------------------------------


def deleteFiles(DIR, fileformat='png'):
    """ Delete files of a given format in a directory (USE CAUTION!) """
    for filename in os.listdir(DIR):
        file_path = os.path.join(DIR, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                if os.path.splitext(file_path)[1] == fileformat:
                    os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

# ----------------------------------------------------------------------------


def dataDisplayConfig(data):
    scale = None
    vmin = None
    vmax = None
    ticks = None
    ticklabels = None
    barlabel = None

    if data == 'n':
        scale = 1E-6
        # scale = 1
        vmin = 0.5
        vmax = 50.
        vmin = 0.1
        vmax = 70.
        vmin = 0.005
        vmax = 70.
        vmin = 0.01
        vmax = 70.
        ticks = [0.5, 1., 5., 10., 50.]
        ticks = [0.1, 0.5, 1., 5., 10., 70.]
        ticks = [0.005, 0.5, 1., 5., 10., 70.]
        ticks = [0.01, 0.1, 1., 10., 70.]
        ticklabels = ['0.5', '1', '5', '10', '50']
        ticklabels = ['0.1', '0.5', '1', '5', '10', '70']
        ticklabels = ['0.005', '0.5', '1', '5', '10', '70']
        ticklabels = ['0.01', '0.1', '1', '10', '70']
        barlabel = r'Water Group Ion Density ($cm^{-3})$'
        barlabel = r'Plasma Density ($cm^{-3})$'
    if data == 'B':
        scale = 1E9
        vmin = 1
        vmax = 10000
        ticks = [1, 10, 100, 1000, 10000]
        # ticklabels = ['10', '100', '1000', '10000']
        ticklabels = [r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$']
        barlabel = r'Magnetic Field (nT)'
    if data == 'vA':
        scale = 1E-3
        vmin = 1.
        vmax = 10000
        vmin = 10.
        vmax = 10000
        ticks = [1, 10, 100, 1000, 10000]
        ticks = [10, 100, 1000, 10000]
        ticklabels = ['1', '10', '100', '1000', '10000']
        ticklabels = ['10', '100', '1000', '10000']
        ticklabels = [r'$10$', r'$10^2$', r'$10^3$', r'$10^4$']
        barlabel = r'Alfven Velocity ($kms^{-1})$'
    if data == 'w':
        scale = 1.
        vmin = 1.
        vmax = 75
        ticks = [1., 15., 30, 45, 60, 75]
        ticklabels = ['1', '15', '30', '45', '60', '75']
        barlabel = r'Period ($min)$'
    return scale, vmin, vmax, ticks, ticklabels, barlabel

# ----------------------------------------------------------------------------


def fancyVarName(var):
    ''' Fancy text for plotting '''
    if var == 'E':
        return r'$E_{\perp}$'
    elif var == 'b':
        return r'$b_{\perp}$'
    elif var == 'b/B':
        return r'$b_{\perp}/B$'
    elif var == 'xi':
        return r'$\xi_{\perp}$'
    else:
        return ''

# ----------------------------------------------------------------------------


def selectSolutionVar(var, solution, SIM, norm=True):
    ''' Select an appropriate variable of the simulation '''
    if var == 'E':
        y = solution.E
    elif var == 'b':
        y = solution.b
    elif var == 'b/B':
        y = solution.b / SIM.B
    elif var == 'xi':
        y = solution.xi
    return normalize(y) if norm else y

# ----------------------------------------------------------------------------


def toPlotCoords(z, SIM):
    ''' Convert simulation coordinates to plot coordinates '''
    if SIM.coords == 'cartesian':
        return z / SIM.units
    if SIM.coords == 'ds':
        return z / SIM.units
    if SIM.coords == 'cos':
        return np.arccos(z) * 180. / np.pi
    if SIM.coords == 'deg':
        return z * 180. / np.pi

# ----------------------------------------------------------------------------


def fromPlotCoords(z, SIM):
    ''' Convert from plot coordinates to simulation coordinates'''
    if SIM.coords == 'cartesian':
        return z * SIM.units
    if SIM.coords == 'ds':
        return z * SIM.units
    if SIM.coords == 'cos':
        return np.cos(z / 180. * np.pi)
    if SIM.coords == 'deg':
        return z / 180. * np.pi

