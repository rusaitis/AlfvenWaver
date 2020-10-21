#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from wavesolver.model import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LogNorm
from matplotlib.colors import NoNorm
from matplotlib.colors import Normalize
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as pe
from wavesolver.helperFunctions import *
# import plotly.graph_objs as go
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from wavesolver.configurations import *
from wavesolver.io import *
# import mayavi.mlab as m
import matplotlib as mpl
from wavesolver.helperFunctions import *
# plt.style.use('PaperDoubleFig.mplstyle')

mpl.rcParams['hatch.linewidth'] = 0.4  # previous svg hatch linewidth
mpl.rcParams['hatch.linewidth'] = 2.0  # previous svg hatch linewidth
#-----------------------------------------------------------------------------


def calculateSecondGradient(solution, SIM, norm=True):
    y = solution.b
    dy = np.gradient(y)
    dtheta = np.gradient(dy)
    return normalize(dtheta) if norm else dtheta
#-----------------------------------------------------------------------------
def axhlines(ys, ax=None, **plot_kwargs):
    """
    Draw horizontal lines across plot
    :param ys: A scalar, list, or 1D array of vertical offsets
    :param ax: The axis (or none to use gca)
    :param plot_kwargs: Keyword arguments to be passed to plot
    :return: The plot object corresponding to the lines.
    """
    if ax is None:
        ax = plt.gca()
    ys = np.array((ys, ) if np.isscalar(ys) else ys, copy=False)
    lims = ax.get_xlim()
    y_points = np.repeat(ys[:, None], repeats=3, axis=1).flatten()
    x_points = np.repeat(np.array(lims + (np.nan, ))[None, :], \
                         repeats=len(ys), axis=0).flatten()
    plot = ax.plot(x_points, y_points, scalex = False, **plot_kwargs)
    return plot

#-----------------------------------------------------------------------------

def axvlines(xs, ax=None, **plot_kwargs):
    """
    Draw vertical lines on plot
    :param xs: A scalar, list, or 1D array of horizontal offsets
    :param ax: The axis (or none to use gca)
    :param plot_kwargs: Keyword arguments to be passed to plot
    :return: The plot object corresponding to the lines.
    """
    if ax is None:
        ax = plt.gca()
    xs = np.array((xs, ) if np.isscalar(xs) else xs, copy=False)
    lims = ax.get_ylim()
    x_points = np.repeat(xs[:, None], repeats=3, axis=1).flatten()
    y_points = np.repeat(np.array(lims + (np.nan, ))[None, :], \
                         repeats=len(xs), axis=0).flatten()
    plot = ax.plot(x_points, y_points, scaley = False, **plot_kwargs)
    return plot

#-----------------------------------------------------------------------------

def vector_plot_pl(tvects,is_vect=True,orig=[0,0,0]):
    """Plot vectors using plotly"""

    if is_vect:
        if not hasattr(orig[0],"__iter__"):
            coords = [[orig,np.sum([orig,v],axis=0)] for v in tvects]
        else:
            coords = [[o,np.sum([o,v],axis=0)] for o,v in zip(orig,tvects)]
    else:
        coords = tvects

    data = []
    for i,c in enumerate(coords):
        X1, Y1, Z1 = zip(c[0])
        X2, Y2, Z2 = zip(c[1])
        vector = go.Scatter3d(x = [X1[0],X2[0]],
                              y = [Y1[0],Y2[0]],
                              z = [Z1[0],Z2[0]],
                              marker = dict(size = [0,5],
                                            color = ['blue'],
                                            line=dict(width=5,
                                                      color='DarkSlateGrey')),
                              name = 'Vector'+str(i+1))
        data.append(vector)

    layout = go.Layout(
             margin = dict(l = 4,
                           r = 4,
                           b = 4,
                           t = 4)
                  )
    fig = go.Figure(data=data,layout=layout)
    fig.show()

    p0 = [0.799319, -3.477045e-01, 0.490093]
    p1 = [0.852512, 9.113778e-16, -0.522708]
    p2 = [0.296422, 9.376042e-01, 0.181748]

    # PLOTLY
    # vector_plot([p0,p1,p2])

    # PYPLOT
    # origin = [0,0,0]
    # X, Y, Z = zip(origin,origin,origin) 
    # U, V, W = zip(p0,p1,p2)

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.quiver(X,Y,Z,U,V,W,arrow_length_ratio=0.01)
    # plt.show()

    #MAYAVI
    # origin = [0,0,0]
    # X, Y, Z = zip(origin,origin,origin) 
    # U, V, W = zip(p0,p1,p2)

    # m.quiver3d(X,Y,Z,U,V,W)

#-----------------------------------------------------------------------------

def plotxAxis(ioconfig, SIM):
    if ioconfig.plotCoord == 'th':
        xlim = (-90, 90.)
        xaxis = -(np.rad2deg(SIM.th) - 90.)
        label = 'Latitude (deg)'
    if ioconfig.plotCoord == 's':
        xlim = None
        xaxis = SIM.z / SIM.units
        label = r'Distance ($R_S$)'
    return xaxis, label, xlim

def loadFigConfig(type, SIM, ioconfig=None):
    if type == 'wL':
        filename = 'wL'
        # title = 'Root search using a shooting method' 
        # ylabel = 'Boundary overshoot at y[n] /\n Arbitrary amplitude'
        xlabel = r'L / [R_S]'
        config = figConfig(filename=filename, xlabel=xlabel)
    if type == "erf":
        # plt.style.use('ggplot')
        filename = 'erf'
        title = 'Root search using a shooting method' 
        ylabel = 'Boundary overshoot at y[n] /\n Arbitrary amplitude'
        xlabel = '$w$' 
        xlim = [SIM.wlim[0], SIM.wlim[1]]
        ylim = [-1.05, 1.05]
        config = figConfig(filename=filename, title=title, xlabel=xlabel, ylabel=ylabel, xlim=xlim, ylim=ylim)
    if type == "waveSolution":
        filename = 'solution'
        title = 'Wave Solutions' 
        ylabel = 'Arbitrary Amplitude'
        xlabel = SIM.zlabel
        ylim = [-1.05,1.05]
        if ioconfig is not None:
            if ioconfig.plotCoord == 'th':
                xlim = [0., 180.]
            if ioconfig.plotCoord == 's':
                xlim = [SIM.zlim[0], SIM.zlim[1]]
        xlim = None
        legend = True
        config = figConfig(filename=filename, title=title, xlabel=xlabel, ylabel=ylabel, xlim=xlim, ylim=ylim, legend=legend)
    if type == "insetDensity":
        filename = 'density'
        title = None 
        ylabel = '$v_A$'
        xlabel = SIM.zlabel
        ls = ['--']
        color = ['black']
        # ylim = [0., SIM.vA1]
        ylim = None
        xlim = [toPlotCoords(SIM.zlim[0], SIM), toPlotCoords(SIM.zlim[1], SIM)]
        inset = True
        patchAlpha = 0.5
        fontsize=6.
        config = figConfig(filename=filename, title=title, xlabel=xlabel, ylabel=ylabel, xlim=xlim, ylim=ylim, inset=inset, color=color, patchAlpha=0.5, fontsize=fontsize)
    if type == "insetwL":
        filename = 'wL'
        title = None 
        ylabel = r'$\omega$'
        xlabel = '$l$'
        ls = ['--']
        color = ['black']
        # ylim = [0., SIM.vA1]
        ylim = None
        # xlim = [SIM.zlim[0]/SIM.units, SIM.zlim[1]/SIM.units]
        xlim = None
        inset = True
        patchAlpha = 0.5
        fontsize=6.
        config = figConfig(filename=filename, title=title, xlabel=xlabel, ylabel=ylabel, xlim=xlim, ylim=ylim, inset=inset, color=color, patchAlpha=0.5, fontsize=fontsize)
    if type == "insetRatios":
        filename = 'ratios'
        title = None 
        ylabel = None
        xlabel = None
        ls = ['--']
        color = ['white']
        # ylim = [0., SIM.vA1]
        ylim = None
        # xlim = [SIM.zlim[0]/SIM.units, SIM.zlim[1]/SIM.units]
        xlim = None
        inset = True
        patchAlpha = 0.5
        fontsize=6.
        config = figConfig(filename=filename, title=title, xlabel=xlabel, ylabel=ylabel, xlim=xlim, ylim=ylim, inset=inset, color=color, patchAlpha=0.5, fontsize=fontsize)
    return config

class figConfig():
    """Input/Ouput Configuration"""
    def __init__(self, filename='plot', size=[8,6], title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, legend=None, id=0, legendLoc='lower right', hline=0, vline=None, save=False, logx=False, logy=False, movie=False, color=['black'], ls=['-'], lw=[1], inset=False, patchAlpha=None, fontsize=None):
        self.filename = filename
        self.size = size
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xlim = xlim
        self.ylim = ylim
        self.legend = legend
        self.legendLoc = legendLoc
        self.hline = hline
        self.vline = vline
        self.id = id
        self.color = color
        self.ls = ls
        self.lw = lw
        self.save = save
        self.logx = logx
        self.logy = logy
        self.movie = movie
        self.inset = inset
        self.patchAlpha = patchAlpha
        self.fontsize = fontsize

    def setup(self, fig, ax):
        # if not self.inset:
            # if self.size is not None: fig.set_size_inches(self.size[0],self.size[1])
        if self.title is not None: ax.set_title(self.title)
        if self.xlabel is not None: ax.set_xlabel(self.xlabel)
        if self.ylabel is not None: ax.set_ylabel(self.ylabel)
        # if self.xlim is not None: ax.set_xlim(self.xlim[0],self.xlim[1])
        if self.ylim is not None: ax.set_ylim(self.ylim[0],self.ylim[1])
        # if self.legend is not None: ax.legend(loc = self.legendLoc, fancybox=True, framealpha=0.5)
        if self.hline is not None: ax.axhline(y=0, lw=1, ls='--', c='black', alpha=0.3)
        if self.vline is not None: ax.axvline(x=0, lw=1, ls='--', c='black', alpha=0.3)
        if self.patchAlpha is not None: 
            ax.patch.set_alpha(self.patchAlpha)
        if self.fontsize is not None: 
            ax.xaxis.label.set_size(self.fontsize)
            ax.yaxis.label.set_size(self.fontsize)
            ax.xaxis.label.set_alpha(self.patchAlpha)
            ax.yaxis.label.set_alpha(self.patchAlpha)
            ax.tick_params(axis = 'both', which = 'major', labelsize = self.fontsize)
    def render(self, fig, ioConfig, hold=False):
        if self.save:
            if ioConfig.path is None: fullpath = '' 
            if ioConfig.path is not None: fullpath = ioConfig.path + '/'
            # if self.id is None: fullpath += self.filename + ioConfig.format
            # if self.id is not None: fullpath += self.filename + '_%03d' % self.id + ioConfig.format
            if ioConfig.id is not None: fullpath += self.filename + ioConfig.name + '_%03d' % ioConfig.id + ioConfig.format
            plt.savefig(fullpath)
            if not hold: plt.close(fig)
        else:
            plt.show()
    def processMovie(self, ioconfig, path=None, res='5136x2880', fps=15, filename='movie.mp4', output='movie.mp4'):
        import os
        if path is not None: fullpath = path + "/" + self.filename 
        if path is None: fullpath = ioconfig.path + ioconfig.sep + self.filename
        os.system("ffmpeg -framerate "+str(fps)+ " -i " + fullpath + "_%03d.png -vcodec mpeg4 -s:v " + res + " -y " + output)
        # os.system("ffmpeg -framerate 15 -i Output/solution_%03d.png -vcodec mpeg4 -s:v 5136x2880 -y movie.mp4")
    def plotInset(self, fig, ax, datax, datay, insetConfig, titles=None):
        # inset_axes = plt.axes([.7, .7, .15, .15], facecolor='w')
        inset_ax = inset_axes(ax,
                        width="20%", # width = 30% of parent_bbox
                        height="15%", # height : 1 inch
                        # height=1., # height : 1 inch
                        loc=1)
        insetConfig.setup(fig, inset_ax)
        n=0

        for y in datay:
            ax = inset_ax
            ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
            if n > 0: 
                # ax.set_yscale('log')
                ax = plt.twinx(inset_ax)
                ax.plot(datax, y, color='red', lw=insetConfig.lw[0], ls='--')
                ax.yaxis.label.set_size(8)
                ax.yaxis.label.set_alpha(0.9)
                ax.yaxis.label.set_color('red')
                # ax.set_ylim(0,1.)
                ax.tick_params(axis='y', labelcolor='red')
                ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
                # ax.set_ylabel(titles[n])
            else:
                ax.plot(datax, y, color=insetConfig.color[0], lw=insetConfig.lw[0], ls=insetConfig.ls[0])
            if titles is not None: 
                title = titles[n]
                ax.set_ylabel(title)
            n += 1
        # inset_ax.set_yscale('log')
        # inset_ax.set_ylim(0,300)
        inset_ax.set_xlim(0,180)
    def plotInset2(self, fig, ax, data, insetConfig, highlight=None):
        # inset_axes = plt.axes([.7, .7, .15, .15], facecolor='w')
        inset_ax = inset_axes(ax,
                        width="30%", # width = 30% of parent_bbox
                        height="25%", # height : 1 inch
                        # height=1., # height : 1 inch
                        loc=2)
        insetConfig.setup(fig, inset_ax)
        plotEigenfrequenciesFieldLineLengthSmall(fig, inset_ax, data[0], data[1], plotExact=True)
        if highlight is not None:
            inset_ax.scatter(highlight[0], highlight[1], s=25, color='black', alpha=0.4)
            inset_ax.axvline(x=highlight[0], color='black', ls='--', alpha=0.4)
        # inset_ax.plot(data[0], data[1], color=insetConfig.color[0], lw=insetConfig.lw[0], ls=insetConfig.ls[0])
    def plotInset3(self, fig, ax, data, insetConfig, highlight=None):
        # inset_axes = plt.axes([.7, .7, .15, .15], facecolor='w')
        ax1 = fig.add_axes([0.02, 0.02, 0.03, 0.06])
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(bottom=False)
        ax1.tick_params(labelbottom=False)
        ax1.tick_params(axis='both', which='major', labelsize=4)
        ax1.bar(0, data[1][0], color='yellow', alpha=0.5, width=3.0, align='center')
        ax1.bar(8, data[1][1], color='yellow', alpha=0.5, width=3.0, align='center')
        insetConfig.setup(fig, ax1)
        ax1.set_ylabel('travel time (min)', fontsize=8)
        # ax2.set_ylim([0,60])
        ax1.text(0, data[1][0], round(data[1][0],1), fontsize=10,  rotation=0, color='white', ha='center', va='bottom', alpha=0.8)
        ax1.text(0, -2, r'$t_{in}$', fontsize=10,  rotation=0, color='white', ha='center', va='bottom')
        ax1.text(8, data[1][1], round(data[1][1],1), fontsize=10,  rotation=0, color='white', ha='center', va='bottom', alpha=0.8)
        ax1.text(8, -2, r'$t_{out}$', fontsize=10,  rotation=0, color='white', ha='center', va='bottom')

def plotSmallAxes(ax, data):
    # inset_axes = plt.axes([.7, .7, .15, .15], facecolor='w')
    color='xkcd:sea green'
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(bottom=False)
    ax.tick_params(labelbottom=False)
    ax.tick_params(axis='both', which='major', labelsize=4)
    ax.bar(0, data[0], color=color, alpha=0.5, width=3.0, align='center')
    ax.bar(4, data[1], color=color, alpha=0.5, width=3.0, align='center')
    # insetConfig.setup(fig, ax)
    ax.set_ylabel('time (min)', fontsize=14)
    # ax2.set_ylim([0,60])
    ax.text(0, data[0], r'$t_{in}=$' + str(round(data[0],1)), fontsize=14,  rotation=0, color='white', ha='center', va='bottom', alpha=0.8)
    # ax.text(0, -2, r'$t_{in}$', fontsize=14,  rotation=0, color='white', ha='center', va='bottom')
    ax.text(4, data[1], r'$t_{out}=$' + str(round(data[1],1)), fontsize=14,  rotation=0, color='white', ha='center', va='bottom', alpha=0.8)
    # ax.text(4, -2, r'$t_{out}$', fontsize=14,  rotation=0, color='white', ha='center', va='bottom')

#-----------------------------------------------------------------------------

def plotOvershoot(parametersAll, boundaryOvershoots, roots, save=False, show=True, id=0, ylim=None, xlim=None, xunits='angFreq', markers=None, path=None, filename='overshoot', markerLabels=None,insetData=None, insetSettings=None):
    

    # ax.plot( errors.t, errors.y, ls='--', marker='o', c='b', ms=3)
    # for root in rootsOptimized:
    #     ax.axvline(x=root, ls='--', c='red')
    # for sol in solutions:
    #     ax.axvline(x=sol.roots[0], ls='--', c='green')
    # ax.scatter( errors.zeros, [0]*len(errors.zeros), c='r', s=20)
    # # ax.set_plot(t2, xs2.y[1,:], 'g-')
    # ax.set_title( 'Shooting Method' )
    # ax.set_xlabel( '$w^2$' )
    # ax.set_ylabel( '$erf$' )
    # # ax.set_legend( ( '%3d points' % n, '%3d points' % n ), loc='lower right' )
    # ax.axhline(y=0, ls='--', c='black')


    fig = plt.figure(figsize=(12, 3))
    ax = fig.add_subplot(111) 
    xlabel = r'Angular frequency $\omega$ / $rads^{-1}$'
    if xunits == 'mins':
        parametersAll = angular2mins(parametersAll)
        roots = angular2mins(roots)
        xlim = [0, max(roots)+10]
        xlabel = r'Period / $min$'
        if markers is not None:
            markers = angular2mins(markers)
    ax.plot(parametersAll, boundaryOvershoots, marker='o', markersize=2, color=clist[2])
    ax.scatter(roots, len(roots)*[0], s=35, c='r')
    ax.set_title( 'Root search using a shooting method' )
    ax.set_ylabel( 'Boundary overshoot at y[n] /\n Arbitrary amplitude' )
    ax.set_xlabel( xlabel )
    ax.axhline(y=0,ls='--', color='black', alpha=0.8)
    # for i in range(0, len(roots)):
        # ax.axvline(x=roots[i],ls='--', color='red', alpha=0.8)
    ax.set_xlim(left=0)
    if ylim is not None: ax.set_ylim([ylim[0],ylim[1]]) 
    if xlim is not None: ax.set_xlim([xlim[0],xlim[1]]) #xlim=[2,150]
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    if markers is not None:
        for i in range(0, len(markers)):
            ax.axvline(x=markers[i],ls='--', lw=2, color=clist[0], alpha=0.8)
    if markerLabels is not None:
        for i in range(0, len(markerLabels)):
            # ax.text(x=markers[i],ls='--', lw=2, color='orange', alpha=0.5)
            # ax.text(x=markers[i], ylim[0]*0.8, markerLabels[i], ha='center',
            #     va='center', transform=ax.transData, fontsize=12, bbox=dict(facecolor='red', alpha=0.5))
            ax.annotate(markerLabels[i], xy=(markers[i], ylim[0]*0.8), xycoords="data",
                  va="center", ha="center", color="w",
                  bbox=dict(boxstyle="round", fc="black", alpha=0.5))
    # ax.legend(loc='lower right')
    ax.grid()
    ax.minorticks_on()
    # minor_ticks = np.arange(0, 101, 5)
    # ax.set_xticks(minor_ticks, minor=True)
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%d"))
    plt.rcParams['xtick.labelsize']=6
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    if insetData is not None:
        # inset_axes = plt.axes([.7, .7, .15, .15], facecolor='w')
        inset_axes = inset_axes(ax,
                        width="20%", # width = 30% of parent_bbox
                        height=1., # height : 1 inch
                        loc=1)
        plotInset(inset_axes, insetData, insetSettings)
    # fig.tight_layout()
    if save:
        if path is None: fullpath = '' 
        if path is not None: fullpath = path + '/'
        if id is None: fullpath += filename + '.png'
        if id is not None: fullpath += filename + '_%03d.png' % id
        plt.savefig(fullpath)
        plt.close(fig)
    if show:
        plt.show()
#-----------------------------------------------------------------------------

def plotSolutions(z_list, solutions, z0, model='uniform', vA0=vA0_Yates, vA1=vA1_Yates, labels=None, colors=None, title='Solutions', save=False, show=True, id=0, path=None, filename='solutions', figsize=[8,6], ylim=None, insetData=None, insetSettings=None):
    fig = plt.figure(figsize=(figsize[0],figsize[1]))
    ax = fig.add_subplot(111) 
    for i in range(0, len(solutions)):
        z = z_list[i]
        solution = solutions[i]
        if labels is not None: label = labels[i]
        if labels is None: label = 'n=%2d mode' % (i+1)
        if colors is not None: color = colors[i]
        if colors is None: color = clist[i%10]
        # ax.plot( z, solution, color=color, ls='-', lw=2, label=label)
        ax.plot( z, solution, ls='-', lw=2, label=label)
    ax.axhline(y=0, ls='--', color='black', lw=1, alpha=0.5)
    ax.axvline(x=-z0, ls='--', color='black', lw=1, alpha=0.5)
    ax.axvline(x=z0, ls='--', color='black', lw=1, alpha=0.5)
    vAratio = vA0/vA1
    sheetAlpha = 0.7 - 0.7*vAratio
    if model == 'uniform':
        ax.axvspan(-z0, z0, facecolor='gray', alpha=sheetAlpha)
    elif model == 'linear':
        slabRes = 40
        for i in range(0, slabRes):
            zbox_width = 2 * z0/slabRes
            zbox0= -z0 + i*zbox_width
            zbox1= -z0 + i*zbox_width + zbox_width
            vA = vA0_model(-z0 + i*zbox_width, z0, vA0=vA0, vA1=vA1, model=model)
            vAratio = vA/vA1
            sheetAlpha = 0.7 - 0.7*vAratio
            ax.axvspan(zbox0, zbox1, facecolor='gray', alpha=sheetAlpha)
    ax.set_title( title )
    ax.set_xlabel( 'length / $R_S$' )
    ax.set_ylabel( 'Arbitrary Amplitude' )
    ax.legend(loc='lower right')
    if ylim is not None: ax.set_ylim([ylim[0],ylim[1]]) 
    # if ylim is None: ax.set_ylim([-1,1]) 
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    if insetData is not None:
        # inset_axes = plt.axes([.7, .7, .15, .15], facecolor='w')
        inset_axes = inset_axes(ax,
                        width="20%", # width = 30% of parent_bbox
                        height=1., # height : 1 inch
                        loc=1)
        plotInset(inset_axes, insetData, insetSettings)
        # plt.tight_layout()
    # inset_axes.plot(z_list[0], [0]*len(z_list[0]))

    if save:
        if path is None: fullpath = '' 
        if path is not None: fullpath = path + '/'
        if id is None: fullpath += filename + '.png'
        if id is not None: fullpath += filename + '_%03d.png' % id
        plt.savefig(fullpath)
        plt.close(fig)
    if show:
        plt.show()

#-----------------------------------------------------------------------------

def plotEigenfrequenciesDensity(n_list, eigenfrequencyListMinuntes, roots, labels=None, colors=None, path=None, save=False, show=True, id=None, filename='eigenfrequenciesDensity'):
    #Plot Eigenfrequencies against Field Length
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111) 
    # ax.set_yscale('log')
    ax.set_xlabel('Density variation strength')
    ax.set_ylabel('Eigenfrequency / mHz', color='black')
    # ax.set_ylim([0.02,1.])
    # ax.set_xlim([0,60])
    # ax2 = ax.twinx()
    color = 'black'
    for i in range(0, len(eigenfrequencyListMinuntes)):
        n = n_list[i]
        eigenfrequencies = eigenfrequencyListMinuntes[i]
        if labels is not None: label = labels[i]
        if labels is None: label = 'n=%2d mode - Numerical Solution' % i
        if colors is not None: color = colors[i]
        if colors is None: color = clist[i]
        ax.scatter(n, eigenfrequencies, label=label, color=color, s=15)
    ax.axhline(y=60,ls='--', color='black', alpha=0.6, label='60 min')
    ax.text(0.8, 60+10, '60 min', fontsize=9, color=color, alpha=0.8)
    ax.axhline(y=30,ls='--', color='black', alpha=0.6, label='30 min')
    ax.text(0.8, 30+5, '30 min', fontsize=9, color=color, alpha=0.8)
    ax.axhline(y=90,ls='--', color='black', alpha=0.6, label='90 min')
    ax.text(0.8, 90+10, '90 min', fontsize=9, color=color, alpha=0.8)
    ax.set_ylim([1/(1E-3)/60., 1/(0.02E-3)/60.])
    ax.set_ylabel('Eigenperiod / mins', color='black')
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.tick_params(axis='y', labelcolor=color)
    ax.legend(loc="lower right")
    fig.tight_layout()
    if save:
        if path is None: fullpath = '' 
        if path is not None: fullpath = path + '/'
        if id is None: fullpath += filename + '.png'
        if id is not None: fullpath += filename + '_%03d.png' % id
        plt.savefig(fullpath)
        plt.close(fig)
    if show:
        plt.show()

#-----------------------------------------------------------------------------

def plotEigenfrequenciesFieldLineLength(ax, xvalues, eigenfrequencyList, plotAgainst='L', plotStyle=None, plotWeight=None, plotAlpha=None, plotMarker='.', markerSize=6, labelPrefix=''):
    # Plot Eigenfrequencies against Field Length

    # ax.get_xaxis().set_visible(False)
    for i in range(0, len(eigenfrequencyList)):
        eigenfrequencies = eigenfrequencyList[i]
        label = labelPrefix + 'm = %d mode' % (i+1)
        plotWeight = 1.5 if plotWeight is None else plotWeight
        plotStyle = '-' if plotStyle is None else plotStyle
        plotAlpha = 1. if plotAlpha is None else plotAlpha
        plotMarker = None if plotMarker is None else plotMarker

        # if labels is not None: label = labels[i]
        # if labels is None: label = 'n=%2d mode - Numerical Solution' % i
        # if colors is not None: color = colors[i]
        # if colors is None: color = clist[i]
        color = clights[i]
        color = cxkcd[i]
        # if 'Day' in labelPrefix: color = cxkcd[i]
        # if 'Day' in labelPrefix: color = cdarks[i]
        # if 'Night' in labelPrefix: color = clights[i]
        
        # ax.plot(xvalues, eigenfrequencies, label=label, ls=plotStyle, color=color, lw=plotWeight, alpha=plotAlpha, marker=plotMarker, markersize=markerSize, markeredgecolor=None, markeredgewidth=None, path_effects=[pe.Stroke(linewidth=2.2, alpha=0.4, foreground='black'), pe.Normal()])
        ax.plot(xvalues, eigenfrequencies, label=label, ls=plotStyle, color=color, lw=plotWeight, alpha=plotAlpha, marker=plotMarker, markersize=markerSize, markeredgecolor=None, markeredgewidth=None)
        # ax.plot(xvalues, eigenfrequencies, label=label, ls=plotStyle, color=color, lw=plotWeight, alpha=plotAlpha)

        # ax.plot(xvalues, eigenfrequencies, label=label, ls=plotStyle, color=color, lw=plotWeight, alpha=plotAlpha, path_effects=[pe.Stroke(linewidth=plotWeight+0.5, alpha=0.1, foreground='black'), pe.Normal()])
        # ax.scatter(xvalues, eigenfrequencies, label=label, color=color, lw=plotWeight, alpha=plotAlpha, marker=plotMarker)

#-----------------------------------------------------------------------------

def plotEigenfrequenciesFieldLineLengthSmall(fig, ax, L, eigenfrequencyList, plotExact=True):
    # Plot Eigenfrequencies against Field Length
    color='black'
    for i in range(0, len(eigenfrequencyList)):
        eigenfrequencies = eigenfrequencyList[i]
        # if labels is not None: label = labels[i]
        # if labels is None: label = 'n=%2d mode - Numerical Solution' % i
        # if colors is not None: color = colors[i]
        # if colors is None: color = clist[i]
        ax.scatter(L, eigenfrequencies, label='m = %d' % (i+1), s=10)

    if plotExact:
        yatesEvenMode = np.loadtxt('even_mode.csv', skiprows=0)
        yatesOddMode = np.loadtxt('odd_mode.csv', skiprows=0)
        ax.plot(yatesEvenMode[:,0], yatesEvenMode[:,1], label='n=1 (even) mode - Yates', color='tab:blue')
        ax.plot(yatesOddMode[:,0], yatesOddMode[:,1], label='n=2 (odd) mode - Yates', color='tab:red')
    ax.set_yscale('log')
    # ax.set_xlabel(r'length /$R_{\it S}$'))
    ax.set_xlabel(r'l /$R_{\it S}$')
    # ax.set_ylabel('Eigenfrequency / mHz', color='black')
    ax.set_ylabel('$f$ / mHz', color='black')
    ax.set_ylim([0.02,1.])
    ax.set_xlim([0,60])
    ax2 = ax.twinx()
    color = 'black'
    ax2.axhline(y=60,ls='--', color='black', alpha=0.6, label='60 min')
    # ax2.text(60/4*3, 60+10, '60 min', fontsize=9, color=color, alpha=0.8)
    ax2.axhline(y=30,ls='--', color='black', alpha=0.6, label='30 min')
    # ax2.text(60/4*3, 30+5, '30 min', fontsize=9, color=color, alpha=0.8)
    ax2.axhline(y=90,ls='--', color='black', alpha=0.6, label='90 min')
    # ax2.text(60/4*3, 90+10, '90 min', fontsize=9, color=color, alpha=0.8)
    ax2.set_ylim([1/(1E-3)/60., 1/(0.02E-3)/60.])
    # ax2.set_ylabel('Eigenperiod / mins', color='black')
    ax2.set_yscale('log')
    ax2.invert_yaxis()
    ax2.tick_params(axis='y', labelcolor=color)
    # ax.legend(loc="lower right")
    fig.tight_layout()
#-----------------------------------------------------------------------------

def processMovie( path=None, filename='solution_', res='5136x2880', fps=15, output='movie.mp4' ):
    import os
    if path is not None: fullpath = path + "/" + filename 
    if path is None: fullpath = filename
    os.system("ffmpeg -framerate "+str(fps)+ " -i " + fullpath + "_%03d.png -vcodec mpeg4 -s:v " + res + " -y " + output)
    # os.system("ffmpeg -framerate 15 -i Output/solution_%03d.png -vcodec mpeg4 -s:v 5136x2880 -y movie.mp4")
    # os.system("ffmpeg -framerate 15 -i Output/solution_%03d.png -vcodec mpeg4 -s:v 5136x2880 -y movie.mp4")

#-----------------------------------------------------------------------------

def drawSaturn(axis, color='black', dim='2D', alpha=0.8):
    if dim == '2D':
        saturn = plt.Circle((0, 0), 1., facecolor=color, fill=True, alpha=0.5, zorder=1E6)
        ringsright = patches.Rectangle((1.0, -0.07), 1.32, 0.14,facecolor=color, alpha=0.5, zorder=1E6)
        ringsleft = patches.Rectangle((-1.0, -0.07),-1.32, 0.14,facecolor=color, alpha=0.5, zorder=1E6)
        axis.add_artist(saturn)
        axis.add_artist(ringsleft)
        axis.add_artist(ringsright)

    if dim == '3D':
        # draw sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = np.cos(u)*np.sin(v)
        y = np.sin(u)*np.sin(v)
        z = np.cos(v)
        axis.plot_wireframe(x, y, z, color=color, alpha=alpha)

        point  = np.array([0, 0, 0])
        normal = np.array([0, 0, 1])
        # point2 = np.array([10, 50, 50])
        d = -point.dot(normal)

        xx, yy = np.meshgrid(range(-15,15), range(-15,15))
        z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

        # z = 0
        axis.plot_surface(xx, yy, z, alpha=0.2, color='white')


#-----------------------------------------------------------------------------
def findInDic(field, config):
    if field in config:
        return config[field]
    else:
        None

def findTitleInfoString(SIM, type="KMAGmin"):
    ETIME = findInDic("ETIME", SIM.config)
    Dp = findInDic("Dp", SIM.config)
    BY_IMF = findInDic("BY_IMF", SIM.config)
    BZ_IMF = findInDic("BZ_IMF", SIM.config)
    PHI = findInDic("PHI", SIM.config)
    TH = findInDic("TH", SIM.config)
    R = findInDic("R", SIM.config)
    if type == "KMAGmin":
        title = '| Dp: %4.3f | BY_IMF: %3.1f | BZ_IMF: %3.1f' \
         %(Dp, BY_IMF, BZ_IMF)
    if type == "KMAGfull":
        datestring = datetime.datetime.fromtimestamp(ETIME).strftime('%Y-%m-%d')
        title = datestring + ' | Dp: %4.3f | BY_IMF: %3.1f | BZ_IMF: %3.1f' \
         %(Dp, BY_IMF, BZ_IMF)
        if PHI is not None: title += ' | PHI: %.1f' % PHI
        if TH is not None: title += ' | TH: %.1f' % TH
    return title

def colorline(
    ax, x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0),
        linewidth=3, linestyle='-', alpha=1.0, zorder=None):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    import matplotlib.collections as mcoll
    import matplotlib.path as mpath

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, linestyle=linestyle, alpha=alpha, zorder=zorder)

    # ax = plt.gca()
    ax.add_collection(lc)

    return lc
def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


def plotConfiguration(SIMS,
                      ax=None,
                      fieldlines=None,
                      plotData=None,
                      mapGrid=False,
                      plotRsheet=False,
                      ioconfig=None):
    if ax is None:
        newFig = True
        fig, ax = plt.subplots(figsize=(12, 6))
        fig.set_facecolor('black')
    else:
        newFig = False

    if plotData in ['h1', 'h2', 'Xi', 'E', 'b', 'bratio']:
        plotPerturbationVar = plotData
        plotPerturbations = True
        plotData = False
    else:
        plotPerturbations = False

    if ioconfig.plotx is not None:
        xlim = ioconfig.plotx
    else:
        xlim = [-20, 20]
    if ioconfig.ploty is not None:
        ylim = ioconfig.ploty
    else:
        ylim = [-14, 14]

    x = np.linspace(xlim[0], xlim[1], 40)
    z = np.linspace(ylim[0], ylim[1], 40)
    Lmax = ioconfig.Lmax
    Lmin = 3

    PLTHAXIS = 0 # x-axis
    # PLTHAXIS = 1 # y-axis
    PLTVAXIS = 2 # z-axis

    # x = np.linspace(0., 10., 100) # To Compare with Persoon
    # z = np.linspace(-4., 4., 100) # To Compare with Persoon
    DataBounds=[[x[0], x[-1]],[z[0], z[-1]]]
    scale, vmin, vmax, ticks, ticklabels, barlabel = dataDisplayConfig(plotData)

    # fig, ax = plt.subplots()
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax = Axes3D(fig)

    # fig.set_size_inches(10,5)
    # fig.set_size_inches(5,4)
    # plt.style.use('ggplot')
    fontcolor = 'black'
    fontcolor = 'white'
    if fontcolor == 'white':
        plt.style.use('dark_background')
        ax.set_facecolor('black')
    else:
        ax.set_facecolor('xkcd:white')
    # ax.grid(False)
    # ax.w_xaxis.pane.fill = False
    # ax.w_yaxis.pane.fill = False
    # ax.w_zaxis.pane.fill = False
    ax.set_aspect('equal')

    # ax.axhline(y=0, lw=1, ls='--', c='black', alpha=0.3)
    ax.grid(color=fontcolor, linestyle=':', linewidth=1., alpha=0.1)
    # ax.grid(color=fontcolor, linestyle=':', linewidth=1., alpha=0)
    # drawSaturn(ax)
    # print("SIMS[0].BFieldModelName ", SIMS[0].BFieldModelName)

    def polar2car2d(r, th):
        x = r * np.cos(np.deg2rad(th))
        y = r * np.sin(np.deg2rad(th))
        return(x, y)

    def drawLatitudeLine(r, th, DataBounds=None):
        th_display = th if th > 0 else -(th+180)
        if DataBounds:
            r = findMaxRinBox(th, DataBounds)
            text_coord_norm = 1.
        else:
            text_coord_norm = 3.
        x, y = polar2car2d(r, th)
        for yflip in (1, -1):
            ha = 'left' if x < 0 else 'right'
            va = 'bottom' if y*yflip < 0 else 'top'
            ax.text(x/text_coord_norm, yflip * y/text_coord_norm, \
                    r'%d$^{\circ}$' % (yflip * th_display), \
                    color=fontcolor, fontsize=10, alpha=0.5, ha=ha, va=va)
            ax.plot([0., x], [0.,  yflip * y], ls=':', lw=1.2, alpha=0.3, color=fontcolor)

    def findMaxRinBox(th, DataBounds):
        number_of_steps = 100.
        step = np.abs(DataBounds[0][1] - DataBounds[0][0]) / number_of_steps
        # print(step)
        rtemp = 0.
        insidebox = True
        while insidebox:
            x, y = polar2car2d(rtemp, th)
            # print(rtemp)
            if x > DataBounds[0][1] or x < DataBounds[0][0] or \
               y > DataBounds[1][1] or y < DataBounds[1][0]:
                insidebox = False
                rtemp -= step
            else:
                rtemp += step
        return rtemp

    if SIMS[0].BFieldModelName != 'uniform':
        drawSaturn(ax, color=fontcolor, dim='2D')
        # drawSaturn(ax, color='white', dim='3D')

        polarR = 25
        TH_day = np.linspace(15, 60, 4)
        TH_night = np.linspace(- 180 + 15, - 180 + 60, 4)
        if xlim[0] < -2:
            TH_all = np.concatenate((TH_day, TH_night), axis=0)
        else:
            TH_all = TH_day
        for th in TH_all:
            drawLatitudeLine(polarR, th, DataBounds)

    if mapGrid:
        n, B, vA = mapConfigurationGrid(x, z, SIM)
        if plotData == 'n':
            Data = n
        if plotData == 'B':
            Data = B
        if plotData == 'vA':
            Data = vA

        # print('min: ', min(n))
        # print('max: ', max(n))
        Data = np.asarray(Data)
        # Data[np.isnan(Data)] = 0
        Data = Data.transpose() * scale
        # print(Data)
        my_cmap = matplotlib.cm.get_cmap('jet') #'BuPu' | 'magma' | 'cividis'
        my_cmap.set_under('black', alpha=0.)
        # my_cmap.set_over('white')
        my_cmap.set_over('black', alpha=0.)

        # img = ax.matshow(Data)
        img = ax.imshow(Data, cmap=my_cmap, extent=[DataBounds[0][0],DataBounds[0][1],DataBounds[1][0],DataBounds[1][1]], interpolation='none', norm=LogNorm(vmin=vmin, vmax=vmax))
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("rig ht", size="5%", pad=0.05)
        cbar = colorbar(img, cax=cax, ticks=ticks) #| fraction=0.046, pad=0.04
        cbar.set_label(r'Alfven Velocity ($kms^{-1}$)', labelpad=-1.)
        cbar.ax.set_yticklabels(ticklabels)  # vertically oriented colorbar
        cbar.ax.tick_params(labelsize=8) 
        # cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation=45)
    
    fieldlineData = True

    if fieldlines is not None:
        fieldline_polygon_x = []
        fieldline_polygon_z = []
        for n, fieldline in enumerate(fieldlines):
            SIM = SIMS[n]

            color = 'black'
            if 'KMAG' in fieldline.name:
                color = 'red'
            my_cmap = matplotlib.cm.get_cmap('jet')
            my_cmap.set_under(alpha=0.)
            my_cmap.set_over(alpha=0.)

            label = r'L = %.0f $R_S$' % fieldline.L
            label = fieldline.name
            # label = None
            if "KMAG" in fieldline.name:
                color = 'red'
            if "dipole" in fieldline.name:
                color = 'blue'

            if plotData == 'n':
                Data = fieldline.n
            if plotData == 'B':
                Data = fieldline.BT
            if plotData == 'vA':
                Data = fieldline.vA
            # if plotData == 'w':
                # Data = fieldline.vA
            #2D
            my_cmap2 = 'gray'
            vmin2 = 1.
            vmax2 = max(fieldline.traceRTP[:,0])
            plotc = fieldline.traceRTP[:,0]

            # ax.plot(fieldline.traceXYZ[:,PLTHAXIS], fieldline.traceXYZ[:,PLTVAXIS],
                # ls='--', lw=1.1, c='white', alpha=0.4, zorder=1000)
            # ax.plot(fieldline.traceXYZ[:,PLTHAXIS], fieldline.traceXYZ[:,PLTVAXIS],
                # ls='-', lw=1.2, c='white', alpha=1., zorder=1E7)
            # ax.plot(fieldline.traceXYZ[:,0], fieldline.traceXYZ[:,2],
                # ls='-', lw=2.5, c='yellow', alpha=0.7, zorder=1000)

            # FIELD LINE
            if size(fieldline.traceXYZ[:,0]) < 1000:
                sparse = 2
            else:
                sparse = 40
            if SIM.L > Lmin:
                colorline(ax, fieldline.traceXYZ[::sparse,0], fieldline.traceXYZ[::sparse,2], z=plotc[::sparse], cmap=my_cmap2, norm=plt.Normalize(vmin2, vmax2),
                    linewidth=1.2, linestyle='-', alpha=.5, zorder=1E1)

            # ax.scatter(fieldline.traceXYZ[:,0], fieldline.traceXYZ[:,2],
                # s=0.5, marker='_', c=plotc, cmap=my_cmap2, norm=LogNorm(vmin=vmin2, vmax=vmax2), alpha=0.7, zorder=10)
            #3D
            # ax.plot(fieldline.traceXYZ[:,0], fieldline.traceXYZ[:,1], fieldline.traceXYZ[:,2],
            #     ls='--', lw=1, c='white', alpha=1.)
                # marker='.', markersize = 1.2 | 

            nsols = size(SIM.solutions)

            if plotPerturbations and SIM.L < Lmax and SIM.L > Lmin:
                if SIM.solutions:
                    for i, solution in enumerate(SIM.solutions):
                        TH = findInDic("TH", SIM.config)
                        # if TH plot_mode =
                        # plot_mode = nsols-1-SIM.id
                        plot_mode = SIM.id
                        # plot_mode = 2
                        if i == plot_mode:
                            if plotPerturbationVar == 'Xi':
                                E = normalize(solution.xi) * 1.5
                            if plotPerturbationVar == 'E':
                                E = normalize(solution.E) * 2.
                            if plotPerturbationVar == 'b':
                                E = normalize(solution.b) * 1.5
                            if plotPerturbationVar == 'bratio':
                                E = normalize(solution.b / SIM.B) * 1.0
                            if plotPerturbationVar == 'h1':
                                E = normalize(SIM.h1) * 1.5
                            if plotPerturbationVar == 'h2':
                                E = normalize(SIM.h1) * 1.5

                            Efun = interp1d(SIM.z, E, kind='cubic',
                                            fill_value="extrapolate")
                            z_eval = np.linspace(SIM.zlim[0], SIM.zlim[1],
                                                 len(fieldline.B))
                            E = Efun(z_eval)

                            B = fieldline.B
                            R = fieldline.traceXYZ
                            Epos1, Epos2 = perpendicularStep(R, B, step=E)
                            # print(len(Epos))
                            # print(Epos)

                            for j in range(len(B)):
                                if j % 2 == 0:
                                    ax.plot([R[j, PLTHAXIS],
                                             Epos1[j, PLTHAXIS]],
                                            [R[j, PLTVAXIS],
                                             Epos1[j, PLTVAXIS]],
                                            lw=2, c=cxkcd[i],
                                            alpha=0.7)
                else:

                    E = normalize(SIM.h1) * 0.4

                    Efun = interp1d(SIM.z, E, kind='cubic', fill_value="extrapolate")
                    z_eval = np.linspace(SIM.zlim[0], SIM.zlim[1], len(fieldline.B))
                    E = Efun(z_eval)

                    B = fieldline.B
                    R = fieldline.traceXYZ
                    Epos1, Epos2 = perpendicularStep(R, B, step=E)
                    # print(len(Epos))
                    # print(Epos)

                    plotc = Data * scale
                    # norm = matplotlib.colors.Normalize(vmin=10.0, vmax=20.0)
                    my_norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
                    #colormap possible values = viridis, jet, spectral
                    # rgba_color = my_cmap(norm(Data[j]),bytes=True) 

                    # ind = np.where((plotc > vmin) & (plotc < vmax))[0]
                    # colors = pl.cm.jet(np.linspace(0,1,n))
                    # colors = my_cmap(plotc, norm=LogNorm(vmin=vmin, vmax=vmax), alpha=1.0)

                    for j in range(len(B)):
                        if j % 2 == 0:
                            rgba_color = my_cmap(my_norm(Data[j]*scale)) 
                            s = ax.plot([R[j,PLTHAXIS], Epos1[j,PLTHAXIS]], [R[j,PLTVAXIS], Epos1[j,PLTVAXIS]], lw=4, c=rgba_color, alpha=0.9)
                            # ax.plot([R[j,0], Epos1[j,0]], [R[j,2], Epos1[j,2]], lw=2, c=colors[j], alpha=0.5)
                            # s = colorline(ax, [R[j,0], Epos1[j,0]], [R[j,2], Epos1[j,2]], z=plotc, linewidth=2, cmap=my_cmap, norm=plt.Normalize(vmin, vmax), alpha=0.5)

                            # colorline(
                # ax, fieldline.traceXYZ[::40,0], fieldline.traceXYZ[::40,2], z=plotc[::40], cmap=my_cmap2, norm=plt.Normalize(vmin2, vmax2),
                # linewidth=1.0, linestyle='-', alpha=0.5, zorder=1E7)


                        #     ax.plot([R[j,0], Epos1[j,0]], [R[j,2], Epos1[j,2]], lw=2, c=plotc[ind], cmap=my_cmap,
                        # norm=LogNorm(vmin=vmin, vmax=vmax), alpha=1.0, zorder=1E6)
            # ax.plot(Epos1[:,0], Epos1[:,2],
             # ls='--', lw=2, c=cxkcd[i], alpha=0.8)
            # ax.fill_between(Epos1[:,0], yfit, yfit + Epos1[:,2],
                # color='gray', alpha=0.2)


            # fieldlineData = False
            if fieldlineData and SIM.L < Lmax:
            # if fieldlineData:
                if plotData:
                    if plotData == 'w':
                        color = 2 * np.pi / (SIM.solutions[0].roots) / 60.
                        plotc = [color] * len(fieldline.traceXYZ[:,PLTHAXIS])
                    else:
                        plotc = Data * scale
                    #2D
                    # Select points only within the legend range
                    ind = np.where((plotc > vmin) & (plotc < vmax))[0]
                    ind = ind[::10]

                    s = ax.scatter(fieldline.traceXYZ[ind,PLTHAXIS], fieldline.traceXYZ[ind,PLTVAXIS],
                        label = label, s=10, c=plotc[ind], cmap=my_cmap,
                        norm=LogNorm(vmin=vmin, vmax=vmax), alpha=1.0, zorder=1E6)
                        # norm=Normalize(vmin=vmin, vmax=vmax))
                    #3D
                    # s = ax.scatter(fieldline.traceXYZ[:,0], fieldline.traceXYZ[:,1], fieldline.traceXYZ[:,2],
                    #     label = label, s=8, c=plotc, cmap=my_cmap,
                    #     norm=LogNorm(vmin=vmin, vmax=vmax))

            # PLOT HARMONIC MODE REGIONS
            # index=np.mod(n,4*2)
            # if np.mod(index,2) > 0:
            #     fieldline_polygon_x.extend(fieldline.traceXYZ[::-1,0])
            #     fieldline_polygon_z.extend(fieldline.traceXYZ[::-1,2])
            #     alpha = 0.8
            #     fill=True
            #     hatch=None
            #     label = None
            #     if index//2 == 3:
            #         alpha = 1.
            #         fill=False
            #         hatch='///'
            #     if n < 8:
            #         label = 'm='+str(index//2 + 1)
            #     polygon = ax.fill(fieldline_polygon_x, fieldline_polygon_z, color=cxkcd[index//2], alpha=alpha, label=label, edgecolor=None, fill=fill, hatch=hatch, linewidth=0.0)
            #     fieldline_polygon_x = []
            #     fieldline_polygon_z = []
            # else:
            #     fieldline_polygon_x.extend(fieldline.traceXYZ[:,0])
            #     fieldline_polygon_z.extend(fieldline.traceXYZ[:,2])



        if plotData:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.03)
            # cbar = colorbar(s, cax=cax, ticks=ticks) #| fraction=0.046, pad=0.04

            my_norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
            my_cmap = matplotlib.cm.get_cmap('jet')
            # my_cmap = plt.cm.get_cmap('Greys')
            my_cmap.set_under('black')
            my_cmap.set_over('black')

            cbar = mpl.colorbar.ColorbarBase(cax, cmap=my_cmap,
                                norm=my_norm,
                                orientation='vertical', ticks=ticks)

            # cbar.set_label(r'Alfven Velocity / $kms^{-1}$', labelpad=+1))
            cbar.set_label(barlabel, labelpad=-1, fontsize=12)
            cbar.ax.set_yticklabels(ticklabels)  # vertically oriented colorbar
            cbar.ax.tick_params(labelsize=10, direction='in', which='both', width=1.5, length=7, color='white', grid_alpha=0.2)
            cbar.ax.tick_params(labelsize=10, direction='in', which='minor', width=0.7, length=4, color='white', grid_alpha=0.2)
            cbar.ax.ticklabel_format(axis='both')

    if plotRsheet:
        for sheet in SIM.sheet:
            psheet = ax.plot(sheet[0], sheet[2], \
            ls='--', lw=1.2, alpha=0.3, color=fontcolor, label='Plasma Sheet')


    ax.set_xlim(DataBounds[0][0],DataBounds[0][1])
    ax.set_ylim(DataBounds[1][0],DataBounds[1][1])
    ax.set_xlabel(r'x ($R_S$)', labelpad= 0, fontsize=12, color=fontcolor)
    ax.set_ylabel(r'z ($R_S$)', labelpad=-2, fontsize=12, color=fontcolor)
    ax.tick_params(axis='x', colors=fontcolor, direction='in')
    ax.tick_params(axis='y', colors=fontcolor, direction='in')

    ax.spines['bottom'].set_color('0.5')
    ax.spines['top'].set_color('0.5')
    ax.spines['right'].set_color('0.5')
    ax.spines['left'].set_color('0.5')


    # Edit the major and minor ticks of the x and y axes
    # ax.xaxis.set_tick_params(which='major', size=8, width=1, direction='in')
    # ax.xaxis.set_tick_params(which='minor', size=7, width=1, direction='in')
    # ax.yaxis.set_tick_params(which='major', size=8, width=1, direction='in')
    # ax.yaxis.set_tick_params(which='minor', size=7, width=1, direction='in')

    ax.xaxis.set_major_locator(MultipleLocator(2.0))
    ax.yaxis.set_major_locator(MultipleLocator(2.0))
    # ax.get_xaxis().set_visible(False)
    # plt.rcParams.update({'font.size': 20})

    # title = findTitleInfoString(SIM, type="KMAGfull")
    # title = r'Regions of $45 - 75$ minute Field Lines Resonances')
    # ax.set_title(title, fontsize=13, color=fontcolor)
    
   # leg1 = ax.legend(loc="lower left", prop={'size': 10})


    custom_lines2 =[Line2D([0], [0], color=cxkcd[0], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[1], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[2], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[3], ls='-', lw=2.5, alpha=1.)]
    # custom_lines2 =[Line2D([0], [0], color=cxkcd[4], ls='-', lw=2.5, alpha=1.),
                    # Line2D([0], [0], color=cxkcd[3], ls='-', lw=2.5, alpha=1.),
                    # Line2D([0], [0], color=cxkcd[2], ls='-', lw=2.5, alpha=1.),
                    # Line2D([0], [0], color=cxkcd[1], ls='-', lw=2.5, alpha=1.),
                    # Line2D([0], [0], color=cxkcd[0], ls='-', lw=2.5, alpha=1.)]
    # leg1 = ax.legend(custom_lines2, ['m=1', 'm=2', 'm=3', 'm=4'], title="Reference", loc="lower left", prop={'size': 10})
    # leg1 = ax.legend(custom_lines2, ['m=1', 'm=2', 'm=3', 'm=4', 'm=5'], title="Reference", loc="lower left", prop={'size': 10})
    # matplotlib version 3.3.0, you can now directly use the keyword argument labelcolor

    # Add Period and L to legend
    labels = []
    custom_lines = []

    if SIM.solutions:
        for n, fieldline in enumerate(fieldlines):
            SIM = SIMS[n]
            eigs_sorted = []
            Llist, lengthlist, LATlist, eigs_sorted = extractEigenfrequencies(np.atleast_1d(SIM))  
            period_string = 'm = %d,\n$L = $%2.1f $R_S$,\n$T = $%3.0f min' % (n+1, Llist[0], 1./(eigs_sorted[n][0]*1E-3)/60.)
            period_string = 'm = %d' % (n+1)
            labels.append(period_string)
            custom_lines.append(Line2D([0], [0], color=cxkcd[n], ls='-', lw=4.5, alpha=1.))

    # leg = ax.legend(custom_lines2, ['m=1', 'm=2', 'm=3', 'm=4'], handlelength=1.2, loc="upper right", borderaxespad=1.2, prop={'size': 13}, framealpha=0.08, ncol=1)
    leg = ax.legend(custom_lines, labels, handlelength=1.1, labelspacing=1.6, loc="upper right", borderaxespad=1.2, prop={'size': 12}, framealpha=0.01, ncol=1)
    for text in leg.get_texts():
        plt.setp(text, color = color)
    # Color the legend text to match the line color
    for handle, label in zip(leg.legendHandles, leg.texts):
        label.set_color(handle.get_color())
    # ax.add_artist(leg)
    # ax.set_axis_off()
    # ax.axis('equal')

    # plt.tight_layout()
    # ax.set(xlim=(0, 10), ylim=(-4,4))
    # ax.set(xlim=(-15, 15), ylim=(-7,7))
    # ax.set(xlim=(-18, 18), ylim=(-7,7))
    # ax.set(xlim=(-20, 20), ylim=(-5,5), zlim=(-5,5))
    # ax.view_init(elev=10., azim=-90)
    # if ioconfig is not None:
    #     if ioconfig.save:
    #         filename = 'fieldlines' + ioconfig.name
    #         fullpath = 'Output/'
    #         fullpath += filename + '_%03d.png' % ioconfig.id
    #         plt.savefig(fullpath)
    #         plt.close(fig)
    #     else:
    #         plt.show()
    #         exit()
    # else:
    #     plt.show()
    #     exit()
    # figconf.render(fig, ioconfig)

    if newFig:
        #  plt.savefig(os.path.join(dirFile,'ProdCountries.pdf'),dpi=300)
        plt.savefig('./Output/figure_field.pdf', dpi=300,
                    transparent=False, bbox_inches='tight')
        # plt.show()
#-----------------------------------------------------------------------------

def mapConfigurationGrid(x, z, SIM):
    xarray = np.atleast_1d(x)
    zarray = np.atleast_1d(z)

    def f(x,z):
        B  = SIM.BFieldModel([x,0.,z], SIM)
        BT = np.linalg.norm(B, axis=-1)
        n  = SIM.densityModel([x, 0., z], SIM)[0]
        vA = calculate_vA(BT, n, SIM.m, SIM)[0]
        # print(n)
        # print(BT)
        # print(vA)
        # exit()
        return[float(n), float(BT), float(vA)]

    m = [[f(x_,z_) for z_ in zarray] for x_ in xarray]
    m = np.asarray(m)
    n  = m[:,:,0]
    BT = m[:,:,1]
    vA = m[:,:,2]

    return np.asarray(n), np.asarray(BT), np.asarray(vA)

#-----------------------------------------------------------------------------

def plotDensitySpace(ax, SIM, param='vA', ioconfig=None, under='black', over='black', alpha=0.7, res=200):
    scale, vmin, vmax, ticks, ticklabels, barlabel = dataDisplayConfig(param)
    if param == 'n': 
        plotvalues = SIM.nfun(SIM.z)
    if param == 'vA':
        plotvalues = SIM.vAfun(SIM.z)

    my_cmap = matplotlib.cm.get_cmap('jet')
    # my_cmap = matplotlib.cm.get_cmap('viridis')
    # my_cmap = plt.cm.get_cmap('Greys')
    my_cmap.set_under(under)
    my_cmap.set_over(over)

    if ioconfig:
        if ioconfig.plotCoord == 'th':
            span = [-(np.rad2deg(SIM.th[0])-90.), -(np.rad2deg(SIM.th[-1])-90.)]

            nfun = interp1d(SIM.th, plotvalues, kind='cubic', fill_value="extrapolate")
            thlist = np.linspace(SIM.th[0], SIM.th[-1], res)
            plotvalues = nfun(thlist) * scale

        if ioconfig.plotCoord == 's':
            span = [SIM.z[0]/SIM.units, SIM.z[-1]/SIM.units]
            plotvalues = plotvalues * scale
    else:
        span = [toPlotCoords(SIM.zlim[0], SIM), toPlotCoords(SIM.zlim[1], SIM)]

    barprops = dict(aspect='auto', cmap=my_cmap, norm=LogNorm(vmin=vmin, vmax=vmax), interpolation='nearest', alpha=alpha)
    im = ax.imshow(plotvalues.reshape((1, -1)), extent=(span[0], span[1], -1, 1), **barprops)
    return im



#-----------------------------------------------------------------------------

def fieldlinePlot3D(fl):
    # 3D FIELD LINE PLOT
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(fl.traceXYZ[:,0], fl.traceXYZ[:,1], fl.traceXYZ[:,2], color='black', label='Field line')
    ax.scatter(fl.traceP1[:,0], fl.traceP1[:,1], fl.traceP1[:,2], color='red', label='P1')
    ax.scatter(fl.traceP2[:,0], fl.traceP2[:,1], fl.traceP2[:,2], color='blue', label='P2')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_aspect(1.0)
    ax.set_xlim(0, 20)
    ax.set_ylim(-10, 10)
    ax.set_zlim(-10, 10)
    plt.legend()
    plt.show()

#-----------------------------------------------------------------------------

def fieldlineScalingFactorPlot(fl, extraData=None, extraDataNames=None):
    fig, ax = plt.subplots()
    fig.set_size_inches(8,6)
    theta = fl.traceRTP[:,1]
    theta = np.rad2deg(theta)
    # print(len(fl.h1_scale))
    # print(len(fl.h2_scale))
    # print(len(fl.h1_scale_th))
    # print(len(fl.h2_scale_th))
    if extraData is not None:
        for i, data in enumerate(extraData):
            ax.plot(theta, normalize(data), ls=':', color='black', label=extraDataNames[i])
    # print('theta', theta)
    # print('fl.h1_scale', fl.h1_scale)
    # print(normalize(fl.h1_scale))
    # print(fl.h2_scale/fl.h2_scale_th)
    # ax.plot(theta, fl.h1_scale, ls=':', color='blue', label='h1 scale factor (model)')
    
    h1 = fl.h1_scale
    h1th = fl.h1_scale_th
    h1 = normalize(h1)
    h1th = normalize(h1th)

    h2 = fl.h2_scale
    h2th = fl.h2_scale_th
    h2 = normalize(h2)
    h2th = normalize(h2th)

    # h1th = np.abs(h1th)
    # h2th = np.abs(h2th)

    ax.plot(theta, h1, ls=':', color='red', lw=2, label='h1 scale factor')
    ax.plot(theta, h1th, ls='--', color='red', lw=2, label='h1 scale factor (theory)')
    # ax.plot(theta, fl.h1_scale_th, ls='--', color='blue', lw=2, label=r'h1 scale factor (theory) ($r\sin\theta$)')
    ax.plot(theta, h2, ls=':', color='blue', label=r'h2 scale factor ($1/(r B\sin\theta) $)')
    ax.plot(theta, h2th, ls='--', color='blue', label=r'h2 scale factor (theory) $)')
    ax.set_ylabel('Normalized Scaling Factors')
    ax.set_xlabel('Colatitude / deg')
    plt.legend()
    plt.show()

#-----------------------------------------------------------------------------

def plotErrorFunctions(series, SIM, ioconfig):
    figconf = loadFigConfig('erf', SIM)
    # insetConfig = loadFigConfig('insetDensity', SIM)
    figconf.save = ioconfig.save
    figconf.movie = ioconfig.movie
    fig, ax = plt.subplots()
    figconf.setup(fig, ax)
    ax.set_xlim(0., np.max(series.t))
    ax.set_ylim(np.min(series.y), np.max(series.y))
    # ax.plot( series.t, series.y, ls='--', marker='o', c='b', ms=3)
    # for root in rootsOptimized:
        # ax.axvline(x=root, ls='--', c='red')
    #toroidal
    ax.axvline(x=0.1201, ls='--', c='red', alpha=0.7, lw=1.5)
    ax.axvline(x=0.3168, ls='--', c='red', alpha=0.7, lw=1.5)
    ax.axvline(x=0.5122, ls='--', c='red', alpha=0.7, lw=1.5)

    #poloidal
    ax.axvline(x=0.08707484, ls='--', c='blue', alpha=0.7, lw=1.5)
    ax.axvline(x=0.3114820, ls='--', c='blue', alpha=0.7, lw=1.5)
    ax.axvline(x=0.51022397, ls='--', c='blue', alpha=0.7, lw=1.5)


    for i, t in enumerate(series.t):
        figconf.id = i
        ax.set_title('Error Root Function \n Iteration #%d' %i)
        ax.axvline(x=t, ls='--', c='black', alpha=0.7, lw=1.)
        ax.scatter( series.t[i], series.y[i], s=12, c='black', alpha=0.8)
        if ioconfig.movie: figconf.render(fig, ioconfig, hold=True)
    # for solution in solutions:
        # ax.axvline(x=solution.roots[0], ls='--', c='green')
    # ax.scatter( series.zeros, [0]*len(series.zeros), c='r', s=20)
    # fig, ax, datax, datay, insetConfig, titles=None)
    # figconf.plotInset(fig, ax, SIM.z / SIM.units, SIM.vA * 1E-3, insetConfig)
    if ioconfig.movie:
        figconf.processMovie(ioconfig, fps=15, output='erf.mp4')
    else:
        figconf.render(fig, ioconfig)

#

def colorAxisLabels(ax, color, axis='y', spine='left', axisColor=None):
    if axisColor:
        ax.spines[spine].set_color(color)
        ax.tick_params(axis=axis, colors=color)
    if axis == 'x':
        ax.xaxis.label.set_color(color)
    if axis == 'y':
        ax.yaxis.label.set_color(color)
#-----------------------------------------------------------------------------

def plotDensity(SIM, ax, ioconfig):
    ax1color = 'xkcd:sea green'
    ax2color = 'xkcd:light red'
    plt.style.use('dark_background')
    ax.set_facecolor('black')
    ax.grid(color='white', linestyle=':', linewidth=1.5, alpha=0.2)

    B = SIM.Bfun(SIM.z) * 1E9
    vA = SIM.vAfun(SIM.z) * 1E-3
    n = SIM.nfun(SIM.z) * 1E-6
    xaxis = -(np.rad2deg(SIM.th) - 90.)
    # ax.set_xlabel('Colatitude (degrees)', fontsize=8)
    ax.set_xlabel('Latitude (degrees)', fontsize=10)
    ax.set_ylabel(r'Density ($cm^{-3}$)', fontsize=10)
    ax.plot( xaxis, n, lw=2., alpha=0.7, ls='-', label = 'Density', color=ax1color)
    colorAxisLabels(ax, ax1color, axis='y', spine='left', axisColor=True)
    colorAxisLabels(ax, 'white', axis='x', spine='bottom', axisColor=True)
    ax.set_yscale('log')
    ax.set_ylim([0.01, 100])
    ax.set_xlim([-90, 90])

    ax2 = ax.twinx()
    # ax.set_yscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim([1E0, 5*1E5])
    ax2.set_ylabel(r'$v_A$ ($kms^{-1}$)', fontsize=10)
    # ax2.plot( xaxis, B, lw=2., alpha=0.7, ls='-', label = 'B strength', color=ax2color)
    ax2.plot( xaxis, vA, lw=2., alpha=0.7, ls='-', label = r'$v_A$ ($kms^{-1}$)', color=ax2color)
    colorAxisLabels(ax2, ax2color, axis='y', spine='right', axisColor=True)   

def plot_cbar(ax, vmin, vmax, label, ticks=None, ticklabels=None, cmap='jet', under=None, over=None):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.03)
    my_norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
    my_cmap = matplotlib.cm.get_cmap('jet')
    if under is not None: my_cmap.set_under('black')
    if under is not None: my_cmap.set_over('black')
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=my_cmap,
                        norm=my_norm,
                        orientation='vertical', ticks=ticks)
    cbar.set_label(label, labelpad=-1, fontsize=14)
    if ticklabels is not None: cbar.ax.set_yticklabels(ticklabels)  # vertically oriented colorbar
    cbar.ax.tick_params(labelsize=10)

def plotSolutionsPub(SIM, fig, ax, ioconfig, plotVar='E', plotBG='n', norm=True, fontcolor='white', bgcolor='black'):
    solutions = SIM.solutions
    
    if bgcolor == 'black':
        plt.style.use('dark_background')
        ax.set_facecolor(bgcolor)
        fig.set_facecolor(bgcolor)

    # ax.title = 'L = %.1f' % (SIM.L)

    if norm: ax.set_ylim((-1., 1.))

    fancyName = fancyVarName(plotVar)
    ax.set_ylabel('Arbitrary Amplitude, ' + fancyName, labelpad=-50, color=fontcolor, fontsize=11, alpha=0.8)

    ax.tick_params(axis='x', colors=fontcolor)
    ax.tick_params(axis='y', colors=fontcolor)
    ax.spines['bottom'].set_color('0.5')
    ax.spines['top'].set_color('0.5')
    ax.spines['right'].set_color('0.5')
    ax.spines['left'].set_color('0.5')
    ax.grid(b=True, which='major', color=fontcolor, alpha=0.12, linestyle='--')
    ax.grid(b=True, which='minor', color=fontcolor, alpha=0.12, linestyle='--')

    if ioconfig.plotCoord == 'th':
        Rmax_index = np.argmax(SIM.rthetaphi[:,0], axis=0)
        th_plasma_sheet = SIM.rthetaphi[Rmax_index, 1]
        ax.axvline(x=(90.-np.rad2deg(th_plasma_sheet)),      lw=1.5, ls='--', c=fontcolor, alpha=0.5)
        ax.axvline(x=(90.-np.rad2deg(SIM.rthetaphi[ 0, 1])), lw=1.5, ls=':', c=fontcolor, alpha=0.35)
        ax.axvline(x=(90.-np.rad2deg(SIM.rthetaphi[-1, 1])), lw=1.5, ls=':', c=fontcolor, alpha=0.35)
        ax.axhline(y=0, lw=1.5, ls=':', c=fontcolor, alpha=0.35)
    
    if plotBG is not None:
        im = plotDensitySpace(ax, SIM, param=plotBG, ioconfig=ioconfig, under=bgcolor, over=bgcolor, alpha=0.7, res=400)
        scale, vmin, vmax, ticks, ticklabels, barlabel = dataDisplayConfig(plotBG)
        # plot_cbar(ax, vmin, vmax, barlabel, ticks=ticks, ticklabels=ticklabels, cmap='jet', under=bgcolor, over=bgcolor)
    
    xaxis, xlabel, xlim = plotxAxis(ioconfig, SIM)
    if xlim: ax.set_xlim(xlim)
    ax.set_xlabel(xlabel, fontsize=12, color=fontcolor)


    eigs_sorted = []
    if SIM.solutions:
        Llist, lengthlist, LATlist, eigs_sorted = extractEigenfrequencies(np.atleast_1d(SIM))  
    # title = 'Eigenperiods (min) | ' + ' | '.join(['T%d = %7.1f' % (n+1, 1./(eigs_sorted[n][0]*1E-3)/60.) for n in range(len(eigs_sorted)) ])

    nsol = 0
    for nsol, solution in enumerate(solutions):
        y = selectSolutionVar(plotVar, solution, SIM, norm=norm)
        period_string = '%3.0f min' % (1./(eigs_sorted[nsol][0]*1E-3)/60.)
        label = solution.name + ', ' + period_string

        # A magnification to make 'b' look clearer
        if plotVar == 'b': y = y * 2

        # Plot each higher mode with a thinner line
        lw = 6.5 - nsol * 1.4
        # Dashed lines
        # ls = '-' if nsol < 2 else (0, (5,0.5))
        ls = '-'
        alpha = 0.9
        # if nsol == 3:
            # dtheta = calculateSecondGradient(solution, SIM, norm=True)
            # vA = SIM.vAfun(SIM.z) * 1E-3
            # vA2 = normalize(np.power(vA, 2))
            # xi = np.divide(1., solution.xi)
            # rc_calc = normalize(np.divide(np.power(vA, 2), solution.xi))
            # rc = np.divide(1, dtheta)
            # rc = dtheta
            # ax.plot(xaxis, rc, lw=2, color=cxkcd[0], label=r'R_C')
            # ax.plot(xaxis, rc_calc, lw=2, color=cxkcd[1], label=r'v_A^2/xi')

        ax.plot(xaxis, y, lw=lw, alpha=alpha, ls=ls, label=label, color=cxkcd[nsol],
            path_effects=[pe.Stroke(linewidth=lw+0.7, foreground='black', alpha=1.), pe.Normal()])

    # Plot legend
    legloc = "upper right" if plotVar == 'E' else "lower right"
    leg = ax.legend(handlelength=1.2, loc=legloc, prop={'size': 10}, framealpha=0.08)

    # Change the color of the legend text
    for text in leg.get_texts():
        plt.setp(text, color = fontcolor)

    # Color the legend text to match the line color
    for handle, label in zip(leg.legendHandles, leg.texts):
        label.set_color(handle.get_color())

    # fig.tight_layout()
    # figconf.render(fig, ioconfig)

def plotFluxTubeParameters(SIM, ax, ioconfig, plotVar='n'):
    color = "black"
    color = "white"

    ax.set_facecolor('black')

    B = SIM.Bfun(SIM.z) * 1E9
    n = SIM.nfun(SIM.z) * 1E-6
    vA = SIM.vAfun(SIM.z) * 1E-3

    solutions = SIM.solutions
    solution = solutions[1]

    # xaxis = -(np.rad2deg(SIM.th) - 90.)
    xaxis = SIM.z / SIM.units

    n = np.power(n, 1./2)
    yaxis1 = normalize(n)
    # yaxis = solution.E
    # yaxis1 = normalize(solution.b)
    # yaxis = normalize(solution.b/SIM.B)
    # yaxis3 = normalize(B)
    # yaxis3 = normalize(np.power(vA, -1.))
    yaxis2 = normalize(np.power(solutions[0].b, 2.) / B)
    yaxis3 = normalize(np.power(solutions[1].b, 2.) / B)
    yaxis4 = normalize(np.power(solutions[2].b, 2.) / B)
    # yaxis2 = solution.y

            # path_effects=[pe.Stroke(linewidth=3.3, foreground='black'), pe.Normal()])
    ax.plot( xaxis, yaxis1, lw=3.5, alpha=0.8, ls='-', label=r'$n^{1/2}$', color=cxkcd[4])
    ax.plot( xaxis, yaxis2, lw=2.5, alpha=0.8, ls='--', label=r'$b_{\perp m=1}^2 / B$', color=cxkcd[0])
    ax.plot( xaxis, yaxis3, lw=2.5, alpha=0.8, ls='--', label=r'$b_{\perp m=2}^2 / B$', color=cxkcd[1])
    ax.plot( xaxis, yaxis4, lw=2.5, alpha=0.8, ls='--', label=r'$b_{\perp m=3}^2 / B$', color=cxkcd[2])
    ax.legend()

    ax.set_xlabel(r'Length ($R_{\it S}$)', color=color, fontsize=14)
    ax.grid(color=color, linestyle=':', linewidth=1., alpha=0.2)
    # ax.set_title(title, fontsize=14, alpha=0.8, color=color)
    ax.spines['bottom'].set_color(color)
    ax.spines['left'].set_color(color)
    ax.spines['right'].set_color(color)
    ax.tick_params(axis='x', which='both', colors=color)
    ax.tick_params(axis='y', colors=color)


def plotEigenfrequenciesSubplots(SIMS, ioconfig):
    # fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.5))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4), sharey=True)
    # fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    # fig = plt.figure(num=None, figsize=(8.5, 4.5), dpi=80, facecolor='black', edgecolor='k')
    # fig = plt.figure(figsize=(7.2, 2), dpi=300, facecolor='black')

    # ax1 = plt.subplot(121)
    # ax2 = plt.subplot(122)
    ax1 = axes[0]
    ax2 = axes[1]

    plotEigenfrequenciesPub(SIMS, ax1, ioconfig)
    plotEigenfrequenciesPub(SIMS, ax2, ioconfig)
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    fig.subplots_adjust(hspace=0)
    # plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    plt.setp([a.get_yticklabels() for a in fig.axes[:-1]], visible=False)

    # scale, vmin, vmax, ticks, ticklabels, barlabel = dataDisplayConfig(plotBG)

    # my_norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
    # my_cmap = matplotlib.cm.get_cmap('jet')
    # if under is not None: my_cmap.set_under('black')
    # if under is not None: my_cmap.set_over('black')
    # cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat], pad=0.005)
    # # plt.colorbar(im, cax=cax, **kw)

    # cbar = mpl.colorbar.ColorbarBase(cax, cmap=my_cmap,
    # # cbar = mpl.colorbar.ColorbarBase(cax, cmap=my_cmap,
    #                     norm=my_norm,
    #                     orientation='vertical', ticks=ticks)
    # cbar.set_label(barlabel, labelpad=-1, fontsize=14)
    # if ticklabels is not None: cbar.ax.set_yticklabels(ticklabels)  # vertically oriented colorbar
    # cbar.ax.tick_params(labelsize=10)

    if ioconfig.path is None: fullpath = '' 
    if ioconfig.path is not None: fullpath = ioconfig.path + '/'
    fullpath += ioconfig.name + '_%03d' % ioconfig.id + ioconfig.format
    plt.savefig('./Output/figure.pdf', dpi=300, transparent=False, bbox_inches='tight')
    # plt.show()
    # plt.savefig(fullpath)

def plotSolutionSubplots(SIMS, fieldlines, ioconfig, plotBG='n', under=None):
    # fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.5))
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5), sharex=True, sharey=True)
    # fig, axes = plt.subplots(2, 1, figsize=(7.2, 10), sharex=True, sharey=True)
    # fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    SIM = SIMS[0]
    # fig = plt.figure(num=None, figsize=(8.5, 4.5), dpi=80, facecolor='black', edgecolor='k')
    # fig = plt.figure(figsize=(7.2, 2), dpi=300, facecolor='black')

    # ax1 = plt.subplot(121)
    # ax2 = plt.subplot(122)
    ax1 = axes[0]
    ax2 = axes[1]

    plotSolutionsPub(SIMS[0], fig, ax1, ioconfig, plotVar='E', plotBG=plotBG)

    plotSolutionsPub(SIMS[0], fig, ax2, ioconfig, plotVar='b', plotBG=plotBG)

    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    scale, vmin, vmax, ticks, ticklabels, barlabel = dataDisplayConfig(plotBG)

    my_norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
    my_cmap = matplotlib.cm.get_cmap('jet')
    if under is not None: my_cmap.set_under('black')
    if under is not None: my_cmap.set_over('black')
    cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat], pad=0.005)
    # plt.colorbar(im, cax=cax, **kw)

    cbar = mpl.colorbar.ColorbarBase(cax, cmap=my_cmap,
    # cbar = mpl.colorbar.ColorbarBase(cax, cmap=my_cmap,
                        norm=my_norm,
                        orientation='vertical', ticks=ticks)
    cbar.set_label(barlabel, labelpad=-1, fontsize=14)
    if ticklabels is not None: cbar.ax.set_yticklabels(ticklabels)  # vertically oriented colorbar
    cbar.ax.tick_params(labelsize=10)

    if ioconfig.path is None: fullpath = '' 
    if ioconfig.path is not None: fullpath = ioconfig.path + '/'
    fullpath += ioconfig.name + '_%03d' % ioconfig.id + ioconfig.format
    plt.savefig('./Output/figure.pdf', dpi=300, transparent=False, bbox_inches='tight')
    plt.show()
    # plt.savefig(fullpath)

def plotConfigSubplots(SIMS, fieldlines, ioconfig, under=None):
    # fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.5))
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 8), sharex=False, sharey=False)
    fig.set_facecolor('black')
    ax1 = axes[0]
    ax2 = axes[1]


    plotConfiguration(SIMS, ax1, fieldlines=fieldlines, plotData='n', ioconfig=ioconfig, plotRsheet=True)
    plotConfiguration(SIMS, ax2, fieldlines=fieldlines, plotData='vA', ioconfig=ioconfig, plotRsheet=True)
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    fig.subplots_adjust(hspace=-0.1)
    # plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    if ioconfig.path is None: fullpath = '' 
    if ioconfig.path is not None: fullpath = ioconfig.path + '/'
    fullpath += ioconfig.name + '_%03d' % ioconfig.id + ioconfig.format
    plt.savefig('figure.pdf', dpi=300, transparent=False, bbox_inches='tight')
    plt.show()
#-----------------------------------------------------------------------------

def plotSubplots(SIMS, fieldlines, ioconfig):
    # fig, axes = plt.subplots(2, 2, figsize=(6,8))
    SIM = SIMS[0]
    plt.close('all')
    # fig = plt.figure(num=None, figsize=(8.5, 4.5), dpi=80, facecolor='black', edgecolor='k')
    fig = plt.figure(num=None, facecolor='black')

    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    plotSolutionsPub(SIMS[0], fig, ax3, ioconfig, plotVar='E')

    plotSolutionsPub(SIMS[0], fig, ax4, ioconfig, plotVar='b')

    plotConfiguration(SIMS, ax1, fieldlines=fieldlines, plotData='bratio', ioconfig=ioconfig)

    plotEigenfrequenciesPub(SIMS, ax2, ioconfig)

    # plotDensity(SIMS[0], ax4, ioconfig)
    B = SIM.Bfun(SIM.z) * 1E9
    n = SIM.nfun(SIM.z) * 1E-6
    vA = SIM.vAfun(SIM.z) * 1E-3
    dz = SIM.z[1] - SIM.z[0]
    T = dz / (vA * 1E3)
    T=np.cumsum(T)
    Tin = 0.
    if any(n > 0.1):
        nindices = np.where(n > 0.1)
        nindices0 = nindices[0][0]
        nindices1 = nindices[0][-1]
        T0 = T[nindices0]
        T1 = T[nindices1]
        Tin = T1 - T0
    Tout = T[-1] - Tin
    ax5 = fig.add_axes([0.1, 0.04, 0.2, 0.08])
    plotSmallAxes(ax5, [Tin/60., Tout/60.])

    ax6 = fig.add_axes([0.35, 0.04, 0.2, 0.08])
    plotDensity(SIM, ax6, ioconfig)

    title = findTitleInfoString(SIMS[0], type="KMAGfull")
    fig.suptitle(title, fontsize=14, color='white', x=0.8, y=0.08)
    fig.tight_layout()
    # plt.subplots_adjust(bottom=0.15, wspace=0.3, hspace=0.2)
    plt.subplots_adjust(bottom=0.2)
    # fig.tight_layout(pad=0)
    # plt.show()
    # plt.savefig('Output/pic.png')
    # exit()
    # figconf.render(fig, ioconfig)
    if ioconfig.path is None: fullpath = '' 
    if ioconfig.path is not None: fullpath = ioconfig.path + '/'
    fullpath += ioconfig.name + '_%03d' % ioconfig.id + ioconfig.format
    plt.savefig(fullpath)


#-----------------------------------------------------------------------------

def plotEigenfrequenciesPub(SIMS,
                            ioconfig,
                            ax=None):
    if ax is None:
        newFig = True
        fig, ax = plt.subplots(figsize=(7.2, 4))
        fig.set_facecolor('black')
        ax.set_facecolor('black')
    else:
        newFig = False

    phi = np.rad2deg(SIMS[0].rthetaphi[0, 2])
    if (phi - 0.) < 10:
        meridian = 'day'
    elif (phi - 180.) < 10:
        meridian = 'night'
    elif (phi - 90.) < 10:
        meridian = 'dusk'
    elif (phi - 270.) < 10:
        meridian = 'dawn'

    # if SIMS[0].rthetaphi[0,2]
    # meridian

    eigfList = []
    eigs_sorted = []
    plotReference = ioconfig.plotRefEigs
    tanh
    Llist = []
    lengthlist = []
    LATlist = []
    SIMS = np.atleast_1d(SIMS)
    SIM = SIMS[0]
    SIM.zlim = [0., 10.]

    if SIM.solutions:
        Llist, lengthlist, LATlist, eigs_sorted = extractEigenfrequencies(SIMS)  

    # print("lengthlist", lengthlist)
    lengthlist = np.asarray(lengthlist) / 2.
    figconf = loadFigConfig('wL', SIM)
    # plt.rcParams.update({'font.size': 20})
    insetConfig = loadFigConfig('insetDensity', SIM)
    figconf.save = ioconfig.save
    figconf.movie = ioconfig.movie
    figconf.id = SIM.id
    figconf.id = ioconfig.id
    # fig, ax = plt.subplots()
    color='black'
    color='white'

    # plt.style.use('ggplot')
    if color == 'black':
        ax.set_facecolor('xkcd:white')
    else:
        # plt.set_facecolor('black')
        ax.set_facecolor('black')
    # fig.set_facecolor('black')
    ax.set_facecolor('black')
    ax.axhline(y=0, lw=1, ls='--', c=color, alpha=0.3)
    
    if ioconfig.plotCoord == 's':
        plotAgainst='L'
        text_xpos = 18.0
        plotX = lengthlist
        plotX = Llist
        ax.set_xlabel(r'L / ($R_{\it S}$)', fontsize=14, color=color)
        # ax.set_xlabel(r'Length ($R_{\it S}$)', color=color, fontsize=14))
        ax.set_xlim([0, 42.])
        ax.set_xlim([0, 30.])
        ax.set_xlim([0, 21.])
    if ioconfig.plotCoord == 'th':
        plotAgainst='latitude'
        text_xpos = 74.80
        plotX = LATlist
        # ax.set_xlabel(r'Colatitude / degrees', fontsize=13))
        ax.set_xlabel(r'Invariant Latitude (degrees)', fontsize=14, labelpad=-0.5)
        ax.set_xlim([64,76.5])

    # Highlight the first set
    # if SIM.solutions:    
        # ax.axvline(x=plotX[0], lw=2, ls='--', c=color, alpha=0.7)

    ax.grid(color=color, linestyle=':', linewidth=1., alpha=0.2)

    # if "Dp" in SIMS[0].config:
        # title = 'Dp: %4.3f, BY_IMF: %3.1f, BZ_IMF: %3.1f' %(SIM.config["Dp"], SIM.config["BY_IMF"], SIM.config["BZ_IMF"])
    # if "PHI" in SIMS[0].config:
        # title += ' PHI: %.1f' % SIMS[0].config["PHI"]
    # if "PHI" in SIMS[0].config["fieldlines"]:
        # title += ' PHI: %.1f' % SIMS[0].config["base"]["PHI"]
    if meridian=='day':
        title = 'Day-side Field Line Eigenmodes'
    if meridian=='night':
        title = 'Night-side Field Line Eigenmodes'
    # title = 'Dusk-side Field Line Eigenmodes'
    # title = 'Dawn-side Field Line Eigenmodes'
    # print(1./(eigs_sorted[0]*1E-3)/60.)
    title = 'Eigenperiods (min) | ' + ' | '.join(['T%d = %7.1f' % (n+1, 1./(eigs_sorted[n][0]*1E-3)/60.) for n in range(len(eigs_sorted)) ])
    # title = 'Eigenperiods (min) | T1 = %6.1f | T2 = %6.1f | T3 = %6.1f | T4 = %6.1f ' %(1./(eigs_sorted[0][0]*1E-3)/60., 1./(eigs_sorted[1][0]*1E-3)/60., 1./(eigs_sorted[2][0]*1E-3)/60.)
    # ax.set_title(title, fontsize=14, alpha=0.8, color=color)

    if SIM.solutions:
        plotEigenfrequenciesFieldLineLength(ax, plotX, eigs_sorted, plotAgainst=plotAgainst, plotWeight=2.6, plotAlpha=0.9, plotMarker='o', markerSize=7)

    # Plot Precomputed Eigenfrequencies
    plotYatesReference = False
    if plotYatesReference:
        yatesEvenMode = np.loadtxt('even_mode.csv', skiprows=0)
        yatesOddMode = np.loadtxt('odd_mode.csv', skiprows=0)
        # ax2 = ax.twiny()
        # ax2.set_xlabel(r'Length /$R_{\it S}$', color='grey')
        ax.plot(yatesEvenMode[:,0], yatesEvenMode[:,1], label='m=1 (Yates et al., 2016)', color='tab:grey', alpha=0.4)
        ax.plot(yatesOddMode[:,0], yatesOddMode[:,1], label='m=2 (Yates et al., 2016)', color='tab:grey', alpha=0.4)

    if plotReference:
        if meridian=='night':
            filename1 = "ToroidalNight"
            filename2 = "PoloidalNight"
        if meridian=='day':
            filename1 = "ToroidalDay"
            filename2 = "PoloidalDay"
        L1, th1, w1 = loadEigenfrequencies(filename=filename1)
        L2, th2, w2 = loadEigenfrequencies(filename=filename2)
        plotEigenfrequenciesFieldLineLength(ax, th1,  w1, plotAgainst=plotAgainst, plotStyle='-', plotAlpha=0.6, plotWeight=2.8, labelPrefix='Toroidal Day | ')
        plotEigenfrequenciesFieldLineLength(ax, th2,  w2, plotAgainst=plotAgainst, plotStyle='--', plotAlpha=0.6, plotWeight=2., labelPrefix='Poloidal Day | ')

    # Fixed for reference
    ax.set_yscale('log')
    ax.set_ylim([0.02,1.])

    ax.spines['bottom'].set_color(color)
    ax.spines['left'].set_color(color)
    ax.spines['right'].set_color(color)
    ax.tick_params(axis='x', which='both', colors=color, size=6, labelsize=12)
    ax.tick_params(axis='y', colors=color, size=6)
    ax.tick_params(axis='y', colors=color, size=4, which='minor', pad=-10)
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(1.))
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))

    ax.set_ylabel('Eigenfrequency (mHz)', color=color, fontsize=14, labelpad=-1)


    ax2 = ax.twinx()
    # text_xpos = max(xvalues)/10.*9.
    # ax2.axhline(y=45,ls=':', lw=1.6,  color='xkcd:teal', alpha=0.9, label='45 min')
    # ax2.text(text_xpos, 45-2, '45 min', fontsize=11, color='xkcd:teal', alpha=0.8)
    # ax2.axhline(y=75,ls=':', lw=1.6,  color='xkcd:teal', alpha=0.9, label='75 min')
    # ax2.text(text_xpos, 75-2, '75 min', fontsize=11, color='xkcd:teal', alpha=0.8)
    # ax2.axhspan(45, 75, alpha=0.1, color='xkcd:teal')    

    line_alpha = 0.6
    ax2.axhline(y=60,ls=':', lw=1.6,  color=color, alpha=line_alpha, label='60 min')
    ax2.text(text_xpos, 60-3, r'$\mathbf{60}$ min', fontsize=12, color=color, alpha=0.8)
    ax2.axhline(y=30,ls=':', lw=1.6,  color=color, alpha=line_alpha, label='30 min')
    ax2.text(text_xpos, 30-1, r'$\mathbf{30}$ min', fontsize=12, color=color, alpha=0.8)
    ax2.axhline(y=90,ls=':', lw=1.6,  color=color, alpha=line_alpha, label='90 min')
    ax2.text(text_xpos, 90-4, r'$\mathbf{90}$ min', fontsize=12, color=color, alpha=0.8)
    ax2.axhline(y=120,ls=':', lw=1.6,  color=color, alpha=line_alpha/2, label='120 min')
    ax2.text(text_xpos, 120-4, r'$\mathbf{2}$ h', fontsize=12, color=color, alpha=0.8)
    # ax2.axhline(y=360,ls=':', lw=1.6,  color=color, alpha=line_alpha/2, label='360 min')
    # ax2.text(text_xpos, 360-4, '6 h', fontsize=12, color=color, alpha=0.8)
    ax2.axhline(y=300,ls=':', lw=1.6,  color=color, alpha=line_alpha/2, label='300 min')
    ax2.text(text_xpos, 300-10, r'$\mathbf{5}$ h', fontsize=12, color=color, alpha=0.8)
    ax2.axhline(y=600,ls=':', lw=1.6,  color=color, alpha=line_alpha/2, label='600 min')
    ax2.text(text_xpos, 600-20, r'$\mathbf{10}$ h', fontsize=12, color=color, alpha=0.8)
    ax2.set_ylim([1/(1E-3)/60., 1/(0.02E-3)/60.])
    ax2.set_ylabel('Eigenperiod (mins)', fontsize=12, color=color)
    ax2.set_yscale('log')
    ax2.invert_yaxis()
    ax2.tick_params(axis='y', colors=color, size=4, direction='in', which='minor')
    ax2.tick_params(axis='y', colors=color, size=8, direction='in', which='major', width=1., pad=-20)
    ax2.set_yticklabels([])
    ax2.tick_params(axis='x', colors=color)
    ax2.spines['right'].set_color(color)
    ax2.spines['bottom'].set_color(color)
    ax2.spines['left'].set_color(color)
    ax2.spines['top'].set_color(color)
    ax.tick_params(axis='y', colors=color)
    ax.tick_params(axis='x', colors=color)
    leg1 = ax.legend(loc="lower left", prop={'size': 10})


    custom_lines2 =[Line2D([0], [0], color=cxkcd[5], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[4], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[3], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[2], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[1], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[0], ls='-', lw=2.5, alpha=1.)]
    custom_lines2 =[Line2D([0], [0], color=cxkcd[3], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[2], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[1], ls='-', lw=2.5, alpha=1.),
                    Line2D([0], [0], color=cxkcd[0], ls='-', lw=2.5, alpha=1.)]
    # leg1 = ax.legend(custom_lines2, ['m=1', 'm=2', 'm=3', 'm=4', 'm=5', 'm=6'], title="Reference", loc="lower left", prop={'size': 10})
    # matplotlib version 3.3.0, you can now directly use the keyword argument labelcolor

    # ax.legend(handlelength=1.2, loc=legloc, prop={'size': 10}, framealpha=0.08)

    leg = ax.legend(custom_lines2, ['m=4', 'm=3', 'm=2', 'm=1'], handlelength=1.2, loc="lower left", prop={'size': 14}, markerscale=2, facecolor='black', framealpha=0.6)
    # leg = ax.legend(custom_lines2, ['m=6', 'm=5', 'm=4', 'm=3', 'm=2', 'm=1'], handlelength=1.2, loc="lower left", prop={'size': 14}, markerscale=2, facecolor='black', framealpha=0.6)
    for text in leg.get_texts():
        plt.setp(text, color = 'white')

    # Color the legend text to match the line color
    for handle, label in zip(leg.legendHandles, leg.texts):
        label.set_color(handle.get_color())

    ax.add_artist(leg)

    # cmap = plt.cm.coolwarm
    # color=cmap(0.)
    custom_lines = [Line2D([0], [0], color=color, ls='-', lw=2., alpha=1.),
                    Line2D([0], [0], color=color, ls='--',  lw=2., alpha=0.5)]
    # leg2 = ax.legend(custom_lines, ['Toroidal Mode', 'Poloidal Mode'], loc="upper center", prop={'size': 11}, framealpha=0.2)
    # for text in leg2.get_texts():
        # plt.setp(text, color = color)

    # ax.xaxis.grid(True, which='minor')
    # ax.grid()
    # ax.grid(which='major', linestyle='--', linewidth='2.', color='black')
    # ax.grid(which='minor', linestyle=':', linewidth='0.2', color='black')
    # ax2.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
    
    # Bottom Second Axis for Equatorial Crossing Distance
    # ax3 = ax.twiny()
    # ax3.set_xlim([64,76.5])

    # # Add some extra space for the second axis at the bottom
    # plt.subplots_adjust(bottom=0.2)

    # # Move twinned axis ticks and label from top to bottom
    # ax3.xaxis.set_ticks_position("bottom")
    # ax3.xaxis.set_label_position("bottom")

    # # Offset the twin axis below the host
    # ax3.spines["bottom"].set_position(("axes", -0.15))

    # # Turn on the frame for the twin axis, but then hide all 
    # # but the bottom spine
    # ax3.set_frame_on(True)
    # ax3.patch.set_visible(False)

    # # as @ali14 pointed out, for python3, use this
    # for sp in ax3.spines.values():
    # # and for python2, use this
    # # for sp in ax3.spines.itervalues():
    #     sp.set_visible(False)
    # ax3.spines["bottom"].set_visible(True)

    # #Night
    # if meridian=='night':
    #     LAT_labels = [74.73, 74.08, 73.35, 71.86, 66.44, 63.59]
    #     L_labels = [30,20,15,10,5,4]
    # #Day
    # if meridian=='day':
    #     LAT_labels = [75.63,74.08,71.97,66.32,63.59]
    #     L_labels = [20,15,10,5 ,4]
    # ax3.set_xticks(LAT_labels)
    # ax3.tick_params(axis='x', colors=color)
    # ax3.spines['bottom'].set_color(color)
    # ax3.set_xticklabels(L_labels)
    # ax3.set_xlabel(r"Equatorial Crossing Distance  ($R_S$)", fontsize=14, color=color)
    # plt.tight_layout()


    # figconf.render(fig, ioconfig)
    if newFig:
        # plt.show()
        plt.savefig('./Output/figure_eigfreq.pdf', dpi=300, transparent=False,
                    bbox_inches='tight', facecolor=fig.get_facecolor(),
                    edgecolor='none')
        plt.close(fig)
#-----------------------------------------------------------------------------

def plotCummingsModes(SIMS, ioconfig):
    # for i, SIM in enumerate(SIMS):
        # Llist.append(SIM.L)
        # LATlist.append(90. - todeg(SIM.rthetaphi[0,1]))
    for SIM in SIMS:
        solutions = SIM.solutions
        # print('coords: ', SIM.coords)
        # print('zlim: ', SIM.zlim)
        roots = []
        if SIM.component == 'poloidal':
            color = 'blue'
        else:
            color = 'red'

        marker = 'x'
        if SIM.coords == 'ds':
            # scalefactor = 0.223
            marker = '+'

        for solution in solutions:
            roots.append(angularFreqConversion(solution.roots, format='f'))
        print('m mode: ', SIM.m)
        # print('component: ', SIM.component)
        # print('roots: ', roots)
        for i, root in enumerate(roots):
            plot(SIM.m, root, color = 'green', ls='--', label='Toroidal h=%s, Numerical' %str(i+1))
            scatter(SIM.m, root, color = color, marker=marker, s=50.)

    # if ncomp == 0: roots_toroidal.append(roots / 2. / np.pi)
    # if ncomp == 1: roots_poloidal.append(roots / 2. / np.pi)
    
    # roots_toroidal_sorted = sortEigenmodes(roots_toroidal, SIM.modes)
    # roots_poloidal_sorted = sortEigenmodes(roots_poloidal, SIM.modes)
    # plot(mList, roots_toroidal_sorted[0], color = 'red', ls='--', label='Toroidal h=1, Numerical')
    # plot(mList, roots_toroidal_sorted[1], color = 'red', ls='--', label='Toroidal h=2, Numerical', alpha=0.8)
    # plot(mList, roots_toroidal_sorted[2], color = 'red', ls='--', label='Toroidal h=3, Numerical', alpha=0.6)
    # plot(mList, roots_toroidal_sorted[3], color = 'red', ls='--', label='Toroidal h=4, Numerical', alpha=0.4)
    # plot(mList, roots_poloidal_sorted[0], color = 'blue',ls='--', label='Toroidal h=1, Numerical')
    # plot(mList, roots_poloidal_sorted[1], color = 'blue',ls='--', label='Toroidal h=2, Numerical', alpha=0.8)
    # plot(mList, roots_poloidal_sorted[2], color = 'blue',ls='--', label='Toroidal h=3, Numerical', alpha=0.6)
    # plot(mList, roots_poloidal_sorted[3], color = 'blue',ls='--', label='Toroidal h=4, Numerical', alpha=0.4)

    # Plot density index m: [0, 1, 2, 3, 4, 5] against Hz (cm-3 density)
    scatter([0, 1, 2, 3, 4, 5], [0.019, 0.019, 0.018, 0.018, 0.017, 0.016], label='Toroidal h=1, Cummings et al. (1969)', color='red', alpha=0.2)
    scatter([0, 1, 2, 3, 4, 5], [0.051, 0.048, 0.045, 0.042, 0.039, 0.035], label='Toroidal h=2, Cummings et al. (1969)', color='red', alpha=0.2)
    scatter([0, 1, 2, 3, 4, 5], [0.082, 0.077, 0.072, 0.066, 0.060, 0.054], label='Toroidal h=3, Cummings et al. (1969)', color='red', alpha=0.2)
    scatter([0, 1, 2, 3, 4, 5], [0.113, 0.106, 0.099, 0.091, 0.082, 0.072], label='Toroidal h=4, Cummings et al. (1969)', color='red', alpha=0.2)

    scatter([0, 1, 2, 3, 4, 5], [0.014, 0.014, 0.013, 0.013, 0.012, 0.012], label='Poloidal h=1, Cummings et al. (1969)', color='blue', alpha=0.2)
    scatter([0, 1, 2, 3, 4, 5], [0.050, 0.047, 0.045, 0.042, 0.038, 0.035], label='Poloidal h=2, Cummings et al. (1969)', color='blue', alpha=0.2)
    scatter([0, 1, 2, 3, 4, 5], [0.082, 0.077, 0.072, 0.066, 0.060, 0.053], label='Poloidal h=3, Cummings et al. (1969)', color='blue', alpha=0.2)
    scatter([0, 1, 2, 3, 4, 5], [0.113, 0.106, 0.099, 0.091, 0.082, 0.072], label='Poloidal h=4, Cummings et al. (1969)', color='blue', alpha=0.2)
    plt.title('Eigenfrequencies of the Uncoupled Toroidal and Poloidal Wave Equations')
    plt.xlabel('Density Index (m)')
    plt.ylabel(r'Frequency, Hz (1 $cm^{-3}$ density)')
    plt.grid(True)
    # plt.legend()
    plt.show()
    exit()

#-----------------------------------------------------------------------------

def plotScaleHeightsAgainstL():
    # PLOT SCALE HEIGHTS AGAINST L
    x = np.linspace(0., 15, 100)
    x = np.linspace(0., 40, 100)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    R = np.asarray([x,y,z]).T
    H = fH_bagenal(np.abs(R.T[0]))
    H2 = 0.047*np.power(np.abs(R.T[0]), 1.8) 
    H2 = 0.040*np.power(np.abs(R.T[0]), 1.8) 
    H2 = fH_persoon(np.abs(R.T[0]))
    plt.figure(figsize=(6,4))
    plt.ticklabel_format(style='plain',useOffset=False, axis='x')
    plt.scatter(H_bagenal[0], H_bagenal[1], label='Bagenal & Delamere [2011]', s=10, color=cxkcd[1])
    # plt.scatter(H_persoon[0], H_persoon[1], label='Persoon et al. [2013]', s=10, color=cxkcd[0])
    plt.plot(x, H, label='Bagenal Extrapolated', color=cxkcd[1], lw=2)
    # plt.plot(x, H2, label='Persoon et al. [2013], $H=0.047 L^{1.5}$', color=cxkcd[0], lw=2)
    plt.axhline(y=5.0, color='black', alpha=0.5, ls='--')
    plt.text(20., 5.0, '5 $R_S$', fontsize=9, color='black', alpha=0.8)
    plt.axhline(y=2.0, color='black', alpha=0.5, ls='--')
    plt.text(20., 2.0, '2 $R_S$', fontsize=9, color='black', alpha=0.8)
    plt.axhline(y=1.0, color='black', alpha=0.5, ls='--')
    plt.text(20., 1.0, '1 $R_S$', fontsize=9, color='black', alpha=0.8)
    plt.xlim(1., 15.)
    plt.xlim(1., 28.)
    plt.title('Scale Height (H)')
    plt.xlabel(r'L ($R_S$)')
    plt.ylabel(r'Scale Height, H ($R_S$)')
    plt.yscale('log')
    plt.legend(fontsize=9)
    plt.show()

#-----------------------------------------------------------------------------

def plotVerticalDensityProfiles():
    # PLOT VERTICAL DENSITY PROFILES FOR DIFFERENT L
    x = np.array([4.0, 6.0, 8.0])
    z = np.linspace(-5., 5., 100)
    y = np.zeros_like(z)
    R = [[np.full_like(z,x_), np.zeros_like(z), z] for x_ in x]
    R = np.asarray(R)
    R1 = R[0].T; R2 = R[1].T; R3 = R[2].T
    n1=[]; n2=[]
    for R0 in R:
        R0 = np.asarray(R0).T
        SIM.densityModelName = 'bagenal'
        n1.append(SIM.densityModel(R0, SIM) * 1E-6)
        SIM.densityModelName = 'persoon'
        n2.append(SIM.densityModel(R0, SIM) * 1E-6)
    plt.figure(figsize=(6,3))
    plt.plot(z, n1[0], label='Bagenal L='+str(x[0]), color=cxkcd[0], lw=2)
    plt.plot(z, n2[0], label='Persoon L='+str(x[0]), color=cxkcd[0], lw=2, ls='--')
    plt.plot(z, n1[1], label='Bagenal L='+str(x[1]), color=cxkcd[1], lw=2)
    plt.plot(z, n2[1], label='Persoon L='+str(x[1]), color=cxkcd[1], lw=2, ls='--')
    plt.plot(z, n1[2], label='Bagenal L='+str(x[2]), color=cxkcd[2], lw=2)
    plt.plot(z, n2[2], label='Persoon L='+str(x[2]), color=cxkcd[2], lw=2, ls='--')
    plt.ylim(1E-2, 100.)
    plt.title('Vertical (z) Density Profiles at different L')
    plt.xlabel(r'Vertical Distance, z ($R_S$)')
    plt.ylabel(r'cm$^{-3}$')
    plt.yscale('log')
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.show()

#-----------------------------------------------------------------------------

def plotVerticalvAProfiles():
    # PLOT VERTICAL vA PROFILES FOR DIFFERENT L
    x = np.array([4.0, 6.0, 8.0])
    z = np.linspace(-5., 5., 100)
    y = np.zeros_like(z)
    R = [[np.full_like(z,x_), np.zeros_like(z), z] for x_ in x]
    R = np.asarray(R)
    R1 = R[0].T; R2 = R[1].T; R3 = R[2].T
    vA1=[]; vA2=[]
    for R0 in R:
        R0 = np.asarray(R0).T

        B = SIM.BFieldModel(R0, SIM)
        BT = np.linalg.norm(B, axis=-1)

        SIM.densityModelName = 'bagenal'
        n1 = SIM.densityModel(R0, SIM)
        vA1.append(calculate_vA(BT, n1, SIM.m0)*1E-3)

        SIM.densityModelName = 'persoon'
        n2 = SIM.densityModel(R0, SIM)
        vA2.append(calculate_vA(BT, n2, SIM.m0)*1E-3)
    plt.figure(figsize=(6,3))
    plt.plot(z, vA1[0], label='Bagenal L='+str(x[0]), color=cxkcd[0], lw=2)
    plt.plot(z, vA2[0], label='Persoon L='+str(x[0]), color=cxkcd[0], lw=2, ls='--')
    plt.plot(z, vA1[1], label='Bagenal L='+str(x[1]), color=cxkcd[1], lw=2)
    plt.plot(z, vA2[1], label='Persoon L='+str(x[1]), color=cxkcd[1], lw=2, ls='--')
    plt.plot(z, vA1[2], label='Bagenal L='+str(x[2]), color=cxkcd[2], lw=2)
    plt.plot(z, vA2[2], label='Persoon L='+str(x[2]), color=cxkcd[2], lw=2, ls='--')
    plt.axhline(y=100, color='black', alpha=0.5, ls='--')
    plt.text(4., 100, '100 $km/s$', fontsize=9, color='black', alpha=0.8)
    plt.axhline(y=1000, color='black', alpha=0.5, ls='--')
    plt.text(4., 1000, '1000 $km/s$', fontsize=9, color='black', alpha=0.8)
    # plt.ylim(1E-2, 100.)
    plt.title('Vertical (z) vA Profiles at different L')
    plt.xlabel(r'Vertical Distance, z ($R_S$)')
    plt.ylabel(r'km/s')
    plt.yscale('log')
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.show()

#-----------------------------------------------------------------------------

def plotEquatorialDensityProfile(SIM):
    # PLOT DENSITY PROFILE AT THE EQUATOR
    x = np.linspace(0., 40., 100)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    R = np.asarray([x,y,z]).T
    # SIM.densityModelName = 'bagenal'
    # n = SIM.densityModel(R, SIM) * 1E-6
    n = fn_bagenal(x)
    # n = fH_bagenal(x)
    # print(SIM.densityModel([22., 0., 0.], SIM) * 1E-6)
    # SIM.densityModelName = 'persoon'
    # n2 = SIM.densityModel(R, SIM) * 1E-6
    # n3 = fn_persoon(x)
    plt.figure(figsize=(6,4))
    plt.ticklabel_format(style='plain',useOffset=False, axis='x')
    plt.scatter(n_bagenal[0], n_bagenal[1], label='Bagenal & Delamere [2011]', s=12, color=cxkcd[1])
    # plt.scatter(neq_persoon[0], neq_persoon[1], label='Persoon et al.[2013]', s=12, color=cxkcd[0])
    plt.plot(x, n, label='Bagenal Extrapolated', color=cxkcd[1], lw=2)
    # plt.plot(x, n2, label='Persoon Extrapolated', color=cxkcd[0], lw=2)
    # plt.plot(x, n3, label='Persoon Extrapolated 2', color=cxkcd[2], lw=2)
    plt.axvline(x=3.0, color='black', alpha=0.5, ls='--')
    plt.axvline(x=20.0, color='black', alpha=0.5, ls='--')
    plt.xlim(0., 13.)
    plt.xlim(0., 30.)
    plt.ylim(0.01, 100.)
    # plt.ylim(0., 100.)
    plt.title('Density Profile at the Magnetic Equator (W+)')
    plt.xlabel(r'Radial Distance ($R_S$)')
    plt.ylabel(r'cm$^{-3}$')
    plt.yscale('log')
    plt.legend(fontsize=9)
    plt.show()
    # ax.get_xaxis().get_major_formatter().set_useOffset(False)
    # ax.get_xaxis().get_major_formatter().set_scientific(False)
    # ax.ticklabel_format(style='plain')
    # ax.ticklabel_format(useOffset=False, style='plain')
    # ax.ticklabel_format(useOffset=False)

# ----------------------------------------------------------------------------
def compareFieldLines():
    print('L: ', fieldlines[0].L)
    print('L2: ', fieldlines2[0].L)
    fieldlineScalingFactorPlot(fieldlines[0])
    plotConfiguration(SIM, fax, ieldlines, saveFig=False, id=0)
    plotConfiguration(SIM, [ax, fieldlines[0], fieldlines2[0]], saveFig=False, id=0)
    exit()

# ----------------------------------------------------------------------------
def plotBasicConfiguration(SIM):
    fig, ax = plt.subplots(figsize=(7.2,5))
    plotEquatorialDensityProfile(SIM)
    # plotScaleHeightsAgainstL()
    fig.set_facecolor('black')
    # plotEigenfrequenciesPub(SIM, ax, ioconfig,
                            # plotReference=True, meridian='night')
    # plotEigenfrequenciesSubplots(SIM, ioconfig)
    plt.savefig('figure.pdf', dpi=300,
                transparent=False, bbox_inches='tight')
    plt.show()

#-----------------------------------------------------------------------------

if __name__ == "__main__":
    from pylab import *