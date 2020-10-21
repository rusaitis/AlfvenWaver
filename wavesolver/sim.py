#!/usr/bin/env python3

import numpy as np
from wavesolver.model import *
from wavesolver.shoot import *
from wavesolver.configurations import *
from scipy import integrate
from scipy.interpolate import interp1d
from copy import deepcopy

# ------------------------------------------------------------------------------
#                                Description
#
# This file contains the class for simulation configurations (sim), and
# some basic operations of these configuration objects (e.g. import).
#
# ----------------------------------------------------------------------------


class sim:
    """ A simulation configuration class. Everything needed for standing
    Alfven calculation is put in an object of this class. """
    def __init__(self,
                 fun=None,
                 BC=[0., 0.],
                 dydz=[1., 0.5],
                 wlim=[0., 1.],
                 zlim=None,
                 planet='saturn',
                 densityModel=None,
                 BFieldModel=None,
                 coords='ds',
                 vA=None,
                 B=None,
                 n=None,
                 z=None,
                 zlabel=None,
                 component='toroidal',
                 m=amu,
                 length=0,
                 xyz=None,
                 rthetaphi=None,
                 th=None,
                 vAfun=None,
                 rfun=None,
                 THfun=None,
                 nfun=None,
                 Bfun=None,
                 dBfun=None,
                 L=None,
                 N=0.,
                 mIndex=1.0,
                 tol=1E-7,
                 maxStep=1.,
                 res=None,
                 wRes=100,
                 rtol=1E-5,
                 atol=1E-6,
                 integrationMethod='RK45',
                 modes=6,
                 units=1.,
                 id=0,
                 h1=None,
                 h2=None,
                 h1fun=None,
                 h2fun=None,
                 dh1=None,
                 dh2=None,
                 dh1fun=None,
                 dh2fun=None,
                 dlnh1B=None,
                 dlnh2B=None,
                 dlnh1Bfun=None,
                 dlnh2Bfun=None,
                 name=None,
                 dipM=dipM_Saturn,
                 densityModelName='yates',
                 BFieldModelName='box',
                 solutions=None,
                 traceFrom='equator',
                 planetRadius=[1., 1., 1.],
                 step=0.1,  # tracing step size
                 hstep=1.0,  # step size for adjacent field line tracing
                 sumh1h2=False,  # EXPERIMENTAL: cumalative sum of h1 and h2
                 traceh1h2=False,  # Trace adjacent field lines for h1 & h2
                 maxIter=1E6,
                 directions=[-1., 1.],
                 config={},
                 flconfig=None,
                 sheet=RSHEET,
                 zmax=None,
                 scalingFactor=None):
        self.fun = fun
        self.BC = BC
        self.dydz = dydz
        self.zlim = zlim
        self.wlim = wlim
        self.res = res
        self.wRes = wRes
        self.planet = planet
        self.densityModel = densityModel
        self.densityModelName = densityModelName
        self.BFieldModel = BFieldModel
        self.BFieldModelName = BFieldModelName
        self.vA = vA
        self.n = n
        self.B = B
        self.z = z
        self.xyz = xyz
        self.rthetaphi = rthetaphi
        self.th = th
        self.vAfun = vAfun
        self.nfun = nfun
        self.rfun = rfun
        self.Bfun = Bfun
        self.dBfun = dBfun
        self.THfun = THfun
        self.h1 = h1
        self.h2 = h2
        self.dh1 = dh1
        self.dh2 = dh2
        self.h1fun = h1fun
        self.h2fun = h2fun
        self.dh1fun = dh1fun
        self.dh2fun = dh2fun
        self.dlnh1B = dlnh1B
        self.dlnh2B = dlnh2B
        self.dlnh1Bfun = dlnh1Bfun
        self.dlnh2Bfun = dlnh2Bfun
        self.L = L
        self.N = N
        self.m = m
        self.mIndex = mIndex
        self.component = component
        self.length = length
        self.units = units
        self.coords = coords
        self.zlabel = zlabel
        self.modes = modes
        self.integrationMethod = integrationMethod
        self.id = id
        self.name = name
        self.tol = tol
        self.rtol = rtol
        self.atol = atol
        self.dipM = dipM
        self.planetRadius = planetRadius
        self.traceFrom = traceFrom
        self.step = step
        self.hstep = hstep
        self.sumh1h2 = sumh1h2
        self.traceh1h2 = traceh1h2
        self.maxIter = maxIter
        self.directions = directions
        self.config = config
        self.flconfig = flconfig
        self.sheet = sheet
        self.zmax = zmax
        self.scalingFactor = scalingFactor

        # SOLUTIONS
        self.solutions = solutions
        if vAfun is not None:
            self.interpolate()
        if coords == 'cartesian':
            self.zlabel = r'length / $R_P$'
        elif coords == 'cos':
            self.zlabel = 'colatitude / deg'
        elif coords == 'deg':
            self.zlabel = 'colatitude / deg'
        elif coords == 'ds':
            self.zlabel = r'distance along the field line / $R_P$'
        else:
            self.zlabel = ''

    def integrate(self):
        """ Integration Method """
        if self.res is not None:
            t_eval = np.linspace(self.zlim[0], self.zlim[1], self.res)
        return integrate.solve_ivp(lambda t, y: self.fun(t, y, self),
                                   (self.zlim[0], self.zlim[1]),
                                   [self.BC[0], self.dydz[0], self.wlim[0]],
                                   method=self.integrationMethod,
                                   dense_output=True,
                                   t_eval=t_eval,
                                   rtol=self.tol * 1E-1,
                                   atol=self.tol * 1E-3)

    def interpolate(self):
        """ Establish interpolation functions for variables used in
        PDE's """
        kargs = {"kind": 'cubic', "fill_value": "extrapolate"}
        self.vAfun = interp1d(self.z, self.vA, **kargs)
        self.rfun = interp1d(self.z, self.rthetaphi[:, 0], **kargs)
        self.THfun = interp1d(self.z, self.rthetaphi[:, 1], **kargs)
        self.nfun = interp1d(self.z, self.n, **kargs)
        self.Bfun = interp1d(self.z, self.B, **kargs)
        # self.dBfun = interp1d(self.z, self.dB, **kargs)
        self.h1fun = interp1d(self.z, self.h1, **kargs)
        self.h2fun = interp1d(self.z, self.h2, **kargs)

        stepsize = self.step * self.units
        dh1 = np.gradient(self.h1, stepsize, edge_order=2)
        dh2 = np.gradient(self.h2, stepsize, edge_order=2)
        lnh1B = np.log(np.power(self.h1, 2.) * self.B)
        lnh2B = np.log(np.power(self.h2, 2.) * self.B)
        dlnh1B = np.gradient(lnh1B, stepsize, edge_order=2)
        dlnh2B = np.gradient(lnh2B, stepsize, edge_order=2)
        self.dh1fun = interp1d(self.z, dh1, **kargs)
        self.dh2fun = interp1d(self.z, dh2, **kargs)
        self.dlnh1Bfun = interp1d(self.z, dlnh1B, **kargs)
        self.dlnh2Bfun = interp1d(self.z, dlnh2B, **kargs)

        z_eval = np.linspace(self.zlim[0], self.zlim[1], self.res)
        self.h1 = self.h1fun(z_eval)
        self.h2 = self.h2fun(z_eval)
        self.B = self.Bfun(z_eval)
        self.th = self.THfun(z_eval)
        self.z = z_eval

    def update(self, conf):
        """ Update the SIM configuration with a new dictionary """
        if "densityModelName" in conf:
            self.densityModelName = conf["densityModelName"]
        if "BFieldModelName" in conf:
            self.BFieldModelName = conf["BFieldModelName"]
        if "coords" in conf:
            self.coords = conf["coords"]
        if "component" in conf:
            self.component = conf["component"]
        if "m" in conf:
            self.m = conf["m"]
        if "modes" in conf:
            self.modes = conf["modes"]
        if "step" in conf:
            self.step = conf["step"]
            # Perp step for scaling factors h1 & h2
            self.hstep = conf["step"] * 20.
        if "hstep" in conf:
            self.hstep = conf["hstep"]
        if "sumh1h2" in conf:
            self.sumh1h2 = conf["sumh1h2"]
        if "traceh1h2" in conf:
            self.traceh1h2 = conf["traceh1h2"]
        if "maxIter" in conf:
            self.maxIter = conf["maxIter"]
        if "traceFrom" in conf:
            self.traceFrom = conf["traceFrom"]
        if "planetRadius" in conf:
            self.planetRadius = conf["planetRadius"]
        if "wlim" in conf:
            self.wlim = conf["wlim"]
        if "res" in conf:
            self.res = conf["res"]
        if "wRes" in conf:
            self.wRes = conf["wRes"]
        if "dydz" in conf:
            self.dydz = conf["dydz"]
        if "BC" in conf:
            self.BC = conf["BC"]
        if "tol" in conf:
            self.tol = conf["tol"]
        if "integrationMethod" in conf:
            self.integrationMethod = conf["integrationMethod"]
        if "fun" in conf:
            self.fun = conf["fun"]
        if "units" in conf:
            self.units = conf["units"]
        if "densityModel" in conf:
            self.densityModel = conf["densityModel"]
        if "BFieldModel" in conf:
            self.BFieldModel = conf["BFieldModel"]
        if "dipM" in conf:
            self.dipM = conf["dipM"]
        if "mIndex" in conf:
            self.mIndex = conf["mIndex"]
        if "mIndex" in conf:
            self.mIndex = conf["mIndex"]
        if "mIndex" in conf:
            self.mIndex = conf["mIndex"]
        if "planet" in conf:
            self.planet = conf["planet"]
        if "R" in conf:
            self.config["R"] = conf["R"]
        if "TH" in conf:
            self.config["TH"] = conf["TH"]
        if "PHI" in conf:
            self.config["PHI"] = conf["PHI"]
        if "ETIME" in conf:
            self.config["ETIME"] = conf["ETIME"]
        if "z0" in conf:
            self.config["z0"] = conf["z0"]
        if "Dp" in conf:
            self.config["Dp"] = conf["Dp"]
        if "name" in conf:
            self.name = conf["name"]
        if "zmax" in conf:
            self.zmax = conf["zmax"]
        if "scalingFactor" in conf:
            self.scalingFactor = conf["scalingFactor"]

    def updateFieldLines(self,
                         configuration=0,
                         NFieldLines=1,
                         L=[20],
                         th=None,
                         phi=[0],
                         component='toroidal'):
        self.flconfig = fieldlineConfigurations(configuration,
                                                NFieldLines=NFieldLines,
                                                L=L,
                                                th=th,
                                                phi=phi,
                                                component=component)


# ----------------------------------------------------------------------------

def loadsim(config):
    """ Load a configuration into a new SIM object """

    # Choose the planet
    planet = config['base']['planet']

    SIM = sim()
    SIM.update(configBasicSolve)

    if planet == 'box':
        SIM.update(Box)
    if planet == 'saturn':
        SIM.update(Saturn)
    if planet == 'saturn_uniform':
        SIM.update(SaturnUniform)
    if planet == 'earth':
        SIM.update(Earth)

    SIM.config = config['base']
    SIM.update(config['base'])
    SIM.config.update(KMAGConfigStandard)
    SIM.flconfig = config['fieldlines']
    return SIM

# ----------------------------------------------------------------------------

# Set up the Simulations


def configureSIMS(SIM_template, fieldlines):
    """  Create a new SIM object for each field line. Return array of SIMS """
    SIMS = []

    for i, fieldline in enumerate(fieldlines):

        SIM = deepcopy(SIM_template)
        SIM.update(fieldline.config)
        SIM.config = fieldline.config
        SIM.L = fieldline.L
        SIM.length = fieldline.length

        if SIM.coords == 'ds':
            SIM.z = fieldline.z
            SIM.zlim = [np.min(SIM.z), np.max(SIM.z)]
        elif SIM.coords == 'cos':
            costh = np.cos(fieldline.traceRTP[:, 1])
            SIM.z = costh
            SIM.zlim = [-pow(1. - 1. / SIM.L, 0.5), pow(1. - 1. / SIM.L, 0.5)]
        elif SIM.coords == 'cartesian':
            SIM.z = fieldline.z
            SIM.zlim = [fieldline.z[0], fieldline.z[-1]]

        SIM.vA = fieldline.vA
        SIM.B = fieldline.BT
        SIM.n = fieldline.n
        SIM.xyz = fieldline.traceXYZ
        SIM.rthetaphi = fieldline.traceRTP
        SIM.h1 = fieldline.h1_scale
        SIM.h2 = fieldline.h2_scale
        SIM.interpolate()
        SIMS.append(SIM)
    return SIMS


# ----------------------------------------------------------------------------

if __name__ == "__main__":
    from pylab import *
