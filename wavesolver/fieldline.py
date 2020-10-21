#!/usr/bin/env python3
import numpy as np
from wavesolver.linalg import *
from wavesolver.helperFunctions import *
from wavesolver.model import *
from KMAGhelper import KmagFunctions
from wavesolver.plot import *

# ----------------------------------------------------------------------------


class fieldline:
    """ fieldline class for trace parameters """
    def __init__(self,
                 L=None,
                 length=None,
                 z=None,
                 vA=None,
                 n=None,
                 BT=None,
                 dBT=None,
                 B=None,
                 traceXYZ=None,
                 traceRTP=None,
                 h1_scale=None,
                 h2_scale=None,
                 dh1_scale=None,
                 dh2_scale=None,
                 ddh1_scale=None,
                 ddh2_scale=None,
                 h1_scale_th=None,
                 h2_scale_th=None,
                 lnh1B=None,
                 lnh2B=None,
                 dlnh1B=None,
                 dlnh2B=None,
                 traceP1=None,
                 traceP2=None,
                 name='fieldline',
                 nIter=0,
                 maxIter=1E6,
                 step=1.,
                 direction=1,
                 BFluxConstant=1.,
                 method='RK4',
                 Rsheet=None,
                 mIndex=None,
                 config=None,
                 SIM=None):
        self.L = L
        self.length = length
        self.z = z if z else []
        self.n = n if n else []
        self.BT = BT if BT else []
        self.dBT = dBT if dBT else []
        self.B = B if B else []
        self.vA = vA if vA else []
        self.name = name
        self.BFluxConstant = None
        self.h1_scale = h1_scale if h1_scale else []
        self.dh1_scale = dh1_scale if dh1_scale else []
        self.ddh1_scale = ddh1_scale if ddh1_scale else []
        self.h2_scale = h2_scale if h2_scale else []
        self.dh2_scale = dh2_scale if dh2_scale else []
        self.ddh2_scale = ddh2_scale if ddh2_scale else []
        self.lnh1B = lnh1B if lnh1B else []
        self.lnh2B = lnh2B if lnh2B else []
        self.dlnh1B = dlnh1B if dlnh1B else []
        self.dlnh2B = dlnh2B if dlnh2B else []
        self.h1_scale_th = h1_scale_th if h1_scale_th else []
        self.h2_scale_th = h2_scale_th if h2_scale_th else []
        self.traceXYZ = traceXYZ if traceXYZ else []
        self.traceRTP = traceRTP if traceRTP else []
        self.traceP1 = traceP1 if traceP1 else []
        self.traceP2 = traceP2 if traceP2 else []
        self.Rsheet = Rsheet if Rsheet else []
        # Tracing Parameters
        self.mIndex = mIndex
        self.nIter = nIter
        self.maxIter = maxIter
        self.step = step
        self.direction = direction
        self.method = method
        self.config = config
        self.SIM = SIM

    def addParameters(self, R=None, B=None):
        """ Add initial parameters to the field line trace """
        BT = None
        # If position is not provided, calculate for the last point
        if R is None:
            R = self.traceXYZ[-1]
            B = self.B[-1]
            BT = self.BT[-1]
        # If a position is provide, first calculate the field
        else:
            B = self.SIM.BFieldModel(R, self.SIM)
            BT = np.sqrt(B.dot(B))
            # self.addPositionAndField(R, B, BT=BT)

        self.n.append(densityModels(R, self.SIM)[0])
        self.vA.append(calculate_vA(BT, self.n[-1], self.SIM.m)[0])
        self.z.append(self.nIter * self.direction)
        h_scales = self.scaleFactors(R, B)
        h_scales_th = self.scaleFactorModel(R, B)
        self.h1_scale.append(h_scales[0])
        self.h2_scale.append(h_scales[1])
        self.h1_scale_th.append(h_scales_th[0])
        self.h2_scale_th.append(h_scales_th[1])

    def parameters(self):
        """ Calculate parameters for the traced field line """
        self.convert2NPArrays()
        if self.SIM.BFieldModelName == 'uniform':
            L_index = np.argmin(self.traceRTP[:, 0], axis=0)
        else:
            L_index = np.argmax(self.traceRTP[:, 0], axis=0)
        L = self.traceRTP[L_index, 0]
        self.L = L
        self.SIM.L = L
        self.Rsheet = self.traceXYZ[L_index]
        arraySecondPartLength = size(self.BT) - (L_index + 1)

        self.z = (np.zeros(size(self.BT))
                  + np.arange(-L_index, arraySecondPartLength + 1))
        self.z = np.asarray(self.z) * self.step * self.SIM.units

        self.length = abs(self.z[-1] - self.z[0]) / self.SIM.units

        R = self.traceXYZ
        B = self.B
        BT = self.BT

        n = densityModels(R, self.SIM, fieldline=self)
        vA = calculate_vA(BT, n, self.SIM.m)

        h_scales = self.scaleFactorModel(R, B)
        h_scales_th = self.scaleFactorModel(R, B)
        # h_scales_th = self.scaleFactors(R, B)

        self.n = np.array(n)
        self.vA = np.array(vA)
        self.h1_scale_th = np.asarray(h_scales_th)[0]
        self.h2_scale_th = np.asarray(h_scales_th)[1]
        self.h1_scale = np.asarray(h_scales)[0]
        self.h2_scale = np.asarray(h_scales)[1]
        self.convert2NPArrays()

    def reverseData(self):
        """ Reverse field line data """
        varList = [self.z, self.n, self.BT, self.B, self.vA, self.traceXYZ,
                   self.traceRTP, self.h1_scale, self.h2_scale,
                   self.h1_scale_th, self.h2_scale_th,
                   self.traceP1, self.traceP2]
        return [i.reverse() for i in varList]

    def convert2NPArrays(self):
        """ Convert to NP arrays for easier manipulation """
        self.z = np.asarray(self.z)
        self.n = np.asarray(self.n)
        self.BT = np.asarray(self.BT)
        self.B = np.asarray(self.B)
        self.vA = np.asarray(self.vA)
        self.traceXYZ = np.asarray(self.traceXYZ)
        self.traceRTP = np.asarray(self.traceRTP)
        self.h1_scale = np.asarray(self.h1_scale)
        self.h2_scale = np.asarray(self.h2_scale)
        self.h1_scale_th = np.asarray(self.h1_scale_th)
        self.h2_scale_th = np.asarray(self.h2_scale_th)
        self.traceP1 = np.asarray(self.traceP1)
        self.traceP2 = np.asarray(self.traceP2)

    def meetTraceConditions(self, R=None):
        """ Field line tracing boundary conditions (e.g. planet) """
        if self.nIter >= self.maxIter:
            return False
        else:
            # If it's a Box Model
            if self.SIM.BFieldModelName == 'uniform':
                if R is None:
                    R = self.traceXYZ[-1]
                R = np.asarray(R)
                # print(R)
                if abs(R[2]) >= self.SIM.zmax:
                    return False
                else:
                    return True
            # If it's a dipole or a more realistic 2D/3D model
            else:
                r = 1.
                if R is None:
                    R = self.traceXYZ[-1]
                R = np.asarray(R)
                r = np.sqrt(R.dot(R))
                if r > 1.:
                    return True
                else:
                    return False

    def addPositionAndField(self, R, B=None, BT=None):
        """ Add position and field data to the field line trace """
        if B is None:
            B = self.SIM.BFieldModel(R, self.SIM)
        if BT is None:
            BT = np.linalg.norm(B, axis=-1)

        if size(R) > 3:
            R = np.transpose(R)
            r = np.linalg.norm(R, axis=0)
            self.traceRTP = car2sph(R)
            self.traceXYZ = R
            self.B = B
            self.BT = BT
        else:
            self.traceRTP.append(car2sph(R))
            self.traceXYZ.append(R)
            self.B.append(B)
            self.BT.append(BT)

# ----------------------------------------------------------------------------

    def traceEuler(self, R=None, direction=None):
        """ Euler method tracing """
        B = self.SIM.BFieldModel(R, self.SIM)
        BT = np.sqrt(B.dot(B))

        if direction is None: 
            direction = self.direction
        dd = self.step * direction
        step = [B[0] / BT * dd, B[1] / BT * dd, B[2] / BT * dd]
        R = [sum(x) for x in zip(R, step)]
        return R, B, BT

# ----------------------------------------------------------------------------

    def traceRK4(self, R=None, direction=None):
        """
        Fourth order Runge-Kutta method tracing
        Solving x' = f(x,t) with x(t[0]) = x0
        """
        if direction is None:
            direction = self.direction
        dd = self.step * self.direction

        def f(R, step):
            B = self.SIM.BFieldModel(R + step, self.SIM)
            return unitVector(B)

        k1 = f(R, 0.0 * dd)
        k2 = f(R, 0.5 * dd * k1)
        k3 = f(R, 0.5 * dd * k2)
        k4 = f(R, 1.0 * dd * k3)

        kavg = unitVector((k1 + 2.0 * (k2 + k3) + k4) / 6.0)

        R = R + 1.0 * dd * kavg
        B = self.SIM.BFieldModel(R, self.SIM)
        BT = np.linalg.norm(B, axis=-1)
        return R, B, BT

# ----------------------------------------------------------------------------
    def newTracePoint(self,
                      R=None,
                      direction=None,
                      parameters=False,
                      collect=True):
        """ Next trace point """
        B, BT = [None, None]

        if R is None:
            R = self.traceXYZ[-1]
            R = np.asarray(R)

        if self.method == 'RK4':
            R, B, BT = self.traceRK4(R=R, direction=direction)
        elif self.method == 'euler':
            R, B, BT = self.traceEuler(R=R, direction=direction)

        if collect is True:
            self.addPositionAndField(R, B=B, BT=BT)
            if parameters is True:
                self.addParameters()
        return R

# ----------------------------------------------------------------------------
    # Experimental
    def scaleFactors(self, R, B, R1=None, R2=None):
        """ Calculate scale factors locally without tracing """
        if self.SIM.BFieldModelName == 'uniform':
            h1 = np.full_like(B, 1.)
            h2 = np.full_like(B, 1.)
        else:
            R = np.asarray(R)
            RTP = car2sph(R).T
            r, TH, PHI = [np.asarray(RTP[i]) for i in range(0,3)]

            step = self.step
            step = 1E10

            if R1 is None:
                R1, R2 = perpendicularStep(R, B, step=step)

            phi_unit = unitVector(np.asarray([-R.T[1],
                                  R.T[0],
                                  np.zeros_like(R.T[0])]).T)

            B1 = self.SIM.BFieldModel(R1, self.SIM)
            B1_unit = unitVector(B1)
            B2 = self.SIM.BFieldModel(R2, self.SIM)
            B2_unit = unitVector(B2)

            perp1_unit, perp2_unit = perpendicularUnitVectors(R, B)
            B = np.asarray(B)
            B_unit = unitVector(B)
            dB1 = np.subtract(B1_unit, B_unit)
            dB2 = np.subtract(B2_unit, B_unit)
            # dB1 = np.subtract(B1, B)
            # dB2 = np.subtract(B2, B)
            # h1 = np.dot(dB1, perp1_unit)
            # h2 = np.dot(dB2, perp2_unit)
            h1 = np.einsum('...i,...i', dB1, perp1_unit)
            h2 = np.einsum('...i,...i', dB2, perp2_unit)
            # np.dot(a,b.T).diagonal()
            # h1_.append(h1)
            # h2_.append(h2)
            # return np.asarray(h1_), np.asarray(h2_)
        return np.asarray(h1), np.asarray(h2)

# ----------------------------------------------------------------------------

    def ScaleFactorsFromSeperation(self, R=None, B=None, R1=None, R2=None):
        """ Calculate scale factors by tracing adjacent field lines """

        def minTraceSeparation(trace1, trace2):
            trace_diff = [np.amin(np.linalg.norm(np.subtract(trace2, R0),
                          axis=1)) for R0 in trace1]
            return np.asarray(trace_diff)

        h1_scale = minTraceSeparation(self.traceXYZ, self.traceP2)
        h2_scale = minTraceSeparation(self.traceXYZ, self.traceP1)

        h1_scale = np.absolute(np.asarray(h1_scale))
        h2_scale = np.absolute(np.asarray(h2_scale))
        h1_scale = np.nan_to_num(h1_scale)
        h2_scale = np.nan_to_num(h2_scale)
        h1_scale *= 1. / np.amax(h1_scale, axis=0)
        # h1_scale *= self.SIM.units
        h1_scale *= self.L * self.SIM.units
        h2_scale *= 1. / np.amax(h2_scale, axis=0)
        h2_scale *= np.max(np.reciprocal(np.multiply(self.BT, h1_scale)))
        h2_scale *= 1E-9

        self.h1_scale = h1_scale
        self.h2_scale = h2_scale
        return h1_scale, h2_scale

# ----------------------------------------------------------------------------

    def scaleFactorModel(self, R, B):
        """ Theoretical values for scale factors for a dipole """
        R = np.asarray(R)
        RTP = car2sph(R).T
        r, TH, PHI = [np.asarray(RTP[i]) for i in range(0, 3)]
        B = np.asarray(B)
        BT = np.linalg.norm(B, axis=-1)
        h1 = np.multiply(r * self.SIM.units, np.sin(TH))  # Toroidal Mode
        h2 = np.reciprocal(np.multiply(BT, h1))  # poloidal
        if self.SIM.BFieldModel == 'uniform':
            h1 = np.full_like(B, 1.)
            h2 = np.full_like(B, 1.)
        return h1, h2
# ----------------------------------------------------------------------------


def traceFieldline(R0, SIM):
    """ Trace a field line from R0 until tracing conditions are met """
    fl = fieldline(step=SIM.step, maxIter=SIM.maxIter, SIM=SIM)

    if SIM.BFieldModelName == 'KMAG':
        SIM.config["step"] = SIM.step
        KmagFunctions.traceKMAG([R0], SIM.config)
        R, B, BT = KmagFunctions.KMAGfieldlines(SIM)
        for n in range(0, len(BT)):
            fl.addPositionAndField(np.asarray(R[n]),
                                   B=np.asarray(B[n]),
                                   BT=BT[n])
        fl.name = 'KMAG'
    elif SIM.BFieldModelName == 'dipole':
        fl.nIter = 0
        fl.addPositionAndField(R0)
        reversedDirections = False
        for ndir, direction in enumerate(SIM.directions):
            fl.direction = direction
            fl.newTracePoint(R=R0)
            while fl.meetTraceConditions():
                fl.nIter += 1
                fl.newTracePoint()
            if len(SIM.directions) > 1 and fl.meetTraceConditions(
                    R=fl.newTracePoint(R=R0,
                                       direction=SIM.directions[1],
                                       collect=False)):
                fl.reverseData()
                reversedDirections = True
        if reversedDirections:
            fl.reverseData()
        fl.name = 'DIP'
    elif SIM.BFieldModelName == 'uniform':
        fl.nIter = 0
        fl.addPositionAndField(R0)
        reversedDirections = False
        for ndir, direction in enumerate(SIM.directions):
            fl.direction = direction
            fl.newTracePoint(R=R0)
            while fl.meetTraceConditions():
                fl.nIter += 1
                fl.newTracePoint()
            if len(SIM.directions) > 1 and fl.meetTraceConditions(
                    R=fl.newTracePoint(R=R0,
                                       direction=SIM.directions[1],
                                       collect=False)):
                fl.reverseData()
                reversedDirections = True
        if reversedDirections:
            fl.reverseData()
        fl.name = 'UNI'

    fl.convert2NPArrays()
    fl.parameters()

    return fl

# ----------------------------------------------------------------------------

def traceFieldlines(SIM):
    """ Trace all field lines configured in the simulation (SIM) """
    configs = SIM.flconfig
    configBase = SIM.config

    fls = []
    for i, config in enumerate(configs):

        SIM.update(config)

        configFull = configBase
        configFull.update(config)
        R = configFull["R"] if "R" in configFull else 1.
        TH = configFull["TH"] if "TH" in configFull else -20.
        PHI = configFull["PHI"] if "PHI" in configFull else 180.

        if configFull["BFieldModelName"] == "uniform":
            traceStart = sph2car([R,
                                  np.deg2rad(TH),
                                  np.deg2rad(PHI)])
            SIM.dipM = DIPOLE_BEQ_f(R) * 1E-9
            SIM.zmax = DIPOLE_Length_f(R) / 2.
            # SIM.dipM = KMAG_BEQ_f(R) * 1E-9
            # SIM.zmax = KMAG_Length_f(R) / 2.
        else:
            traceStart = sph2car([R + SIM.step,
                                  np.deg2rad(TH),
                                  np.deg2rad(PHI)])

        if "name" in config:
            fieldlineName = config["name"]
        else:
            fieldlineName = r'TH = %.0f | PHI = %.0f$^0$' % (TH, PHI)

        fl = traceFieldline(traceStart, SIM)
        fl.name = fieldlineName
        fl.config = configFull
        print('FINISHED TRACING #%d. | L: %3.2f | Length: %3.2f'
              % (i, fl.L, fl.length))

        if SIM.sumh1h2:
            fl.h1_scale_th = normalize(abs(np.cumsum(fl.h1_scale_th)))
            fl.h2_scale_th = normalize(abs(np.cumsum(fl.h2_scale_th)))

        if SIM.traceh1h2:
            L_index = np.argmax(fl.traceRTP[:, 0], axis=0)
            psheet_start = fl.traceXYZ[L_index, :]
            B = fl.B[L_index]
            perpStep = perpendicularStep(psheet_start, B, step=SIM.hstep)
            P1 = traceFieldline(perpStep[0], SIM)
            P2 = traceFieldline(perpStep[1], SIM)
            fl.traceP1 = P1.traceXYZ
            fl.traceP2 = P2.traceXYZ
            fl.ScaleFactorsFromSeperation()

        # fieldlineScalingFactorPlot(fl)
        # exit()
        fls.append(fl)
    return fls

# ----------------------------------------------------------------------------


def findPlasmaSheet(fieldlines):
    """ FIND AND PLOT PLASMA SHEET COORDINATES IN KMAG """
    Rsheet = []
    for fl in fieldlines:
        Rsheet.append(fl.Rsheet)
    Rsheet = np.asarray(Rsheet)
    dayside = Rsheet[Rsheet[:, 0] > 0, :]
    nightside = np.flip(Rsheet[Rsheet[:, 0] < 0, :], axis=0)
    Rsheet = np.concatenate((nightside, dayside), axis=0)
    Rsheets = []
    if nightside.size:
        Rsheets.append(nightside.T)
    if dayside.size:
        Rsheets.append(dayside.T)
    Rsheet = np.asarray(Rsheet).T
    return Rsheets
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from pylab import *
