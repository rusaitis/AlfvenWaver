#!/usr/bin/env python3

import numpy as np
from scipy import optimize
from wavesolver.model import *
from wavesolver.linalg import *
from wavesolver.helperFunctions import *
# ----------------------------------------------------------------------------
#
# Shooting method and root finding procedures.
#
# ----------------------------------------------------------------------------


class timeseries:
    """ A timeseries class to store eigenfunctions """
    def __init__(self, y=None, t=None, dy=None, name=None, N=None,
                 zeros=None, zeroIntervals=None, roots=None, mode=None,
                 A=None, erf=None, errors=None, b=None, E=None, xi=None):
        self.y = y if y else []
        self.t = t if t else []
        self.dy = dy if dy else []
        self.name = name if name else ''
        self.N = N if N else 0
        self.zeros = zeros if zeros else []
        self.zeroIntervals = zeroIntervals if zeroIntervals else []
        self.roots = roots if roots else []
        self.signFlipIndices = None
        self.signFlips = None
        self.A = A
        self.mode = mode if mode else 0
        self.erf = erf if erf else []
        self.errors = errors if errors else []
        self.b = b
        self.E = E
        self.xi = xi

    def normalize(self, A=1.0):
        """ Normalize eigenfunction """
        Ay = np.amax(np.abs(self.y), axis=0)
        if Ay > 0:
            self.y = np.divide(self.y * A, Ay)
        if len(self.dy) > 0:
            Ady = np.amax(np.abs(self.dy), axis=0)
            if Ady > 0:
                self.dy = np.divide(self.dy * A, Ady)

    def normalizeFields(self, A=1.0):
        """ Normalize Electric filed, Magnetic field, and displacement """
        self.amplitudes()
        self.E = np.divide(self.E * A, self.A[0])
        self.b = np.divide(self.b * A, self.A[1])
        self.xi = np.divide(self.xi * A, self.A[2])

    def amplitudes(self):
        """ Find the Amplitudes of E, b, and xi """
        AE = np.amax(np.abs(self.E), axis=0)
        Ab = np.amax(np.abs(self.b), axis=0)
        Axi = np.amax(np.abs(self.xi), axis=0)
        self.A = [AE, Ab, Axi]

    def scale(self, A=1., normfirsts=False):
        """ Rescale the eigenfunctions by the respective amplitudes """
        # Scale each eigenfunctions individually
        if normfirsts:
            self.E *= A[0]
            self.b *= A[1]
            self.xi *= A[2]
        # Scale each eigenfunction by the same amplitude
        else:
            self.E *= A[0]
            self.b *= A[0]
            self.xi *= A[0]

    def findSignFlips(self):
        """ Find the indices and number of sign flips in the eigf """
        self.signFlipIndices = np.where(np.sign(self.y[:-1])
                                        != np.sign(self.y[1:]))[0] + 1
        self.signFlips = len(self.signFlipIndices)

    def findModeNumber(self):
        """ Deduce the mode number by the num of sign flips in the eigf """
        dySignFlips = np.where(np.sign(self.dy[:-1])
                               != np.sign(self.dy[1:]))[0] + 1
        self.mode = len(dySignFlips)
        self.name = "m=%d" % self.mode

    def calculate_erf(self, w1, w2, res, SIM):
        """ Calculate the error of the solution at the boundary """
        self.t = np.linspace(w1, w2, res)
        erf = []
        for i, w in enumerate(self.t):
            SIM.wlim = [w, w]
            sol = SIM.integrate()
            erf.append(sol.y[0, -1] - SIM.BC[1])
        self.y = np.asarray(erf)
        self.findSignFlips()

    def calculate_E_and_b(self, SIM):
        """ Compute electric field, magnetic field, and displacement """
        self.xi = displacement(self, SIM)
        self.E = electricPerturbation(self, SIM)
        self.b = magneticPerturbation(self, SIM)

    def optimizeRoots(self, SIM, VERBOSE=True):
        """ Root optimization (experimental) """
        def funScan(w, SIM):
            SIM.wlim = [w, w]
            sol = SIM.integrate()
            return sol.y[0, -1] - SIM.BC[1]
        sols = []

        for w in self.zeros:
            sol = optimize.fsolve(funScan, w, args=(SIM))
            sols.append(sol[0])
        return sols

    def shootSolutions(self, SIM, VERBOSE=True, normalize=False,
                       normfirst=True, normfirsts=False, scaleFactor=None):
        """ Compute the eigenfunctions for each boundary-satisfying
        eigenfrequency (zeros of the boundary error function) """
        solutions = []
        Ascale = 1.
        for i in range(0, len(self.zeros)):
            # Configure the eigenfrequency range for the solution
            SIM.wlim = self.zeroIntervals[i]
            # Apply a shooting method
            solution = shoot(SIM, VERBOSE=VERBOSE)
            # Calculate fields and displacement
            solution.calculate_E_and_b(SIM)

            # Normalization / Scaling
            if scaleFactor is not None:
                solution.scale(scaleFactor)
            if normalize:
                solution.normalizeFields(1.)
            solution.amplitudes()
            if normfirst:
                if i == 0:
                    Ascale = np.asarray(solution.A)
                solution.scale(1. / Ascale, normfirsts)
            solution.findModeNumber()
            solutions.append(solution)
        return solutions

    def addLOG(self, *seriesv):
        """ Keep a log of the solution refinements """
        if self.t is None:
            self.t = []
        if self.y is None:
            self.y = []
        for series in seriesv:
            self.t = np.concatenate((self.t, series.t), axis=0)
            self.y = np.concatenate((self.y, series.y), axis=0)

# ----------------------------------------------------------------------------


def rootErrorFunction(SIM, VERBOSE=True, plotLOG=False):
    """ A semi-intelligent root search algorithm that starts with an initial
    range of frequencies given in SIM configuration. It works by finding the
    number of sign flips in the root error function, and redefining the span
    of frequency search to find a required number of roots. """
    errors = timeseries()
    errors2 = timeseries()
    LOG = timeseries()

    def refineSearch(w1, w2, signFlips, modes):
        """ Refines the range of frequencies to search """
        safebuffer = 1.
        if signFlips > 0:
            rootInterval = (w2 - w1) / signFlips
            w1 = w1
            w2 = w2 + (modes + safebuffer - signFlips) * rootInterval
        else:
            w1 = w1
            w2 = w2 + (w2 - w1)
        return w1, w2

    w1 = SIM.config["wlim"][0]
    w2 = SIM.config["wlim"][1]
    res = 20
    max_number_of_refines = 100

    # ROOT SEARCH DIAGNOSTICS
    if plotLOG:
        from wavesolver.plot import io, plotErrorFunctions
        ioconfig = io()
        ioconfig.save = False
    # ------------------------------------------------------------------------
    for i in range(max_number_of_refines):
        errors.calculate_erf(w1, w2, res, SIM)

        if VERBOSE:
            print('w1: %6.5f | w2: %6.5f | res: %d | Sign Flips: %d'
                  % (w1, w2, res, errors.signFlips))

        if plotLOG:
            plotErrorFunctions(errors, SIM, ioconfig)

        errors2.calculate_erf(w1, w2, 2 * res, SIM)
        if errors2.signFlips > errors.signFlips and \
           errors2.signFlips < 2. * SIM.modes:
            res *= 2  # The w root number is changing, let's refine
        else:
            if plotLOG:
                plotErrorFunctions(errors2, SIM, ioconfig)
            if (errors2.signFlips > 2. * SIM.modes) or \
               (errors2.signFlips < SIM.modes):
                w1, w2 = refineSearch(w1, w2, errors2.signFlips, SIM.modes)
            else:
                modeTruncate = np.min([errors2.signFlips, SIM.modes])
                signFlipIndices = errors2.signFlipIndices[:modeTruncate]
                errors2.zeros = errors2.t[signFlipIndices]
                errors2.zeroIntervals = (
                    np.stack((errors2.t[signFlipIndices[:] - 1],
                              errors2.t[signFlipIndices[:]]),
                             axis=-1))
                break
        LOG.addLOG(errors, errors2)
    errors2.normalize(1.)
    if VERBOSE:
        print("Found %d roots." % (len(errors2.zeros)))
    return errors2, LOG

# ----------------------------------------------------------------------------


def rk4(f, x0, t):
    """Fourth-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.

    USAGE:
        x = rk4(f, x0, t)
        y = rk4( f, [a,z1,s1], t)
        w1 = y[-1,0]

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a NumPy array.  In this
                case f must return a NumPy array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or NumPy array
                if a system of equations is being solved.
        t     - list or NumPy array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - NumPy array containing solution values corresponding to each
                entry in t array.  If a system is being solved, x will be
                an array of arrays.
    """

    n = len(t)
    x = np.array([x0] * n)
    for i in range(n - 1):
        h = t[i + 1] - t[i]
        k1 = h * f(x[i], t[i])
        k2 = h * f(x[i] + 0.5 * k1, t[i] + 0.5 * h)
        k3 = h * f(x[i] + 0.5 * k2, t[i] + 0.5 * h)
        k4 = h * f(x[i] + k3, t[i + 1])
        x[i + 1] = x[i] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0

    return x

# ----------------------------------------------------------------------------


def shoot(SIM, VERBOSE=True, max_iter=10):
    """ Shooting Method based on the tangent method """

    soln = timeseries()

    def printlog(i, w, error):
        print("%2d: w = %6.5f, error = %3.2e" % (i, w / 2. / np.pi, error))

    w1, w2 = (SIM.wlim[0], SIM.wlim[1])

    SIM.wlim = [w1, w1]
    sol = SIM.integrate()
    erf1 = sol.y[0, -1]

    if VERBOSE:
        printlog(0, w1, SIM.BC[1] - erf1)

    for i in range(max_iter):

        SIM.wlim = [w2, w2]
        sol = SIM.integrate()
        erf2 = sol.y[0, -1]

        if VERBOSE:
            printlog(i, w2, SIM.BC[1] - erf2)

        nonzeroCount = np.count_nonzero(sol.y[:, 0])

        # Break the search for a root if the tolerance is met
        if (abs(SIM.BC[1] - erf2) < SIM.tol
                and nonzeroCount > 0
                and not np.isnan(sol.y[:, 0]).any()):
            break

        # Tangent method for the refined root
        w1, w2 = (w2,
                  w2 + (w2 - w1) / (erf2 - erf1) * (SIM.BC[1] - erf2))
        erf1 = erf2

    if abs(SIM.BC[1] - erf2) >= SIM.tol:
        soln.erf.append("Maximum number of iterations (%d) exceeded" %
                        max_iter)
        if VERBOSE:
            print(soln.erf[-1])

    soln.y = sol.y[0, :]
    soln.dy = sol.y[1, :]
    soln.t = sol.t
    soln.roots = w2
    soln.erf = SIM.BC[1] - erf2
    return soln

# ----------------------------------------------------------------------------


if __name__ == "__main__":
    from pylab import *
