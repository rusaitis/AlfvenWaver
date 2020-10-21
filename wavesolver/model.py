#!/usr/bin/env python3

import numpy as np
from scipy import optimize
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from wavesolver.linalg import *
from wavesolver.helperFunctions import *
import matplotlib.pyplot as plt
from pylab import *

# ----------------------------------------------------------------------------
#                         __  __           _      _
#                        |  \/  | ___   __| | ___| |
#                        | |\/| |/ _ \ / _` |/ _ \ |
#                        | |  | | (_) | (_| |  __/ |
#                        |_|  |_|\___/ \__,_|\___|_|
#
# This file contains physical parameters, density models, analytical
# solutions, and a dipole field model.
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
#                               CONSTANTS
# ----------------------------------------------------------------------------
c = 299792458.  # ms-1
mu0 = 4. * np.pi * 1E-7  # Permeability, Hm-1
RS = 60268E3  # m
RE = 6356E3  # m polar
RE = 6378E3  # m equator
dipM_Earth = 3.12E-5
dipM_EarthSI = 8E15
dipM_Saturn = 20E-6
amu = 1.660539E-27  # Atomic mass unit
mp = 1.6726219E-27  # Mass of a proton
# ----------------------------------------------------------------------------
#                            SIM PARAMETERS
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
#                           Bagenal et al. (2011)
# ----------------------------------------------------------------------------
# Bagenal, F., & Delamere, P. A. (2011). Flow of mass and energy in the
# magnetospheres of Jupiter and Saturn. Journal of Geophysical Research:
# Space Physics, 116(A5). https://doi.org/10.1029/2010JA016294
# The densities at Saturn are derived from Cassini CAPS data (<5 RS) by
# Sittler et al. [2008] and (>5 RS) by Thomsen et al. [2010].
# Digitized using the WebPlotDigitizer.
# https://automeris.io/WebPlotDigitizer/
#
# Figure 2. Water group (W+) densities, marked by the dashed blue line.
# First column: Radial Distance
# Second column: Water-Group Density (cm^-3)
n_bagenal = [[3.01, 3.13, 3.24, 3.35, 3.45, 3.52, 3.58, 3.67,
3.78, 3.96, 4.14, 4.34, 4.54, 4.74, 4.92, 5.10, 5.32, 5.49, 5.69,
5.90, 6.09, 6.28, 6.48, 6.64, 6.91, 7.06, 7.26, 7.41, 7.62, 7.81,
7.97, 8.26, 8.54, 8.81, 9.19, 9.56, 9.88, 10.22, 10.57, 10.90, 11.24,
11.50, 11.69, 11.95, 12.29, 12.61, 12.78, 13.15, 13.56, 13.90, 14.34,
14.62, 15.00, 15.51, 16.04, 16.82, 17.30, 17.89, 18.66, 19.19, 19.68],
[17.56,21.54,26.05,31.97,39.81,52.56,63.56,77.99,86.40,86.40,
85.14,82.69,81.49,79.14,76.86,71.44,65.44,59.08,53.33,48.14,42.83,
38.10,34.90,31.97,30.16,25.68,20.92,18.08,15.17,12.73,10.53,9.23,
8.21,6.50,5.61,4.71,3.84,3.08,2.48,1.96,1.62,1.34,1.09,0.93,0.75,
0.62,0.51,0.43,0.39,0.33,0.28,0.24,0.20,0.17,0.13,0.12,0.11,0.10,
0.09,0.08,0.07]]
# Figure 5.
# First column: Radial Distance (R_S)
# Second column: Saturn Scale Height for W+ (water group) (R_S)
H_bagenal = [[2.99, 3.16, 3.33, 3.55, 3.78, 3.98, 4.31, 4.52,
4.75,5.05, 5.34, 5.70, 6.12, 6.47, 7.00, 7.46, 8.03, 8.86, 9.59,
10.22, 11.00, 11.72, 12.49, 13.09, 13.84, 14.55, 15.18, 15.96,
16.52, 17.28, 18.07, 18.80, 19.77],
[0.22, 0.24, 0.26, 0.27, 0.31, 0.31, 0.39, 0.48, 0.56, 0.64,
0.74, 0.87, 1.04, 1.24, 1.47, 1.62, 1.81, 1.97, 2.08, 2.23, 2.40,
2.55, 2.82, 2.99, 3.30, 3.62, 3.84, 4.11, 4.40, 4.64, 4.78, 4.78,
4.76]]

# Interpolation functions for Plasma Scale Height and plasma density
fn_bagenal = interp1d(n_bagenal[0], n_bagenal[1],
                      bounds_error=False, fill_value=0.07)
# fill_value=(np.min(b_lines), np.max(b_lines))
# fn_bagenal = interp1d(n_bagenal[0], n_bagenal[1], \
            # kind='linear', bounds_error=False, fill_value=0.0)
# fn_bagenal = CubicSpline(n_bagenal[0], n_bagenal[1], extrapolate=True)
# fn_bagenal = InterpolatedUnivariateSpline(n_bagenal[0], n_bagenal[1], bbox=[0.,30.])
# fn_bagenal = interpolate.splrep(n_bagenal[0], n_bagenal[1])
# fn_bagenal = interpolate.splev(n_bagenal[0], n_bagenal[1])
# fn_bagenal =  UnivariateSpline(n_bagenal[0], n_bagenal[1], k=1, ext=3)

fH_bagenal = CubicSpline(H_bagenal[0], H_bagenal[1], extrapolate=True)
# fH_bagenal = interp1d(H_bagenal[0], H_bagenal[1], \
            # kind='cubic', bounds_error=False)
# fH_bagenal = UnivariateSpline(H_bagenal[0], H_bagenal[1], k=3)
# z = np.polyfit(n_bagenal[0], n_bagenal[1], 50)
# fn_bagenal = np.poly1d(z)
# ----------------------------------------------------------------------------
#                           Yates et al. (2016)
# ----------------------------------------------------------------------------
vA0_Yates = 105E3  # ms-1 Yates et al. (2016)
vA1_Yates = 3500E3  # ms-1 Yates et al. (2016)
BT_Yates = 2E-9  # T, eq. plasma, Arridge et al. [2011], water ions only
n0_Yates = 1.5E4  # m-3 eq. plasma, Arridge et al. [2011], water ions only
z0_RS_Yates = 2.  # RS
L_RS_Yates = 60.
n0_Yates = 9659.3954
n1_Yates = 8.693455
BT = BT_Yates  # 2nT
# Computation of the plasma densities based on the Alfven velocities given
# in the (Yates et al., 2016) paper.
# vA = np.divide(B, np.sqrt(mu0 * n * m))
# (2E-9)^2/((105E3)^2 * 4.*pi*1E-7*18*1.660539E-27)
# n0 = 9659.395401578
# (2E-9)^2/((3500E3)^2 * 4.*pi*1E-7*18*1.660539E-27)
# n1 = 8.693455861

# ----------------------------------------------------------------------------
#                           Persoon et al. (2013)
# ----------------------------------------------------------------------------
m_Persoon = 4.0
n_Persoon = 4.8
n0_Persoon = 72.0E6  # m-3
# m0_Persoon = 18*1.661E-27
m0_Persoon = 18 * amu
# m0_Persoon = 13 * amu # (Cramm et al., 1998)
# m0_Persoon = 23 * amu
R0_Persoon = 4.6  # RS

# Persoon, A. M., Gurnett, D. A., Leisner, J. S., Kurth, W. S., Groene,
# J. B., & Faden, J. B. (2013). The plasma density distribution in the
# inner region of Saturn’s magnetosphere: SATURN’S INNER MAGNETOSPHERIC
# DENSITIES. Journal of Geophysical Research: Space Physics, 118(6),
# 2970–2974. https://doi.org/10.1002/jgra.50182
# Digitized using the WebPlotDigitizer.
# https://automeris.io/WebPlotDigitizer/
#
# Figure 2. Plasma scale height values, marked by the black circles.
# First column: L-Shell
# Second column: H, Plasma Scale Height (R_S)
H_persoon = [[2.710, 2.777, 2.844, 2.900, 2.927, 2.962, 2.988, 3.014, 3.041,
3.068, 3.104, 3.181, 3.261, 3.319, 3.356, 3.393, 3.423, 3.453, 3.488, 3.534,
3.577, 3.619, 3.660, 3.700, 3.782, 3.867, 3.917, 3.952, 3.987, 4.022, 4.061,
4.107, 4.164, 4.212, 4.267, 4.323, 4.376, 4.420, 4.479, 4.534, 4.604, 4.685,
4.756, 4.857, 4.942, 5.023, 5.089, 5.128, 5.167, 5.210, 5.258, 5.312, 5.404,
5.508, 5.637, 5.729, 5.801, 5.887, 6.044, 6.152, 6.221, 6.284, 6.352, 6.410,
6.465, 6.510, 6.552, 6.586, 6.640, 6.700, 6.838, 6.905, 7.065, 7.237, 7.414,
7.602, 7.722, 7.799, 7.851, 7.868, 7.970, 7.993, 8.021, 8.093, 8.209, 8.320,
8.431, 8.491, 8.549, 8.628, 8.703, 8.850, 8.899, 8.945, 8.985, 9.027, 9.061,
9.100, 9.204, 9.256, 9.272, 9.374, 9.422, 9.445, 9.470],
[0.353, 0.345, 0.338, 0.331, 0.323, 0.314, 0.306, 0.299, 0.292, 0.285, 0.278,
0.277, 0.275, 0.277, 0.286, 0.295, 0.302, 0.310, 0.319, 0.328, 0.337, 0.346,
0.354, 0.362, 0.355, 0.350, 0.362, 0.372, 0.383, 0.393, 0.405, 0.417, 0.427,
0.438, 0.449, 0.462, 0.475, 0.487, 0.501, 0.514, 0.528, 0.545, 0.560, 0.575,
0.590, 0.607, 0.623, 0.640, 0.657, 0.676, 0.699, 0.718, 0.738, 0.756, 0.773,
0.791, 0.812, 0.836, 0.853, 0.877, 0.905, 0.930, 0.958, 0.988, 1.019, 1.046,
1.072, 1.098, 1.127, 1.154, 1.132, 1.129, 1.153, 1.176, 1.203, 1.212, 1.201,
1.167, 1.140, 1.117, 1.186, 1.152, 1.225, 1.274, 1.316, 1.351, 1.385, 1.423,
1.476, 1.526, 1.556, 1.517, 1.466, 1.417, 1.370, 1.328, 1.293, 1.256, 1.211,
1.173, 1.153, 1.246, 1.293, 1.209, 1.333]]
# Figure 3(b). Equatorial density measurements averaged in non-overlapping
# L-shell bins of 0.2 RS.
# First column: R, Saturnian Radii (R_S)
# Second column: Equatorial Elctron Density, n_e_eq (cm^-3)
neq_persoon = [[2.50,2.70,2.90,3.10,3.30,3.50,3.70,3.89,
                4.09,4.30,4.50,4.70,4.90,5.09,5.28,5.48,5.68,5.89,6.10,6.29,
                6.49,6.69,6.88,7.09,7.29,7.50,7.69,7.87,8.07,8.28,8.49,8.71,
                8.95,9.20,9.37,9.48,9.67,9.90,10.07],
               [17.27,14.99,21.85,31.83,36.09,45.25,50.88,
                60.66,65.98,71.16,73.00,72.43,68.35,62.38,55.06,47.01,42.19,
                36.93,32.87,29.26,25.82,22.41,19.62,17.03,15.03,13.27,12.00,
                10.68,9.59,8.60,7.85,7.02,6.56,6.18,6.12,5.22,4.49,4.35,3.97]]
# hp1=np.asarray(H_persoon[1])
# fH_persoon = interp1d(H_persoon[0], H_persoon[1], \
#             kind='linear', bounds_error=False, fill_value=0.0)
fH_persoon = CubicSpline(H_persoon[0], H_persoon[1], extrapolate=True)
# fH_persoon = UnivariateSpline(H_persoon[0], H_persoon[1], k=3, ext=3)
fn_persoon = CubicSpline(neq_persoon[0], neq_persoon[1], extrapolate=True)
# ----------------------------------------------------------------------------
#                          Cummings et al. (1969)
# ----------------------------------------------------------------------------
n0_Cummings = 1.E6  # m^-3
m0_Cummings = mp
BE_Cummings = 3.12E-5  # w/m^2

# ----------------------------------------------------------------------------
#                    Position of the center of the plasma sheet
#                      for nominal conditions at Saturn.
#                          (turning point of B_r)
#                         Used only for plotting.
# ----------------------------------------------------------------------------

RSHEET=[[ -2.03 ,  -2.07 ,  -2.111,  -2.155,  -2.2  ,  -2.247,  -2.296,\
        -2.348,  -2.401,  -2.458,  -2.516,  -2.578,  -2.643,  -2.71 ,\
        -2.781,  -2.856,  -2.934,  -3.017,  -3.104,  -3.195,  -3.292,\
        -3.394,  -3.503,  -3.617,  -3.739,  -3.869,  -4.008,  -4.155,\
        -4.314,  -4.484,  -4.667,  -4.866,  -5.082,  -5.319,  -5.58 ,\
        -5.87 ,  -6.197,  -6.569,  -7.   ,  -7.505,  -8.105,  -8.821,\
        -9.678, -10.72 , -12.043, -13.871, -16.59 , -20.323, -26.801,\
       -40.802],\
[-0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00,\
       -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00,\
       -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00,\
       -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00, -0.00e+00,\
       -1.00e-03, -1.00e-03, -1.00e-03, -1.00e-03, -1.00e-03, -1.00e-03,\
       -1.00e-03, -1.00e-03, -2.00e-03, -2.00e-03, -3.00e-03, -3.00e-03,\
       -4.00e-03, -5.00e-03, -6.00e-03, -8.00e-03, -1.10e-02, -1.50e-02,\
       -2.20e-02, -3.10e-02, -4.70e-02, -7.80e-02, -1.40e-01, -2.61e-01,\
       -6.70e-01, -1.78e+00],\
[ 0.039,  0.045,  0.063,  0.044,  0.057,  0.063,  0.062,  0.053,\
        0.068,  0.046,  0.067,  0.062,  0.041,  0.065,  0.062,  0.045,\
        0.062,  0.046,  0.035,  0.07 ,  0.062,  0.071,  0.038,  0.074,\
        0.069,  0.064,  0.02 ,  0.069,  0.04 ,  0.046,  0.059,  0.039,\
        0.041,  0.016, -0.012, -0.038, -0.037, -0.07 , -0.056, -0.103,\
       -0.121, -0.171, -0.256, -0.342, -0.479, -0.754, -1.2  , -2.068,\
       -4.016, -9.763]]

# ----------------------------------------------------------------------------
# Comparison of KMAG and dipole field line length and equatorial strength
# ----------------------------------------------------------------------------
# DIPOLE_L = [2.0, 3.4, 4.9, 6.4, 7.9, 9.3, 10.8, 12.3, 13.7, 15.2, 16.7,
#             18.2, 19.6, 21.1, 22.6, 24.1, 25.5, 27.0, 28.5, 30.0]
# DIPOLE_BEQ = [2356, 461, 161, 74, 40.0, 24.0, 15.5, 10.6, 7.5, 5.5,
#               4.2, 3.2, 2.6, 2.1, 1.7, 1.4, 1.1, 1.0, 0.8, 0.73]
# DIPOLE_Length = [3.46, 7.58, 11.66, 15.74, 19.82, 23.88, 27.96, 32.02,
#                  36.1, 40.16, 44.22, 48.30, 52.36, 56.44, 60.5, 64.58,
#                  68.64, 72.7, 76.78, 80.84]
# DIPOLE_BEQ_f = CubicSpline(DIPOLE_L, DIPOLE_BEQ, extrapolate=True)
# DIPOLE_Length_f = CubicSpline(DIPOLE_L, DIPOLE_Length, extrapolate=True)
# KMAG_L = [2.0, 3.2, 4.5, 5.8, 7.1, 8.4, 9.7, 11.0, 12.3, 13., 15.1, 16.5,
#           17.9, 19.3, 20.6]
# KMAG_BEQ = [2488, 571, 205, 90, 44.5, 25.4, 16.7, 11.7, 8.5, 6.5, 5.5, 5.06,
#             4.9, 4.9, 5.0]
# KMAG_Length = [3.5, 7.1, 10.5, 13.9, 17.0, 20.0, 23.0, 25.9, 28.9, 31.9,
#                35.1, 38.3, 41.8, 45.6, 49.7]
# KMAG_BEQ_f = CubicSpline(KMAG_L, KMAG_BEQ, extrapolate=True)
# KMAG_Length_f = CubicSpline(KMAG_L, KMAG_Length, extrapolate=True)


# ----------------------------------------------------------------------------


def eigenfunctions(z,
                   w,
                   mode=1,
                   l=L_RS_Yates / 2 * RS,
                   z0=z0_RS_Yates * RS,
                   A0=vA0_Yates,
                   A1=vA1_Yates,
                   C=1.):
    """ Compute the analytical eigenfunctions (Southwood & Kivelson, 1986)
    for a parallel inhomogeneity """
    soln = []
    if mode % 2 != 0:  # m=1 Even Mode
        for zval in z:
            if abs(zval) <= z0:
                soln.append(C * np.sin(w * (l - z0) / A1)
                              * np.cos(w * zval / A0))
            else:
                soln.append(C * np.sin(w * (l - abs(zval)) / A1)
                              * np.cos(w * z0 / A0))
    elif mode % 2 == 0:  # m=2 Odd Mode
        for zval in z:
            if abs(zval) <= z0:
                eigvalue = (C * np.sin(w * (l - z0) / A1)
                              * np.sin(w * abs(zval) / A0))
            else:
                eigvalue = (C * np.sin(w * (l - abs(zval)) / A1)
                              * np.sin(w * z0 / A0))
            if zval < 0:
                eigvalue *= -1
            soln.append(eigvalue)
    return soln
# ----------------------------------------------------------------------------


def f_analytical_root(wt, vA0, vA1, z0, l, mode=1):
    """ Compute the analytical roots (Southwood & Kivelson, 1986) """
    if mode % 2 != 0:  # m=1 (Even)
        LHS = (1. / vA0) * np.tan(wt * z0 / vA0)
        RHS = (1. / vA1) / np.tan(wt * (l - z0) / vA1)
    if mode % 2 == 0:  # m=2 (Odd)
        LHS = 1. / (vA0 * np.tan(wt * z0 / vA0))
        RHS = 1. / (-vA1 * np.tan(wt * (l - z0) / vA1))
    return LHS - RHS
# ----------------------------------------------------------------------------


def analyticalSolutions(z,
                        vA0,
                        vA1,
                        z0,
                        l,
                        rootEstimates=None,
                        mode=1,
                        verbose=False):
    """ Compute the analytical solutions (Southwood & Kivelson, 1986) """
    if rootEstimate is None:
        rootEstimate = [0.001E-3, 0.1E-3, 1E-3, 4E-3]
    roots = []
    solutions = []
    for m in range(1, mode + 1):
        # FIND ANALYTICAL SOLUTION ROOTS
        root = optimize.newton(f_analytical_root, rootEstimates[m - 1],
                               args=(vA0, vA1, z0, l, m))
        # CALCULATE ANALYTICAL EIGENFUNCTIONS
        solution = eigenfunctions(z, root, mode=m, l=l,
                                  z0=z0, A0=vA0, A1=vA1, C=1.)
        # SHOW RESULTS IN CONSOLE
        if verbose:
            print('Solution f(m = %2d) = %.3f mHz (w = %.3f mrads-1)'
                  % (m, root / 2 / np.pi * 1E3, root * 1E3))
        roots.append(root)
        solutions.append(solution)
    return solutions, roots

# ----------------------------------------------------------------------------


def searchAnalyticalSolutions(z, vA0, vA1, z0, l,
                              rootMax=0.1E-2,
                              parameterBins=100,
                              mode=1,
                              nrootlim=3,
                              verbose=False):
    """ Look for roots in the PDE's of (Southwood & Kivelson, 1986) """
    roots = []
    solutions = []
    errors = []
    parameterBounds = [0., rootMax]
    parameterSpace = np.linspace(parameterBounds[0], parameterBounds[1],
                                 parameterBins)
    error1 = f_analytical_root(parameterSpace[0], vA0, vA1, z0, l, mode=mode)
    errors.append(error1)
    derivative1 = 0
    for i in range(1, len(parameterSpace)):
        # Break if all requested roots are found
        if len(roots) > nrootlim:
            break
        error2 = f_analytical_root(parameterSpace[i], vA0, vA1, z0, l,
                                   mode=mode)
        derivative2 = ((error2 - error1)
                       / (parameterSpace[i] - parameterSpace[i - 1]))
        if (error2) * (error1) < 0:
            if abs(derivative2 / derivative1) > 0.5 and \
               abs(derivative2 / derivative1) < 1.5:
                roots.append(0.5 * (parameterSpace[i]
                                    + parameterSpace[i - 1]))
        errors.append(error2)
        error1 = error2
        derivative1 = derivative2

    return roots, parameterSpace, errors
# ----------------------------------------------------------------------------


def fde(z, y, SIM):
    """ Function for the Alfven (standing) Wave PDE

    DEFINITIONS:
    y[0] = y
    y[1] = y'
    y[2] = w (frequency)

    The first column is simply the first derivative, y'
    The second column is the expression for y'', in terms of y', y, and z.
    i.e.

    (y', dy_factor * y' - y_factor * y (w / vA)^2)

    where dy_factor and y_factor are determined above, and depend on the
    geometry of the field line.
    """

    # Singer et al. (1981) - solution along the field line (distance = s)
    if SIM.coords == 'ds':
        yfactor = 1.
        # dyfactor = 0.
        if SIM.component == 'toroidal':
            dyfactor = SIM.dlnh1Bfun(z)
        if SIM.component == 'poloidal':
            dyfactor = SIM.dlnh2Bfun(z)
    # Cummings et al. (1969) - dipole solution as a function of cos(theta)
    elif SIM.coords == 'cos':
        r0 = SIM.L * SIM.units
        yfactor = np.power(r0, 2.) * (1. + 3. * np.power(z, 2.))
        if SIM.component == 'toroidal':
            dyfactor = 0.
        if SIM.component == 'poloidal':
            dyfactor = - 6. * z / (1. + 3. * np.power(z, 2.))
    # Magnetospheric box solution (uniform field)
    elif SIM.coords == 'cartesian':
        dyfactor = 0.
        yfactor = 1.

    return np.array([y[1],
                     - dyfactor * y[1]
                     - yfactor * y[0] * np.power(y[2] / SIM.vAfun(z), 2.),
                     0.0])

# ----------------------------------------------------------------------------


def displacement(solution, SIM):
    """ Return displacement from the solution """
    if SIM.coords == 'cos':
        return solution.y
    elif SIM.coords == 'ds':
        if SIM.component == 'toroidal':
            h = SIM.h1
        if SIM.component == 'poloidal':
            h = SIM.h2
        return solution.y * h
    else:
        return solution.y

# ----------------------------------------------------------------------------


def electricPerturbation(solution, SIM):
    """ Return Electric field perturbation """
    if SIM.coords == 'cos':
        z = SIM.z
        if SIM.component == 'toroidal':
            E_nu = np.divide(np.power(1. + 3. * np.power(z, 2.), 0.5),
                             SIM.L ** 2
                             * np.power(1. - np.power(z, 2.), 3. / 2.))
            return np.multiply(solution.xi, E_nu)
        if SIM.component == 'poloidal':
            E_phi = np.divide(1., SIM.L * np.power(1. - np.power(z, 2.),
                                                   3. / 2.))
            return np.multiply(solution.xi, E_phi)
    if SIM.coords == 'ds':
        return -solution.roots * solution.xi * SIM.B
    else:
        h = 1.0
        return np.multiply(solution.y, h)

# ----------------------------------------------------------------------------


def magneticPerturbation(solution, SIM):
    """ Return Magnetic Field Pertubation """
    if SIM.coords == 'cos':
        # z = np.cos(xaxis*np.pi/180.)
        z = SIM.z
        w = solution.roots
        if SIM.component == 'poroidal':
            B_nu = np.divide(1., w * np.power(SIM.L, 2.)
                   * np.power(1. + 3 * np.power(z, 2.), 1. / 2.)
                   * np.power(1. - np.power(z, 2.), 3. / 2.))
            return np.multiply(solution.dy, B_nu)
        if SIM.component == 'toroidal':
            B_phi = np.divide(-1., w * np.power(SIM.L, 3.)
                    * np.power(1. - np.power(z, 2.), 3. / 2.))
            return np.multiply(solution.dy, B_phi)
    if SIM.coords == 'ds':
        if SIM.component == 'toroidal':
            h = SIM.h1
        if SIM.component == 'poloidal':
            h = SIM.h2
        return h * SIM.B * solution.dy
    else:
        B = np.divide(-1., w * np.power(SIM.L, 2.))
        return np.multiply(np.gradient(solution), B)

# ----------------------------------------------------------------------------


def uniformField(R, SIM):
    """ Return a grid with values for a unform field """
    M = SIM.dipM
    # M = BT_Yates
    R = np.asarray(R)
    if np.size(R) > 3:
        R = np.transpose(R)
    # X,Y,Z = [np.atleast_1d(R[i]) for i in range(0,3)]

    r = np.linalg.norm(R, axis=0)
    X, Y, Z = [R[i] for i in range(0, 3)]
    return np.asarray([np.zeros_like(X),
                       np.zeros_like(Y),
                       np.full_like(Z, M)])

# ----------------------------------------------------------------------------


def dipField(R, SIM, external=False):
    """ Return the magnitude of a dipole field """
    M = SIM.dipM
    R = np.asarray(R)
    if np.size(R) > 3:
        R = np.transpose(R)
    r = np.linalg.norm(R, axis=0)
    X, Y, Z = [R[i] for i in range(0, 3)]

    BX = np.divide(3. * M * X * Z, np.power(r, 5.))
    BY = np.divide(3. * M * Y * Z, np.power(r, 5.))
    BZ = np.divide(M * (3. * np.power(Z, 2.) - np.power(r, 2.)),
                   np.power(r, 5.))

    if external:
        sigma = 4.  # RS
        A_ext = CURRENT_SHEET_STRENGTH * 250E-9
        BX_ext = (A_ext
                  * (Z / np.power(sigma, 2))
                  * (1. / (sigma * np.sqrt(2. * np.pi)))
                  * np.exp(- 0.5 * np.power((Z - 0.) / sigma, 2.)))
        # BX_ext = 0
        BY_ext = np.zeros_like(BX_ext)
        BZ_ext = np.zeros_like(BX_ext)

        BX = BX + BX_ext
        BY = BY + BY_ext
        BZ = BZ + BZ_ext

    B = np.asarray([BX, BY, BZ])
    return B.T

# ----------------------------------------------------------------------------


def calculate_vA(B, n, m, SIM=None, maxv=c):
    """ Return the Alfven velocity for a given B, n, and m """
    n = np.nan_to_num(n)
    with np.errstate(divide='ignore'):
        vA = np.divide(np.atleast_1d(B), np.sqrt(mu0 * np.atleast_1d(n) * m))

    # Relativistic correction
    vA = np.sqrt(np.divide(np.power(vA, 2) * c ** 2,
                           np.power(vA, 2) + c ** 2))
    vA = np.nan_to_num(vA, nan=c, posinf=c, neginf=c)

    # vA[vA > c] = c
    # vA[vA < 0] = c
    return vA

# ----------------------------------------------------------------------------


def densityModels(R, SIM, fieldline=None):
    """ Return the plasma density at a given location R """
    R = np.asarray(R)

    if np.size(R) > 3:
        R = np.transpose(R)

    r = np.linalg.norm(R, axis=0)
    rho = np.linalg.norm(R[0:2], axis=0)
    X, Y, Z = [np.atleast_1d(R[i]) for i in range(0, 3)]

    # ------------------------------------------------------------------------
    # Boxcar Function Plasma Sheet, Yates et al. [2017]
    # ------------------------------------------------------------------------
    if SIM.densityModelName == 'yates':
        # L = np.full_like(r, fieldline.L)
        L = fieldline.L
        neq = fn_bagenal(L) * 1E6
        H = fH_bagenal(L)
        z0 = H * 0.6
        n0 = neq
        # n0 = n0_Yates
        # z0 = z0_RS_Yates
        n = [n0 if (z < z0 * 2 and z > -z0 * 2) else n1_Yates for z in Z]
        # n = [SIM.dipM if (z < z0_RS_Yates 
               # and z > -z0_RS_Yates) else n1_Yates for z in Z]
        n = np.asarray(n)
        n[np.abs(X) < 2.] = n1_Yates

    # ------------------------------------------------------------------------
    # Boxcar Function Plasma Sheet, Yates et al. [2017], decreasing with x
    # ------------------------------------------------------------------------
    elif SIM.densityModelName == 'yates_decreasing':
        n = [n0_Yates * np.exp(-np.power(np.abs(x) - 10.0, 2.)
                               / (20.)) for x in X]
        n = np.asarray(n)
        n[(Z > -z0_RS_Yates) & (Z < z0_RS_Yates)] = n1_Yates
        n[np.abs(X) < 2.] = n1_Yates
        n[n < n1_Yates] = n1_Yates
    # ------------------------------------------------------------------------
    # Persoon et al. [2005]
    # ------------------------------------------------------------------------
    elif SIM.densityModelName == 'persoon':
        magLat = np.arctan2(Z, rho)
        L = 0.
        if SIM.BFieldModelName == 'dipole':
            L = np.divide(r, np.power(np.cos(magLat), 2.))
        else:
            L = SIM.L
        # Density model equation (approximation)
        # neq = np.divide(n0_Persoon * 2. , \
            # ( np.power(R0_Persoon / L, m_Persoon) + \
            # np.power(L / R0_Persoon, n_Persoon)))
        # Density data interpolation
        neq = fn_persoon(L) * 1e6
        # H = 0.047*np.power(rho, 1.8)  # [Persoon et al., 2006]
        # H = 0.047*np.power(rho, 1.5)  # [Persoon et al., 2013]
        # H = 0.047*np.power(L, 1.5)  # [Persoon et al., 2013]
        H = fH_persoon(L)
        ne = neq * np.exp(-np.power(L, 2.)
                          * (1. - np.power(np.cos(magLat), 6.))
                          / 3.0 / np.power(H, 2.))
        n = ne
    # -----------------------------------------------------------------------
    # Bagenal & Delamere [2011]
    # ------------------------------------------------------------------------
    elif SIM.densityModelName == 'bagenal':
        n = []

        def f_s(r, th, L):
            """ Calculate distance along a dipole field line """
            # s = np.abs(L/4.0*(2.*th + np.sin(2.*th)))
            # Split the field line by latitude into many straight segments
            th_range = np.linspace(0., np.abs(th), 50)
            s = 0.
            r_0 = L
            th_0 = 0.
            # A rather slow method, just for double checking
            for th_ in th_range:
                r_ = L * np.power(np.cos(th_), 2.)  # Dipole field line
                R_0 = np.asarray([r_0 * np.cos(th_0), 0., r_0 * np.sin(th_0)])
                R_ = np.asarray([r_ * np.cos(th_), 0., r_ * np.sin(th_)])
                ds = np.linalg.norm(R_0 - R_, axis=0)
                s += ds
                r_0 = np.linalg.norm(R_, axis=0)
                th_0 = th_
            return s

        if SIM.BFieldModelName == 'dipole':
            magLat = np.arctan2(Z, rho)
            L = np.divide(r, np.power(np.cos(magLat), 2.))
            L = np.atleast_1d(L)
            r = np.atleast_1d(r)
            magLat = np.atleast_1d(magLat)
            s = [f_s(r[i], magLat[i], L[i]) for i in range(0, np.size(r))]
            s = np.asarray(s)
        else:
            L = np.full_like(r, fieldline.L)
            s = np.abs(fieldline.z / SIM.units)

        neq = fn_bagenal(L) * 1E6
        H = fH_bagenal(L)

        # Scaling adjustment to test sensitivity of the results
        if SIM.scalingFactor is not None:
            neq *= SIM.scalingFactor
            # H *= SIM.scalingFactor

        n = neq * np.exp(-np.power(s / H, 2.))
    # ------------------------------------------------------------------------
    # Cummmings & Coleman [1969]
    # ------------------------------------------------------------------------
    elif SIM.densityModelName == 'cummings':

        n = n0_Cummings * np.power(SIM.L / r, SIM.mIndex)

    # If density model name wasn't recognized, choose n = 1 m-3
    else:
        np.full_like(X, 1.0)

    return np.asarray(n)

# ----------------------------------------------------------------------------


if __name__ == "__main__":
    from pylab import *
