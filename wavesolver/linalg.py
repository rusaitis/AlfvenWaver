#!/usr/bin/env python3
import numpy as np

# ----------------------------------------------------------------------------
# Useful Linear Algebra Functions
# ----------------------------------------------------------------------------

def normalize(y, A=1.0):
    """ Normalize an array (or provide an amplitude to normalize with) """
    norm = np.amax(np.abs(y))
    if norm > 0:
        return A * np.asarray(y) / norm
    else:
        return np.asarray(y)

# ----------------------------------------------------------------------------


def normalize2d(y, A=1.0):
    """ Normalize a 2D array """
    y = np.atleast_2d(y)

    norm = [np.amax(np.abs(yy)) for yy in y]
    norm = np.asarray(norm)
    norm = np.atleast_1d(norm)
    norm[norm > 0] = 1. / norm[norm > 0]

    return A * (np.asarray(y).T * norm).T

# ----------------------------------------------------------------------------


def car2sph(R):
    """ Cartesian to spherical coordinates """
    R = np.asarray(R)

    if np.size(R) > 3:
        R = np.transpose(R)

    r = np.linalg.norm(R, axis=0)
    rho = np.linalg.norm(R[0:2], axis=0)
    X, Y, Z = [R[0], R[1], R[2]]
    R_polar = [r,
               np.arctan2(rho, Z),
               np.arctan2(Y, X)]

    return np.asarray(R_polar).T

# ----------------------------------------------------------------------------


def sph2car(R):
    """ Spherical to cartesian coordinates """
    R = np.asarray(R)
    if np.size(R) > 3:
        R = np.transpose(R)
    X, Y, Z = [R[0], R[1], R[2]]
    R_cartesian = [X * np.sin(Y) * np.cos(Z),
                   X * np.sin(Y) * np.sin(Z),
                   X * np.cos(Y)]
    return np.asarray(R_cartesian).T
# ----------------------------------------------------------------------------


def unitVector(R):
    """ Return a unit vector """
    R = np.asarray(R)
    norm = np.linalg.norm(R, axis=-1)
    norm = np.atleast_1d(norm)
    norm[norm > 0] = 1. / norm[norm > 0]
    return (np.asarray(R).T * norm).T

# ----------------------------------------------------------------------------


def perpendicularUnitVectors(R, B):
    """ Return two perpendicular unit vectors to field vector B at position R.

    INPUT: R (position vector)
           B (magnetic field vector)

    The azimuthal vector, phi, is computed using R.

    OUTPUT: PERP1_UNIT = B_unit X phi_unit
            PERP2_UNIT = B_unit X perp1_unit
    """

    R = np.asarray(R)
    B = np.asarray(B)

    # BT = np.linalg.norm(B, axis=-1)
    B_unit = unitVector(B)
    # R_unit = unitVector(R)
    phi_unit = unitVector(np.asarray([-R.T[1],
                                      R.T[0],
                                      np.zeros_like(R.T[0])]).T)
    perp1 = np.cross(B_unit, phi_unit)
    perp1_unit = unitVector(perp1)
    perp2 = np.cross(B_unit, perp1_unit)
    perp2_unit = unitVector(perp2)
    return perp1_unit, perp2_unit

# ----------------------------------------------------------------------------


def perpendicularStep(R, B, step=1.):
    """ Returns positions of two perpendicular steps to field B at position R.
    The step size can be specified (default = 1.)
    """
    R = np.asarray(R)
    perp1_unit, perp2_unit = perpendicularUnitVectors(R, B)
    step = np.asarray(step)
    R1 = np.add(R, np.multiply(perp1_unit.T, step).T)
    R2 = np.add(R, np.multiply(perp2_unit.T, step).T)
    return R1, R2

# ----------------------------------------------------------------------------


if __name__ == "__main__":
    from pylab import *

    ###### TEST: unitVector(B) and perpendicularUnitVectors(R, B)
    # R = [[-1., 0., 0.],[1., 0., 0.], [2., 0., 0.], [3., 0., 0.]]
    # R = [-1., 0., 0.]
    # B = [[0., 0., 1.],[2., 0., 2.], [0., 0., 1.], [-1., -1., 1.]]
    # B = [0., 0., 2.]
    # R = np.asarray(R)
    # B = np.asarray(B)
    # uv = unitVector(B)
    # p1, p2 = perpendicularUnitVectors(R, B)
    # print('p1: ', p1)
    # print('p2: ', p2)

    ###### TEST: normalize(y)
    # y = [[10., 5., 0., 0., 5.], [0., 0., 0., 0., 0.]]
    # y = [0., 0., 0., 5., 2.5]
    # print('y: ', y)
    # y = normalize(y, A = 2.0)
    # print('y_norm: ', y)

    ###### TEST: angular2mins(w)
    # w = [10., 20., 30., 0., 50.]
    # tmin = angular2mins(w)
    # print('tmin: ', tmin)

    ###### TEST: sph2car(R)
    # R = [[1, 0., 0.], [1., np.pi/2., np.pi/4.]]
    # R = [1, 0., 0.]
    # R = sph2car(R)
    # print('R: ', R)