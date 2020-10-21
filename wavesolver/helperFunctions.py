#!/usr/bin/env python3

import numpy as np

# ----------------------------------------------------------------------------


def sortEigenmodes(eigfList, modes):
    """ Sort eigenfrequency modes for easier plotting """
    eigs_sorted = []
    for m in range(1, modes + 1):
        m_eigs = []
        for eigs in eigfList:
            if len(eigs) >= m:
                m_eigs.append(eigs[m - 1])
            else:
                m_eigs.append(0.)
        eigs_sorted.append(m_eigs)
    return eigs_sorted

# ----------------------------------------------------------------------------


def extractEigenfrequencies(SIMS, save=False):
    """ Extract and export eigenfrequencies """
    L = []
    lengths = []
    LAT = []
    eigs = []
    eigsMins = []

    for i, SIM in enumerate(SIMS):
        L.append(SIM.L)
        lengths.append(SIM.length)
        LAT.append(90. - np.rad2deg(SIM.rthetaphi[0, 1]))
        roots = []
        for solution in SIM.solutions:
            roots.append(solution.roots)
        roots = np.asarray(roots)
        eigs.append(angularFreqConversion(roots, format='fmHz'))
        eigsMins.append(angularFreqConversion(roots, format='Tmin'))
    eigs = sortEigenmodes(eigs, SIM.modes)

    if save:
        labels = ['L (RS)', 'LAT (degrees)']
        labelsModes = ['m=%d' % i for i in range(0, len(eigs))]
        labels = np.hstack((labels, labelsModes))
        labels = ['{:^16s}'.format(label) for label in labels]
        vstack = np.vstack((L, LAT, eigs))
        if SIM.name is not None:
            filename = SIM.name + '_eigenmodes.txt'
        else:
            filename = 'eigenmodes.txt'
        np.savetxt('Output/' + filename,
                   vstack.T,
                   fmt="%15.8f",
                   header="".join(labels))

    return L, lengths, LAT, eigs
# ----------------------------------------------------------------------------


def loadEigenfrequencies(filename=None):
    """ Load eigenfrequencies from an external file """

    if filename is not None:
        filename = filename + '_eigenmodes.txt'
    else:
        filename = 'eigenmodes.txt'

    data = np.loadtxt('Output/' + filename, skiprows=1).T
    L = data[0]
    LAT = data[1]
    LAT0 = data[1]
    eigs = data[2:]
    # L = np.asarray(L)
    # LAT = np.asarray(LAT)
    # eigs = np.asarray(eigs)
    LAT_LIM = 74.8
    if 'Night' in filename:
        LAT = LAT[LAT0 < LAT_LIM]
        L = L[LAT0 < LAT_LIM]
        eigs = eigs[:, LAT0 < LAT_LIM]
    return L, LAT, eigs

# ----------------------------------------------------------------------------


def angularFreqConversion(w, format='fmHz'):
    if format == 'f':
        return w / 2. / np.pi
    if format == 'fmHz':
        return w / 2. / np.pi * 1E3
    if format == 'T':
        return np.divide(1, w / 2. / np.pi)
    if format == 'Tmin':
        return np.divide(1, w / 2. / np.pi) / 60.

# ----------------------------------------------------------------------------


def angular2freq(angfreq, units=1E3):
    """ Angular frequency to frequency conversion (in specifient units) """
    angfreq = np.asarray(angfreq)
    freq = np.zeros_like(angfreq)
    freq[angfreq > 0] = angfreq[angfreq > 0] / 2. / np.pi * units
    return freq
# ----------------------------------------------------------------------------


def angular2mins(angfreq, time_units=60.):
    """ Angular frequency to period in minutes """
    angfreq = np.asarray(angfreq)
    Tmin = np.zeros_like(angfreq)
    Tmin[angfreq > 0] = 1. / (angfreq[angfreq > 0] / 2. / np.pi) / time_units
    return Tmin
    # return 1./(angfreq/2./np.pi)/60.
    # eigenfrequencyListMinuntes.append(1/(roots[1]/2/np.pi)/60.)

# ----------------------------------------------------------------------------


def toArray(A, datatype="float"):
    # Get the dimension of t and make sure that t is an n-element vector
    if type(A) != np.ndarray:
        if type(A) == list:
            A = np.array(A)
        else:
            if datatype == "float":
                A = np.array([float(A)])
            if datatype == 'int':
                A = np.array([int(A)])
            if datatype == 'str':
                A = np.array([str(A)])
    return A

# ----------------------------------------------------------------------------


if __name__ == "__main__":
    from pylab import *
