#!/usr/bin/env python3
# ----------------------------------------------------------------------------
#          _  ____  __    _    ____ _          _
#         | |/ /  \/  |  / \  / ___| |__   ___| |_ __   ___ _ __
#         | ' /| |\/| | / _ \| |  _| '_ \ / _ \ | '_ \ / _ \ '__|
#         | . \| |  | |/ ___ \ |_| | | | |  __/ | |_) |  __/ |
#         |_|\_\_|  |_/_/   \_\____|_| |_|\___|_| .__/ \___|_|
#                                               |_|
#
# Written by:                             Liutauras (Leo) Rusaitis
#                                         10-13-2020
#
#                                         Space Physics PhD Student,
#                                         Earth, Planetary, and Space Sciences,
#                                         University of California, Los Angeles
#                                         GitHub: https://github.com/rusaitis
#                                         Contact: rusaitis@ucla.edu
# ----------------------------------------------------------------------------
#
# Useful functions for interfacing with the KMAG magnetic field model and
# tracing procedures.
#
# ----------------------------------------------------------------------------

import numpy as np
import os

def saveKmagParameters(config,
                       COMMENT='KMAG',
                       filename='./KMAGhelper/kmag_params.txt'):
    """ Save KMAG and tracing parameters """
    ETIME = config["ETIME"]
    BY_IMF = config["BY_IMF"]
    BZ_IMF = config["BZ_IMF"]
    Dp = config["Dp"]
    STEP = config["step"]
    IN_COORD = config["IN_COORD"]
    OUT_COORD = config["OUT_COORD"]
    OUT_CARSPH = config["OUT_CARSPH"]
    TRACE = config["TRACE"]

    CARSPH = OUT_CARSPH
    PARAMS_fields = ['ETime(sec)', 'BY_IMF', 'BZ_IMF', 'Dp', 'IN_COORD',
                     'OUT_COORD', 'COMMENT', 'TRACE(Y/N)', 'CARSPH', 'STEP']
    PARAMS = [ETIME, BY_IMF, BZ_IMF, Dp, IN_COORD,
              OUT_COORD, COMMENT, TRACE, CARSPH, STEP]
    PARAMS = np.vstack([PARAMS_fields, PARAMS])
    np.savetxt(filename, PARAMS, fmt='%14s', delimiter=' ')

# ----------------------------------------------------------------------------


def traceKMAG(traceStartList, config=None):
    """ Write out the trace starting points and call the KMAG tracing code """
    ETIME = config["ETIME"]
    DATA_ROWS = []  # DATA FOR KMAG INPUT
    DATA_FIELDS = ['Etime(sec)', 'X(RS)', 'Y(RS)', 'Z(RS)']
    for R in traceStartList:
        X = R[0]
        Y = R[1]
        Z = R[2]
        DATA_ROWS.append([ETIME, X, Y, Z])

    DATA_SAVE = np.vstack([DATA_FIELDS, DATA_ROWS])
    np.savetxt('./KMAGhelper/kmag_input.txt',
               DATA_SAVE,
               delimiter=' ',
               fmt='%14s')

    if config:
        i = 0
        config["COMMENT"] = "{0:0=3d}".format(i)
        saveKmagParameters(config)

    # -w to inhibit Fortran2018 warnings
    os.system("gfortran -c ./KMAGhelper/KMAGtracer.f "
              + "-w -o ./KMAGhelper/KMAGtracer.o")
    os.system("gfortran -c ./KMAG2012/KMAG2012.f "
              + "-w -o ./KMAG2012/KMAG2012.o")
    os.system("gfortran -o ./KMAGhelper/KMAG "
              + "./KMAGhelper/KMAGtracer.o "
              + "./KMAG2012/KMAG2012.o")
    os.system("./KMAGhelper/KMAG")

# ----------------------------------------------------------------------------


def KMAGfieldlines(SIM, ioconfig=None):
    """ Read in the field lines generated in the tracing procedure """
    path_output = ioconfig.path if ioconfig is not None else 'Output'
    for filename in os.listdir(path_output):
        if ((os.path.splitext(filename)[1][1:] == 'txt')
            and ('eigenmodes' not in filename)):
            DATA = np.loadtxt(path_output + '/' + filename,
                              skiprows=1, dtype='U')
            DESCRIPTIONS = DATA[:, -1]
            e_ind = np.where(DESCRIPTIONS == 'EVENT')[0]

            DATA_START = DATA[:e_ind[1]]
            DATA_END = np.flip(DATA[e_ind[1]:])
            DATA_END = DATA[e_ind[1]:]
            DATA_END = DATA_END[::-1]
            DATA = np.concatenate((DATA_END, DATA_START), axis=0)

            TIME = DATA[:, 0].astype(str)
            X = DATA[:, 1].astype(float)
            Y = DATA[:, 2].astype(float)
            Z = DATA[:, 3].astype(float)
            BX = DATA[:, 4].astype(float) * 1E-9  # From nT
            BY = DATA[:, 5].astype(float) * 1E-9  # From nT
            BZ = DATA[:, 6].astype(float) * 1E-9  # From nT
            BT = DATA[:, 7].astype(float) * 1E-9  # From nT
            R = DATA[:, 1:4].astype(float)
            B = DATA[:, 4:7].astype(float)
    return R, B, BT