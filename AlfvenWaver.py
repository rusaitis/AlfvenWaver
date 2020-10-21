#!/usr/bin/env python3
# ------------------------------------------------------------------------------
#           _    _  __              __        __
#          / \  | |/ _|_   _____ _ _\ \      / /_ ___   _____ _ __
#         / _ \ | | |_\ \ / / _ \ '_ \ \ /\ / / _` \ \ / / _ \ '__|
#        / ___ \| |  _|\ V /  __/ | | \ V  V / (_| |\ V /  __/ |
#       /_/   \_\_|_|   \_/ \___|_| |_|\_/\_/ \__,_| \_/ \___|_|
#
# Written by:                             Liutauras (Leo) Rusaitis
#                                         10-13-2020
#
#                                         Space Physics PhD Student,
#                                         Earth, Planetary, and Space Sciences,
#                                         University of California, Los Angeles
#                                         GitHub: https://github.com/rusaitis
#                                         Contact: rusaitis@ucla.edu
# ------------------------------------------------------------------------------
#                                Description
#
# AlfvenWaver is a command-line interfaced program package that is aimed at
# computing shear Alven wave resonances along planetary magnetic field lines.
#
# Currently, it is using a package called WaveSolver, which is solely
# focused on solving Alfven standing modes along magnetic field lines.
#
#
# ------------------------------------------------------------------------------
import argparse
from wavesolver import *
from wavesolver.model import *
from wavesolver.sim import *
from wavesolver.plot import *
from wavesolver.fieldline import *
from wavesolver.io import *
from KMAGhelper.KmagFunctions import *
from wavesolver.configurations import *

if __name__ == "__main__":

    # --------------------------------------------------------------------------
    # Command line aguments
    parser = argparse.ArgumentParser(description='Process the waves.')
    parser.add_argument('--calc', dest='calculate', action='store_true',
                        help='Compute the Alfven waves')
    parser.add_argument('--sol', dest='plotSolutions',
                        action='store_true', help='Plot the eigenfunctions')
    parser.add_argument('--field', dest='plotFieldLines',
                        action='store_true', help='Plot the field lines')
    parser.add_argument('--eig', dest='plotEigenfrequencies',
                        action='store_true', help='Plot the eigenfrequencies')
    parser.add_argument('--errors', dest='plotErrorFunc',
                        action='store_true', help='Plot the error function')
    parser.add_argument('--mov', dest='movie',
                        action='store_true', help='Generate a movie')
    parser.add_argument('--save', dest='save',
                        action='store_true', help='Save the figures')
    parser.add_argument('--v', dest='verbose',
                        action='store_true', help='Verbose mode')
    parser.add_argument('--dipole', dest='dipole',
                        action='store_true', help='Use a dipole field')
    parser.add_argument('--uniform', dest='uniform',
                        action='store_true', help='Use a uniform field')
    parser.add_argument('--plotRefEigs', dest='plotRefEigs',
                        action='store_true', help='Plot reference '
                                                  + 'eigenfrequencies')
    parser.add_argument('--conf', dest='configuration', type=int,
                        default=0,
                        help='Field line configuration number')
    parser.add_argument('--plotvar', dest='plotvar', type=str,
                        default='vA',
                        choices=['Xi', 'E', 'b', 'bratio', 'vA', 'n'],
                        help='Variable to plot along the field lines')
    parser.add_argument('--plotCoord', dest='plotCoord', type=str,
                        default='th', choices=['th', 's'],
                        help='Field line coordinate to plot against: '
                             + 'latitude (th) or distance along the field (s)')
    parser.add_argument('--component', dest='component', type=str,
                        default='toroidal', choices=['toroidal', 'poloidal'],
                        help='Toroidal or poloidal component for solution')
    parser.add_argument('--modes', dest='modes', type=int,
                        default=1,
                        help='Number of modes to solve')
    parser.add_argument('--L', dest='L', type=int, nargs='*',
                        # default=[20],
                        help=('Trace a field line from an equatorial crossing '
                              + 'distance, L (in planet radii). (Multiple '
                              + 'values accepted)'))
    parser.add_argument('--th', dest='th', type=int, nargs='*',
                        # default=[18],
                        help=('Trace a field line from a colatitude of '
                              + 'theta (in degrees) staring in the northern '
                              + 'hemisphere (multiple values accepted)'))
    parser.add_argument('--phi', dest='phi', type=int, nargs='*',
                        # default=[0],
                        help=('Trace a field line from the azimuth of '
                              + 'phi (in degrees) staring from the sunward '
                              + 'direction (CCW). (Multiple values accepted)'))
    parser.add_argument('--fast', dest='fast', action='store_true',
                        help='Enable a faster calculation (less accurate)')
    parser.add_argument('--nfl', dest='nfl', type=int, default=1,
                        help='Number of field lines in latitude / radial range')
    parser.add_argument('--plotx', dest='plotx', type=int,
                        default=[-20, 20], nargs=2,
                        help='X limits of the plot')
    parser.add_argument('--ploty', dest='ploty', type=int,
                        default=[-14, 14], nargs=2,
                        help='Z limits of the plot')
    parser.add_argument('--Lmax', dest='Lmax', type=int,
                        default=21,
                        help=('Furthest field line for which to plot'
                              + ' field line parameters, measured by equatorial'
                              + ' crossing distance, L (in planet radii).'))
    args = parser.parse_args()
    # --------------------------------------------------------------------------

    # Set up the output folder and other plotting parameters
    ioconfig = io(path='./Output',
                  plotCoord=args.plotCoord,  # 's' | 'th'
                  save=False,
                  plotx=args.plotx,
                  ploty=args.ploty,
                  Lmax=args.Lmax,
                  plotRefEigs=args.plotRefEigs)

    # Delete old *.png files in the output directory
    deleteFiles(ioconfig.path, fileformat='.pdf')

    # Load Saturn configuration (browse configurations.py for more options)
    SIM = loadsim(configSaturnNominal)

    if args.dipole:
        SIM.BFieldModelName = 'dipole'

    SIM.updateFieldLines(configuration=args.configuration,
                         NFieldLines=args.nfl,
                         L=args.L,
                         th=args.th,
                         phi=args.phi,
                         component=args.component)

    if args.fast:
        SIM.update(configFastSolve)

    modes_to_calculate = 1  # Number of modes to calculate
    save_eigenfrequencies = False  # Do not export eigenfrequencies by default

    # Batches are used for changing one or more parameters between similar sets
    batches = batchNominal  # DEFAULT: no batches (nominal)
    SIMSbatch = []
    fieldlinesbatch = []

    # For testing basic configuration (e.g. equatorial density)
    # plotBasicConfiguration(SIM)

    # --------------------------------------------------------------------------
    for batch_n in range(0, len(batches)):
        # Update the simulation parameters with the new batch parameters
        SIM.update(batches[batch_n])

        if len(batches) > 1:
            print('Batch id: ', batch_n)

        # Update output configuration
        ioconfig.id = batch_n
        ioconfig.name = SIM.name
        # COMMENT = "{0:0=3d}".format(batch_n)

        # Trace field lines
        fieldlines = traceFieldlines(SIM)

        # Find the plasma sheet
        SIM.sheet = findPlasmaSheet(fieldlines)

        # Store the field line information in simulation (SIM) objects
        SIMS = configureSIMS(SIM, fieldlines)

        # For testing purposes: plotting the scalling factors, h1 and h2
        # fieldlineScalingFactorPlot(fieldlines[0])
        # exit()

        # Compute the standing Alfven waves
        if args.calculate:
            # Each simulation (SIM) can have multiple field lines configured
            # -> compute the standing waves for each individually
            for n, SIM in enumerate(SIMS):
                SIM.id = n
                SIM.modes = args.modes
                print('Calculating field line #%3d | L = %3.1f'
                      % (n, SIM.L))

                # Compute the errors for different values of w (ang freq)
                errors, LOG = rootErrorFunction(SIM,
                                                VERBOSE=args.verbose,
                                                plotLOG=False
                                                )

                # Compute the solutions using a shooting method
                # NOTE: eigenfunctions are normalized by default
                solutions = errors.shootSolutions(SIM,
                                                  VERBOSE=args.verbose,
                                                  normalize=True,
                                                  normfirst=True,
                                                  normfirsts=True,
                                                  scaleFactor=None
                                                  )
                SIM.solutions = solutions

                # Plot error function (useful if solutions are not found)
                if args.plotErrorFunc:
                    plotErrorFunctions(LOG,
                                       SIM,
                                       ioconfig
                                       )
                # Plot various parameters along the fiel (for testing):
                # plotFluxTubeParameters(SIM, ax, ioconfig, plotVar='n')

                # Plot the eigenfunctions
                if args.plotSolutions:
                    plotSolutionSubplots(SIMS,
                                         fieldlines,
                                         ioconfig,
                                         plotBG='vA'
                                         )
                    # plotSolutionsPub(SIM, fig, ax, ioconfig, plotVar='b')
                    # plotSubplots(SIMS, fieldlines, ioconfig)

        # Plot field lines and other data like density or Alfven velocity
        if args.plotFieldLines:
            plotConfiguration(SIMS,
                              fieldlines=fieldlines,
                              plotData=args.plotvar,
                              ioconfig=ioconfig,
                              plotRsheet=True
                              )
        # Export eigenfrequencies into a file
        if save_eigenfrequencies:
            extractEigenfrequencies(SIMS, save=True)

        # Plot eigenfrequencies versus inv lat or equatorial crossing distance
        if args.plotEigenfrequencies:
            plotEigenfrequenciesPub(SIMS, ioconfig)

        # Collect the simulations from this batch
        SIMSbatch.extend(SIMS)
        fieldlinesbatch.extend(fieldlines)

    # Process a movie from figures (useful for batches)
    if args.movie:
        filenames = ['solution', 'wL', 'fieldlines', 'kmag']
        # os.path.join(dirFile, ('movie_' + filnames[3] + '.pdf'))
        processMovie(path=ioconfig.path,
                     filename=filenames[3],
                     res='5136x2880',
                     fps=15,
                     output=ioconfig.path + '/movie_' + filenames[3] + '.mp4'
                     )
