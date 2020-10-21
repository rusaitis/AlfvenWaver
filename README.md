# KMAG

## Introduction

AlfvenWaver is a command-line interfaced program package that is aimed at
computing shear Alven wave resonances along planetary magnetic field lines.

Currently, it is using a package called WaveSolver, which is solely
focused on solving Alfven standing modes along magnetic field lines, and a
realistic Saturn magnetic field model, KMAG.

## Installation

Download or clone the whole program directory to anywhere on your computer. Make
sure you have Fortran compilers and Python 3 installed, and configured in your
PATH.

Tested using GNU Fortran (Homebrew GCC 10.2.0) and Python 3.7.6 running on
Mac OS 10.15.7 and ZSH on iTerm2 (Build 3.3.12).

The directories assume UNIX style path seperators, thus, it might not work on
Windows. A newer release will address this. 

## Usage
Pass -h or --help argument in the terminal to get a help with the command-line
driven interface to the program.

E.g. 'python3 AlfvenWaver.py --help'

Optional arguments:
  -h, --help            show this help message and exit
  --calc                Compute the Alfven waves
  --sol                 Plot the eigenfunctions
  --field               Plot the field lines
  --eig                 Plot the eigenfrequencies
  --errors              Plot the error function
  --mov                 Generate a movie
  --save                Save the figures
  --v                   Verbose mode
  --dipole              Use a dipole field
  --uniform             Use a uniform field
  --plotRefEigs         Plot reference eigenfrequencies
  --conf CONFIGURATION  Field line configuration number
  --plotvar {Xi,E,b,bratio,vA,n}
                        Variable to plot along the field lines
  --plotCoord {th,s}    Field line coordinate to plot against: latitude (th)
                        or distance along the field (s)
  --component {toroidal,poloidal}
                        Toroidal or poloidal component for solution
  --modes MODES         Number of modes to solve
  --L [L [L ...]]       Trace a field line from an equatorial crossing
                        distance, L (in planet radii). (Multiple values
                        accepted)
  --th [TH [TH ...]]    Trace a field line from a colatitude of theta (in
                        degrees) staring in the northern hemisphere (multiple
                        values accepted)
  --phi [PHI [PHI ...]]
                        Trace a field line from the azimuth of phi (in
                        degrees) staring from the sunward direction (CCW).
                        (Multiple values accepted)
  --fast                Enable a faster calculation (less accurate)
  --nfl NFL             Number of field lines in latitude / radial range
  --plotx PLOTX PLOTX   X limits of the plot
  --ploty PLOTY PLOTY   Z limits of the plot
  --Lmax LMAX           Furthest field line for which to plot field line
                        parameters, measured by equatorial crossing distance,
                        L (in planet radii).

## Program Directory Tree

The program includes the KMAG magnetic field model for Saturn (Khurana et al.,
2016), put seperately in the directory called KMAG2012. Python routines that
handle the interactions with the KMAG model are contained in the KMAGhelper
directory.

The standing Alfven wave model is contained within the wavesolver package, and
the calls to the package are handled through the main program file 
AlfvenWaver.py. It is assumed that the KMAG2012 and the KMAGhelper files stay
in the same relative directories.

The program save figures in the Output directory. For example, an eigenfrequency
plot will be created in it as "figure_eigfreq.pdf".

PROGRAM TREE 
* AlfvenWaver.py
  "The main program, mostly handling calls to the wavesolver package"
* KMAG2012
    * KMAG2012.f
      "KMAG Saturn Magnetic Field Model"
    * README.md
* KMAGhelper
  * KMAGtracer.f
    "Field line tracer for the KMAG model, processed by KmagFunction.py"
  * KmagFunctions.py
    "Python functions to communicate with the KMAG model"
* Output
  * figure_eigfreq.pdf
  "An example figure produced after running the program"
* README.md
* wavesolver
  * configurations.py
    "Basic field line, PDE solver, and model parameters"
  * fieldline.py
    "Field line class and associated functions to store field line data"
  * helperFunctions.py
    "Miscellaneous general-use functions"
  * io.py
    "Input / Output functions"
  * linalg.py
    "Linear Algebra"
  * model.py
    "Plasma density models, PDE functions, and dipole field model"
  * plot.py
    "Plotting routines"
  * shoot.py
    "Shooting Method PDE solver"
  * sim.py
    "Simulation class for storing all relevant parameters for calculations"

## Example

In order to trace a single field line from an equatorial position of 10 RS
(in KSM): 

```bash
python3 AlfvenWaver.py --field --L 10
```

The --field flag tells the program to create a figure showing the field line in
the Output directory.

If we would like to calculate the eigenfrequencies of the Alfven standing modes
for the first four modes, we can type:

```bash
python3 AlfvenWaver.py --calc --fast --L 10 --modes 4 --sol
```

or, if we want to see the eigenfunctions, change the --eig flag to --sol:

```bash
python3 AlfvenWaver.py --calc --fast --L 10 --modes 4 --sol
```

In order to have more accurate caclculations, remove the --fast flag.

You can also map field lines using the invariant latitude. For example, to map
a field line from a colatitude of 17 degrees in the northern hemisphere, at
azimuthal angle of 0 degrees (i.e. looking sunward), type:

```bash
python3 AlfvenWaver.py --field --th 17 --phi 0
```

We can map a field line at 17 degrees colatitude both in the day-side and the
night-side meridian by simple extending the list of azimuthal angles, i.e.

```bash
python3 AlfvenWaver.py --field --th 17 --phi 0 180
```

## Instructions for Reproducing Figures in Rusaitis et al. (2020, in-review)

In order to reproduce the figures as given in Rusaitis et al. (2020), run the
following (bash) shell script from the terminal in the main code folder. 

Note: add a --fast flag to compute field lines and eigenfrequencies more 
quickly (sucrifices accuracy).

### Figure 1a.
Color plot of the plasma density (Bagenal & Delamere, 2011) in Saturn's magnetic
field model (Khurana et al., 2016).
```bash
python3 AlfvenWaver.py --field --L -24 24 --nfl 60 --plotx -20 20 --ploty -7 7 
--plotvar n
```

### Figure 1b.
Color plot of the Alfven velocity (Bagenal & Delamere, 2011) in Saturn's 
magnetic field model (Khurana et al., 2016).
```bash
python3 AlfvenWaver.py --field --L -24 24 --nfl 60 --plotx -20 20 --ploty -7 7 
--plotvar vA
```

### Figure 2a and Figure 2b.
The wave electric field, E, for the first four harmonics of a field line with
an equatorial crossing distance of 19 RS, each normalized to a maximum arbitrary
amplitude of 1. 
```bash
python3 AlfvenWaver.py --calc --sol --L 19 --modes 4 --fast
```

### Figure 2c.
Field displacement, Xi, normalized to a maximum arbitrary amplitude of 1, for
the first four harmonics, shown separately along four different field lines with
equatorial crossing distances of 6, 11, 14, and 19 RS.
```bash
python3 AlfvenWaver.py --calc --field --plotvar Xi  --modes 4 --fast 
--L 6 11 14 19 --plotx -2 28  --ploty -7 7
```

### Figure 3a.
Eigenfrequencies for the first 6 toroidal (solid) and poloidal (dashed) modes
for varying invariant latitudes of field lines in the day side of Saturn’s
magnetosphere. The two figures are combined in the paper.
```bash
python3 AlfvenWaver.py --calc --eig --fast --th 15 30 --nfl 20 --modes 6 
--component toroidal --phi 0

python3 AlfvenWaver.py --calc --eig --fast --th 15 30 --nfl 20 --modes 6 
--component poloidal --phi 0
```

### Figure 3b.
Eigenfrequencies for the first 6 toroidal (solid) and poloidal (dashed) modes
for varying invariant latitudes of field lines in the night side of Saturn’s
magnetosphere. The two figures are combined in the paper.
```bash
python3 AlfvenWaver.py --calc --eig --fast --th 16 30 --nfl 20 --modes 6 
--component toroidal --phi 180

python3 AlfvenWaver.py --calc --eig --fast --th 16 30 --nfl 20 --modes 6 
--component poloidal --phi 180
```

### Figure 4a and Figure 4b.
A comparison of field lines mapped from the same equatorial positions for a
dipole (a) and a realistic field model (b). 
```bash
python3 AlfvenWaver.py --field --L 4 20 --nfl 6 --plotx -2 21 --ploty -8 8 --dip
python3 AlfvenWaver.py --field --L 4 20 --nfl 6 --plotx -2 21 --ploty -8 8
```

### Figure 4c and Figure 4d.
The respective field line eigenfrequencies of the first four modes for the
dipole (c) and a realistic field model (d).
```bash
python3 AlfvenWaver.py --calc --fast --eig --L 4 20 --nfl 6 --modes 4 
--plotCoord=s --dip

python3 AlfvenWaver.py --calc --fast --eig --L 4 20 --nfl 6 --modes 4 
--plotCoord=s
```

## References:
Khurana, Krishan K. (2020, October 12). KMAG - Kronian Magnetic Field Model 
(Version 1.0). Zenodo. http://doi.org/10.5281/zenodo.4080294

Khurana, K. K., Arridge, C. S., Schwarzl, H., & Dougherty, M. K. (2006). A model
of Saturn’s magnetospheric field based on latest Cassini observations. In AGU
Spring Meeting Abstracts (Vol. 2007, pp. P44A-01).

## Version
AlfvenWaver v1.0 10-21-2020.

## License
[MIT](https://choosealicense.com/licenses/mit/)
Liutauras Rusaitis
