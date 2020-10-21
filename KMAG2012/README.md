# KMAG

## Introduction

KMAG is a global magnetic field of the Saturn’s field (Khurana et al., 2006) that includes modules that specify (1) the internal spherical harmonic field, (2) the ring current and the magnetotail current system, (3) the field from the radial current system which reinforces corotation on the outflowing plasma, (4) the shielding fields from an axially symmetric but hinged magnetopause and (5) the interconnection magnetic field from the solar wind IMF. The current sheet field is based on models of disk-shaped current sheets by Tsyganenko and Peredo (1994). The tilt, and hinging of the current sheet is based on the general deformation technique of Tsyganenko (1998, 2002a, 2002b). The shielding fields are derived using a model of the magnetopause by Arridge et al. (2006) constructed from magnetopause crossings observed by both Cassini and Voyager. Saturn’s internal field uses magnetic moments derived from the Cassini Saturn Orbital Insertion (SOI) measurements (Dougherty, 2005). KMAG is based on a previous Euler-potential model of Jupter's magnetospheric field (Khurana, 1997).

## Installation

The model is fully self-contained in a single Fortran file. Call the relevant
subroutines from your custom Fortran code.

Tested using GNU Fortran (Homebrew GCC 10.2.0) 10.2.0.

## Usage
The main subroutine that computes the magnetic field strength is
called KMAG(). It requires the following parameters to be passed through:

### INPUT:
 * TIME (REAL) - a double precision variable denoting number of seconds from 
   an epoch
 * EPOCH (STRING) - a five letter character which can be either 'ctime',
   'etime' or 'J2000'/upper or lower case
 * R, THETA, PHI (REAL) -  position of the field in a right-handed System III
   coordinate system (THETA and PHI are in radians)
 * BY_IMF, BZ_IMF (REAL) - Y and Z component of the interplanetary magnetic
   field in right-handed KSM coordinate system.
### OUTPUT:
 * RLT (REAL) - local time
 * BRM, BTM, BPM (REAL) - model magnetic field in a spherical SYSTEM III 
   coordinate system.

The supported coordinate systems are:
* S3C System III Cartesian (right-handed)
* KSO Kronian-Sun-Orbital
* KSM Kronian-Sun-Magnetic (So far same as KSO)
* DIP Dipole (cartesian)

### Example

To call the main subroutine from your custom Fortran code, you can simply 
write:

```fortran
CALL KMAG(TIME,EPOCH,R,THETA,PHI,RLT,BY_IMF,BZ_IMF,BRM,BTM,BPM,Dp)
```

In order to compile and run the model from within Python, you can use the
'os' module. For example, if your custom Fortran code is called
'KMAG_pipe.f', you can pass the following commands to compile and run the
code, assuming 'KMAG2012.f' is in the same directory, and you have a
gfortran compiler installed:

```python
import os
# -w to inhibit Fortran2018 warnings
os.system("gfortran -c KMAG_pipe.f -w")
os.system("gfortran -c KMAG2012.f -w")
os.system("gfortran -o KMAG KMAG_pipe.o KMAG2012.o")
os.system("./KMAG")
```

## References:
Arridge, C. S., Achilleos, N., Dougherty, M. K., Khurana, K. K., & Russell, C. T. (2006). Modeling the size and shape of Saturn’s magnetopause with variable dynamic pressure. Journal of Geophysical Research, 111(A11), A11227. https://doi.org/10.1029/2005JA011574

Dougherty, M. K. (2005). Cassini magnetometer observations during Saturn orbit insertion. Science, 307(5713), 1266–1270. https://doi.org/10.1126/science.1106098

Khurana, K. K., Arridge, C. S., Schwarzl, H., & Dougherty, M. K. (2006). A model of Saturn’s magnetospheric field based on latest Cassini observations. In AGU Spring Meeting Abstracts (Vol. 2007, pp. P44A-01).

Khurana, Krishan K. (1997). Euler potential models of Jupiter’s magnetospheric field. Journal of Geophysical Research: Space Physics, 102(A6), 11295–11306. https://doi.org/10.1029/97JA00563

Tsyganenko, N. A. (1998). Modeling of twisted/warped magnetospheric configurations using the general deformation method. Journal of Geophysical Research: Space Physics, 103(A10), 23551–23563. https://doi.org/10.1029/98JA02292

Tsyganenko, N. A. (2002a). A model of the near magnetosphere with a dawn-dusk asymmetry 1. Mathematical structure. Journal of Geophysical Research: Space Physics, 107(A8), SMP 12-1-SMP 12-15. https://doi.org/10.1029/2001JA000219

Tsyganenko, N. A. (2002b). A model of the near magnetosphere with a dawn-dusk asymmetry 2. Parameterization and fitting to observations. Journal of Geophysical Research: Space Physics, 107(A8), SMP 10-1-SMP 10-17. https://doi.org/10.1029/2001JA000220

Tsyganenko, N. A., & Peredo, M. (1994). Analytical models of the magnetic field of disk-shaped current sheets. Journal of Geophysical Research, 99(A1), 199. https://doi.org/10.1029/93JA02768


## Version
KMAG v1.0 10-13-2020.

## License
[MIT](https://choosealicense.com/licenses/mit/)
