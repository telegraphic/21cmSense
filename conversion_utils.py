"""
# Conversion_utils.py

Conversion utilities, and other miscellaneous methods & shared variables that are useful
across multiple scripts.
"""

import aipy as a
import numpy as np
import sys

#=========================COSMOLOGY/BINNING FUNCTIONS=========================

class LinePrint():
    """
    Print things to stdout on one line dynamically
    """
    def __init__(self,data):
        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()


#Convert frequency (GHz) to redshift for 21cm line.
def f2z(fq):
    F21 = 1.42040575177
    return (F21 / fq - 1)

#Convert frequency (GHz) to redshift for 21cm line.
def z2f(z):
    F21 = 1.42040575177
    return (F21 / (z + 1))

#Multiply by this to convert an angle on the sky to a transverse distance in Mpc/h at redshift z
def dL_dth(z):
    """[h^-1 Mpc]/radian, from Furlanetto et al. (2006)"""
    return 1.9 * (1./a.const.arcmin) * ((1+z) / 10.)**.2

#Multiply by this to convert a bandwidth in GHz to a line of sight distance in Mpc/h at redshift z
def dL_df(z, omega_m=0.266):
    """[h^-1 Mpc]/GHz, from Furlanetto et al. (2006)"""
    return (1.7 / 0.1) * ((1+z) / 10.)**.5 * (omega_m/0.15)**-0.5 * 1e3

#Multiply by this to convert a baseline length in wavelengths (at the frequency corresponding to redshift z) into a tranverse k mode in h/Mpc at redshift z
def dk_du(z):
    """2pi * [h Mpc^-1] / [wavelengths], valid for u >> 1."""
    return 2*np.pi / dL_dth(z) # from du = 1/dth, which derives from du = d(sin(th)) using the small-angle approx

#Multiply by this to convert eta (FT of freq.; in 1/GHz) to line of sight k mode in h/Mpc at redshift z
def dk_deta(z):
    """2pi * [h Mpc^-1] / [GHz^-1]"""
    return 2*np.pi / dL_df(z)

#scalar conversion between observing and cosmological coordinates
def X2Y(z):
    """[h^-3 Mpc^3] / [str * GHz]"""
    return dL_dth(z)**2 * dL_df(z)

#A function used for binning
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def fname_to_z(z_filename):
    fparts = z_filename.split("_")
    redshift =  float(fparts[1][1:])
    return redshift

def in_range(z, lower, upper):
    return z >= lower and z <= upper

#============================SIMPLE GRIDDING FUNCTION=======================

def beamgridder(xcen,ycen,size):
    crds = np.mgrid[0:size,0:size]
    cen = size/2 - 0.5 # correction for centering
    xcen += cen
    ycen = -1*ycen + cen
    beam = np.zeros((size,size))
    if round(ycen) > size - 1 or round(xcen) > size - 1 or ycen < 0. or xcen <0.:
        return beam
    else:
        beam[int(round(ycen)), int(round(xcen))] = 1. #single pixel gridder
        return beam

def closest(arr, val):
     return np.argmin(np.abs(arr - val))

def T_sky(f):
    """ Compute sky temperature in mK for freq in GHz. """
    Tsky = 60e3 * (3e8 / (f * 1e9)) ** 2.55
    return Tsky

def B_cosmo(f, f0=.150, B0=0.008):
    """ Compute cosmological bandwidth for given frequency.

    Keeps fractional bandwidth fixed, extrapolated from a base frequency & bandwidth.

    Parameters:
        f:  freq in GHz
        f0: base frequency, defaults to 150 MHz.
        B0: base bandwidth, defaults to 0.008
    """
    B = f / f0 * B0
    return B


# Gridding scale factors
WL_100MHZ = 2.99            # Wavelength at 100 MHz
N_SCALE   = 0.5             # Scale factor for bin size