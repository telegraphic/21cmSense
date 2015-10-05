import aipy as a
import numpy as np

#=========================COSMOLOGY/BINNING FUNCTIONS=========================

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
    '''[h^-1 Mpc]/radian, from Furlanetto et al. (2006)'''
    return 1.9 * (1./a.const.arcmin) * ((1+z) / 10.)**.2

#Multiply by this to convert a bandwidth in GHz to a line of sight distance in Mpc/h at redshift z
def dL_df(z, omega_m=0.266):
    '''[h^-1 Mpc]/GHz, from Furlanetto et al. (2006)'''
    return (1.7 / 0.1) * ((1+z) / 10.)**.5 * (omega_m/0.15)**-0.5 * 1e3

#Multiply by this to convert a baseline length in wavelengths (at the frequency corresponding to redshift z) into a tranverse k mode in h/Mpc at redshift z
def dk_du(z):
    '''2pi * [h Mpc^-1] / [wavelengths], valid for u >> 1.'''
    return 2*np.pi / dL_dth(z) # from du = 1/dth, which derives from du = d(sin(th)) using the small-angle approx

#Multiply by this to convert eta (FT of freq.; in 1/GHz) to line of sight k mode in h/Mpc at redshift z
def dk_deta(z):
    '''2pi * [h Mpc^-1] / [GHz^-1]'''
    return 2*np.pi / dL_df(z)

#scalar conversion between observing and cosmological coordinates
def X2Y(z):
    '''[h^-3 Mpc^3] / [str * GHz]'''
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
        beam[round(ycen),round(xcen)] = 1. #single pixel gridder
        return beam

def closest(arr, val):
     return np.argmin(np.abs(arr - val))
