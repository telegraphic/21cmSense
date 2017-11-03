#! /usr/bin/env python
"""
# generate_uv_coverage.py

Creates an array file for use by calc_sense.py.  The main product is the uv coverage produced by the array during
the time it takes the sky to drift through the primary beam; other array parameters are also saved.  Array specific
information comes from an aipy cal file. If opts.track is set, produces the uv coverage for the length specified
instead of that set by the primary beam.
"""

import aipy as a, numpy as n
import optparse, sys
import hickle as hkl

from conversion_utils import *
from generate_array import *
from antenna import lwa_antenna


def generate_uv_coverage(aa, freq, bl_max=None, track=None):
    """ Generate UV coverage for a given antenna array.

    Parameters:
        aa (AntennaArray): Antenna array object, made with generate_array
        bl_max (None / float): Maximum baseline length (filter out long baselines that aren't sensitive)
        track (None / float): Length of observation, assuming that it is a tracked observation (default is drift scan).
    """
    #==============================READ ARRAY PARAMETERS=========================

    #load cal file
    n_ants = len(aa)
    prms = aa.get_arr_params()
    if track:
        obs_duration=60.*track
        name = prms['name'] + 'track_%.1fhr' % track
    else:
        obs_duration = 12.0 * 60  # A dipole-type antenna sees the entire sky for 12 hours
        name = prms['name'] + 'drift'

    lwa = lwa_antenna.LwaBeamPattern()
    lwa.generate(freq * 1e3)
    dish_size_in_lambda = lwa.estimate_dish_size_in_lambda()

    #==========================FIDUCIAL OBSERVATION PARAMETERS===================

    #while poor form to hard code these to arbitrary values, they have very little effect on the end result

    #observing time
    t_int = 60. #how many seconds a single visibility has integrated
    cen_jd = 2454600.90911
    start_jd = cen_jd - (1./24)*((obs_duration/t_int)/2)
    end_jd = cen_jd + (1./24)*(((obs_duration-1)/t_int)/2)
    times = n.arange(start_jd,end_jd,(1./24/t_int))
    print 'Observation duration (julian days): %.3f' % (end_jd - start_jd)

    #================================MAIN CODE===================================

    cnt = 0
    uvbins = {}

    cat = a.src.get_catalog('z')                        # create zenith source object
    aa.set_jultime(cen_jd)
    obs_lst = aa.sidereal_time()
    obs_zen = a.phs.RadioFixedBody(obs_lst, aa.lat)
    obs_zen.compute(aa)                                 # observation is phased to zenith of the
                                                        # center time of the drift

    #find redundant baselines
    bl_len_max = 0.
    for i in xrange(n_ants):
        LinePrint('working on antenna %i of %i' % (i, len(aa)))
        for j in xrange(n_ants):
            if i == j:
                continue #no autocorrelations
            u, v, w = aa.gen_uvw(i, j, src=obs_zen)
            bl_len = n.sqrt(u ** 2 + v ** 2)
            if bl_len > bl_len_max:
                bl_len_max = bl_len
            uvbin = '%.1f,%.1f' % (u, v)
            cnt +=1
            if not uvbins.has_key(uvbin):
                uvbins[uvbin] = ['%i,%i' % (i, j)]
            else:
                uvbins[uvbin].append('%i,%i' % (i, j))
    print '\nThere are %i baseline types' % len(uvbins.keys())

    print 'The longest baseline is %.2f meters' % bl_len_max
    if bl_max:
        bl_len_max = bl_max
    print 'The longest baseline being included is %.2f m' % bl_len_max

    # grid each baseline type into uv plane
    dim = int(n.round(bl_len_max / dish_size_in_lambda) + 1)
    uvsum, quadsum = n.zeros((dim, dim)), n.zeros((dim, dim))  # quadsum adds all non-instantaneously-redundant
                                                               # baselines incoherently
    for cnt, uvbin in enumerate(uvbins):
        LinePrint('working on %i of %i uvbins' % (cnt+1, len(uvbins)))
        uvplane = n.zeros((dim,dim))
        for t in times:
            aa.set_jultime(t)
            lst = aa.sidereal_time()
            obs_zen.compute(aa)
            bl = uvbins[uvbin][0]
            nbls = len(uvbins[uvbin])
            i, j = bl.split(',')
            i, j = int(i), int(j)
            u, v, w = aa.gen_uvw(i, j, src=obs_zen) / dish_size_in_lambda
            _beam = beamgridder(xcen=u, ycen=v, size=dim)
            uvplane += nbls * _beam
            uvsum += nbls * _beam
        quadsum += uvplane ** 2
    print("")

    quadsum **= .5

    uv_coverage = {
        'uv_coverage': uvsum,
        'uv_coverage_pess': quadsum,
        'name': name,
        'obs_duration': obs_duration,
        'Trx': prms['Trx'],
        't_int': t_int,
        'freq_ghz': freq,
        'dish_size_in_lambda': dish_size_in_lambda
    }

    return uv_coverage

def save_uv_coverage(uv_coverage, directory_out):
    """ Save pre-computed uv-coverage to file. """
    fout = os.path.join(directory_out, "%s_arrayfile.hkl" % uv_coverage['name'])
    hkl.dump(uv_coverage, fout)

def load_uv_coverage(filename, directory="uv_coverages"):
    uv = hkl.load(os.path.join(directory, filename))
    return uv
    
def plot_uv_coverage(uv_coverage, directory_out="figures"):
    """ Plot UV coverage """
    fout = os.path.join(directory_out, "%s_uv_coverage.pdf" % uv_coverage['name'])
    uvc = uv_coverage['uv_coverage']
    ds  = uv_coverage['dish_size_in_lambda']
    plt.imshow(uvc,
               interpolation='nearest',
               aspect='auto',
               extent=(-uvc.shape[0] / 2 * ds, uvc.shape[0] / 2 * ds, -uvc.shape[1] / 2 * ds, uvc.shape[1] / 2 * ds)
               )
    plt.xlabel("U ($\lambda$)")
    plt.ylabel("V ($\lambda$)")
    plt.colorbar()
    plt.savefig(fout)
    plt.show()

if __name__ == "__main__":
    o = optparse.OptionParser()
    o.set_usage('generate_uv_coverage.py [calfile]')
    #a.scripting.add_standard_options(o, cal=True)
    o.add_option('--track', dest='track', default=None, type=float,
        help="If set, calculate sensitivity for a tracked observation of this duration in hours; otherwise, "
             "calculate for a drift scan.")
    o.add_option('--bl_max', dest='bl_max', default=None, type=float,
        help="Set the maximum baseline (in meters) to include in the uv plane.  Use to exclude outriggers with "
             "little EoR sensitivity to speed up calculation.")
    opts, args = o.parse_args(sys.argv[1:])

    ant_pos = generate_hexagon(7, 10, 1.0)
    plot_antenna_array(ant_pos)
    save_antenna_array(ant_pos, "brawl127.txt")

    ovro = ('37.240391', '-118.281667', 1184)

    #Set other array parameters here
    params = {
        'name': 'BRAWL',
        'Trx': 500 * 1e3  # receiver temp in mK, T_sky is taken care of later
    }

    aa = generate_aa(ovro, ant_pos, params)
    uv = generate_uv_coverage(aa, bl_max=opts.bl_max, track=opts.track)
    save_uv_coverage(uv, "uv_coverages")
    plot_uv_coverage(uv)