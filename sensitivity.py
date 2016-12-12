#! /usr/bin/env python
"""
# sensitivity.py

Calculates the expected sensitivity of a 21cm experiment to a given 21cm power spectrum.
Requires as input an array .npz file created with generate_uv_coverage.py.
"""

import numpy as n, optparse, sys
from scipy import interpolate
import os
from conversion_utils import *
import hickle as hkl
from antenna import lwa_antenna

def calc_sense(arr_file, ps_model, opts, print_report=True):
    """ Calculate sensitivity of an antenna array.

    Args:
        arr_file (str): filename of input array, as made with generate_uv_coverage.py. Should
            end with .npz.
        opts (dict): A dictionary of observation options / parameters. Keys are:
            'model' (str): The model of the foreground wedge to use.  Three options are 'pess' (all k modes inside
                horizon + buffer are excluded, and all baselines are added incoherently), 'mod' (all k modes inside
                horizon + buffer are excluded, but all baselines within a uv pixel are added coherently), and 'opt'
                (all modes k modes inside the primary field of view are excluded).
                See Pober et al. 2014 for more details.
            'buff' (float): The size of the additive buffer outside the horizon to exclude in the pessimistic and
                moderate models.
            'eor' (str): Filename for the model epoch of reionization power spectrum.  The code is built to handle
                output power spectra from 21cmFAST.
            'ndays' (int): The total number of days observed. The default is 180, which is the maximum a particular R.A.
                can be observed in one year if one only observes at night.
                The total observing time is ndays*n_per_day.
            'n_per_day' (float): The number of good observing hours per day.  This corresponds to the size of a
                low-foreground region in right ascension for a drift scanning instrument.  The total observing time is
                ndays*n_per_day.  Default is 6.  If simulating a tracked scan, n_per_day should be a multiple of the
                length of the track (i.e. for two three-hour tracks per day, n_per_day should be 6).
            'bwidth' (float): Cosmological bandwidth in GHz.  Note this is not the total instrument bandwidth, but the
                redshift range that can be considered co-eval.  Default is 0.008 (8 MHz).
            'nchan' (int): Integer number of channels across cosmological bandwidth.  Defaults to 82, which is
                equivalent to 1024 channels over 100 MHz of bandwidth.  Sets maximum k_parallel that can be probed,
                but little to no overall effect on sensitivity.
            'no_ns' (bool): Remove pure north/south baselines (u=0) from the sensitivity calculation.
                These baselines can potentially have higher systematics, so excluding them represents a
                conservative choice.
            'z' (float): Redshift of input power spectra model. Observation frequency is computed from this value.
        print_report (bool): Flag to suppress (False) or show (True) printing of sensitivity info.

    Returns:
        outdict (dict): A dictionary containing computed sensitivity and other information. Keys are:
            'snr' (float): computed signal-to-noise ratio,
            'f_mhz' (float): frequency of observation in MHz,
            'z' (float): redshift of observation corresponding to f_mhz,
            'T_sky' (float): Sky temperature in K, at f_mhz
    """

    # ====================OBSERVATION/COSMOLOGY PARAMETER VALUES====================

    # Load in data from array file; see generate_uv_coverage.py for definitions of the parameters
    array = hkl.load(arr_file)
    name = array['name']
    obs_duration = array['obs_duration']
    Trx = array['Trx']
    t_int = array['t_int']
    if opts['model'] == 'pess':
        uv_coverage = array['uv_coverage_pess']
    else:
        uv_coverage = array['uv_coverage']

    h = 0.7
    B = opts['bwidth']
    z = opts['z']
    freq = z2f(z)                           # Frequency is now computed from input redshift
    wl_meters = 2.99e8 / (freq * 1e9)
    if opts['scale_bandwidth']:
        B = B_cosmo(freq, f0=.100, B0=B)        # Scale bandwidth with frequency
    else:
        B = B
    #dish_size_in_lambda = dish_size_in_lambda*(freq/.150)
    # Dish size is no longer frequency scaled!

    nchan = opts['nchan']
    kpls = dk_deta(z) * n.fft.fftfreq(nchan, B / nchan)

    Tsky = 60e3 * (3e8 / (freq * 1e9)) ** 2.55  # sky temperature in mK
    n_lstbins = opts['n_per_day'] * 60. / obs_duration

    # Precalculate bm and Tsys
    Tsys = Tsky + Trx
    lwa = lwa_antenna.LwaBeamPattern()
    lwa.generate(freq * 1e3)
    bm_eff = lwa.compute_bm_eff()
    dish_size_in_lambda = lwa.estimate_dish_size_in_lambda()

    if print_report:
        print "Redshift of model:        %2.2f" % z
        print "Frequency:                %2.2f MHz" % (freq * 1e3)
        print "Wavelength:               %2.2f m"   % wl_meters
        print "Sky temperature:          %2.1f K" % (Tsky / 1e3)
        print "Cosmological bandwidth:   %2.2f MHz" % (B * 1e3)
        print "Min, max k_parallel:      %2.4f, %2.4f" % (np.min(np.abs(kpls)), np.max(kpls))
        print "Beam eff solid angle:     %2.4fsr" % bm_eff
        print "Eff. dish size in lambda: %2.2f Lambda" % dish_size_in_lambda
        print "UV coverage max:          %2.2f" % np.max(uv_coverage)

    #print "BM: %2.2f    Tsys: %2.2fK" % (bm_eff, Tsys / 1e3)


    #===============================EOR MODEL===================================

    #You can change this to have any model you want, as long as mk, mpk and p21 are returned

    #This is a dimensionless power spectrum, i.e., Delta^2
    #modelfile = opts['eor']
    #model = n.loadtxt(modelfile)
    mk, mpk = ps_model[:, 0] / h, ps_model[:, 1]  #k, Delta^2(k)
    #note that we're converting from Mpc to h/Mpc
    #interpolation function for the EoR model

    p21 = interpolate.interp1d(mk, mpk, kind='linear')

    #=================================MAIN CODE===================================

    #set up blank arrays/dictionaries
    kprs = []
    #sense will include sample variance, Tsense will be Thermal only
    sense, Tsense = {}, {}

    uv_coverage *= t_int
    SIZE = uv_coverage.shape[0]

    # Cut unnecessary data out of uv coverage: auto-correlations & half of uv plane
    # (which is not statistically independent for real sky)
    uv_coverage[SIZE / 2, SIZE / 2] = 0.
    uv_coverage[:, :SIZE / 2] = 0.
    uv_coverage[SIZE / 2:, SIZE / 2] = 0.
    if opts['no_ns']: uv_coverage[:, SIZE / 2] = 0.

    #loop over uv_coverage to calculate k_pr
    nonzero = n.where(uv_coverage > 0)

    for iu, iv in zip(nonzero[1], nonzero[0]):

        # convert u, v back into m, and then back into wavelength
        u, v = (iu - SIZE / 2) * dish_size_in_lambda, (iv - SIZE / 2) * dish_size_in_lambda
        umag = n.sqrt(u ** 2 + v ** 2)
        kpr = umag * dk_du(z)
        kprs.append(kpr)
        #calculate horizon limit for baseline of length umag
        if opts['model'] in ['mod', 'pess']:
            hor = dk_deta(z) * umag / freq + opts['buff']
        elif opts['model'] in ['opt']:
            hor = dk_deta(z) * (umag / freq) #* n.sin(first_null / 2)
        elif opts['model'] in ['trott']:
            hor = 0.1
        elif opts['model'] in ['none']:
            hor = 0.00
        else:
            print '%s is not a valid foreground model; Aborting...' % opts['model'];
            sys.exit()
        if not sense.has_key(kpr):
            sense[kpr] = n.zeros_like(kpls)
            Tsense[kpr] = n.zeros_like(kpls)
        for i, kpl in enumerate(kpls):
            #exclude k_parallel modes contaminated by foregrounds
            if n.abs(kpl) < hor:
                continue
            k = n.sqrt(kpl ** 2 + kpr ** 2)
            if k < n.min(mk):
                continue
            #don't include values beyond the interpolation range (no sensitivity anyway)
            if k > n.max(mk):
                continue

            delta21 = p21(k)

            scalar = X2Y(z) * bm_eff * B * k ** 3 / (2 * n.pi ** 2)
            tot_integration = uv_coverage[iv, iu] * opts['ndays']
            Trms = Tsys / n.sqrt(2 * (B * 1e9) * tot_integration)
            #add errors in inverse quadrature
            sense[kpr][i] += 1. / (scalar * Trms ** 2 + delta21) ** 2
            Tsense[kpr][i] += 1. / (scalar * Trms ** 2) ** 2

    #bin the result in 1D
    delta = dk_deta(z) * (1. / B)  #default bin size is given by bandwidth
    #kmag = n.arange(delta, n.max(mk), delta)
    kmag = n.arange(n.min(mk), n.max(mk), delta)

    #print "delta: %s max: %s (min: %s)" % (delta, n.max(mk), n.min(mk))

    kprs = n.array(kprs)
    sense1d = n.zeros_like(kmag)
    Tsense1d = n.zeros_like(kmag)
    for ind, kpr in enumerate(sense.keys()):
        #errors were added in inverse quadrature, now need to invert and take square root
        # to have error bars; also divide errors by number of indep. fields
        sense[kpr] = sense[kpr] ** -.5 / n.sqrt(n_lstbins)
        Tsense[kpr] = Tsense[kpr] ** -.5 / n.sqrt(n_lstbins)
        for i, kpl in enumerate(kpls):
            k = n.sqrt(kpl ** 2 + kpr ** 2)
            if k > n.max(mk): continue
            #add errors in inverse quadrature for further binning
            sense1d[find_nearest(kmag, k)] += 1. / sense[kpr][i] ** 2
            Tsense1d[find_nearest(kmag, k)] += 1. / Tsense[kpr][i] ** 2

    #invert errors and take square root again for final answer
    for ind, kbin in enumerate(sense1d):
        sense1d[ind] = kbin ** -.5
        Tsense1d[ind] = Tsense1d[ind] ** -.5

    #calculate significance with least-squares fit of amplitude
    A = p21(kmag)
    M = p21(kmag)

    sense_dict = {
        'name': name,
        'model': opts['model'],
        'freq': freq,
        'ks': kmag,
        'errs': sense1d,
        'T_errs': Tsense1d,
        'p21': A
    }

    ps_model = os.path.basename(opts['ps_model'])
    hkl.dump(sense_dict, 'sensitivities/%s_%s_%.3f_%s.hkl' % (name, opts['model'], freq, ps_model))



    wA, wM = A * (1. / sense1d), M * (1. / sense1d)
    wA, wM = n.matrix(wA).T, n.matrix(wM).T
    amp = (wA.T * wA).I * (wA.T * wM)
    #errorbars
    Y = n.float(amp) * wA
    dY = wM - Y
    s2 = (len(wM) - 1) ** -1 * (dY.T * dY)
    X = n.matrix(wA).T * n.matrix(wA)
    err = n.sqrt((1. / n.float(X)))

    snr = float(amp / err)
    #print 'total snr = %s' % snr

    outdict = {'z': z,
               'snr': snr,
               'f_mhz': freq,
               'T_sky': Tsky
    }

    if print_report:
        print "SNR:               %.4f" % snr

    return outdict


if __name__ == "__main__":
    o = optparse.OptionParser()
    o.set_usage('calc_sense.py [options] *.npz')
    o.set_description(__doc__)
    o.add_option('-m', '--model', dest='model', default='mod',
                 help="The model of the foreground wedge to use.  Three options are 'pess' (all k modes inside horizon "
                      "+ buffer are excluded, and all baselines are added incoherently), 'mod' (all k modes inside "
                      "horizon + buffer are excluded, but all baselines within a uv pixel are added coherently), and "
                      "'opt' (all k modes inside the primary field of view are excluded).  "
                      "See Pober et al. 2014 for more details.")
    o.add_option('-b', '--buff', dest='buff', default=0.1, type=float,
                 help="The size of the additive buffer outside the horizon to exclude in the pessimistic "
                      "and moderate models.")
    o.add_option('--eor', dest='eor', default='ps_no_halos_nf0.521457_z9.50_useTs0_zetaX-1.0e+00_200_400Mpc_v2',
                 help="The model epoch of reionization power spectrum.  The code is built to handle output power "
                      "spectra from 21cmFAST.")
    o.add_option('--ndays', dest='ndays', default=180., type=float,
                 help="The total number of days observed.  The default is 180, which is the maximum a particular R.A. "
                      "can be observed in one year if one only observes at night.  The total observing time is n"
                      "days*n_per_day.")
    o.add_option('--n_per_day', dest='n_per_day', default=6., type=float,
                 help="The number of good observing hours per day.  This corresponds to the size of a low-foreground "
                      "region in right ascension for a drift scanning instrument.  "
                      "The total observing time is ndays*n_per_day.  Default is 6.  If simulating a tracked scan,"
                      " n_per_day should be a multiple of the length of the track (i.e. for two three-hour tracks "
                      "per day, n_per_day should be 6).")
    o.add_option('--bwidth', dest='bwidth', default=0.008, type=float,
                 help="Cosmological bandwidth in GHz.  Note this is not the total instrument bandwidth, but the "
                      "redshift range that can be considered co-eval.  Default is 0.008 (8 MHz).")
    o.add_option('--nchan', dest='nchan', default=82, type=int,
                 help="Integer number of channels across cosmological bandwidth.  Defaults to 82, which is equivalent "
                      "to 1024 channels over 100 MHz of bandwidth.  Sets maximum k_parallel that can be probed, but "
                      "little to no overall effect on sensitivity.")
    o.add_option('--no_ns', dest='no_ns', action='store_true',
                 help="Remove pure north/south baselines (u=0) from the sensitivity calculation.  "
                      "These baselines can potentially have higher systematics, so excluding them represents a "
                      "conservative choice.")
    o.add_option("--redshift", dest="z", default=20, type=float,
                 help="Redshift of power spectra model. An observational frequency is computed from this value.")

    opts, args = o.parse_args(sys.argv[1:])
    arr_file = args[0]

    calc_sense(arr_file, opts)

