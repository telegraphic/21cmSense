#! /usr/bin/env python
"""
Calculates the expected 2D thermal sensitivity of a 21cm experiment.  Since it is thermal sensitivity only,
no power spectrum model is required.  Requires as input an array .npz file created with generate_uv_coverage.py.
"""

import numpy as n, optparse, sys
import pylab as p
from conversion_utils import *

import hickle as hkl

def calc_sense2d(opts, args):

    #Load in data from array file; see generate_uv_coverage.py for definitions of the parameters
    array = hkl.load(args[0])
    name = array['name']
    obs_duration = array['obs_duration']
    dish_size_in_lambda = array['dish_size_in_lambda']
    Trx = array['Trx']
    t_int = array['t_int']
    uv_coverage = array['uv_coverage']

    h = 0.7
    B = .008        # largest bandwidth allowed by "cosmological evolution", i.e., the maximum line of sight
                    # volume over which the universe can be considered co-eval
    z = f2z(opts.freq)

    #dish_size_in_lambda = dish_size_in_lambda*(opts.freq/.150) # linear frequency evolution, relative to 150 MHz
    first_null = 1.22/dish_size_in_lambda           # for an airy disk, even though beam model is Gaussian
    bm = 1.13*(2.35*(0.45/dish_size_in_lambda))**2
    nchan = opts.nchan
    kpls = dk_deta(z) * n.fft.fftfreq(nchan, B/nchan)
    print kpls

    Tsky = 60e3 * (3e8/(opts.freq*1e9))**2.55  # sky temperature in mK
    n_lstbins = opts.n_per_day*60./obs_duration

    #=================================MAIN CODE===================================

    #set up blank arrays/dictionaries
    Tsense = {}

    p.imshow(uv_coverage, interpolation='nearest')
    p.colorbar()
    p.show()

    uv_coverage *= t_int
    SIZE = uv_coverage.shape[0]

    # Cut unnecessary data out of uv coverage: auto-correlations & half of uv plane
    # (which is not statistically independent for real sky)
    uv_coverage[SIZE/2,SIZE/2] = 0.
    uv_coverage[:,:SIZE/2] = 0.
    uv_coverage[SIZE/2:,SIZE/2] = 0.
    if opts.no_ns: uv_coverage[:,SIZE/2] = 0.

    #loop over uv_coverage to calculate k_pr
    nonzero = n.where(uv_coverage > 0)
    for iu,iv in zip(nonzero[1], nonzero[0]):
       u, v = (iu - SIZE/2) * dish_size_in_lambda, (iv - SIZE/2) * dish_size_in_lambda
       umag = n.sqrt(u**2 + v**2)
       kpr = umag * dk_du(z)
       if not Tsense.has_key(kpr):
           Tsense[kpr] = n.zeros_like(kpls)
       for i, kpl in enumerate(kpls):
           k = n.sqrt(kpl**2 + kpr**2)
           tot_integration = uv_coverage[iv,iu] * opts.ndays
           Tsys = Tsky + Trx
           bm2 = bm/2. # beam^2 term calculated for Gaussian; see Parsons et al. 2014
           bm_eff = bm**2 / bm2 # this can obviously be reduced; it isn't for clarity
           scalar = X2Y(z) * bm_eff * B  * k**3 / (2*n.pi**2)
           Trms = Tsys / n.sqrt(2*(B*1e9)*tot_integration)
           # add errors in inverse quadrature
           Tsense[kpr][i] += 1./(scalar*Trms**2)**2

    print Tsense

    # bin in annuli of k_perp
    deltak_perp = dk_du(z) * dish_size_in_lambda #bin on beam scale
    kprs = n.arange(0,0.3,deltak_perp) #binning range
    Tsense2d = n.zeros((len(kprs),len(kpls)))
    cen = len(kpls)/2
    for kpr in Tsense.keys():
        ind = find_nearest(kprs,kpr)
        Tsense2d[ind] += n.append(Tsense[kpr][cen:], Tsense[kpr][:cen])

    max_idx = 0
    kpls = n.append(kpls[cen:],kpls[:cen])
    for ind, kpr in enumerate(kprs):
        if np.average(Tsense2d[ind]) >= 1e-20:
            max_idx += 1
            print "%2.3f %2.3e" % (kpr, np.average(Tsense2d[ind]))
            Tsense2d[ind] = Tsense2d[ind]**-.5 / n.sqrt(n_lstbins)
        else:
            continue

    Tsense2d = Tsense2d[:max_idx]
    kprs     = kprs[:max_idx]


    #fold over k-parallel
    Tsense2d = Tsense2d[:,cen:][:,::-1]
    Tsense2d[:,:-1] /= n.sqrt(2)
    kpls = kpls[cen:]

    #save results to output hkl
    output = {
        'kprs' : kprs,
        'kpls': kpls,
        'T_errs': Tsense2d
    }
    hkl.dump(output, 'sensitivities/%s_%.3f_2dsense.hkl' % (name,opts.freq))

    if opts.plot_lin:
        p.imshow(Tsense2d.T, aspect='auto', interpolation='nearest',
                 extent=[kprs[0],kprs[-1],kpls[0],kpls[-1]])#,vmin=3,vmax=4)
    else:
        p.imshow(10 * n.log10(Tsense2d.T), aspect='auto', interpolation='nearest',
                 extent=[kprs[0],kprs[-1],kpls[0],kpls[-1]])#,vmin=3,vmax=4)
    p.colorbar()

    #plot horizon
    #umag = n.arange(0.,200.)
    #hor = C.pspec.dk_deta(z) * umag/opts.freq
    #kprs_hor = umag*C.pspec.dk_du(z)

    #p.plot(kprs_hor,hor,color='k')

    #p.ylim(0.,1.0)
    #p.xlim(0.,0.2)

    p.ylabel(r'$k_{||}$')
    p.xlabel(r'$k_{\perp}$')

    p.show()

if __name__ == "__main__":

    o = optparse.OptionParser()
    o.set_usage('calc_sense.py [options] *.npz')
    o.set_description(__doc__)
    o.add_option('-f', '--freq', dest='freq', default=.135, type=float,
        help="The center frequency of the observation in GHz.  If you change from the default, be sure to use a "
             "sensible power spectrum model from that redshift.  Note that many values in the code are "
             "calculated relative to .150 GHz and are not affected by changing this value.")
    o.add_option('--ndays', dest='ndays', default=180., type=float,
        help="The total number of days observed.  The default is 180, which is the maximum a particular R.A. can "
             "be observed in one year if one only observes at night.  The total observing time is ndays*n_per_day.")
    o.add_option('--n_per_day', dest='n_per_day', default=6., type=float,
        help="The number of good observing hours per day.  This corresponds to the size of a low-foreground r"
             "egion in right ascension for a drift scanning instrument.  The total observing time is ndays*n_per_day.  "
             "Default is 6.  If simulating a tracked scan, n_per_day should be a multiple of the length of the "
             "track (i.e. for two three-hour tracks per day, n_per_day should be 6).")
    o.add_option('--nchan', dest='nchan', default=82, type=int,
        help="Integer number of channels across cosmological bandwidth of 8 MHz.  Defaults to 82, which is equivalent"
             " to 1024 channels over 100 MHz of bandwidth.  Sets maximum k_parallel that can be probed, but "
             "little to no overall effect on sensitivity.")
    o.add_option('--no_ns', dest='no_ns', action='store_true',
        help="Remove pure north/south baselines (u=0) from the sensitivity calculation.  These baselines can "
             "potentially have higher systematics, so excluding them represents a conservative choice.")
    o.add_option('--linplot', dest='plot_lin', action='store_true',
        help="Plot in linear scale. Defaults to False (dB scale).")
    opts, args = o.parse_args(sys.argv[1:])

    calc_sense2d(opts, args)

