"""
# run_batch.py

Master script used for running batch simulations for computing sensitivity at different redshifts.

## Workflow
  1) Edit generate_models.py to create power spectra models in correct format for 21cmsense
  2) Edit leda127.py to set the parameters of the antenna array to be used
  3) Edit this script to set observational parameters (obs_opts), and set runtime parameters
  4) Run this master script.

## Directory structure
   anastasia/     input power spectra, as provided by Anastasia
   models/        generated power spectra in 21cmfast format
   output/        computed outputs (not used currently)
   tex/           latex files for memos and notes
   figures/       output directory for generated figures
"""

import numpy as np
import pylab as plt
import os

from sensitivity import calc_sense
from conversion_utils import *

# Runtime parameters
arr_module = 'leda127'                              # Name of antenna array python module file
arr_file = 'leda127drift_arrayfile.npz'             # Name of generated array file from mk_array_file.py
redshift_range = (15, 25)                           # Set redshift range of interest (ignore other PS model files)
regenerate_models = False                           # Regenerate PS models via generate_models.py
remake_array = True                                 # Regenerate array file via mk_array_file.py

if regenerate_models:
    os.system("generate_models.py")

if remake_array:
    os.system("python mk_array_file.py -C %s" % arr_module)


zs, snrs_mod, snrs_opt, snrs_pess = [], [], [], []
for modelfile in os.listdir('models'):
    if modelfile.endswith('.txt') and in_range(fname_to_z(modelfile), *redshift_range):
        print modelfile
        redshift = fname_to_z(modelfile)

        freq = z2f(redshift)

        obs_opts = {'model': 'mod',
                'buff': 0.1,
                'eor': 'models/%s' % modelfile,
                'ndays': 180,
                'n_per_day': 6,
                'bwidth': 0.008,
                'nchan': 82,
                'no_ns': False,
                'z': redshift}

        zs.append(redshift)

        obs_opts['model'] = 'pess'
        sense_dict = calc_sense(arr_file, obs_opts, print_report=False)
        snrs_pess.append(sense_dict['snr'])

        obs_opts['model'] = 'mod'
        sense_dict = calc_sense(arr_file, obs_opts, print_report=False)
        snrs_mod.append(sense_dict['snr'])

        obs_opts['model'] = 'opt'
        sense_dict = calc_sense(arr_file, obs_opts, print_report=True)
        snrs_opt.append(sense_dict['snr'])

        print("")
        #time.sleep(1)

# Output data as a table
snr_table = np.column_stack((zs, snrs_pess, snrs_mod, snrs_opt))
np.savetxt('output/snr-detection.txt', snr_table, header='# Redshift SNR_PESS SNR_MOD SNR_OPT')
np.savetxt('tex/snr-detection.txt', snr_table[:82:4], delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')

# Create pretty plots of SNR vs redshift
plt.figure("SNR")
plt.plot(zs, snrs_pess, c='#00cc00', label='Conservative')
plt.plot(zs, snrs_opt, c='#0000cc', label="Optimistic")
plt.xlabel("Redshift (z)")
plt.ylabel("Detection SNR")
plt.axvline(x=f2z(0.04), c='#666666', ls='dashed')
plt.axvline(x=f2z(0.0875), c='#666666', ls='dashed')
plt.axvline(x=f2z(0.108), c='#666666', ls='dashed')
plt.axvspan(f2z(0.0875), f2z(0.108), color='#333333', alpha=0.5, lw=0)
plt.xlim(10, 30)
plt.minorticks_on()
plt.legend(frameon=False)
plt.savefig('figures/detection-snr.pdf')

# Same plot, in frequency space instead of redshift
plt.figure("SNR FREQUENCY")
fs = z2f(np.array(zs)) * 1e3
plt.plot(fs, snrs_pess, c='#00cc00', label='Conservative')
plt.plot(fs, snrs_opt, c='#0000cc', label="Optimistic")
plt.xlabel("Frequency (MHz)")
plt.ylabel("Detection SNR")
plt.axvspan(87.5, 108, color='#333333', alpha=0.5, lw=0)
plt.minorticks_on()
plt.legend(frameon=False)
plt.savefig('figures/detection-snr-freq.pdf')
plt.show()