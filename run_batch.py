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

from conversion_utils import *
from generate_array import *
from generate_uv_coverage import *
from generate_models import *
from sensitivity import calc_sense


# Runtime parameters
ps_model_dir = "models/fialkov_h2"
arr_file     = 'brawl127.txt'                               # Name of antenna array python module file
uv_file      = 'uv_coverages/BRAWLdrift_arrayfile.hkl'      # Name of generated array file from generate_uv_coverage.py
redshift_range = (15, 30)                                   # Set redshift range of interest
n_redshifts    = 15                                        # Number of redshift values to calculate within range

ovro = ('37.240391', '-118.281667', 1184)

regenerate_array       = True                               # Regenerate array geometry
regenerate_uv_coverage = True                               # Regenerate array file via generate_uv_coverage.py
run_calc_sense         = True                                # Run sensitivity calculation
save_snr_table         = False


z, k, ps = load_21cm_model(ps_model_dir)
print z.shape, k.shape, ps.shape
plot_power_spectrum(z, k, ps)


if regenerate_array:
    ant_pos = generate_hexagon(7, 6, 1.0)
    save_antenna_array(ant_pos, "brawl127.txt")
    plot_antenna_array(ant_pos)
else:
    ant_pos = load_antenna_array(arr_file)

if regenerate_uv_coverage:
    #Set other array parameters here
    uv_params = {
        'name': 'BRAWL',
        'Trx': 400 * 1e3  # receiver temp in mK, T_sky is taken care of later
    }

    aa = generate_aa(ovro, ant_pos, uv_params)
    uv = generate_uv_coverage(aa)
    save_uv_coverage(uv, "uv_coverages")
    plot_uv_coverage(uv)

if run_calc_sense:
    ps_z, ps_k, ps_model = load_21cm_model(ps_model_dir)
    z_targets = np.linspace(redshift_range[0], redshift_range[1], n_redshifts)
    zs, snrs_s0, snrs_s1, snrs_s2, snrs_h0, snrs_h1, snrs_h2 = [], [], [], [], [], [], []
    for z in z_targets:
            freq = z2f(z)

            obs_opts = {
                    'model': 'trott',
                    'buff': 0.1,
                    'ndays': 180,
                    'n_per_day': 12,
                    'bwidth': 0.024,
                    'nchan': 128,
                    'no_ns': False,
                    'ps_model': ps_model_dir,
                    'z': z}

            zs.append(z)

            ps_z, ps_k, ps_model = load_21cm_model("models/fialkov_h01")
            ps_model_z = slice_model_along_z(z, ps_z, ps_k, ps_model)
            obs_opts['ps_model'] = "models/fialkov_h01"
            sense_dict = calc_sense(uv_file, ps_model_z, obs_opts, print_report=False)
            snrs_h0.append(sense_dict['snr'])

            ps_z, ps_k, ps_model = load_21cm_model("models/fialkov_h1")
            ps_model_z = slice_model_along_z(z, ps_z, ps_k, ps_model)
            obs_opts['ps_model'] = "models/fialkov_h1"
            sense_dict = calc_sense(uv_file, ps_model_z, obs_opts, print_report=False)
            snrs_h1.append(sense_dict['snr'])

            ps_z, ps_k, ps_model = load_21cm_model("models/fialkov_h2")
            ps_model_z = slice_model_along_z(z, ps_z, ps_k, ps_model)
            obs_opts['ps_model'] = "models/fialkov_h2"
            sense_dict = calc_sense(uv_file, ps_model_z, obs_opts, print_report=False)
            snrs_h2.append(sense_dict['snr'])

            ps_z, ps_k, ps_model = load_21cm_model("models/fialkov_s01")
            ps_model_z = slice_model_along_z(z, ps_z, ps_k, ps_model)
            obs_opts['ps_model'] = "models/fialkov_s01"
            sense_dict = calc_sense(uv_file, ps_model_z, obs_opts, print_report=False)
            snrs_s0.append(sense_dict['snr'])

            ps_z, ps_k, ps_model = load_21cm_model("models/fialkov_s1")
            ps_model_z = slice_model_along_z(z, ps_z, ps_k, ps_model)
            obs_opts['ps_model'] = "models/fialkov_s1"
            sense_dict = calc_sense(uv_file, ps_model_z, obs_opts, print_report=False)
            snrs_s1.append(sense_dict['snr'])

            ps_z, ps_k, ps_model = load_21cm_model("models/fialkov_s2")
            ps_model_z = slice_model_along_z(z, ps_z, ps_k, ps_model)
            obs_opts['ps_model'] = "models/fialkov_s2"
            sense_dict = calc_sense(uv_file, ps_model_z, obs_opts, print_report=True)
            snrs_s2.append(sense_dict['snr'])

            print("")
            #time.sleep(1)

    # Output data as a table
    if save_snr_table:
        snr_table = np.column_stack((zs, snrs_1, snrs_6, snrs_12))
        np.savetxt('output/snr-detection-int-time.txt', snr_table, header='# Redshift SNR_180 SNR_1080 SNR_2160')
        np.savetxt('tex/snr-detection-int-time.txt', snr_table[:82:4], delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')

    # Create pretty plots of SNR vs redshift
    plt.figure("SNR")
    plt.plot(zs, snrs_h0, c='#00cc00', label='Weak heating, hard x-ray')
    plt.plot(zs, snrs_h1, c='#0000cc', label="Normal heating, hard x-ray")
    plt.plot(zs, snrs_h2, c='#cc0000', label="Strong heating, hard x-ray")
    plt.plot(zs, snrs_s0, c='#00cc00', ls="dashed", label='Weak heating, soft x-ray')
    plt.plot(zs, snrs_s1, c='#0000cc', ls="dashed", label="Normal heating, soft x-ray")
    plt.plot(zs, snrs_s2, c='#cc0000', ls="dashed", label="Strong heating, soft x-ray")
    plt.xlabel("Redshift (z)")
    plt.ylabel("Detection SNR")
    plt.axvline(x=f2z(0.04), c='#666666', ls='dashed')
    plt.axvline(x=f2z(0.0875), c='#666666', ls='dashed')
    plt.axvline(x=f2z(0.108), c='#666666', ls='dashed')
    plt.axvspan(f2z(0.0875), f2z(0.108), color='#333333', alpha=0.5, lw=0)
    plt.xlim(redshift_range[0], redshift_range[1])
    plt.minorticks_on()
    plt.legend(frameon=False)
    plt.savefig('figures/detection-snr-int-time.pdf')

    # Same plot, in frequency space instead of redshift
    plt.figure("SNR FREQUENCY")
    fs = z2f(np.array(zs)) * 1e3
    plt.plot(fs, snrs_h0, c='#00cc00', label='Weak heating')
    plt.plot(fs, snrs_h1, c='#0000cc', label="Normal heating")
    plt.plot(fs, snrs_h2, c='#cc0000', label="Strong heating")
    plt.plot(fs, snrs_s0, c='#00cc00', ls="dashed", label='Weak heating')
    plt.plot(fs, snrs_s1, c='#0000cc', ls="dashed", label="Normal heating")
    plt.plot(fs, snrs_s2, c='#cc0000', ls="dashed", label="Strong heating")
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Detection SNR")
    plt.axvspan(87.5, 108, color='#333333', alpha=0.5, lw=0)
    plt.minorticks_on()
    plt.legend(frameon=False)
    plt.savefig('figures/detection-snr-freq-int-time.pdf')
    plt.show()
