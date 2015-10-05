"""
# generate_array.py

Utilities for generating power spectra models & displaying them.
"""

import numpy as np
import pylab as plt
import os
from conversion_utils import *

def load_model(directory):
    """ Load a model power spectrum from a given directory

    Parameters:
        directory: name of directory that holds the power spectrum model. Requires three text files:
                z.txt, which holds redshift information
                K.txt, which holds k-mode numbers
                PS.txt, the power spectrum, with dimensions z and k

    Returns:
        (z, k, ps): tuple of z, k and power spectrum data
    """
    z   = np.genfromtxt(os.path.join(directory, 'z.txt'))
    k   = np.genfromtxt(os.path.join(directory, 'K.txt'))
    ps = np.genfromtxt(os.path.join(directory, 'PS.txt'))
    return z, k, ps

def plot_power_spectrum(z, k, ps):
    plt.imshow(ps, extent=(k[0], k[-1], z[0], z[-1]), aspect='auto', cmap='jet')
    plt.colorbar(label='mK$^2$')
    plt.xlabel('k (1/Mpc)')
    plt.ylabel('Redshift (z)')
    plt.savefig('figures/ps-model.pdf')
    plt.show()

if __name__ == "__main__":
    z, k, ps = load_model('fialkov')

    plot_power_spectrum(z, k, ps)

    z_targets = np.arange(10, 40, 0.25)

    for z_target in z_targets:
        z_idx    = closest(z, z_target)

        model_out = np.column_stack((k, ps[z_idx], np.zeros_like(k)))
        np.savetxt('models/fialkov_z%2.2f_psh.txt' % z_target, model_out, delimiter='\t')

        plt.figure('Model k vs mK')
        plt.plot(model_out[:, 0], model_out[:, 1], c='#333333')
        plt.xlabel('k (1/Mpc)')
        plt.ylabel('mK$^2$')
        plt.title('z = 20.5')
        plt.savefig('figures/fialkov_z%2.2f_psh.pdf' % z_target)
        plt.clf()

