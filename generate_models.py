"""
# generate_array.py

Utilities for generating power spectra models & displaying them.
"""

import numpy as np
import pylab as plt
import os
from conversion_utils import *

def load_21cm_model(directory):
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

def plot_power_spectrum(z, k, ps, show=True, savefig=False):
    plt.imshow(np.log(ps), extent=(k[0], k[-1], z[0], z[-1]), aspect='auto', cmap='inferno')
    plt.colorbar(label='log mK$^2$')
    plt.xlabel('k (1/Mpc)')
    plt.ylabel('Redshift (z)')
    if savefig:
        plt.savefig('figures/ps-model.pdf')
    if show:
        plt.show()

def slice_model_along_z(z_slice, z, k, ps):
    """ Slice a 2D (z vs k) power spectrum model along z, for a given redshift z.

    Arguments:
        z_slice (float): Redshift value along which to slice. Closest value of z will be returned.
        z (np.array): redshift values for power spectrum
        k (np.array): wave number values for power spectrum
        ps (np.array): 2D power model, in terms of z vs k

    Returns:
        model_out (np.array): Output model for closest slice along z-axis. Output is given in the
            same format as 21cmfast code, i.e. three columns, (k, ps, 0) (no idea what the values
            in the last column actually mean).
    TODO: Apply interpolation along axes instead of just picking closest value.
    """
    z_idx  = closest(z, z_slice)
    model_out = np.column_stack((k, ps[z_idx], np.zeros_like(k)))
    return model_out

if __name__ == "__main__":
    plt.figure("Fialkov models", figsize=(10, 8))

    z, k, ps = load_21cm_model('models/fialkov_h01')
    plt.subplot(3,2,2)
    plt.title("Weak heating, hard x-ray")
    plot_power_spectrum(z, k, ps, show=False)

    z, k, ps = load_21cm_model('models/fialkov_h1')
    plt.subplot(3,2,4)
    plt.title("Normal heating, hard x-ray")
    plot_power_spectrum(z, k, ps, show=False)

    z, k, ps = load_21cm_model('models/fialkov_h2')
    plt.subplot(3,2,6)
    plt.title("Strong heating, hard x-ray")
    plot_power_spectrum(z, k, ps, show=False)

    z, k, ps = load_21cm_model('models/fialkov_s01')
    plt.subplot(3,2,1)
    plt.title("Weak heating, soft x-ray")
    plot_power_spectrum(z, k, ps, show=False)

    z, k, ps = load_21cm_model('models/fialkov_s1')
    plt.subplot(3,2,3)
    plt.title("Normal heating, soft x-ray")
    plot_power_spectrum(z, k, ps, show=False)

    z, k, ps = load_21cm_model('models/fialkov_s2')
    plt.subplot(3,2,5)
    plt.title("Strong heating, soft x-ray")
    plot_power_spectrum(z, k, ps, show=False)

    plt.tight_layout()
    plt.savefig("figures/fialkov_models.pdf")
    plt.show()