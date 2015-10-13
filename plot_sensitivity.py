"""
# plot_sensivity.py

Function for loading sensitivty output files (.hkl) and plotting them.
"""
import pylab as plt
import hickle as hkl

def plot_sensitivity_file(filename, show=True):
    d = hkl.load(filename)
    kmag = d['ks']
    p21 = d['p21']
    sense1d = d['errs']

    print p21

    plt.subplot(211)
    plt.plot(kmag, p21)
    plt.xlabel("$k_{mag}$")
    plt.ylabel("P21($k_{mag}$)")
    plt.subplot(212)
    plt.plot(kmag, sense1d)
    plt.xlabel("$k_{mag}$")
    plt.ylabel("Sensivitity (mK$^2$)")

    if show:
        plt.show()

if __name__ == "__main__":

    plot_sensitivity_file("sensitivities/BRAWLdrift_none_0.081_fialkov_s01.hkl", show=False)
    plot_sensitivity_file("sensitivities/BRAWLdrift_none_0.081_fialkov_s1.hkl", show=False)
    plot_sensitivity_file("sensitivities/BRAWLdrift_none_0.081_fialkov_s2.hkl", show=True)
