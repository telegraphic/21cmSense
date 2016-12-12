"""
# generate_array.py

Utilities for generating antenna array geometries. Use to
"""

import aipy as a
import numpy as np
import pylab as plt
import os

class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.array_params = {}
    def get_ant_params(self, ant_prms={'*':'*'}):
        prms = a.fit.AntennaArray.get_params(self, ant_prms)
        for k in ant_prms:
            top_pos = np.dot(self._eq2zen, self[int(k)].pos)
            if ant_prms[k] == '*':
                prms[k].update({'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]})
            else:
                for val in ant_prms[k]:
                    if   val == 'top_x': prms[k]['top_x'] = top_pos[0]
                    elif val == 'top_y': prms[k]['top_y'] = top_pos[1]
                    elif val == 'top_z': prms[k]['top_z'] = top_pos[2]
        return prms

    def set_ant_params(self, prms):
        changed = a.fit.AntennaArray.set_params(self, prms)
        for i, ant in enumerate(self):
            ant_changed = False
            top_pos = np.dot(self._eq2zen, ant.pos)
            try:
                top_pos[0] = prms[str(i)]['top_x']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[1] = prms[str(i)]['top_y']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[2] = prms[str(i)]['top_z']
                ant_changed = True
            except(KeyError): pass
            if ant_changed: ant.pos = np.dot(np.linalg.inv(self._eq2zen), top_pos)
            changed |= ant_changed
        return changed

    def get_arr_params(self):
        return self.array_params

    def set_arr_params(self, prms):
        for param in prms:
            self.array_params[param] = prms[param]
        return self.array_params

    def view(self):
        ant_pos = []
        for ant in self.ants:
            ant_pos.append(ant.pos)
        plot_antenna_array(ant_pos)


def generate_hexagon(n_side, l, scale_factor):
    """ Generate antennas configured in an hexagonal array.

    Parameters:
        nside (int): Number of antennas on each side of the hexagon
        l (float):   Spacing between antennas in meters
        scale_fac (float): scale factor between rows of antennas

    Returns:
        antpos: a numpy array with shape (n_ant, 3), giving x,y,z
            coords of antennas in meters
    """
    dl = l * scale_factor    #close packed hex
    antpos = []
    cen_y, cen_z = 0, 0
    for row in np.arange(n_side):
        for cen_x in np.arange((2 * n_side - 1) - row):
            dx = row / 2.0
            antpos.append(((cen_x + dx) * l, row * dl, cen_z))
            if row != 0:
                antpos.append(((cen_x + dx) * l, -row * dl, cen_z))
    print "Antenna array has %i antennas" % len(antpos)
    return np.array(antpos)


def plot_antenna_array(antpos):
    """ Plot an antenna array (2D X-Y positions). """
    antpos = np.array(antpos)

    plt.plot(antpos[:,0], antpos[:, 1], 'o', c='#333333')
    plt.xlabel('x-pos [m]')
    plt.ylabel('y-pos [m]')
    #plt.xlim(np.min(antpos[:, 0]) - 2, np.max(antpos[:, 0]) + 2)
    #plt.ylim(np.min(antpos[:, 1]) - 2, np.max(antpos[:, 1]) + 2)
    plt.minorticks_on()
    plt.savefig("figures/antenna-positions.pdf")
    plt.show()

def save_antenna_array(antpos, filename_out, directory_out="array_geometries"):
    """ Save an antenna array geometry to a CSV file.

    Args:
        antpos: Numpy array of antenna xyz positions, shape (n_ant, 3)
        filename_out: Name of output file, e.g. lwa_ovro.txt
        directory_out: Output directory. Defaults to array_geometries/
    """
    fout = os.path.join(directory_out, filename_out)
    np.savetxt(fout, antpos, header="antpos_x antpos_y antpos_z")

def load_antenna_array(filename, directory="array_geometries"):
    """Load an antenna array from file."""
    fin = os.path.join(directory, filename)
    ant_pos = np.genfromtxt(fin)
    return ant_pos

def generate_aa(arr_pos, ant_pos, params):
    """Return the AntennaArray to be used for simulation."""
    antennas = []

    n_ants = len(ant_pos)
    for i in range(n_ants):
        beam = a.fit.Beam(np.array([1.0]))     # Beam is computed at 100 MHz. Beam object is used for UVW calcs
                                               # Which are returned in wavelength (I think)
        antennas.append(a.fit.Antenna(0, 0, 0, beam))
    aa = AntennaArray(arr_pos, antennas)
    p = {}
    for i in range(n_ants):
        top_pos = ant_pos[i]
        p[str(i)] = {'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]}
    aa.set_ant_params(p)
    aa.set_arr_params(params)
    return aa


if __name__ == "__main__":

    ant_pos = generate_hexagon(7, 10, 1.0)
    #plot_antenna_array(antpos)
    save_antenna_array(ant_pos, "brawl127.txt")

    ovro = ('37.240391', '-118.281667', 1184)

    #Set other array parameters here
    params = {
        'name': 'BRAWL',
        'Trx': 500 * 1e3  # receiver temp in mK, T_sky is taken care of later
    }

    aa = generate_aa(ovro, ant_pos, params)
    aa.view()
