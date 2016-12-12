import aipy as a, numpy as n, os

class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.array_params = {}
    def get_ant_params(self, ant_prms={'*':'*'}):
        prms = a.fit.AntennaArray.get_params(self, ant_prms)
        for k in ant_prms:
            top_pos = n.dot(self._eq2zen, self[int(k)].pos)
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
            top_pos = n.dot(self._eq2zen, ant.pos)
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
            if ant_changed: ant.pos = n.dot(n.linalg.inv(self._eq2zen), top_pos)
            changed |= ant_changed
        return changed 
    def get_arr_params(self):
        return self.array_params
    def set_arr_params(self, prms):
        for param in prms:
            self.array_params[param] = prms[param]
            if param == 'dish_size_in_lambda':
                FWHM = 2.35*(0.45/prms[param]) #radians
                self.array_params['obs_duration'] = 60.*FWHM / (15.*a.const.deg)# minutes it takes the sky to drift through beam FWHM
            if param == 'antpos':
                bl_lens = n.sum(n.array(prms[param])**2,axis=1)**.5
        return self.array_params

#===========================ARRAY SPECIFIC PARAMETERS==========================

#Set antenna positions here; for regular arrays like Hera we can use an algorithm; otherwise antpos should just be a list of [x,y,z] coords in light-nanoseconds
nside = 7. #hex number
L = 10.0
scalefac = 1.0 # 1400. / 1212
dL = L * scalefac    #close packed hex
antpos = []
cen_y, cen_z = 0, 0
for row in n.arange(nside):
    for cen_x in n.arange((2*nside-1)-row):
        dx = row/2
        antpos.append(((cen_x + dx)*L, row*dL, cen_z))
        if row != 0:
            antpos.append(((cen_x + dx)*L, -row*dL, cen_z))

import pylab as plt
import numpy as np

antpos = np.array(antpos)
plt.plot(antpos[:,0], antpos[:, 1], 'o', c='#333333')
plt.xlabel('x-pos [m]')
plt.ylabel('y-pos [m]')
plt.xlim(np.min(antpos[:, 0]) - 2, np.max(antpos[:, 0]) + 2)
plt.ylim(np.min(antpos[:, 1]) - 2, np.max(antpos[:, 1]) + 2)
plt.minorticks_on()
plt.savefig("figures/antenna-positions.pdf")
plt.show()


#Set other array parameters here
prms = {
    'name': os.path.basename(__file__)[:-3],    #remove .py from filename
    'loc': ('37.240391', '-118.281667', 1184),  # Owens Valley
    'antpos': antpos,
    'beam': a.fit.Beam2DGaussian,
    'dish_size_in_lambda': 2.0,
    'Trx': 500 * 1e3  # receiver temp in mK, T_sky is taken care of later
}

#=======================END ARRAY SPECIFIC PARAMETERS==========================

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    print freqs
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        beam = prms['beam'](freqs, xwidth=n.pi/8, ywidth=n.pi/8)
        # as it stands, the size of the beam as defined here is not actually used anywhere in this package,
        # # but is a necessary parameter for the aipy Beam2DGaussian object
        antennas.append(a.fit.Antenna(0, 0, 0, beam))
    aa = AntennaArray(prms['loc'], antennas)
    p = {}
    for i in range(nants):
        top_pos = prms['antpos'][i]
        p[str(i)] = {'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]}
    aa.set_ant_params(p)
    aa.set_arr_params(prms) 
    return aa

def get_catalog(*args, **kwargs): return a.src.get_catalog(*args, **kwargs)
