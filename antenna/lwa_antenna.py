import numpy as np
import pylab as plt
import h5py

from scipy.interpolate import interp1d, RectBivariateSpline

H5_DATA = "antenna/antenna_data.h5"

class LwaBeamPattern(object):
    """ Beam pattern class for the LWA antenna """
    def __init__(self):
        self.h5 = h5py.File(H5_DATA)
        self.data   = self.h5["data"][:]
        self.theta  = self.h5["theta"][:]
        self.phi    = self.h5["phi"][:]
        self.freqs  = self.h5["frequency"][:]
        
        self.freq_interpolator = interp1d(self.freqs, self.data, axis=0, kind='cubic')
        
        self.generated_beam_data  = None
        self.generated_beam_freq  = None
        self.generated_beam_theta = None
        self.generated_beam_phi   = None
        
    def generate(self, freq, n_theta=None, n_phi=None, linear=False):
        """ Generate a beam for the LWA antenna at a given frequency 
        
        Parameters
        ----------
        freq: float
            Frequency in MHz for beam pattern
        n_theta: int
            Number of points to generate on theta axis (0 to 360 deg). Optional.
        n_phi: int
            Number of points to generate on phi axis (-90 to 90 deg). Optional.
        linear: bool
            Output beam pattern in linear units. Defaults to False (dB).
        
        Notes
        -----
        This uses EM simulations done in 4NEC2, every 5MHz, 1 degree resolution,
        from 20 MHz to 90 MHz. File used:
        `LWA_DUAL_0-90DEG-DRIVE_300x300x1cmGRID_Symbols.nec`
        
        Interpolation (cubic) is done over frequency, to generate a 1 degree 
        resolution pattern at the new frequency. If more points are required 
        (i.e. higher res) then interpolation is done once again to achieve this,
        this time using RectBivariateSpline().
        """
        beam = self.freq_interpolator(freq)
        
        if n_theta is not None and n_phi is not None:
            theta_new = np.linspace(-90, 90, n_theta)
            phi_new   = np.linspace(0, 360, n_phi)
            
            # Do fit in linear space
            beam_lin = 10**(beam / 10.0)
            b_interp = RectBivariateSpline(self.phi, self.theta, beam_lin)
            beam     = 10 * np.log10(b_interp(phi_new, theta_new))
        else:
            theta_new, phi_new = self.theta, self.phi
            
        self.generated_beam_data  = beam
        self.generated_beam_theta = theta_new
        self.generated_beam_phi   = phi_new
        self.generated_beam_freq  = freq
        
        return beam
    
    def view(self, show=True):
        """ View generated beam """
        
        if self.generated_beam_data is None:
            raise RuntimeError("Beam pattern not generated yet. Run generate() first.")

        plt.figure(figsize=(8,4))
        plt.subplot(121)
        tmin, tmax = self.generated_beam_theta[0], self.generated_beam_theta[-1]
        pmin, pmax = self.generated_beam_phi[0], self.generated_beam_phi[-1]
        plt.imshow(10**(self.generated_beam_data / 10), extent=(tmin, tmax, pmin, pmax), aspect='auto')
        plt.xlabel("Theta [deg]")
        plt.ylabel("Phi [deg]")
        #plt.colorbar(orientation='horizontal')
        
        plt.subplot(122)
        beam_slice = self.generated_beam_data[self.generated_beam_data.shape[0]/2]
        print self.generated_beam_phi.shape, beam_slice.shape
        plt.plot(self.generated_beam_theta, beam_slice, c='#333333')
        plt.xlabel("Theta [deg]")
        plt.ylabel("Normalized gain [dB]")
        plt.xlim(-91, 91)
        plt.ylim(-30, 3)
        plt.minorticks_on()
        
        plt.tight_layout()
        if show:
            plt.show()

    def view2(self, show=True):
        """ View generated beam """
        
        if self.generated_beam_data is None:
            raise RuntimeError("Beam pattern not generated yet. Run generate() first.")

        plt.figure(figsize=(8,4))
        plt.subplot(121)
        tmin, tmax = self.generated_beam_theta[0], self.generated_beam_theta[-1]
        pmin, pmax = self.generated_beam_phi[0], self.generated_beam_phi[-1]
        plt.imshow(10**(self.generated_beam_data / 10), extent=(tmin, tmax, pmin, pmax), aspect='auto')
        plt.xlabel("Theta [deg]")
        plt.ylabel("Phi [deg]")
        #plt.colorbar(orientation='horizontal')
        
        plt.subplot(122)
        beam_slice = self.generated_beam_data[self.generated_beam_data.shape[0]/2]
        beam_slice2 = 10*np.log10((10.**(beam_slice / 10.))**2)
        beam_slice = 10*np.log10(10.0**(beam_slice / 10.))
        print self.generated_beam_phi.shape, beam_slice.shape
        plt.plot(self.generated_beam_theta, beam_slice, c='#333333')
        plt.plot(self.generated_beam_theta, beam_slice2, c='#cc0000')
        plt.xlabel("Theta [deg]")
        plt.ylabel("Normalized gain [dB]")
        plt.xlim(-91, 91)
        #plt.ylim(-30, 3)
        plt.minorticks_on()
        
        plt.tight_layout()
        if show:
            plt.show()
    
    def compute_bm_eff(self):
        """ Compute the omega prime (effective beam) from Pober et. al. 2013
        
        Notes:
            Computes the (solid angle of the power beam) squared 
            divided by the solid angle of (the sqaure of the power beam).
        
        """
        beam_lin = 10.0**(self.generated_beam_data / 10.0)
        beam_lin2 = beam_lin**2
        theta = np.deg2rad(self.generated_beam_theta)  + np.pi/2
        beam_weighted = beam_lin * np.sin(theta)
        beam_weighted2 = beam_lin2 * np.sin(theta)
        bm  = np.average(beam_weighted) / np.max(beam_weighted) 
        bm2 = np.average(beam_weighted2) / np.max(beam_weighted2) 
        bm_eff = bm**2 / bm2
        
        #print bm, bm2, bm_eff
        return bm_eff

    def compute_fwhm(self):
        pass

if __name__ == "__main__":
    lwa = LwaBeamPattern()
    lwa.generate(60)
    #lwa.view2()
    lwa.bm_eff()

        
        