    #! /usr/bin/env python3
# coding: utf-8
"""
Module with functions for identifying the best tradeoff for the detection of accreting planets with the SPHERE+ medium-resolution IFU
"""
__author__ = 'M. Bonnefoy, R. Gratton'

__all__ = ['band_def','throughput','pu','thermal_background','background_noise','noise_budget','dist_modulus']


import numpy as np

def band_def(name):
    """ Return tabulated values for the wavelength, zero points, Strehl ratio and band widths to expect for observations in different bands

    Inputs:
        name: string. Name of the observation band. Can be either V, J,  or H band

    Outputs:
        tupple giving the central wavelength in microns, zero point in ph/s/A/cm^2, Strehl ratio in percent, and band width in micron for the given band.
    """
    if name == 'V':
        WAVELENGTH=0.65 #microns
        ZP=1050.  #ph/s/A/cm^2
        STREHL=0.5 #percents
        EW=1. #width of the band in microns
    elif name == 'J':
        WAVELENGTH=1.25
        ZP=177.4
        STREHL=0.75
        EW=0.1
    elif name == 'H':
        WAVELENGTH=1.6
        ZP=98.8
        STREHL=0.85
        EW=0.1
    else:
        raise ValueError

    return (WAVELENGTH, ZP, STREHL, EW)

def throughput(mode, lenslet=0.7, fibre=0.23, telescope=0.7, sphere=0.6,spectrograph=0.7, detector=0.9):
    """ Computes the full throughput assuming either a lenslet or fiber-based IFU spectrograph
    Inputs:
        mode: string. Spectrograph design. Can be either "lenslet" or "fiber"
        lenslet: float. Throughput of the lenslets
        fibre: float. Throughput of the fibers (including injection)
        telescope: float. telescope throughput (percent)
        sphere: float. SPHERE instrument throughput (comon-path + IR or VIS arm)
        spectrograph: float. spectrograph throughput
        detector: float. detector throughput

    Outputs:
        float value giving the throughput"""

    if mode =='lenslet':
        return lenslet*telescope*sphere*spectrograph*detector
    elif mode == 'fiber':
        return fibre*telescope*sphere*spectrograph*detector
    else:
        raise ValueError

def pu(fwhm=2.5, pacd=3,pald=1):
    """ Computes the  number of pixel units used on the detector to form a given voxel
    Inputs:
        :fwhm: float. Size of the fwhm in spaxel size unit
        :pacs: float. pixels accross dispersion
        :pald: float. pixels along dispersion
    Outputs:
        float value giving the number of detector pixels integrated translating to a voxel"""

    return pacd*pald*(fwhm**2)

def thermal_background(spectemp):
    """Computes the flux level due to the thermal background
    Inputs:
        :spectemp: float.  Temperature of the instrument in Kelvin
    Outputs:
        float value giving the thermal background flux level in e-/s/pixel"""

    return 5.5*10**(0.044*(spectemp-293.))

def background_noise(bandname, intime, spectemp, fwhm=2.5, pacd=3,pald=1):
    """Computes the background noise
    Inputs:
        :bandname: string. Name of the observation band. Can be either V, J,  or H band
        :intime: float. total integration time in second
        :spectemp: float. Temperature of the instrument in Kelvin
        :fwhm: float. Size of the fwhm in spaxel size unit
        :pacs: float. pixels accross dispersion
        :pald: float. pixels along dispersion
    Outputs:
        float value giving the background noise level in electrons unit"""


    #pixels used for the integration
    pix_used=pu(fwhm=fwhm, pacd=pacd, pald=pald)

    #thermal_noise per spectral channel
    if bandname == 'V':
        tn_psc=0.0
    elif bandname == 'J' or  bandname == 'H':
        tb=thermal_background(spectemp)
        tn_psc=tb*np.sqrt(pix_used)
    else:
        raise ValueError

    #integrated thermal noise
    tn=np.sqrt(intime*pix_used)*tn_psc

    return tn

def noise_budget(bandname, magobj, contrast, EW_emissionline, central_wl, specres, DIT, NDIT, spectemp, ronlevel, mode="lenslet", nchannels=35, fwhm=2.5, pacd=3, pald=1, throughput_lenslet=0.7,throughput_fiber=0.23):
    """ Computes the  noise sources for a given observation with the IFU os SPHERE+
    Inputs:
        :bandname: string. Name of the observation band. Can be either V, J,  or H band
        :magobj: float. apparent magnitude of the companion
        :contrast: float. planet-to-star contrast in magnitudes
        :EW_emissionline: float. Equivalent width in Angstroms of the emission line for an putative accreting object
        :central_wl: float. Central wavelength in microns.
        :specres: float.  resolution  R of the spectrograph
        :DIT: float. individual data integration time
        :NDIT: float. number of integrations
        :spectemp: float. Temperature of the instrument in Kelvin
        :ronlevel: float. read-out-noise level in e-
        :mode: string. Spectrograph design. Can be either "lenslet" or "fiber"
        :nchanels: float. Number of spectral channels
        :fwhm: float. Size of the fwhm in spaxel size unit
        :pacs: float. pixels accross dispersion
        :pald: float. pixels along dispersion
        :throughput_lenslet: float. Throughput of the lenslets
        :throughput_fiber: float. Throughput of the fibers (including injection)
    Outputs:
        tuple giving the companion signal, companion's photon noise, star's photon noise, thermal noise, ron,  and calibration noise in electron units"""

    TEL_DIAMETER=800 #telescope diameter in cm
    CENTRAL_OBSTRUCTION = 0.9 #central obstruction

    #collecting area
    telcol_area=np.pi*(TEL_DIAMETER/2.)**2*CENTRAL_OBSTRUCTION #collecting area. cm**2

    #basic numerical values for a given observing band.
    wavelength, zp, strehl, ew=band_def(bandname)

    #channel widths in angstroms
    channel_width=wavelength/(2*specres)*1e4

    #photons per second detected for a given target
    tp=throughput(mode, lenslet=throughput_lenslet, fibre=throughput_fiber)
    pd=tp*telcol_area*channel_width*zp*strehl #photon detected for Vega (zero point)
    pdt=pd*10**(-0.4*magobj) #photon detected for a given target

    int_time=DIT*NDIT

    #companion photon noise
    signal_comp=int_time*pdt*EW_emissionline/(channel_width*3.)*ew #signal in the HI line
    photon_noise=np.sqrt(signal_comp)

    #stellar photon noise
    psfpos=4e-4*(1.-strehl)
    signal = int_time*pdt*psfpos*10**(0.4*contrast)
    photon_noise_star=np.sqrt(signal)

    #thermal noise
    thermal_noise=background_noise(bandname, int_time, spectemp)

    #ron
    pix_used=pu(fwhm=fwhm, pacd=pacd,pald=pald)
    ron_comp=ronlevel*np.sqrt(NDIT*pix_used)

    #calibration_noise
    calnoise=np.sqrt(np.nanmax([photon_noise_star,signal])**2+ron_comp**2+thermal_noise**2)/np.sqrt(nchannels)
    snr=signal_comp/np.sqrt(photon_noise**2+thermal_noise**2+ron_comp**2+photon_noise_star**2+calnoise**2)

    return (snr, signal_comp, photon_noise, photon_noise_star, thermal_noise, ron_comp, calnoise)


def dist_modulus(distance):
    return 5-5*np.log10(distance)
