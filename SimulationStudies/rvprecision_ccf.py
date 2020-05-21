import os
from glob import glob
import numpy as np
import pandas as pd
from astropy.constants import c
from scipy import interpolate

from expres.rv.ccf import (generate_ccf_mask, order_wise_ccfs, coadd_ccfs, fit_rv)

harps_mask = '/expres/lib/Data/xcor_masks/HARPS/G2.mas'
def ccf_rv(wavelength, flux, errs=None,
           vspace=None, Vrange=1e7,
           mask='/expres/lib/Data/xcor_masks/ESPRESSO/ESPRESSO_G2.fits',
           vwidth=None, platescale=1., window='cos',
           V0=0):
    wave = np.array([wavelength.copy()])
    spec = np.array([flux.copy()])
    if errs is None:
        errs = np.sqrt(spec.copy()) # spec is continuum normalized though?
    else:
        errs = np.array([errs.copy()])
    
    sig_v, masks = generate_ccf_mask(mask, wave, vwidth=vwidth)
    
    if vspace is not None:
        sig_v = float(vspace)
        
    V = V0 + np.arange(-Vrange, Vrange, sig_v / platescale)
    orders, ccfs, ccf_errs = order_wise_ccfs(masks, V, wave, spec, errs, window=window, ccor=False)
    ccf = ccfs[0]
    e_ccf = ccf_errs[0]
    try:
        v, delta_v = fit_rv(V, ccf, e_ccf, fit_opts={'xtol':1e-10, 'gtol':1e-10,
                                                'ftol':1e-10,'epsfcn':1e-5,
                                                'factor': 1.}, debug=True)
    except (ValueError, RuntimeError) as e:
        print("Fitting Error")
        v, delta_v, chi2 = -9999999999, -9999999999, -1
    
    # Correct δv by platescale
    delta_v = delta_v * np.sqrt(platescale)
    
    return v/100.0

def hermitegauss1(x, mu, sig):
    z = (x - mu)/sig
    return 2*z*np.exp(-(z**2)/2)/np.sqrt(sig*2*np.sqrt(np.pi))


def lowerSNR(flux, sn):
    orig_mean = np.mean(flux)
    scaled_spec = np.array(flux)*(sn**2)/orig_mean
    return np.array(np.random.poisson(lam=scaled_spec, size=len(flux))*(orig_mean/(sn**2)))

temp = pd.read_csv("51peg/217014smoothtemp.csv")
ospec = pd.read_csv("51peg/217014_190605.1077ctd.csv")
ospec = ospec.iloc[np.where((ospec.Wavelength.values >= 4850) & (ospec.Wavelength.values <= 6800))[0]]


wvl = np.array(list(ospec.Wavelength.values))
tempfunc = interpolate.interp1d(temp.Wavelength.values, temp.Flux.values, kind='cubic', fill_value='extrapolate')
tflx2 = tempfunc(wvl)
offset = ccf_rv(wvl, tflx2, mask=harps_mask)
print(offset)

def estRV(twvl, tflx, snr):
    flx = lowerSNR(tflx, snr)
    rvest = ccf_rv(twvl, flx, mask=harps_mask)
    return rvest-offset

rv = 0.01
doppfact = 1 + rv/299792458
ds_spec = interpolate.interp1d(doppfact*temp.Wavelength.values, temp.Flux.values, kind='cubic', fill_value='extrapolate')
trueflx = ds_spec(wvl)

RVs = []
for k in range(20):
    rvhats = [estRV(wvl, trueflx, 250) for i in range(100)]
    print(np.sqrt(np.mean((np.array(rvhats) - rv)**2)))
    RVs.append(rvhats)

RVs = np.array(RVs).flatten()
df = pd.DataFrame({"RV": RVs})
df.to_csv("ccf_sn250_rv0.01.csv")

