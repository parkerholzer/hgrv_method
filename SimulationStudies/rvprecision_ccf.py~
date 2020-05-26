import os
from glob import glob
import numpy as np
import pandas as pd
from multiprocessing import Pool
from astropy.constants import c
from scipy import interpolate
#import statsmodels.api as sm
#import h5py

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
    #ccf, e_ccf = coadd_ccfs(ccfs, ccf_errs, orders, blue=34, red=74) # defautls from config file
    
    #ccf, e_ccf = OrderAveragedCCF(masks, V, wave, spec, errs, colour_correct=False)
    # Single-"order" version
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
    
    #return V, ccf, e_ccf, v, delta_v, chi2
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

#keep3 = np.where((temp.Wavelength.values >= 4850) & (temp.Wavelength.values <= 6800))[0]
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

#rvhats = []
#for i in range(100):
#    rvhats.append(estRV(wvl, trueflx, 1000))

#def mapfunc(i):
#    rvs = []
#    for i in range(10):
#        rvs.append(estRV(wvl, trueflx, 200))
#    return rvs

#with Pool(processes=10) as pool:
#    rvhats = np.array(pool.map(mapfunc, range(100)))

#rvhats = []
#for i in range(100):
#    p = Process(target=mapfunc, args=(i,))
#    rvhats.append(p)
#    p.start()

#for p in rvhats:
#    p.join()

def mapfunc(i):
    return estRV(wvl, trueflx, 250)

RVs = []
for k in range(20):
    #rvhats = []
    #for i in range(100):
    #    rvhats.append(estRV(wvl, trueflx, 250))
    #rvhats = np.array(rvhats).flatten()
    #with Pool(processes=10) as pool:
    #    rvhats = pool.map(mapfunc, range(100))
    rvhats = [estRV(wvl, trueflx, 250) for i in range(100)]
    print(np.sqrt(np.mean((np.array(rvhats) - rv)**2)))
    RVs.append(rvhats)

RVs = np.array(RVs).flatten()
df = pd.DataFrame({"RV": RVs})
df.to_csv("ccf_sn250_rv0.01.csv")






#spec = h5py.File('chris/res-1000-1years_long_id10.h5', 'r')
#temp = pd.DataFrame({"Wavelength": np.array(spec['lambdas']), "Flux": np.array(spec['quiet']/np.max(spec['quiet']))})

#def estRV(twvl, tflx, snr):
#    flx = lowerSNR(tflx, snr)
#    rvest = ccf_rv(twvl, flx, mask=harps_mask)
#    return rvest

#keep = np.where((temp.Wavelength.values > 4800) & (temp.Flux.values > 0.00001))[0]
#offset = ccf_rv(temp.Wavelength.values[keep], temp.Flux.values[keep], mask=harps_mask)
#tempwvl = np.copy(temp.Wavelength.values[keep])
#tempflx = np.copy(temp.Flux.values[keep])
#rv = 0.01                                                            
#owvl = np.copy(tempwvl)*(1 + rv/c.value)                     
#flx = np.copy(tempflx)                                            
#ds_spec = interpolate.interp1d(owvl, flx, kind='cubic', fill_value='extrapolate')                                                                     
#tflx = ds_spec(tempwvl)

#def mapfunc(i):
#    return estRV(tempwvl, tflx, 1000)-offset

#with Pool(processes=20) as pool:
#rvs = []
#for j in range(100):
#    rvs.append(mapfunc(j))

#print(np.sqrt(np.mean((np.array(rvs) - rv)**2)))


#t=100
#oflx = np.array(spec['active'][100,:]/np.max(spec['active'][100,:]))
#keep = np.where((temp.Wavelength.values > 4800) & (temp.Flux.values > 0.00001))[0]

#offset = ccf_rv(temp.Wavelength.values[keep], temp.Flux.values[keep], mask=harps_mask)
#print("Offset = ", offset)

#def estRV(twvl, tflx, snr):
#    flx = lowerSNR(tflx, snr)
#    return ccf_rv(twvl, flx, mask=harps_mask)

#rv = []
#np.random.seed(7)
#for i in range(100):
#    np.random.seed(i)
#    t = np.random.randint(0, 2*365, 1)
#    owvl = np.array(spec['lambdas'])
#    oflx = np.array(spec['active'][t,:]/np.max(spec['active'][t,:]))
#    tempflx = np.array(spec['quiet']/np.max(spec['quiet']))
#    keep = np.where((owvl > 4800) & (tempflx > 0.00001))[0]
#    rv.append(estRV(owvl[keep], oflx[0,keep], 10000) - offset)

#with Pool(processes = 4) as pool:
#    rv = np.array(pool.map(mapfunc, range(100)))

#print(np.sqrt(np.mean(np.array(rv)**2)))

#temp = pd.read_csv("smoothednso_expres.csv")

#ftrs = pd.read_csv("GoodFeatures.csv")
#hgrv = np.zeros(len(temp.Wavelength.values))
#keep2 = []
#for i in range(ftrs.shape[0]):
#    keep = np.where((temp.Wavelength.values >= ftrs.Gauss_mu.values[i] - 1) & (temp.Wavelength.values <= ftrs.Gauss_mu.values[i] + 1))[0]
#    keep2 = keep2 + list(np.where((temp.Wavelength.values >= ftrs.Wv_lbounds.values[i]) & (temp.Wavelength.values <= ftrs.Wv_ubounds.values[i]))[0])
#    coef = -np.sqrt(np.sqrt(np.pi))*ftrs.Gauss_mu.values[i]*ftrs.Gauss_amp.values[i]/(299792458*np.sqrt(2*ftrs.Gauss_sig.values[i]))
#    hgrv[keep] = hgrv[keep] + coef*hermitegauss1(temp.Wavelength.values[keep], ftrs.Gauss_mu.values[i], ftrs.Gauss_sig.values[i])

#keep2 = np.unique(keep2)

#def estRV(twvl, tflx, snr):
#    flx = lowerSNR(tflx, snr)
#    rvest_ccf = ccf_rv(twvl, flx, mask=harps_mask)
#    diff = flx[keep2] - temp.Flux.values[keep2]
#    mdl = sm.WLS(diff, hgrv[keep2], weights = 1./tflx[keep2]).fit()
#    rvest_hgrv = mdl.params[0]
#    return np.array([rvest_ccf, rvest_hgrv])



#print('Initializing Arrays')
#p = np.linspace(-2, 1, 4)
#s = np.array([100, 300, 500, 700, 900]) #np.linspace(100, 250, 16)

#print("Let's Run Through Some RV Shifts")
#for j in range(4):
#    rv = 10**(p[j])
#    wvl = np.copy(temp.Wavelength.values)*(1 + rv/c.value)
#    flx = np.copy(temp.Flux.values)
#    ds_spec = interpolate.interp1d(wvl, flx, kind='cubic', fill_value='extrapolate')
#    tflx = ds_spec(temp.Wavelength.values)
#    for i in range(5):
#        # We're going to save each RV/SNR run independently
#        # Skip RV/SNR combinations that have already been run
#        if os.path.isfile(f"./Scratch/vals2000_rv{j}_snr{i}.csv"):
#            continue
#        
#        # Define SNR and new SNR-Specific function
#        snr = s[i]
#        pbar = tqdm(total=500, desc=f"RV ({j+1}/4), SNR ({i+1}/5)")
#        #pbar = tqdm(total=128, desc=f"RV ({j+1}/4), SNR ({i+1}/5)")s
#        def mapfunc(k):
#            rvest = estRV(temp.Wavelength.values, tflx, snr)
#            pbar.update() # update progress bar
#            return rvest
#        
#        # Iterate through many simulations
#        with Pool(processes=4) as pool:
#            vhat = np.array(pool.map(mapfunc, np.arange(2000)))
#            #vhat = np.array(pool.map(mapfunc, np.arange(500)))
#        
#        # Save values for this specific RV/SNR combo
#        mse = np.mean((vhat[:,0] -rv)**2)
#        bias = np.mean(vhat[:,0])-rv
#        var = np.mean((vhat[:,0]-np.mean(vhat[:,0]))**2)
#        with open(f"./Scratch/ccfvals2000_rv{j}_snr{i}.csv","w+") as file:
#            np.savetxt(file,[mse,bias,var],delimiter=',')
#        mse = np.mean((vhat[:,1] -rv)**2)
#        bias = np.mean(vhat[:,1])-rv
#        var = np.mean((vhat[:,1]-np.mean(vhat[:,1]))**2)
#        with open(f"./Scratch/hgrvvals2000_rv{j}_snr{i}.csv","w+") as file:
#            np.savetxt(file,[mse,bias,var],delimiter=',')
#        pbar.close() # close progress bar

## Make Data Frames
#MSE = pd.DataFrame({"SNR": s})
#BIAS = pd.DataFrame({"SNR": s})
#VAR = pd.DataFrame({"SNR": s})
#for pk in p:
#    MSE["%.1f"%pk] = np.zeros(5)
#    BIAS["%.1f"%pk] = np.zeros(5)
#    VAR["%.1f"%pk] = np.zeros(5)

## Load in information from saved files
#for j in range(4):
#    for i in range(5):
#        vals = np.loadtxt(f"./Scratch/ccfvals2000_rv{j}_snr{i}.csv",delimiter=',')
#        # Find wanted values and write to data frames
#        MSE.iloc[i,j], BIAS.iloc[i,j], VAR.iloc[i,j] = vals

## Save dataframes
#MSE.to_csv("ccf_mse.csv")
#BIAS.to_csv("ccf_bias.csv")
#VAR.to_csv("ccf_var.csv")

## Make Data Frames
#MSE = pd.DataFrame({"SNR": s})
#BIAS = pd.DataFrame({"SNR": s})
#VAR = pd.DataFrame({"SNR": s})
#for pk in p:
#    MSE["%.1f"%pk] = np.zeros(16)
#    BIAS["%.1f"%pk] = np.zeros(16)
#    VAR["%.1f"%pk] = np.zeros(16)

## Load in information from saved files
#for j in range(4):
#    for i in range(5):
#        vals = np.loadtxt(f"./Scratch/hgrvvals2000_rv{j}_snr{i}.csv",delimiter=',')
#        # Find wanted values and write to data frames
#        MSE.iloc[i,j], BIAS.iloc[i,j], VAR.iloc[i,j] = vals

## Save dataframes
#MSE.to_csv("hgrv_mse.csv")
#BIAS.to_csv("hgrv_bias.csv")
#VAR.to_csv("hgrv_var.csv")

