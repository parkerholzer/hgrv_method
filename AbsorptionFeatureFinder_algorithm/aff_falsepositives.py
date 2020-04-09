import pandas as pd
import numpy as np
from spectra_functions import lowerSNR, findabsorptionfeatures
from multiprocessing import Pool

l = open("nso_spectrum.txt").readlines()
spec1 = pd.read_csv("nso_spectrum.txt", skiprows=17, skipfooter = 1864179-1040309 ,delim_whitespace=True, engine='python', header=None)
spec2 = pd.read_csv("nso_spectrum.txt", skiprows=1040310, skipfooter = 1864179-1286410 ,delim_whitespace=True, engine='python', header=None)
spec3 = pd.read_csv("nso_spectrum.txt", skiprows=1286411,delim_whitespace=True, engine='python', header=None)

spec1.columns = ["A", "B", "C", "D", "E", "F"]
spec2.columns = ["A", "B", "C", "D", "E", "F"]
spec3.columns = ["A", "B", "C", "D"]

Wvl = 1e8/np.hstack((np.hstack((spec1.A.values, spec2.A.values)), spec3.A.values))

keep = np.where((Wvl > 5000) & (Wvl < 6000))[0]
wvl = np.array(list(Wvl[keep]))
trueflx = np.ones(len(keep))

def linedepths(i):
    np.random.seed(i)
    flx = lowerSNR(trueflx, 500)
    wvbounds, minwvs, minflxs, maxflxs = findabsorptionfeatures(wvl, flx, pix_range=25, minlinedepth=0.0, alpha=0.05, gamma=0.01)
    return np.array(maxflxs - minflxs)

rslt = np.hstack(np.array(Pool(processes=20).map(linedepths, np.arange(20))))
df = pd.DataFrame({"LineDepths": rslt})
df.to_csv("aff_falsepos.csv")

