# hgrv_2019paper
The datasets, algorithms, and code associated with Holzer et al. (in prep).

There are 5 Jupyter notebooks with Python 3 code that are designed to give interested readers a more hands on understanding of each portion of the paper. Algorithm 1 described in Section 2 of the paper is programmed in as a function in the file "abfeature_functions.py", and the reader can practice using it on the median spectrum for 51 Pegasi in the notebook "Find Absorption Features.ipynb".

The reader can make more sense of the mathematics of a Doppler-shifted Gaussian as described in Section 3 with the two files "Hermite-Gaussian coefficients math.ipynb" and "Multiple Features Test.ipynb". The first gives the coefficient solutions of Lemma 2 as well as the standardized approximation error whose limit is given in Theorem 1. The second gives the reader an opportunity to test what radial velocities and parameters for a Gaussian absorption feature allow for the problem to be reduced to simple linear regression.

The file "Gaussian Fit Parameters.ipynb" gives the details of fitting a Gaussian to each absorption feature detected in the median spectrum of 51 Peg. These features are given in the "Features.csv" file. Features that are well fitted with a Gaussian, or a sum of two Gaussians, are described with their fit parameters in "GoodFeatures.csv".

Since the raw datasets are proprietary, we are unable to give the full spectrum at each time instance for 51 Peg. We do however provide the derived radial velocities in "Derived_RVs_masked.csv" together with the time of observation (in units of Days after 5/28/2018). We also provide the radial velocities derived with the traditional CCF method in "ccf_rvs.txt". The fitted orbital parameters for each set are derived in "RV Regression.ipynb" and compared with each other.
