library(parallel)
temp = read.csv("smoothednso_expres.csv")

hermitegauss1 = function(x, mu, sig){
  z = (x - mu)/sig
  return(2*z*exp(-(z^2)/2)/sqrt(sig*2*sqrt(pi)))
}

lowerSNR = function(flx, sn){
  orig_mean = mean(flx)
  scaled_spec = flx*(sn^2)/orig_mean
  return(rpois(length(flx), scaled_spec)*(orig_mean/(sn^2)))
}

ftrs = read.csv("GoodFeatures.csv")
hgrv = rep(0, dim(temp)[1])
keep2 = c()
for(i in 1:(dim(ftrs)[1])){
  keep = which((temp$Wavelength >= ftrs$Gauss_mu[i] - 1) & 
                 (temp$Wavelength <= ftrs$Gauss_mu[i] + 1))
  keep2 = c(keep2, which((temp$Wavelength >= ftrs$Wv_lbounds[i]) & 
                  (temp$Wavelength <= ftrs$Wv_ubounds[i])))
  coef = -sqrt(sqrt(pi))*ftrs$Gauss_mu[i]*ftrs$Gauss_amp[i]/(299792458*sqrt(2*ftrs$Gauss_sig[i]))
  hgrv[keep] = hgrv[keep] + coef*hermitegauss1(temp$Wavelength[keep], ftrs$Gauss_mu[i], ftrs$Gauss_sig[i])
}
keep2 = unique(keep2)


estRV = function(ty, snr){
  flx = lowerSNR(temp$Flux[keep2], snr)
  return(as.numeric(lm(ty - flx ~ 0 + hgrv[keep2], weights = 1/ty)$coefficients))
}

p = seq(-2,1,length.out = 16)
s = seq(100,250,length.out = 16)
MSE <- BIAS <- VAR <- as.data.frame(matrix(NA, nrow = 16, ncol = 16))
names(MSE) <- names(BIAS) <- names(VAR) <- as.character(p)
row.names(MSE) <- row.names(BIAS) <- row.names(VAR) <- as.character(s)

for(i in 1:16){
  snr = s[i]
  for(j in 1:16){
    rv = 10^(p[j])
    wvl = temp$Wavelength[keep2]*(1 - rv/299792458)
    tflx = spline(temp$Wavelength, temp$Flux, xout = wvl)
    vhat = unlist(mclapply(1:2000, function(k) estRV(tflx$y, snr), mc.cores=20))
    MSE[i,j] = mean((vhat - rv)^2)
    BIAS[i,j] = mean(vhat) - rv
    VAR[i,j] = mean((vhat - mean(vhat))^2)
  }
}

write.csv(MSE, file = "hgrv_mse.csv")
write.csv(BIAS, file = "hgrv_bias.csv")
write.csv(VAR, file = "hgrv_var.csv")
