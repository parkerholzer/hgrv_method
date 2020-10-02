library(locfit)
library(parallel)
temp = read.csv("smoothednso_expres.csv")
r1 = which((temp$Wavelength > 5240) & (temp$Wavelength < 5245))
#r1 = which((temp$Wavelength > 5230.8) & (temp$Wavelength < 5245.8))

lowerSNR = function(flx, sn){
  orig_mean = mean(flx)
  scaled_spec = flx*(sn^2)/orig_mean
  return(rpois(length(flx), scaled_spec)*(orig_mean/(sn^2)))
}

tempmse = function(cnt, snr){
  flx = sapply(1:cnt, function(i) lowerSNR(temp$Flux[r1], sn=snr))
  #rms_raw = sqrt(mean((temp$Flux[r1] - flx[,1])^2))
  flx = c(flx)
  wvl = c(sapply(10*sin(runif(cnt, max = 2*pi)), 
                 function(v) temp$Wavelength[r1]*(1 - v/299792458)))
  
  amin = length(which(wvl <= wvl[1] + 0.017))/length(wvl)
  amax = length(which(wvl <= wvl[1] + 0.15))/length(wvl)
  alphas = seq(amin, amax, length.out = 20)
  gcvs = gcvplot(flx ~ wvl, deg=2, alpha=alphas, kern='gauss')
  #plot(gcvs$alpha, gcvs$values)
  bestalpha = gcvs$alpha[which.min(gcvs$values)]
  mdl = locfit(flx ~ wvl, deg=2, alpha=bestalpha, kern='gauss')
  
  #w = which(temp$Wavelength[r1] >= 5243.8 & 
  #            temp$Wavelength[r1] <= 5244.2)
  #plot(temp$Wavelength[r1][w], temp$Flux[r1][w], type = 'l')
  #lines(temp$Wavelength[r1][w], predict(mdl, temp$Wavelength[r1][w]), lty=2, col=2)
  #plot(temp$Wavelength[r1][w], temp$Flux[r1][w] - predict(mdl, temp$Wavelength[r1][w]))
  
  rms = sqrt(mean((temp$Flux[r1] - predict(mdl, temp$Wavelength[r1]))^2))
  return(rms)
}


df <- df2 <- as.data.frame(matrix(NA, nrow=16, ncol=16))
rownames(df) <- rownames(df2) <- as.character(2*c(1:16) - 1)
names(df) <- names(df2) <- as.character(10*c(1:16) + 90)

smry = function(x){
   return(list(mean(x), sd(x)))
}

for(i in 1:(dim(df)[1])){
  output = mclapply(as.integer(names(df)), function(snr) smry(sapply(1:50, function(k) tempmse(as.integer(rownames(df)[i]), snr))), mc.cores = 16)
  for(j in 1:(dim(df)[2])){
    df[i,j] = output[[j]][[1]]
    df2[i,j] = output[[j]][[2]]
  }
}

write.csv(df, file="hgrv_temp_est_rmsmean.csv")
write.csv(df2, file="hgrv_temp_est_rmsstd.csv")

##########################################################
#Does the RMS change much when we consider only intervals of
#absorption features?
library(rvmethod)
ftrs = findabsorptionfeatures(temp$Wavelength[r1], temp$Flux[r1], 
                              minlinedepth = 0.015)
plot(temp$Wavelength[r1], temp$Flux[r1], col=2, type='l')
for(i in 1:(length(ftrs$wvbounds))){
  lines(ftrs$wvbounds[[i]], c(1,1))}

r2 = unique(unlist(sapply(1:(length(ftrs$wvbounds)), 
       function(j) which(temp$Wavelength[r1] >= ftrs$wvbounds[[j]][1] &
                           temp$Wavelength[r1] <= ftrs$wvbounds[[j]][2]))))
plot(temp$Wavelength[r1], temp$Flux[r1])
points(temp$Wavelength[r1][r2], temp$Flux[r1][r2], col=2)
tempmse2 = function(cnt, snr){
  flx = sapply(1:cnt, function(i) lowerSNR(temp$Flux[r1], sn=snr))
  #rms_raw = sqrt(mean((temp$Flux[r1] - flx[,1])^2))
  flx = c(flx)
  wvl = c(sapply(10*sin(runif(cnt, max = 2*pi)), 
                 function(v) temp$Wavelength[r1]*(1 - v/299792458)))
  
  amin = length(which(wvl <= wvl[1] + 0.017))/length(wvl)
  amax = length(which(wvl <= wvl[1] + 0.15))/length(wvl)
  alphas = seq(amin, amax, length.out = 20)
  gcvs = gcvplot(flx ~ wvl, deg=2, alpha=alphas, kern='gauss')
  #plot(gcvs$alpha, gcvs$values)
  bestalpha = gcvs$alpha[which.min(gcvs$values)]
  mdl = locfit(flx ~ wvl, deg=2, alpha=bestalpha, kern='gauss')
  
  #w = which(temp$Wavelength[r1] >= 5243.8 & 
  #            temp$Wavelength[r1] <= 5244.2)
  #plot(temp$Wavelength[r1][w], temp$Flux[r1][w], type = 'l')
  #lines(temp$Wavelength[r1][w], predict(mdl, temp$Wavelength[r1][w]), lty=2, col=2)
  #plot(temp$Wavelength[r1][w], temp$Flux[r1][w] - predict(mdl, temp$Wavelength[r1][w]))
  
  rms = sqrt(mean((temp$Flux[r1][r2] - predict(mdl, temp$Wavelength[r1][r2]))^2))
  return(rms)
}

snr=200
nspec = 20

smry = function(x){
  return(list(mean(x), sd(x)))
}
smry(unlist(mclapply(1:50, function(k) tempmse(nspec, snr), 
                     mc.cores=8)))
smry(unlist(mclapply(1:50, function(k) tempmse2(nspec, snr), 
                     mc.cores=8)))

#################################################################
#How does the nonparametric smoothing and cadence affect the
#estimated RV?

temphat = function(cnt, snr){
  rvs = 10*sin(runif(cnt, max = 2*pi))
  rvs = rvs - mean(rvs)
  flx = sapply(1:cnt, function(i) lowerSNR(temp$Flux[r1], sn=snr))
  #rms_raw = sqrt(mean((temp$Flux[r1] - flx[,1])^2))
  flx = c(flx)
  wvl = c(sapply(rvs, 
                 function(v) temp$Wavelength[r1]*(1 - v/299792458)))
  
  amin = length(which(wvl <= wvl[1] + 0.017))/length(wvl)
  amax = length(which(wvl <= wvl[1] + 0.15))/length(wvl)
  alphas = seq(amin, amax, length.out = 20)
  gcvs = gcvplot(flx ~ wvl, deg=2, alpha=alphas, kern='gauss')
  #plot(gcvs$alpha, gcvs$values)
  bestalpha = gcvs$alpha[which.min(gcvs$values)]
  mdl = locfit(flx ~ wvl, deg=2, alpha=bestalpha, kern='gauss')
  
  
  #w = which(temp$Wavelength[r1] >= 5243.8 & 
  #            temp$Wavelength[r1] <= 5244.2)
  #plot(temp$Wavelength[r1][w], temp$Flux[r1][w], type = 'l')
  #lines(temp$Wavelength[r1][w], predict(mdl, temp$Wavelength[r1][w]), lty=2, col=2)
  #plot(temp$Wavelength[r1][w], temp$Flux[r1][w] - predict(mdl, temp$Wavelength[r1][w]))
  
  #rms = sqrt(mean((temp$Flux[r1] - predict(mdl, temp$Wavelength[r1]))^2))
  return(temp$Flux[r1] - predict(mdl, temp$Wavelength[r1]))
}



snr=80000
nspec=20
estspec = as.data.frame(mclapply(1:50, function(k) temphat(nspec, snr), 
                   mc.cores=8))

hgvar = rep(0, length(r1))
gapp = Gaussfit(temp$Wavelength[r1], temp$Flux[r1], ftrs=ftrs)
for(j in 1:(dim(gapp$parameters)[1])){
  keep = which(temp$Wavelength[r1] >= gapp$parameters$Wv_lbounds[j]-1 &
                 temp$Wavelength[r1] <= gapp$parameters$Wv_ubounds[j]+1)
  coef = sqrt(sqrt(pi))*gapp$parameters$Gauss_amp[j]*gapp$parameters$Gauss_mu[j]/(299792458*sqrt(2*gapp$parameters$Gauss_sig[j]))
  hgvar[keep] = hgvar[keep] + coef*HG1(temp$Wavelength[r1][keep], gapp$parameters$Gauss_mu[j], gapp$parameters$Gauss_sig[j])
}
plot(temp$Wavelength[r1], hgvar)


w = which(temp$Wavelength[r1] > 5243.7 & temp$Wavelength[r1] < 5244.2)
plot(temp$Wavelength[r1][w], temp$Flux[r1][w])
par(new=TRUE)
plot(temp$Wavelength[r1][w], apply(as.matrix(estspec), 1, mean)[w], col=2,
     axes=FALSE, xlab = "", ylab = "")

rv1 = rv2 = rv3 = rv4 = rep(NA, 50)
obsrv = 10
for(i in 1:50){
  flx = lowerSNR(temp$Flux[r1], sn=300)
  wvl = temp$Wavelength[r1]*(1 - obsrv/299792458)
  diff = flx - spline(temp$Wavelength[r1], 
                      temp$Flux[r1] - as.numeric(estspec[,i]), 
                      method='fmm', xout=wvl)$y
  smr = summary(lm(diff ~ 0 + hgvar))
  rv1[i] = smr$coefficients[1,1]
  diff2 = flx - spline(temp$Wavelength[r1], temp$Flux[r1], 
                       method='fmm', xout=wvl)$y
  smr = summary(lm(diff2 ~ 0 + hgvar))
  rv2[i] = smr$coefficients[1,1]
  wvl = temp$Wavelength[r1]*(1 + obsrv/299792458)
  diff = flx - spline(temp$Wavelength[r1], 
                      temp$Flux[r1] - as.numeric(estspec[,i]), 
                      method='fmm', xout=wvl)$y
  smr = summary(lm(diff ~ 0 + hgvar))
  rv3[i] = smr$coefficients[1,1]
  diff2 = flx - spline(temp$Wavelength[r1], temp$Flux[r1], 
                       method='fmm', xout=wvl)$y
  smr = summary(lm(diff2 ~ 0 + hgvar))
  rv4[i] = smr$coefficients[1,1]
}

#For each cadence, is the difference in estimated RV approximately
#20? (Note: 10 - (-10) = 20) In other words, does the cadence and 
#nonparametric smoothing affect estimates of the relative RV?
hist(rv1 - rv3 - 20, main = "Estimated Template")
hist(rv2 - rv4 - 20, main = "True Template")

#standard deviation across different cadences
sd(rv1, na.rm=T)
sd(rv2, na.rm=T)




