library(locfit)
library(parallel)
temp = read.csv("smoothednso_expres.csv")
r1 = which((temp$Wavelength > 5240.8) & (temp$Wavelength < 5245.8))

lowerSNR = function(flx, sn){
  orig_mean = mean(flx)
  scaled_spec = flx*(sn^2)/orig_mean
  return(rpois(length(flx), scaled_spec)*(orig_mean/(sn^2)))
}

tempmse = function(cnt, snr){
  flx = sapply(1:cnt, function(i) lowerSNR(temp$Flux[r1], sn=snr))
  rms_raw = sqrt(mean((temp$Flux[r1] - flx[,1])^2))
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
  
  #plot(temp$Wavelength[r1], temp$Flux[r1], type = 'l')
  #lines(temp$Wavelength[r1], predict(mdl, temp$Wavelength[r1]), lty=2, col=2)
  #plot(temp$Wavelength[r1], temp$Flux[r1] - predict(mdl, temp$Wavelength[r1]))
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

