library(copula)
library(sn)
library(ks)
library(quantmod)
library(MASS)
library(fitdistrplus)
rm(list=ls())
getSymbols("NVDA",src='yahoo',auto.assign = TRUE,from='2013-09-29',to='2023-09-29')
getSymbols("AMD",src='yahoo',auto.assign = TRUE,from='2013-09-29',to='2023-09-29')
nvd_r=monthlyReturn(NVDA$NVDA.Adjusted,type='log')
amd_r=monthlyReturn(AMD$AMD.Adjusted,type='log')
nvd_r=coredata(nvd_r)
amd_r=coredata(amd_r)
n=length(nvd_r)
par(mfrow=c(2,1))
acf(amd_r,main='AMD Monthly Log Returns')
acf(nvd_r,main='NVIDIA Monthly Log Returns ')
##Fit The following Distributions: Normal, Skewed normal,T,skewed T

fit_normaln=fitdistr(nvd_r,densfun='normal')
fit_normala=fitdistr(amd_r,densfun='normal')


fit_tn=fitdistr(nvd_r,densfun='t',upper=30)
fit_ta=fitdistr(amd_r,densfun='t',upper=50)

fit_skt_n=st.mple(y=nvd_r)
fit_skt_a=st.mple(y=amd_r)


fit_skn_n=sn.mple(y=nvd_r)
fit_skn_a=sn.mple(y=amd_r)
fit_ta

##Check AIC For Each of them 

-2*fit_normaln$loglik+2*(2)
-2*fit_normala$loglik+2*(2)

-2*fit_tn$loglik+2*(3)
-2*fit_ta$loglik+2*(3)

-2*fit_skn_n$logL+2*(3)
-2*fit_skn_a$logL+2*(3)

-2*fit_skt_n$logL+2*(4)
-2*fit_skt_a$logL+2*(4)

## Skewed Normal & Normal
est1=fit_normala$estimate
est2=fit_skn_n$cp

u1=pnorm(amd_r,est1)
u2=psn(nvd_r,xi=est2['mean'],omega=est2['s.d.'],alpha=est2['gamma1'])
##U hat Combinations
uhat=cbind(u1,u2)

###Copula
Ct=fitCopula(copula=tCopula(dim=2),data=uhat,method="ml")## T Copula
Ct@estimate
Ct_logL=loglikCopula(param=Ct@estimate,u=uhat,copula=tCopula(dim=2));#compute loglikelihood function
-2*Ct_logL+2*length(Ct@estimate)#AIC
Cfr=fitCopula(copula=frankCopula(dim=2),data=uhat,method="ml")#fit frank copula
Cfr@estimate;
Cfr_logL=loglikCopula(param=Cfr@estimate,u=uhat,copula=frankCopula(dim=2));
-2*Cfr_logL+2*length(Cfr@estimate);
Ccl=fitCopula(copula=claytonCopula(dim=2),data=uhat,method="ml");#fit clayton copula
Ccl@estimate
Ccl_logL=loglikCopula(param=Ccl@estimate,u=uhat,copula=claytonCopula(dim=2));
-2*Ccl_logL+2*length(Ccl@estimate);#compute AIC
Cg=fitCopula(copula=normalCopula(dim=2),data=uhat,method="ml");#fit Gaussian copula
Cg_logL=loglikCopula(param=Cg@estimate,u=uhat,copula=normalCopula(dim=2))
-2*Cg_logL+2*length(Cg@estimate);#compute AIC
##Simulation For Frank Copula
n=length(amd_r)
Simu_U=rCopula(n,frankCopula(dim=2,Cfr@estimate))
summary(Simu_U)
plot(Simu_U[,1],Simu_U[,2],main='Simulated Data-Uniform Distribution',xlab='Simulated AMD',ylab='Simulated NVIDIA')
Simu_X1=qnorm(Simu_U[,1],est1)
Simu_X2=qsn(Simu_U[,2],dp=est2)
par(mfrow=c(1,1))
plot(Simu_X1,Simu_X2,main='Simulated Plot of AMD & NVIDIA',xlab='Simulated_AMD',ylab='Simulated_NVIDIA')
plot(coredata(amd_r),coredata(nvd_r),main='Actual Plot of AMD & NVIDIA',xlab="AMD",ylab='NVIDIA')

