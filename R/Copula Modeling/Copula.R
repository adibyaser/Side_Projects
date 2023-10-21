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
n=length(amd_r)
nvd_r=coredata(nvd_r)
amd_r=coredata(amd_r)
##Fit The following Distributions: Normal, Skewed normal,T,skewed T

fit_normaln=fitdistr(nvd_r,densfun='normal')
fit_normala=fitdistr(amd_r,densfun='normal')


fit_tn=fitdistr(nvd_r,densfun='t',lower=0)
fit_ta=fitdistr(amd_r,densfun='t',lower=0)

fit_skt_n=st.mple(rep(1,length(nvd_r)),nvd_r)
fit_skt_a=st.mple(rep(1,length(amd_r)),amd_r)


fit_skn_n=sn.mple(rep(1,length(nvd_r)),nvd_r)
fit_skn_a=sn.mple(rep(1,length(amd_r)),amd_r)

fit_skn_n
##Check AIC For Each of them 

-2*fit_normaln$loglik+2*(2)
-2*fit_normala$loglik+2*(2)

-2*fit_tn$loglik+2*(3)
-2*fit_ta$loglik+2*(3)

-2*fit_skn_n$logL+2*(3)
-2*fit_skn_a$logL+2*(3)

-2*fit_skt_n$logL+2*(4)
-2*fit_skt_a$logL+2*(4)
## Skewed Normal
est1=fit_normala$estimate
est2=fit_skt_n$cp

u1=pnorm(amd_r,est1)
u2=psn(nvd_r,dp=est2)
##U hat Combinations
uhat=cbind(u1,u2)

## Skewed N

Ct=fitCopula(copula=tCopula(dim=2),data=uhat,method="ml")## T Copula
Ct@estimate
Ct_logL=loglikCopula(param=Ct@estimate,u=uhat,copula=tCopula(dim=2));#compute loglikelihood function
Ct_logL
-2*Ct_logL+2*length(Ct@estimate)#AIC
Cfr=fitCopula(copula=frankCopula(dim=2),data=uhat,method="ml")#fit frank copula
Cfr@estimate;
Cfr_logL=loglikCopula(param=Cfr@estimate,u=uhat,copula=frankCopula(dim=2));
Cfr_logL
-2*Cfr_logL+2*length(Cfr@estimate);
Ccl=fitCopula(copula=claytonCopula(dim=2),data=uhat,method="ml");#fit clayton copula
Ccl@estimate
Ccl_logL=loglikCopula(param=Ccl@estimate,u=uhat,copula=claytonCopula(dim=2));
Ccl_logL
-2*Ccl_logL+2*length(Ccl@estimate);#compute AIC

Cgmbl=fitCopula(copula=gumbelCopula(dim=2),data=uhat,method="ml");#fit gumbel copula
Cgmbl@estimate
Cgmbl_logL=loglikCopula(param=Cgmbl@estimate,u=uhat,copula=gumbelCopula(dim=2))
Cgmbl_logL
-2*Cgmbl_logL+2*length(Cgmbl@estimate);#compute AIC
##Simulation For Gumbel Copula
Simu_U=rCopula(n,gumbelCopula(dim=2,Cgmbl@estimate))
summary(Simu_U)
plot(Simu_U[,1],Simu_U[,2])
Simu_X1=qnorm(Simu_U[,1],est1)
Simu_X2=qsn(Simu_U[,2],dp=est2)
par(mfrow=c(2,1))
plot(Simu_X1,Simu_X2,main='Simulated Plot of AMD & NVIDIA',xlab='Simu_AMD',ylab='Simu_NVIDIA')
plot(coredata(amd_r),coredata(nvd_r),main='Actual Plot of AMD & NVIDIA')

