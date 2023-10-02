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
n=nrow(amd_r)
##Fit The following Distributions: Normal, Skewed normal,T,skewed T

fit_normaln=fitdistr(nvd_r,densfun='normal')
fit_normala=fitdistr(amd_r,densfun='normal')


fit_tn=fitdistr(nvd_r,densfun='t',lower=0)
fit_ta=fitdistr(amd_r,densfun='t',lower=0)

fit_skt_n=st.mple(matrix(1,n,1),y=nvd_r,dp=c(mean(nvd_r),sd(nvd_r),0,10),control=list(maxit=100))
fit_skt_a=st.mple(matrix(1,n,1),y=amd_r,dp=c(mean(amd_r),sd(amd_r),0,10),control=list(maxit=100))


fit_skn_n=sn.mple(matrix(1,n,1),y=nvd_r,cp=c(xi=mean(nvd_r),omega=sd(nvd_r),alpha=1))
fit_skn_a=sn.mple(matrix(1,n,1),y=amd_r,cp=c(xi=mean(amd_r),omega=sd(amd_r),alpha=1))

est1=fit_skt_a$dp
est2=fit_skt_n$dp
u1=numeric(length(amd_r))
for(i in seq_along(amd_r)) {
  u1[i] <- pst(amd_r[i], dp = est1)}
u1
u1=xts(u1,order.by = index(amd_r))


u2=numeric(length(nvd_r))
for(i in seq_along(nvd_r)) {
  u2[i] <- pst(nvd_r[i], dp = est2)}
u2=xts(u2,order.by = index(nvd_r))

## Transforming Marginals to Uniform Distributions 
##est1=fit_skt_a$dp
##est2=fit_skt_n$dp
##u1=pst(amd_r,dp=est1)
##u2=pst(nvd_r,dp=est2)

##2)Skewed Normal

est3=fit_skn_a$cp
est4=fit_skn_n$cp
u3=numeric(length(amd_r))
for(i in seq_along(amd_r)) {
  u3[i] <- psn(amd_r[i], dp = est3)}
u3
u3=xts(u3,order.by = index(amd_r))

u4=numeric(length(nvd_r))
for(i in seq_along(nvd_r)) {
  u4[i] <- psn(nvd_r[i], dp = est4)}
u4
u4=xts(u4,order.by = index(nvd_r))

##3Normal

est5=fit_normala$estimate
est6=fit_normaln$estimate

u5=pnorm(amd_r,est5)
u6=pnorm(nvd_r,est6)
u5
u6
##4 T
est7=fit_ta$estimate
est8=fit_tn$estimate
est7['df']

u7=pt(amd_r,df=est7['df'])
u8=pt(nvd_r,df=est8['df'])


##U hat Combinations
u1=coredata(u1)
u2=coredata(u2)
u3=coredata(u3)
u4=coredata(u4)
u5=coredata(u5)
u6=coredata(u6)
u7=coredata(u7)
u8=coredata(u8)
uhat=cbind(u1,u2)
uhat2=cbind(u3,u4)
uhat3=cbind(u5,u6)
uhat4=cbind(u7,u8)

plot(coredata(amd_r),coredata(nvd_r))


tau=as.numeric(cor.test(u1,u2,method="kendall")$estimate)
omega=sin(tau*pi/2)

### Individual Pair Analysis
## Skewed T

Ct=fitCopula(copula=tCopula(dim=2),data=uhat,method="ml")
Ct@estimate
Ct_logL=loglikCopula(param=Ct@estimate,u=uhat,copula=tCopula(dim=2));#compute loglikelihood function
Ct_logL
-2*Ct_logL+2*length(Ct@estimate)#AIC
Cfr=fitCopula(copula=frankCopula(1,dim=2),data=uhat,method="ml")#fit frank copula
Cfr@estimate;
Cfr_logL=loglikCopula(param=Cfr@estimate,u=uhat,copula=frankCopula(dim=2));
Cfr_logL
-2*Cfr_logL+2*length(Cfr@estimate);
Ccl=fitCopula(copula=claytonCopula(1,dim=2),data=uhat,method="ml");#fit clayton copula
Ccl@estimate
Ccl_logL=loglikCopula(param=Ccl@estimate,u=uhat,copula=claytonCopula(dim=2));
Ccl_logL
-2*Ccl_logL+2*length(Ccl@estimate);#compute AIC

##Simulation
Simu_U=rCopula(n,tCopula(dim=2,Ct@estimate[1],df=Ct@estimate[2]))
head(Simu_U)
summary(Simu_U)
plot(Simu_U[,1],Simu_U[,2])
Simu_X1=qst(Simu_U[,1],dp=est1)
Simu_X2=qst(Simu_U[,2],dp=est2)
plot(Simu_X1,Simu_X2)
plot(coredata(amd_r),coredata(nvd_r))

### Skewed N
Ct=fitCopula(copula=tCopula(dim=2),data=uhat2,method="ml")
Ct@estimate
Ct_logL=loglikCopula(param=Ct@estimate,u=uhat2,copula=tCopula(dim=2));#compute loglikelihood function
Ct_logL
-2*Ct_logL+2*length(Ct@estimate)#AIC
Cfr=fitCopula(copula=frankCopula(1,dim=2),data=uhat2,method="ml")#fit frank copula
Cfr@estimate;
Cfr_logL=loglikCopula(param=Cfr@estimate,u=uhat2,copula=frankCopula(dim=2));
Cfr_logL
-2*Cfr_logL+2*length(Cfr@estimate);
Ccl=fitCopula(copula=claytonCopula(1,dim=2),data=uhat2,method="ml");#fit clayton copula
Ccl@estimate
Ccl_logL=loglikCopula(param=Ccl@estimate,u=uhat2,copula=claytonCopula(dim=2));
Ccl_logL
-2*Ccl_logL+2*length(Ccl@estimate);#compute AIC
##Simulation
Simu_U=rCopula(n,tCopula(dim=2,Ct@estimate[1],df=Ct@estimate[2]))
head(Simu_U)
summary(Simu_U)
plot(Simu_U[,1],Simu_U[,2])
Simu_X1=qsn(Simu_U[,1],xi=est3[1],omega=est3[2],alpha=est3[3])
Simu_X2=qsn(Simu_U[,2],xi=est4[1],omega=est4[2],alpha=est4[3])
plot(Simu_X1,Simu_X2)
plot(coredata(amd_r),coredata(nvd_r))
omega

uhat3
##Normal
Ct=fitCopula(copula=tCopula(dim=2),data=uhat3,method="ml",start=c(omega,0))
Ct@estimate
Ct_logL=loglikCopula(param=Ct@estimate,u=uhat3,copula=tCopula(dim=2));#compute loglikelihood function
Ct_logL
-2*Ct_logL+2*length(Ct@estimate)#AIC
Cfr=fitCopula(copula=frankCopula(1,dim=2),data=uhat3,method="ml")#fit frank copula
Cfr@estimate;
Cfr_logL=loglikCopula(param=Cfr@estimate,u=uhat3,copula=frankCopula(dim=2));
Cfr_logL
-2*Cfr_logL+2*length(Cfr@estimate);
Ccl=fitCopula(copula=claytonCopula(1,dim=2),data=uhat3,method="ml");#fit clayton copula
Ccl@estimate
Ccl_logL=loglikCopula(param=Ccl@estimate,u=uhat3,copula=claytonCopula(dim=2));
Ccl_logL
-2*Ccl_logL+2*length(Ccl@estimate);#compute AIC

##Simulation
Simu_U=rCopula(n,frankCopula(dim = 2,Cfr@estimate))
head(Simu_U)
summary(Simu_U)
plot(Simu_U[,1],Simu_U[,2])
Simu_X1=qnorm(Simu_U[,1])
Simu_X2=qnorm(Simu_U[,2])
plot(Simu_X1,Simu_X2)
plot(coredata(amd_r),coredata(nvd_r))
omega

##T
tau=as.numeric(cor.test(u7,u8,method="kendall")$estimate);## How the density looks like for each of them 
omega=sin(tau*pi/2)
Ct=fitCopula(copula=tCopula(dim=2),data=uhat4,method="ml",start=c(-0.1,8))
print(Ct@estimate)
Ct_logL=loglikCopula(param=Ct@estimate,u=uhat4,copula=tCopula(dim=2));#compute loglikelihood function
Ct_logL
-2*Ct_logL+2*length(Ct@estimate)#AIC
Cfr=fitCopula(copula=frankCopula(1,dim=2),data=uhat4,method="ml")#fit frank copula
Cfr@estimate;
Cfr_logL=loglikCopula(param=Cfr@estimate,u=uhat4,copula=frankCopula(dim=2));
Cfr_logL
-2*Cfr_logL+2*length(Cfr@estimate);
Ccl=fitCopula(copula=claytonCopula(1,dim=2),data=uhat4,method="ml");#fit clayton copula
Ccl@estimate
Ccl_logL=loglikCopula(param=Ccl@estimate,u=uhat4,copula=claytonCopula(dim=2));
Ccl_logL
-2*Ccl_logL+2*length(Ccl@estimate);#compute AIC

Cgauss=fitCopula(copula=normalCopula(dim=2),data=uhat4,method="ml",start=c(omega));#fit Gaussian copula
Cgauss@estimate;
Cgauss_logL=loglikCopula(param=Cgauss@estimate,u=uhat4,copula=normalCopula(dim=2));
Cgauss_logL
-2*Cgauss_logL+2*length(Cgauss@estimate);#compute AIC


Simu_U=rCopula(n,frankCopula(dim = 2,Cfr@estimate))
head(Simu_U)
summary(Simu_U)
plot(Simu_U[,1],Simu_U[,2])
Simu_X1=qt(Simu_U[,1],df=est7[3])
Simu_X2=qt(Simu_U[,2],df=est8[3])
plot(Simu_X1,Simu_X2)
plot(coredata(amd_r),coredata(nvd_r))

