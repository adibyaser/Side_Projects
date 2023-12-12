rm(list=ls())
library(copula)
library(tseries)
library(sn)
library("Ecdat")
library("fGarch")
library("evir")
library("forecast")
library(ks)
library(quantmod)
library(MASS)
library(fitdistrplus)
getSymbols("SPXL",src='yahoo',auto.assign = TRUE,from='2013-09-29',to='2023-09-29')
getSymbols("AMD",src='yahoo',auto.assign = TRUE,from='2013-09-29',to='2023-09-29')
spxl_r=dailyReturn(SPXL$SPXL.Adjusted,type='log')
amd_r=dailyReturn(AMD$AMD.Adjusted,type='log')
spxl_r=coredata(spxl_r)
amd_r=coredata(amd_r)
plot(amd_r)
points(spxl_r)
## Preliminary Stats
plot(spxl_r)
plot(amd_r)
par(mfrow=c(1,1))
acf(spxl_r,main='SPXL Daily Returns')
acf(amd_r, main='AMD Daily Returns')
summary(amd_r)
summary(spxl_r)
##Let's Try seeing the residuals with AR(1) and validate the justification of time dependence & volatility clustering
fitAR1_s=arima(spxl_r,order=c(0,0,2),include.mean=TRUE)#fit ARIMA model
fitAR1_s
res_s=residuals(fitAR1_s)
acf(res_s,main='Residual MA(2)') ## Not white noise. Therefore, it is more important to capture the time dependencies
acf(res_s^2,main='Residual Sq MA(2)') ## Numerous spikes at variety of lags
##AR(1)+GARCH(1,1)
fit_s=garchFit(formula=~arma(1,0)+garch(1,1),data=spxl_r,cond.dist="norm",trace=FALSE)
fit_a=garchFit(formula=~arma(1,0)+garch(1,1),data=amd_r,cond.dist="norm",trace=FALSE)
fit_a
res_a=residuals(fit_a,standardize=TRUE) ##White noise
res_s=residuals(fit_s,standardize=TRUE)
fit_ta=fitdistr(res_a,'t',upper=10)
fit_ts=fitdistr(res_s,'t',upper=10)
fit_ts
par(mfrow=c(1,1))
acf(res_s,main='SPXL Residual AR(1)-GARCH (1,1)')
acf(res_s^2, main='SPXL Residual Sq AR(1)-GARCH(1,1)')
acf(res_a,main='AMD Residual AR(1)-GARCH(1,1)')
acf(res_a^2,main='AMD Residual Sq AR(1)-GARCH(1,1)')
## Copula modeling

est1=fit_ta$estimate
est2=fit_ts$estimate

u1=pt(res_a,df=est1['df'])
u2=pt(res_s,df=est2['df'])

##Combining uniform data
uhat=cbind(u1,u2)

##Copula
Ct=fitCopula(copula=tCopula(dim=2),data=uhat,method="ml")## T Copula
Ct@estimate
Ct_logL=loglikCopula(param=Ct@estimate,u=uhat,copula=tCopula(dim=2));#compute loglikelihood function
Ct_logL
-2*Ct_logL+2*length(Ct@estimate)#AIC
Cg=fitCopula(copula=normalCopula(dim=2),data=uhat,method="ml")#fit Gaussian copula
Cg@estimate;
Cg_logL=loglikCopula(param=Cg@estimate,u=uhat,copula=normalCopula(dim=2));
Cg_logL
-2*Cg_logL+2*length(Cg@estimate);
Ccl=fitCopula(copula=claytonCopula(dim=2),data=uhat,method="ml");#fit clayton copula
Ccl@estimate
Ccl_logL=loglikCopula(param=Ccl@estimate,u=uhat,copula=claytonCopula(dim=2));
Ccl_logL
-2*Ccl_logL+2*length(Ccl@estimate);#compute AIC

Cgmbl=fitCopula(copula=gumbelCopula(dim=2),data=uhat,method="ml");#fit gumbel copula
Cgmbl@estimate
Cgmbl
Cgmbl_logL=loglikCopula(param=Cgmbl@estimate,u=uhat,copula=gumbelCopula(dim=2))
Cgmbl_logL
-2*Cgmbl_logL+2*length(Cgmbl@estimate);

##So T Copula is the best for joint distribution
Ct@estimate
rho=Ct@estimate[1]
correlation_matrix=matrix(c(1,rho,rho,1),nrow=2)
correlation_matrix
confidence_level <- 0.99
n=length(amd_r)
Simu_U=rCopula(n,tCopula(dim=2,Ct@estimate[1],df=Ct@estimate[2]))
Simu_X1=qt(Simu_U[,1],df=est1['df'])
Simu_X2=qt(Simu_U[,2],df=est2['df'])
plot(Simu_X1,Simu_X2,xlab='AMD',ylab='SPXL',main='T-Copula 
     Scatter Plot-Simulated AMD & SPXL Residuals ',col='red')
points(res_a,res_s,col='blue',pch=4)
summary(Simu_U)

forecasted_s=predict(fit_s,n.ahead=1)
forecasted_a=predict(fit_a,n.ahead=1)

forecasted_volatility_a=forecasted_a$standardDeviation
forecasted_volatility_s=forecasted_s$standardDeviation

forecasted_s$meanForecast
forecasted_a$meanForecast## Scaling 
scaled_residuals_x1=Simu_X1*forecasted_volatility_a
scaled_residuals_x2=Simu_X2*forecasted_volatility_s

plot(scaled_residuals_x1,scaled_residuals_x2)

plot(spxl_r)
plot(amd_r)
weight_values=seq(0.1,0.9,by=0.1)
n_simulations=10000
VaR_estimates=numeric(length(weight_values))
VaR_estimates

for (i in seq_along(weight_values)){
  weight=weight_values[i]
  portfolio_returns=weight*(forecasted_a$meanForecast+scaled_residuals_x1)+(1-weight)*(forecasted_s$meanForecast+scaled_residuals_x2)
  VaR_estimates[i]=quantile(portfolio_returns,1-confidence_level)
}
VaR_estimates
##Expected Shortfall
confidence_level=0.975
ES_estimates=numeric(length(weight_values))
VaR_estimates=numeric(length(weight_values))
for (i in seq_along(weight_values)){
  weight=weight_values[i]
  portfolio_returns=weight*(forecasted_a$meanForecast+scaled_residuals_x1)+(1-weight)*(forecasted_s$meanForecast+scaled_residuals_x2)
  VaR_estimates[i]=quantile(portfolio_returns,1-confidence_level)
  ES_estimates[i]=mean(portfolio_returns[portfolio_returns<VaR_estimates[i]])}

ES_estimates
