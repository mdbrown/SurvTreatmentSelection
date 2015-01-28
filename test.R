library(survival)

source("functions.R")

mydat <- SIM.data.singleMarker(nn=500,
                               b.t = log(1),
                               b.y = log(2), 
                               b.yt = log(1.5), 
                               cens.lam = .05, time.max = 8)
#table(mydat$di)/250


coxfit <- coxph(Surv(xi, di) ~A*Y, data = mydat, ties = "breslow")


mydat.a0 <- mydat; mydat.a0$A = 0;
mydat.a1 <- mydat; mydat.a1$A = 1; 

sfit.a0 <- survfit(coxfit, newdata = mydat.a0)
sfit.a1 <- survfit(coxfit, newdata = mydat.a1)
sfit <- survfit(coxfit, newdata = mydat)

for(t0 in 1:8){
#t0 = 2

risk.a0 <- 1-summary(sfit.a0, times = t0)$surv
risk.a1 <- 1-summary(sfit.a1, times = t0)$surv

#pred.curves 
F <- ecdf(mydat$Y)(mydat$Y)
ooo <- order(F)

plot(F[ooo], risk.a0[ooo], type = "l", lty = 2, xlim = c(0,1), ylim = c(0,1))
lines(F[ooo], risk.a1[ooo])
title(t0)
}


summary(coxfit)

### calculate summary measures 
#five year prediction time 
t0 = 5

risk.a0 <- 1-summary(sfit.a0, times = t0)$surv
risk.a1 <- 1-summary(sfit.a1, times = t0)$surv
risk <- 1-summary(sfit, times = t0)$surv

####
# prevalence 
###

#censoring weights 
w.cens <- get.censoring.weights(ti = t0, stime = mydat$xi, status = mydat$di)


#prevalence

#mod
rho.mod <- mean(risk); rho.mod
#emp
rho.emp <- sum(w.cens*I(mydat$xi<t0)*mydat$di)/sum(w.cens); rho.emp



