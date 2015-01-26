library(survival)

source("functions.R")

mydat <- SIM.data.singleMarker(nn=250,
                               b.t = log(4),
                               b.y = log(3), 
                               b.yt = log(1.5), 
                               cens.lam = .05, time.max = 8)
table(mydat$di)/250

fit <- survfit(Surv(xi, di)~1, data = mydat, type = "kaplan-meier")
plot(fit)
summary(fit)


coxfit <- coxph(Surv(xi, di) ~T*Y, data = mydat)
summary(coxfit)



