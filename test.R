library(survival)

source("functions.R")

mydat <- SIM.data.singleMarker(nn=502,
                               b.t = log(1),
                               b.y = log(2), 
                               b.yt = log(1.25), 
                               cens.lam = .05, time.max = 8)

st <- surv.trtsel(time = mydat$xi, event = mydat$di, marker = mydat$Y, trt = mydat$A, 
                  predict.time= 5)

plot(st, plot.type = "risk")
plot(st, plot.type = "treatment effect") 

evaluate(st)


#small.simulation to check if estimates are consistent
set.seed(5)
S = 100
res <- data.frame(matrix(nrow = S*2, ncol  = 6)); names(res) = c('theta', 'B.neg', 'B.pos', 'Pmneg', 'rho', 'estimate')
for(s in 1:S){
  mydat <- SIM.data.singleMarker(nn=1000,
                                 b.t = log(1),
                                 b.y = log(2), 
                                 b.yt = log(1.25), 
                                 cens.lam = .05, time.max = 8)
  
  st <- surv.trtsel(time = mydat$xi, event = mydat$di, marker = mydat$Y, trt = mydat$A, 
                    predict.time= 5)
  ind <- (s*2-1):(s*2)
  res[ind,] <-  evaluate(st)
  
 # if(s%%10 == 0 ) {
    cat( paste((s/S)*100, 'percent complete\n'))
  #}
}

long.res <- melt(res, id.vars = "estimate" )
ggplot(long.res, aes(value, colour = factor(estimate), fill = factor(estimate))) + 
  geom_density(alpha = .5) + facet_grid(.~variable, scales = "free_x")

aggregate(res, by = list(res$estimate), FUN = 'mean', na.rm = TRUE)

#Group.1      theta      B.neg      B.pos   Pmneg       rho estimate
#1       1 0.02940426 0.05715534 0.03303418     NaN 0.4268553        1
#2       2 0.02910461 0.05554046 0.03374389 0.51749 0.4268509        2




##################################################################
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

trt.effect <- risk.a0 - risk.a1

#pred.curves 
F <- ecdf(mydat$Y)(mydat$Y)

ooo <- order(F)
plot(F[ooo], trt.effect[ooo], type = "l", ylim = c(-1, 1)); abline(h=0, lty = 2)
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

trt.effect <- risk.a0 - risk.a1 
marker.neg <- trt.effect <= 0 

####
# prevalence 
###

#censoring weights 
w.cens <- get.censoring.weights(ti = t0, stime = mydat$xi, status = mydat$di)


#prevalence
###########
 #mod
rho.mod <- mean(risk); rho.mod
 #emp
rho.emp <- sum(w.cens*I(mydat$xi<t0)*mydat$di)/sum(w.cens); rho.emp

#proportion marker negative
###########################

p.marker.neg <- mean(marker.neg)

#Average benefit of no treatment among marker negatives
#Bneg
########################################################
 #mod
Bneg.mod <- ifelse(sum(marker.neg) >0, -mean(trt.effect[marker.neg]), 0)
Bneg.mod
#emp
num <- (w.cens*I(mydat$xi<t0)*mydat$di*I(marker.neg)); den <- (w.cens*I(marker.neg))
Bneg.emp <- sum(num[mydat$A==1])/sum(den[mydat$A==1]) - sum(num[mydat$A==0])/sum(den[mydat$A==0])
Bneg.emp

#Average benefit of treatment among marker positives
#Bneg
########################################################
#mod
Bpos.mod <- ifelse(sum(!marker.neg) >0, mean(trt.effect[!marker.neg]), 0)
Bpos.mod
#emp
num <- (w.cens*I(mydat$xi<t0)*I(!marker.neg)); den <- (w.cens*I(!marker.neg))
Bpos.emp <- sum(num[mydat$A==0])/sum(den[mydat$A==0]) - sum(num[mydat$A==1])/sum(den[mydat$A==1])
Bpos.emp


#Decrease in event rate under marker based treatment
#theta
####################################################
Theta.mod <- Bneg.mod*p.marker.neg; Theta.mod
Theta.emp <- Bneg.emp*p.marker.neg; Theta.emp











