##function to calculate measures

surv.trtsel <- function(time, event, marker, trt, predict.time){ 
  
  coxfit <- coxph(Surv(time, event)~marker*trt)
  
  
  sfit.a0 <- survfit(coxfit, newdata = data.frame('time' = time, 'event' = event, 'marker' = marker, 'trt' = 0))
  sfit.a1 <- survfit(coxfit, newdata = data.frame('time' = time, 'event' = event, 'marker' = marker, 'trt' = 1))

  risk.a0 <- c( 1-summary(sfit.a0, times = predict.time)$surv)
  risk.a1 <- c(1-summary(sfit.a1, times = predict.time)$surv)

  w.cens <- get.censoring.weights(ti = predict.time, 
                                  stime = time, 
                                  status = event)
  
  out <- list("derived.data" = data.frame('time' = time, 
                                   'event' = event, 
                                   'marker' = marker, 
                                   'trt' = trt, 
                                   'risk.trt' = risk.a0, 
                                   'risk.no.trt' = risk.a1, 
                                   'trt.effect' = risk.a0 - risk.a1, 
                                   'marker.neg' = risk.a0 <= risk.a1, 
                                   'w.cens' = w.cens), 
       "prediction.time" = predict.time)
  class(out) <- "surv.trtsel"
  out
  
}



surv.trtsel.mult <- function(coxfit, data, predict.time){ 
  
  #coxfit <- coxph(Surv(time, event)~marker*trt)
  data.a0 <- data
  data.a0$trt <- 0
  data.a1 <- data
  data.a1$trt <- 1
  

  sfit.a0 <- survfit(coxfit, newdata = data.a0)
  sfit.a1 <- survfit(coxfit, newdata = data.a1)
  
  risk.a0 <- c( 1-summary(sfit.a0, times = predict.time)$surv)
  risk.a1 <- c(1-summary(sfit.a1, times = predict.time)$surv)
  
  w.cens <- get.censoring.weights(ti = predict.time, 
                                  stime = data$time, 
                                  status = data$event)
  marker.names <- names(coxfit$coefficients)[is.element(names(coxfit$coefficients), names(data))]
  out <- list("derived.data" = data.frame('marker' = data[,marker.names],
                                          'time' = data$time, 
                                          'event' = data$event,
                                          'trt' = data$trt, 
                                          'risk.trt' = risk.a0, 
                                          'risk.no.trt' = risk.a1, 
                                          'trt.effect' = risk.a0 - risk.a1, 
                                          'marker.neg' = risk.a0 <= risk.a1, 
                                          'w.cens' = w.cens), 
              "prediction.time" = predict.time)
  
  class(out) <- "surv.trtsel"
  out
  
}


## evaluation method for surv.trtsel
evaluate<- function(x, ...){
  UseMethod('evaluate')
}

evaluate.surv.trtsel <- function(x, ...){

  t0 <-  x$prediction.time
  
  
  
  #prevalence
  ###########
  attach(x$derived.data)
  #mod
  rho.mod <- mean(risk.no.trt*(1-trt) + risk.trt*trt);
  #emp
  rho.emp <- sum(w.cens*I(time < t0))/sum(w.cens);
  
  #proportion marker negative
  ###########################
  
  p.marker.neg <- mean(marker.neg)
  
  #Average benefit of no treatment among marker negatives
  #Bneg
  ########################################################
  #mod
  
  Bneg.mod <- ifelse(sum(marker.neg) >0, -sum(trt.effect[marker.neg]*w.cens[marker.neg])/sum(marker.neg*w.cens), 0)
  if(sum(marker.neg)==0) print("Bneg.mod is set to zero")
  
  #emp
  num <- (w.cens*I(time < t0)*I(marker.neg)); den <- (w.cens*I(marker.neg))
  
  if(sum(den[trt==0])==0 | sum(den[trt==1]) ==0) {
    Bneg.emp <- 0
    print("Bneg.emp is set to zero")
  }else{
    Bneg.emp <- sum(num[trt==1])/sum(den[trt==1]) - sum(num[trt==0])/sum(den[trt==0])
  }
  
  #Average benefit of treatment among marker positives
  #Bneg
  ########################################################
  #mod
  marker.pos <- !marker.neg
  Bpos.mod <- ifelse(sum(marker.pos) >0, sum(trt.effect[marker.pos]*w.cens[marker.pos])/sum(marker.pos*w.cens), 0)
  if(sum(marker.pos)==0) print("Bpos.mod is set to zero")
  
  #emp

  num <- (w.cens*I(time < t0)*marker.pos); den <- (w.cens*marker.pos)
  if(sum(den[trt==0])==0 | sum(den[trt==1]) ==0) {
    Bpos.emp <- 0
    print("Bpos.emp is set to zero")
  }else{
    Bpos.emp <- sum(num[trt==0])/sum(den[trt==0]) - sum(num[trt==1])/sum(den[trt==1])
  }
  detach(x$derived.data)
  
  #Decrease in event rate under marker based treatment
  #theta
  ####################################################
  Theta.mod <- Bneg.mod*p.marker.neg; 
  Theta.emp <- Bneg.emp*p.marker.neg;
  
  out <- data.frame('theta' = c(Theta.mod, Theta.emp), 
                    'B.neg' = c(Bneg.mod,  Bneg.emp), 
                    'B.pos' = c(Bpos.mod, Bpos.emp), 
                    'Pmneg' = c(p.marker.neg, NA), 
                    'rho'   = c(rho.mod, rho.emp), 
                    'estimate'= factor(c('model-based', 'empirical')))
  out
  
}

## get summary measures 

## plot method for 'surv.trtsel'


plot.surv.trtsel <- function(x, plot.type = c("risk", "treatment effect"), ...){
  require(ggplot2)
  plot.type = match.arg(plot.type)
  
  wi <- x$derived.data$w.cens
  
  if(plot.type == "risk"){

    n <- nrow(x$derived.data)
    
    F.x <- ecdf(x$derived.data$marker)(x$derived.data$marker[wi > 0])
    
    dat <- data.frame('risk' = c(x$derived.data$risk.trt[wi > 0], x$derived.data$risk.no.trt[wi > 0]),
                      'treatment' = factor(rep(c('Treatment', 'No Treatment'), c(sum(wi >0), sum(wi > 0)))), 
                      'F.x' = rep(F.x, 2))
    p <- ggplot(dat, aes(F.x, risk, linetype = treatment, color = treatment)) + 
      geom_step(size = 1.1) + theme(legend.title = element_blank(), text = element_text(size = 18)) + 
      xlab("% below marker value") + ylab("risk given marker")
  
  }else if(plot.type == "treatment effect"){
    
    F.delta <- ecdf(x$derived.data$trt.effect)(x$derived.data$trt.effect[wi>0])
    
    dat <- data.frame('trt.effect' = x$derived.data$trt.effect[wi > 0], 'F.delta' = F.delta)

    p <- ggplot(dat, aes(F.delta, trt.effect))+ 
      geom_hline(yintercept = 0, size = .5, linetype = 1, colour = "grey35")  + 
      geom_step(size = 1.1) + 
      xlab("% below treatment effect") + ylab("treatment effect")  +
      geom_hline(yintercept = sum(dat$trt.effect*wi[wi>0])/sum(wi), size = .5, linetype = 3)
    
  }
  print(p)
  invisible(list('dat' = dat, 'p' = p))
}



##function to simulate survival data

SIM.data.singleMarker <- 
  function(nn,
           mu = 0, 
           Sigma = 1, 
           b.a = log(1),
           b.y = log(2), 
           b.ya = log(1.25),
           lam0 = .1, 
           cens.lam = 0, 
           time.max = NULL)
  {
   
    Y <- rnorm(nn, mu, Sigma)
    T <- rbinom(nn, size = 1,prob = .5)
    mu.i <- Y*b.y + b.t*T + b.yt*T*Y  
    
    #true survival time
    r.ti <- log(-log(runif(nn)))
    ti <-  -mu.i + r.ti
    ti <- exp(ti)/lam0
    
    
    
    #time.max is the followup time. 
    ci = rep(time.max, nn)
    
    
    if(cens.lam > 0){
      ci = rexp(nn, rate = cens.lam)     
    }
    
    
    ci = pmin(ci, time.max)
    
    
    #observed marker is the min of ti and ci        
    xi <- pmin(ti, ci)
    # failure indicator
    di <- ifelse( ti == xi, 1, 0)
    
    #xi is the min of ti and ci
    #di is the indicator for failure, 1 for failure, 0 for censored
    #Y is the marker values
    #browser()
    result <- as.data.frame(cbind(xi, di, T, Y)) 
    names(result) = c( "xi", "di", "A",  "Y")
    return(result)
  }


##################
## Calculate 'true' parameter values
##################

know.the.truth <- 
  function(f_y = function(y){dnorm(y, mean = 0, sd = 1)}, 
           b.a = log(1), #treatment main effect
           b.y = log(2), #marker main effect
           b.ya = log(1.25), #marker treatment interaction
           t0 = 5,  #prediction time
           lam0 = .1)
  {
    
    P.neg <- integrate(f_y, lower = -b.a/b.ya, upper = Inf)$value
    
    #S(t | y, a = 1)*f(y) 
    sy.1 <- function(y) exp(-lam0*t0*exp(b.a + b.y*y+b.ya*y))*f_y(y)
    #S(t | y, a = 0)*f(y)
    sy.0 <- function(y) exp(-lam0*t0*exp(b.y*y))*f_y(y)
    
    #multiply both integrals by .5 = Pr(A = 0) = Pr(A = 1)
    rho <- 1 - (0.5*integrate(sy.1, lower = -Inf, upper = Inf)$value + 
                0.5*integrate(sy.0, lower = -Inf, upper = Inf)$value)
    
   
    
    B.neg <- (integrate(sy.0, lower = -b.a/b.ya, upper = Inf)$value -  integrate(sy.1, lower = -b.a/b.ya, upper = Inf)$value)/P.neg 
    B.pos <- (integrate(sy.1, lower = -Inf, upper = -b.a/b.ya)$value -  integrate(sy.0, lower =-Inf, upper = -b.a/b.ya)$value)/(1-P.neg) 

    Theta <- B.neg*P.neg
    return(list("P.neg" = P.neg, "rho" = rho, "B.neg" = B.neg, "B.pos" = B.pos, "Theta" = Theta))
  }

##################
## subroutines
##################


get.censoring.weights <- function(ti, stime, status){

  ix.ctrl <- stime >= ti
  ix.cens <- stime < ti & status == 0 
  new.times.case.ctrl <- stime
  

    # this saves some time. survival curve doesn't need to be estimated past ti.
    ix <- stime > (ti + 0.1) # the magic number is necessary, otherwise the KM estimate at dc$ti is off.
    stime[ix] <- (ti + 0.1)
    
    cc  <- survfit(Surv(stime, status == 0) ~ 1,
                   se.fit = FALSE, type = 'kaplan-meier')
    new.times.case.ctrl[ix.ctrl] <- ti
    
  new.times.case.ctrl.sorted.incr <- sort(new.times.case.ctrl, decreasing = FALSE, method = 'shell')
  recover.original.order <- rank(new.times.case.ctrl)
  
  cens.weights <- summary(cc, times = new.times.case.ctrl.sorted.incr)$surv[recover.original.order]
  cens.weights <- 1/cens.weights
  cens.weights[ix.cens] <- 0 
  return(cens.weights)
}
