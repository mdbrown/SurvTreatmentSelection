##function to calculate measures

surv.trtsel(cox.object, data){
  
  
  #need risk estimates, treatment effect 
  out <- data.frame(data)
}





##function to simulate survival data

SIM.data.singleMarker <- 
  function(nn,
           mu = 0, 
           Sigma = 1, 
           b.y = log(3),
           b.t = log(2), 
           b.yt = log(1.5),
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
    names(result) = c( "xi", "di", "T",  "Y")
    return(result)
  }


##################
## subroutines
##################


get.censoring.weights <- function(dc, data, weights = NULL, new.times, 
                                  ix.ctrl, time.scale = 'time'){
  
  new.times.case.ctrl <- new.times
  
  if(time.scale == 'time'){
    # this saves some time. survival curve doesn't need to be estimated past ti.
    ix <- data$time > (dc$ti + 0.1) # the magic number is necessary, otherwise the KM estimate at dc$ti is off.
    data$time[ix] <- (dc$ti + 0.1)
    
    cc  <- survfit(Surv(time, status == 0) ~ 1,
                   data = data, weights = weights, se.fit = F, type = 'kaplan-meier')
    new.times.case.ctrl[ix.ctrl] <- dc$ti
    
  }else{ # time.scale == t.star
    ix <- data$t.star > (dc$pred.time + 0.1) 
    data$t.star[ix] <- (dc$pred.time + 0.1)
    
    cc  <- survfit(Surv(t.star, status == 0) ~ 1,
                   data = data, weights = weights, se.fit = F, type = 'kaplan-meier')
    new.times.case.ctrl[ix.ctrl] <- dc$pred.time
  }
  
  new.times.case.ctrl.sorted.incr <- sort(new.times.case.ctrl, decreasing = F, method = 'shell')
  recover.original.order <- rank(new.times.case.ctrl)
  
  cens.weights <- summary(cc, times = new.times.case.ctrl.sorted.incr)$surv[recover.original.order]
  
  return(1/cens.weights)
}
