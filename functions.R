##function to calculate measures






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
    names(result) = c( "xi", "di", "A",  "Y")
    return(result)
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
