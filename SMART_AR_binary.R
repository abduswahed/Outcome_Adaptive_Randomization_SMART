############################################################################################
### Apply SMART-AR design for binary outcome based on our two-stage SMART design setting ###
###              Three possible treatments in each stage (A1,A2,A3)                      ###
############################################################################################


# This code is the modified version of SMART-AR design (Cheung et al., 2015) for binary outcome  (Section 5.3.4 and Web Appendix G);
# This code is modified based on the R code provided in the supplementary file of Cheung et al., 2015;
# Reproduce simulation under scenario S8, n=300, burn-in proportion=0.5, and Q-function is correctly specified in Table 2 in the manuscript.


##############################################
################ Functions  ##################
##############################################

# input a1, resp, a2: numeric
# return Q2=prob(y2=1|y1,a1,a2) by transferring logistic regression
Q2 = function(a1,resp,a2,b0,b1,b2,b3) {
  if (resp==0) {
    reg = (1-resp)*(b0 + b1*a1 + b2*a2 + b3*a1*a2)
    return(1/(exp(-reg)+1))
  } else if (resp==1) {
    return(1)
  }
}

# incorrectly specified model
Q2_ic = function(a1,resp,a2,b0,b1) {
  if (resp==0) {
    reg = (1-resp)*(b0 + b1*a1*a2)
    return(1/(exp(-reg)+1))
  } else if (resp==1) {
    return(1)
  }
}

# input a1, resp: numeric
Q2max = function(a1,resp,b0,b1,b2,b3) {
  Q21 = Q2(a1,resp,1,b0,b1,b2,b3)
  Q22 = Q2(a1,resp,2,b0,b1,b2,b3)
  Q23 = Q2(a1,resp,3,b0,b1,b2,b3)
  if (resp==0) {
    # for stage 1 non-responder
    if (a1==1) {
      if (Q22 >= Q23) {
        return ( c(2,Q22) )
      } else {
        return ( c(3,Q23) )
      }
    } else if (a1==2) {
      if (Q21 >= Q23) {
        return( c(1,Q21) )
      } else {
        return( c(3,Q23) )
      }
    } else {
      if (Q21 >= Q22) {
        return( c(1,Q21) )
      } else {
        return ( c(2,Q22) )
      }
    }
  } else if (resp==1) {
    # stage 1 responder => Q21=Q22=Q23=1
    return( c(NA, 1))
  }
  
}

# incorrect specified model
Q2max_ic = function(a1,resp,b0,b1) {
  Q21 = Q2_ic(a1,resp,1,b0,b1)
  Q22 = Q2_ic(a1,resp,2,b0,b1)
  Q23 = Q2_ic(a1,resp,3,b0,b1)
  if (resp==0) {
    # for stage 1 non-responder
    if (a1==1) {
      if (Q22 >= Q23) {
        return ( c(2,Q22) )
      } else {
        return ( c(3,Q23) )
      }
    } else if (a1==2) {
      if (Q21 >= Q23) {
        return( c(1,Q21) )
      } else {
        return( c(3,Q23) )
      }
    } else {
      if (Q21 >= Q22) {
        return( c(1,Q21) )
      } else {
        return ( c(2,Q22) )
      }
    }
  } else if (resp==1) {
    # stage 1 responder => Q21=Q22=Q23=1
    return( c(NA, 1))
  }
}


getRandProb = function(foo, base=2, truemd) {
  
  resp = foo[,2]   # response in stage 1
  a1 = foo[,1]   # a1
  a2 = foo[,3]   # a2
  z1 = (1-resp)  # 1-r1
  z2 = a1*(1-resp)     # a1*(1-r1)
  z3 = a2*(1-resp)     # a2*(1-r1)
  z4 = a1*a2*(1-resp)  # a1*a2*(1-r1)
  
  y = foo[,4]
  # number of patients used to adaptation
  n = length(y) 
  
  # if correctly specified the model as the true model
  if (truemd==T) {
    # no intercept
    fit = glm(y~z1+z2+z3+z4-1, family = "binomial")
    coef = fit$coefficients
    coef[which(is.na(coef))] = 0
    b0 <- coef[1]
    b1 <- coef[2]
    b2 <- coef[3]
    b3 <- coef[4]
    
    # pseudo outcome: phat <- max_{a2\inS2}Q2, yhat~bernoulli(phat)
    phat = rep(NA,n)  
    yhat = rep(NA,n)  
    for (i in 1:n) {
      # Q2max = function(a1,resp,b0,b1,b2,b3)
      phat[i] = Q2max(a1[i],resp[i],b0,b1,b2,b3)[2]
      yhat[i] = rbinom(1,1,phat[i])
    }
  } else if (truemd==F) {
    # not correctly specified the true Q2 model
    fit = glm(y~z1+z4-1, family = "binomial")
    coef = fit$coefficients
    coef[which(is.na(coef))] = 0
    b0 <- coef[1]
    b1 <- coef[2]
    
    
    # pseudo outcome: phat <- max_{a2\inS2}Q2, yhat~bernoulli(phat)
    phat = rep(NA,n)  
    yhat = rep(NA,n)  
    for (i in 1:n) {
      phat[i] = Q2max_ic(a1[i],resp[i],b0,b1)[2]
      yhat[i] = rbinom(1,1,phat[i])
    }
  }
  
  fit1 = glm(yhat~a1, family = "binomial")
  q1hat = c(mean(fit1$fitted.values[a1==1]), mean(fit1$fitted.values[a1==2]), mean(fit1$fitted.values[a1==3]))
  d1 = q1hat - min(q1hat)
  d1 = d1/(sum(fit1$residuals^2)/fit1$df.residual)
  d1 = pmin(d1, log(100000)/log(base))
  PI1 = base^d1 / sum(base^d1)
  PI2 = matrix(rep(NA, 15), nrow=3)
  colnames(PI2) = c("a1","resp","prob(a2=1)","prob(a2=2)", "prob(a2=3)")
  iter = 1
  
  if (truemd==T) {
    for (rcur in 1) {
      for (acur in c(1,2,3)) { 
        # force only non-responders in stage 1 enter stage 2
        if (acur==1) {
          q2hat = c(Q2(acur,0,2,b0,b1,b2,b3),
                    Q2(acur,0,3,b0,b1,b2,b3))
          d2 = q2hat - min(q2hat)
          d2 = d2/(sum(fit$residuals^2)/fit$df.residual)
          d2 = pmin(d2, log(100000)/log(base))
          prand2 = base^d2 / sum(base^d2)
          pd2 = c(0, prand2[1], prand2[2])
        } else if (acur==2) {
          q2hat = c(Q2(acur,0,1,b0,b1,b2,b3),
                    Q2(acur,0,3,b0,b1,b2,b3))
          d2 = q2hat - min(q2hat)
          d2 = d2/(sum(fit$residuals^2)/fit$df.residual)
          d2 = pmin(d2, log(100000)/log(base))
          prand2 = base^d2 / sum(base^d2)
          pd2 = c(prand2[1],0, prand2[2])
        } else if (acur==3) {
          q2hat = c(Q2(acur,0,1,b0,b1,b2,b3),
                    Q2(acur,0,2,b0,b1,b2,b3))
          d2 = q2hat - min(q2hat)
          d2 = d2/(sum(fit$residuals^2)/fit$df.residual)
          d2 = pmin(d2, log(100000)/log(base))
          prand2 = base^d2 / sum(base^d2)
          pd2 = c(prand2[1], prand2[2],0)
        }
        
        PI2[iter,] = c(acur, rcur, pd2)
        iter = iter + 1
        
      }
    }
    
  } else if (truemd==F) {
    for (rcur in 1) {
      for (acur in c(1,2,3)) { 
        # force only non-responders in stage 1 enter stage 2
        if (acur==1) {
          q2hat = c(Q2_ic(acur,0,2,b0,b1),
                    Q2_ic(acur,0,3,b0,b1))
          d2 = q2hat - min(q2hat)
          d2 = d2/(sum(fit$residuals^2)/fit$df.residual)
          d2 = pmin(d2, log(100000)/log(base))
          prand2 = base^d2 / sum(base^d2)
          pd2 = c(0, prand2[1], prand2[2])
        } else if (acur==2) {
          q2hat = c(Q2_ic(acur,0,1,b0,b1),
                    Q2_ic(acur,0,3,b0,b1))
          d2 = q2hat - min(q2hat)
          d2 = d2/(sum(fit$residuals^2)/fit$df.residual)
          d2 = pmin(d2, log(100000)/log(base))
          prand2 = base^d2 / sum(base^d2)
          pd2 = c(prand2[1],0, prand2[2])
        } else if (acur==3) {
          q2hat = c(Q2_ic(acur,0,1,b0,b1),
                    Q2_ic(acur,0,2,b0,b1))
          d2 = q2hat - min(q2hat)
          d2 = d2/(sum(fit$residuals^2)/fit$df.residual)
          d2 = pmin(d2, log(100000)/log(base))
          prand2 = base^d2 / sum(base^d2)
          pd2 = c(prand2[1], prand2[2],0)
        }
        
        PI2[iter,] = c(acur, rcur, pd2)
        iter = iter + 1
        
      }
    }
    
  }
  
  return(list(PI1=PI1, PI2=PI2))
}

learning = function(foo, truemd) {
  resp = foo[,2]   # response in stage 1
  a1 = foo[,1]   # a1
  a2 = foo[,3]   # a2
  z1 = (1-resp)  # 1-r1
  z2 = a1*(1-resp)     # a1*(1-r1)
  z3 = a2*(1-resp)     # a2*(1-r1)
  z4 = a1*a2*(1-resp)  # a1*a2*(1-r1)
  
  
  y = foo[,4]
  # number of patients used to adaptation
  n = length(y) 
  
  if (truemd==T) {
    fit = glm(y~z1+z2+z3+z4-1, family = "binomial")
    coef = fit$coefficients
    coef[which(is.na(coef))] = 0 
    b0 <- coef[1]
    b1 <- coef[2]
    b2 <- coef[3]
    b3 <- coef[4]
    
    # pseudo outcome: phat <- max_{a2\inS2}Q2, yhat~bernoulli(phat)
    phat = rep(NA,n)  
    yhat = rep(NA,n)  
    for (i in 1:n) {
      # Q2max = function(a1,resp,b0,b1,b2,b3)
      phat[i] = Q2max(a1[i],resp[i],b0,b1,b2,b3)[2]
      yhat[i] = rbinom(1,1,phat[i])
    }
  } else if (truemd==F) {
    fit = glm(y~z1+z4-1, family = "binomial")
    coef = fit$coefficients
    coef[which(is.na(coef))] = 0 
    b0 <- coef[1]
    b1 <- coef[2]
    
    
    # pseudo outcome: phat <- max_{a2\inS2}Q2, yhat~bernoulli(phat)
    phat = rep(NA,n)  
    yhat = rep(NA,n)  
    for (i in 1:n) {
      phat[i] = Q2max_ic(a1[i],resp[i],b0,b1)[2]
      yhat[i] = rbinom(1,1,phat[i])
    }
  }
  
  
  fit1 = glm(yhat~a1, family = "binomial")
  q1hat = c(mean(fit1$fitted.values[a1==1]), mean(fit1$fitted.values[a1==2]), mean(fit1$fitted.values[a1==3]))  # values correspond to a1=1, a1=2, a1=3
  
  # optimal DTR
  if (truemd==T) {
    odtr = rep(NA,2)
    if (q1hat[1] >= q1hat[2] & q1hat[1] >= q1hat[3])  {
      odtr[1] = 1
      odtr[2] = Q2max(1,0, b0,b1,b2,b3)[1]
    } else if (q1hat[2] >= q1hat[1] & q1hat[2] >= q1hat[3]) {
      odtr[1] = 2
      odtr[2] = Q2max(2,0, b0,b1,b2,b3)[1]
    } else { 
      odtr[1] = 3
      odtr[2] = Q2max(3,0, b0,b1,b2,b3)[1]
    }
  } else if (truemd==F) {
    odtr = rep(NA,2)
    if (q1hat[1] >= q1hat[2] & q1hat[1] >= q1hat[3])  {
      odtr[1] = 1
      odtr[2] = Q2max_ic(1,0, b0,b1)[1]
    } else if (q1hat[2] >= q1hat[1] & q1hat[2] >= q1hat[3]) {
      odtr[1] = 2
      odtr[2] = Q2max_ic(2,0, b0,b1)[1]
    } else { 
      odtr[1] = 3
      odtr[2] = Q2max_ic(3,0, b0,b1)[1]
    }
  }
  
  
  return(list(odtr=odtr, fit=summary(fit), fit1=summary(fit1), q1hat=q1hat))
}


smart_ar <- function(nsim, n, p.burn, theta, p_R, truemd, S){
  # initial randomization probability
  ## stage 1 initial randomization probability
  pi11 = 1/3
  pi12 = 1/3
  pi13 = 1/3
  ## stage 2 initial randomization probability
  pi21 = c(0, 1/2, 1/2)
  pi22 = c(1/2, 0, 1/2)
  pi23 = c(1/2, 1/2, 0)
  
  # intermediate response (True response rate in stage 1) R|A1 (Y1^A1)~Bernoulli(p1), R|A2 (Y1^A2)~Bernoulli(p2), R|A3 (Y1^A3)~Bernoulli(p3)
  p1 = p_R[1]
  p2 = p_R[2]
  p3 = p_R[3]
  
  # true model for stage 1 non-responders and responder --> true stage 2 response rate for different A1, A2, A3
  Q2sat = function(a1,resp,a2, theta0 = theta) {
    dim(theta0) = c(4,1)
    x = cbind((1-resp[resp==0 & !is.na(resp)]), a1[resp==0 & !is.na(resp)]*(1-resp[resp==0 & !is.na(resp)]), a2[resp==0 & !is.na(resp)]*(1-resp[resp==0 & !is.na(resp)]), a1[resp==0 & !is.na(resp)]*a2[resp==0 & !is.na(resp)]*(1-resp[resp==0 & !is.na(resp)] ))
    Q2 = rep(NA,length(a1))
    Q2[resp==0 & !is.na(resp)] = 1/(1+exp(-(x %*% theta0)))
    Q2[resp==1 & !is.na(resp)]=1
    return(Q2)
  }
  
  V0 = matrix(rep(NA,3*6),nrow=6)
  V0[,1] = c(rep(1,2),rep(2,2),rep(3,2))
  V0[,2] = c(2,3,1,3,1,2)
  for (i in 1:6) {
    d1 = V0[i,1]
    d2 = V0[i,2]
    if (d1==1) w1 = p1 # stage 1 response rate for a1=1
    if (d1==2) w1 = p2 # stage 1 response rate for a1=2
    if (d1==3) w1 = p3 # stage 1 response rate for a1=3
    w0 = 1 - w1 # stage 1 non-response prob
    V0[i,3] = w1 * Q2sat(d1,1,d2) + w0 * Q2sat(d1,0,d2)
  }
  Vest = V0[,3]
  pos = which(Vest==max(Vest))
  dtrBest = V0[pos,1:2]
  
  
  colnames(V0) = c("d1","d2(0)","Value")
  V0row = rep("",nrow(V0))
  V0row[pos] = "OPTIMAL"
  pos = which(Vest==min(Vest))
  V0row[pos] = "WORST"
  rownames(V0) = V0row
  print(V0,digits=3)
  
  # sample size
  n <- n
  # burn-in sample
  Nmin = n1 = round(n*p.burn,0)
  # AR sample
  n2 = n - n1
  base = 10
  tau = 0.75
  tau = n1 * tau^{1/(base-1)}   
 
  
  # estimated optimal DTR
  DTRhat = matrix(rep(NA,nsim*2),nrow=nsim)
  # true ORR for estimated optimal DTR
  val = rep(NA,nsim)
  # final outcome (binary)
  Y = matrix(rep(NA,nsim*n),nrow=nsim)
  # number of patient treated in each DTR
  Nd = matrix(rep(NA,nsim*6), nrow=nsim)
  
  # create N Monte Carlo samples
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    ### Sequential experiment (a1, a2)
    a1 = a2 = R = y = rep(NA,n)
    
    for (i in 1:n) {
      if (i <= n1)  {   
        # AR does not begin until there are n1 complete observations
        # equal randomization
        # stage 1
        a1[i] = sample(1:3, size = 1, replace = TRUE, prob = c(1/3,1/3,1/3))
        if (a1[i]==1)  R[i] = rbinom(1,1,p1)
        else if (a1[i]==2) R[i] = rbinom(1,1,p2)
        else R[i] = rbinom(1,1,p3)
        # stage 2 (add restriction: only define for non-responders in stage 1)
        if (a1[i]==1 & R[i]==0) {
          a2[i] = sample(c(2,3), size = 1, replace = TRUE, prob = c(1/2,1/2))
        } else if (a1[i]==2 & R[i]==0) {
          a2[i] = sample(c(1,3), size = 1, replace = TRUE, prob = c(1/2,1/2))
        } 
        else if (a1[i]==3 & R[i]==0) {
          a2[i] = sample(c(1,2), size = 1, replace = TRUE, prob = c(1/2,1/2))
        } 
        else {
          a2[i] = NA
        } 
      }
      else {
        # Begin AR
        
        # define binary final outcome Y2~Bernoull(Q2sat(a1,R,a2))
        y = rbinom(n, 1, Q2sat(a1,R,a2)) 
        ycomp = y[1:(i-1)]
        a1comp = a1[1:(i-1)]
        a2comp = a2[1:(i-1)]
        Rcomp = R[1:(i-1)]
        foo = cbind(a1comp,Rcomp,a2comp,ycomp)
        
        # calculate the adaptive randomization probabilities
        pfoo = getRandProb(foo,base=base, truemd=truemd)
        
        # weight
        w0 = (tau/(i-1))^{base-1}
        w1 = 1-w0
        
        # stage 1 adaptive randomization probabilities for A1, A2, A3
        pi11hat = pfoo$PI1[1]
        pi12hat = pfoo$PI1[2]
        pi13hat = pfoo$PI1[3]
        rho11 = exp(w0*log(pi11) + w1*log(pi11hat))
        rho12 = exp(w0*log(pi12) + w1*log(pi12hat))
        rho13 = exp(w0*log(pi13) + w1*log(pi13hat))
        pi11til = rho11  / (rho11 + rho12 + rho13)	
        pi12til = rho12  / (rho11 + rho12 + rho13)	
        pi13til = rho13  / (rho11 + rho12 + rho13)	
        
        # stage 2 adaptive randomization probabilities for A1, A2, A3
        pi21hat = pfoo$PI2[,3] 
        pi22hat = pfoo$PI2[,4] 
        pi23hat = pfoo$PI2[,5] 
        
        rho21 = exp(w0*log(pi21) + w1*log(pi21hat))
        rho22 = exp(w0*log(pi22) + w1*log(pi22hat))
        rho23 = exp(w0*log(pi23) + w1*log(pi23hat))
        rho21[1]=rho22[2]=rho23[3]=0
        pi21til = rho21 / (rho21+rho22+rho23)
        pi22til = rho22 / (rho21+rho22+rho23)
        pi23til = rho23 / (rho21+rho22+rho23)
        
        # stage 1
        
        a1[i] = sample(1:3, size = 1, replace = TRUE, prob = c(pi11til,pi12til,pi13til))
        if (a1[i]==1)  {
          R[i] = rbinom(1,1,p1)
        } else if (a1[i]==2) {
          R[i] = rbinom(1,1,p2)
        } else {
          R[i] = rbinom(1,1,p3)
        } 
        # stage 2
        if (a1[i]==1 & R[i]==0) {
          a2[i] = sample(c(2,3), size = 1, replace = TRUE, prob = c(pi22til[1],pi23til[1]))
        } else if (a1[i]==2 & R[i]==0) {
          a2[i] = sample(c(1,3), size = 1, replace = TRUE, prob = c(pi21til[2],pi23til[2]))
        } else if (a1[i]==3 & R[i]==0) {
          a2[i] = sample(c(1,2), size = 1, replace = TRUE, prob = c(pi21til[3],pi22til[3]))
        } else {
          a2[i] = NA
        }
        
        
      }
    }
    
    # number of patient treated in each DTR
    d12 = sum(a1==1 & R==1) + sum(a1==1 & R==0 & a2==2)
    d13 = sum(a1==1 & R==1) + sum(a1==1 & R==0 & a2==3)
    d21 = sum(a1==2 & R==1) + sum(a1==2 & R==0 & a2==1)
    d23 = sum(a1==2 & R==1) + sum(a1==2 & R==0 & a2==3)
    d31 = sum(a1==3 & R==1) + sum(a1==3 & R==0 & a2==1)
    d32 = sum(a1==3 & R==1) + sum(a1==3 & R==0 & a2==2)
    Nd[r,] = c(d12,d13,d21,d23,d31,d32)
    
    ### Final fit after complete follow-up of all subjects
    y = rbinom(n, 1, Q2sat(a1,R,a2)) 
    foo = cbind(a1,R,a2,y)
    fit1 = learning(foo, truemd=truemd)
    
    Y[r,] = y
    # optimal DTR estimated from the model
    DTRhat[r,] = dtrhat = fit1$odtr
    pos = which(V0[,1]==dtrhat[1] & V0[,2]==dtrhat[2])
    # the true ORR for estimated optimal DTR from the model
    val[r] = V0[pos,3]
    
    # progress indicator
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  result = list(V0=V0, Y=Y, DTRhat=DTRhat, val=val,  Nd=Nd)
  saveRDS(result, file=paste0("./output/ar_",n,"_",p.burn,"_",truemd,"_",S, ".Rds"))
  
  
}

#############################################################################################
# Simulation: S8, Q-function correctly specified, sample size=300, burn-in proportion=0.5 ##
#############################################################################################
set.seed(1027)
smart_ar(nsim = 10000, n=300, p.burn = 0.5, theta = c(-1, -0.3, 0.5, -0.2), p_R=c(0.5, 0.35,0.2), truemd = T, S="S8")

# first run the above smart_ar above, then import the data produced by it
df <- readRDS("./output/ar_300_0.5_TRUE_S8.Rds")

V0 <- df$V0
Y <- df$Y
DTRhat <- df$DTRhat
Nd <- df$Nd

# true optimal and worst DTRs
v4 = V0[,3]
opt = which(v4==max(v4))
dmax = V0[opt,1:2]
vmax = V0[opt,3]
wor = which(v4==min(v4))
dmin = V0[wor,1:2]
vmin = V0[wor,3]


# power: proportion of identifying the optimal DTR, calculate under the alternative
pow <- function(x) {
  x[1]==dmax[1] & x[2]==dmax[2]
}
power_hat = mean(apply(DTRhat,1,pow))
cat("Power of identifying the true optimal DTR\n")
print(power_hat)


# type I error: calculate under null scenario
typeI <- function(x) {
  if (x[1]==1 & x[2]==2) {
    y=1
  } else if (x[1]==1 & x[2]==3) {
    y=2
  } else if (x[1]==2 & x[2]==1) {
    y=3
  } else if (x[1]==2 & x[2]==3) {
    y=4
  } else if (x[1]==3 & x[2]==1) {
    y=5
  } else if (x[1]==3 & x[2]==2) {
    y=6
  }
  return(y==sample(1:6,1,replace = T))
}
typeI_hat = mean(apply(DTRhat,1,typeI))
cat("Type I error of identifying the true optimal DTR\n")
print(typeI_hat)


cat("Distribution of the number of response per trial\n")
patientval = apply(Y,1,sum)
print(mean(patientval))

cat("Distribution of the number of patients treated in the true optimal and worst DTR\n")
nopt_fun <- function(x) {
  x[opt]
}
nwor_fun <- function(x) {
  x[wor]
}
nopt <- apply(Nd, 1, nopt_fun)
nwor <- apply(Nd, 1, nwor_fun)
print(mean(nopt))
print(mean(nwor))




