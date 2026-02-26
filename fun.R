##########################################################################################################
# Date generation, Monte Carlo replication, Evaluation functions                                         #
# for the standard SMART design, RA-SMART design (Wang et al.; 2022), and proposed GO-SMART designs      #
# Two-stage design                                                                                       #
# Assume 3 potential treatments, T=(A1, A2, A3)                                                          #
##########################################################################################################

# This file includes 4 functions:
# 1) function: "data_gene_c": date generation for the standard SMART design, and the proposed GO-SMART designs;
# 2) function: "data_gene_ra": date generation for the RA-SMART design (Wang et al., 2022);
# 3) function: "monte": Monte Carlo replication for the designs, and estimate the overall response rate using 4 different estimators on Section 4 of the manuscript;
# 4) function: "eva": Evaluate function for operating characteristics, including bias, variance ratio, coverage probability, number of patients treated in each DTR, number of patients response in each DTR, number of patients response in the trial, type I error and power;


# Input variables for function "data_gene_c":
## n: sample size; e.g. n=500
## pi1.A1, pi1.A2, pi1.A3: stage 1 true response rate; e.g. pi1.A1=0.5; pi1.A2=0.35; pi1.A3=0.2
## pi2.A1A2, pi2.A1A3, pi2.A2A1, pi2.A2A3, pi2.A3A1, pi2.A3A2: stage 2 true response rate; e.g. pi2.A1A2=0.3; pi2.A1A3=0.4; pi2.A2A1=0.35; pi2.A2A3=0.2; pi2.A3A1=0.25; pi2.A3A2=0.1
## p0.burn: stage 1 burn-in sample proportion; e.g. p0.burn=0.25
## p1.burn: stage 2 burn-in sample proportion; e.g. p1.burn=0.5
## AR: GO-SMART AR-1 or AR-2 design; e.g. AR="AR1"
## e: randomization constraint parameter, restrict the randomization probability to the interval [e,1-e]; e.g. e=c(0.1, 0.1, 0.1, 0.1, 0.1)
## c: tuning parameter that controls the degree of dependence of the OAR probabilities on the empirical response rates; e.g. c="n.N" (means c=i/n in the manuscript notation)


data_gene_c <- function(n,pi1.A1,pi1.A2,pi1.A3,pi2.A1A2,pi2.A1A3,pi2.A2A1,pi2.A2A3,pi2.A3A1,pi2.A3A2,p0.burn,p1.burn,AR,e,c) {
  
    
    e1 <- e[1]
    e2 <- e[2]
    e3 <- e[3]
    e4 <- e[4]
    e5 <- e[5]
    
    ##################################################################################################
    ###################################### Stage I ##################################################
    ##################################################################################################
    
    
    # stage I potential outcome: Y_1^A1, Y_1^A2, Y_1^A2
    Y1.A1 <- rbinom(n, 1, pi1.A1)
    Y1.A2 <- rbinom(n, 1, pi1.A2)
    Y1.A3 <- rbinom(n, 1, pi1.A3)
    
    # number(n0.burn) and proportion(p0.burn): first stage burn-in sample
    n0.burn <- round(n*p0.burn, 0)
    
    # treatment assignment A in stage I
    A_1 <- rep(NA, n)
    I1.A1 <- rep(NA, n)
    I1.A2 <- rep(NA, n)
    I1.A3 <- rep(NA, n)
    
    ## i=1,...,n0.burn equal assign
    A_1_n0.burn <- sample(1:3, size = n0.burn, replace = TRUE, prob = c(1/3,1/3,1/3))
    A_1[1:n0.burn] <- A_1_n0.burn
    I1.A1[1:n0.burn]<- ifelse(A_1_n0.burn==1, 1, 0)
    I1.A2[1:n0.burn]<- ifelse(A_1_n0.burn==2, 1, 0)
    I1.A3[1:n0.burn]<- ifelse(A_1_n0.burn==3, 1, 0)
    
    
    # stage I randomization probabilities
    P1.A1 <- rep(NA, n)
    P1.A2 <- rep(NA, n)
    P1.A3 <- rep(NA, n)
    P1.A1[1:n0.burn] <- rep(1/3, n0.burn)
    P1.A2[1:n0.burn] <- rep(1/3, n0.burn)
    P1.A3[1:n0.burn] <- rep(1/3, n0.burn)
    
    
    # observed outcome Y1 in stage I
    Y1 <- rep(NA, n)
    Y1[1:n0.burn] <- Y1.A1[1:n0.burn]*I1.A1[1:n0.burn]+Y1.A2[1:n0.burn]*I1.A2[1:n0.burn]+Y1.A3[1:n0.burn]*I1.A3[1:n0.burn]
    
    ## i=n0+1,...,n, update the randomization probability using all previous data
    if(n0.burn<n) {
      for (i in (n0.burn+1):n) { 
          # P(Y1.A1=1)
          Y1.A1hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])
          # P(Y1.A2=1)
          Y1.A2hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])
          # P(Y1.A3=1)
          Y1.A3hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])
        
        if (c=="n.N") {
          c1 <- i/n
        } else if (c=="n.2N") {
          c1 <- i/(2*n)
        } else {
          c1 <- c[1]
        }
        
        r1.A1 <- ((Y1.A1hat)^c1)/((Y1.A1hat)^c1+(Y1.A2hat)^c1+(Y1.A3hat)^c1)
        r1.A2 <- ((Y1.A2hat)^c1)/((Y1.A1hat)^c1+(Y1.A2hat)^c1+(Y1.A3hat)^c1)
        r1.A3 <- ((Y1.A3hat)^c1)/((Y1.A1hat)^c1+(Y1.A2hat)^c1+(Y1.A3hat)^c1)
        
        # AR(c,e)
        # P1.A1[i]
        if (is.na(r1.A1)){
          P1.A1[i] <- 1/3
        } else if (r1.A1 >= 0 & r1.A1 < e1) {
          P1.A1[i] <- e1
        } else if (r1.A1 > (1-e1) & r1.A1 <= 1) {
          P1.A1[i] <- 1-e1
        } else if (e1<= r1.A1 & r1.A1 <= 1-e1) {
          P1.A1[i] <- r1.A1
        }
        # P1.A2[i]
        if (is.na(r1.A2)){
          P1.A2[i] <- 1/3
        } else if (r1.A2 >= 0 & r1.A2 < e1) {
          P1.A2[i] <- e1
        } else if (r1.A2 > (1-e1) & r1.A2 <= 1) {
          P1.A2[i] <- 1-e1
        } else if (e1<= r1.A2 & r1.A2 <= 1-e1) {
          P1.A2[i] <- r1.A2
        }
        # P1.A3[i]
        if (is.na(r1.A3)){
          P1.A3[i] <- 1/3
        } else if (r1.A3 >= 0 & r1.A3 < e1) {
          P1.A3[i] <- e1
        } else if (r1.A3 > (1-e1) & r1.A3 <= 1) {
          P1.A3[i] <- 1-e1
        } else if (e1<= r1.A3 & r1.A3 <= 1-e1) {
          P1.A3[i] <- r1.A3
        }
        
        
        # check if the prob bound was used!
        P1.old <- c(P1.A1[i], P1.A2[i], P1.A3[i])
        P1.new <- rep(NA,3)
        # if e1=0.1 used, update the probability
        if(e1 %in% P1.old) {
          P1.new[which(e1==P1.old)] <- P1.old[which(e1==P1.old)]
          P1.new[-which(e1==P1.old)] <- (1-sum(P1.old[which(e1==P1.old)]))*P1.old[-which(e1==P1.old)]/sum(P1.old[-which(e1==P1.old)])    
          
          P1.A1[i] <- P1.new[1]
          P1.A2[i] <- P1.new[2]
          P1.A3[i] <- P1.new[3]
        } 
        
        A_1[i] <- sample(1:3, size = 1, replace = TRUE, prob = c(P1.A1[i],P1.A2[i],P1.A3[i]))
        I1.A1[i]<- ifelse(A_1[i]==1, 1, 0)
        I1.A2[i]<- ifelse(A_1[i]==2, 1, 0)
        I1.A3[i]<- ifelse(A_1[i]==3, 1, 0)
        Y1[i] <- Y1.A1[i]*I1.A1[i]+Y1.A2[i]*I1.A2[i]+Y1.A3[i]*I1.A3[i]
        
      }
    }
    
    
    ##################################################################################################
    ###################################### Stage II ##################################################
    ##################################################################################################
    
    # stage II potential outcome: Y_2^A1A2,Y_2^A1A3, Y_2^A2A1,Y_2^A2A3, Y_2^A3A1,Y_2^A3A2
    Y2.A1A2.raw <- rbinom(n, 1, pi2.A1A2)
    Y2.A1A3.raw <- rbinom(n, 1, pi2.A1A3)
    Y2.A2A1.raw <- rbinom(n, 1, pi2.A2A1)
    Y2.A2A3.raw <- rbinom(n, 1, pi2.A2A3)
    Y2.A3A1.raw <- rbinom(n, 1, pi2.A3A1)
    Y2.A3A2.raw <- rbinom(n, 1, pi2.A3A2)
    
    # if no response in stage I, Y_2^AjAl~Bernouli(pi_2^AjAl)
    Y2.A1A2 <- ifelse(Y1.A1==0, Y2.A1A2.raw, NA)
    Y2.A1A3 <- ifelse(Y1.A1==0, Y2.A1A3.raw, NA)
    Y2.A2A1 <- ifelse(Y1.A2==0, Y2.A2A1.raw, NA)
    Y2.A2A3 <- ifelse(Y1.A2==0, Y2.A2A3.raw, NA)
    Y2.A3A1 <- ifelse(Y1.A3==0, Y2.A3A1.raw, NA)
    Y2.A3A2 <- ifelse(Y1.A3==0, Y2.A3A2.raw, NA)
    
    
    # stage II randomization probabilities
    P2.A1A2 <- rep(NA, n)
    P2.A1A3 <- rep(NA, n)
    P2.A2A1 <- rep(NA, n)
    P2.A2A3 <- rep(NA, n)
    P2.A3A1 <- rep(NA, n)
    P2.A3A2 <- rep(NA, n)
    
    # treatment assignment A in stage II
    A_2 <- rep(NA, n)
    I2.A1 <- rep(NA, n)
    I2.A2 <- rep(NA, n)
    I2.A3 <- rep(NA, n)
    
    # observed outcome Y2 in stage II
    Y2 <- rep(NA, n)
    
    # overall observed outcome Y
    Y <- rep(NA, n)
    
    # i=1,...,n0: for those who did not respond in stage I, equally randomized to the rest treatments
    if(n0.burn<=n) {
      for (i in 1:n0.burn) {
        if (I1.A1[i]==1 & Y1.A1[i]==0) {
          P2.A1A2[i] <- 1/2
          P2.A1A3[i] <- 1-P2.A1A2[i]
          A_2[i] <- sample(2:3, size = 1, replace = TRUE, prob = c(P2.A1A2[i], P2.A1A3[i]))
        } else if (I1.A2[i]==1 & Y1.A2[i]==0) {
          P2.A2A1[i] <- 1/2
          P2.A2A3[i] <- 1-P2.A2A1[i]
          A_2[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P2.A2A1[i], P2.A2A3[i]))
        } else if (I1.A3[i]==1 & Y1.A3[i]==0) {
          P2.A3A1[i] <- 1/2
          P2.A3A2[i] <- 1-P2.A3A1[i]
          A_2[i] <- sample(1:2, size = 1, replace = TRUE, prob = c(P2.A3A1[i], P2.A3A2[i]))
        }
        
        I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
        I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
        I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
        
        
        
        Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],
                                        Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],
                                        Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i]), na.rm = T), NA)
        
        Y[i] <- ifelse(Y1[i]==1, Y1[i], Y2[i])
      }
      
      
    }
    
    # number(n1.burn) and proportion(p1.burn): second stage burn-in sample
    n1.burn <- round(n*p1.burn, 0)
    
    
    if(n0.burn<n1.burn) {
      # i=n0+1,...,n1: update the randomization probability using stage I observed response rate
      for (i in (n0.burn+1):n1.burn) {

        # P(Y1.A1=1)
        Y1.A1hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])
        # P(Y1.A2=1)
        Y1.A2hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])
        # P(Y1.A3=1)
        Y1.A3hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])
        
        
        if (c=="n.N") {
          c2 <- i/n
        } else if (c=="n.2N") {
          c2 <- i/(2*n)
        } else {
          c2 <- c[2]
        }
        
        r2.A1A2 <- ((Y1.A2hat)^c2)/((Y1.A2hat)^c2+(Y1.A3hat)^c2)
        r2.A1A3 <- 1-r2.A1A2
        r2.A2A1 <- ((Y1.A1hat)^c2)/((Y1.A1hat)^c2+(Y1.A3hat)^c2)
        r2.A2A3 <- 1-r2.A2A1
        r2.A3A1 <- ((Y1.A1hat)^c2)/((Y1.A1hat)^c2+(Y1.A2hat)^c2)
        r2.A3A2 <- 1-r2.A3A1
        
        # AR(c,e)
        
        if (I1.A1[i]==1 & Y1.A1[i]==0) {
          # P2.A1A2, P2.A1A3
          if (is.na(r2.A1A2)){
            P2.A1A2[i] <- 1/2
          } else if (r2.A1A2 >= 0 & r2.A1A2 < e2) {
            P2.A1A2[i] <- e2
          } else if (r2.A1A2 > (1-e2) & r2.A1A2 <= 1) {
            P2.A1A2[i] <- 1-e2
          } else if (e2<= r2.A1A2 & r2.A1A2 <= 1-e2) {
            P2.A1A2[i] <- r2.A1A2
          }
          P2.A1A3[i] <- 1-P2.A1A2[i]
          A_2[i] <- sample(2:3, size = 1, replace = TRUE, prob = c(P2.A1A2[i], P2.A1A3[i]))
          
        } else if (I1.A2[i]==1 & Y1.A2[i]==0) {
          # P2.A2A1, P2.A2A3
          if (is.na(r2.A2A1)){
            P2.A2A1[i] <- 1/2
          } else if (r2.A2A1 >= 0 & r2.A2A1 < e2) {
            P2.A2A1[i] <- e2
          } else if (r2.A2A1 > (1-e2) & r2.A2A1 <= 1) {
            P2.A2A1[i] <- 1-e2
          } else if (e2<= r2.A2A1 & r2.A2A1 <= 1-e2) {
            P2.A2A1[i] <- r2.A2A1
          }
          P2.A2A3[i] <- 1-P2.A2A1[i]
          A_2[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P2.A2A1[i], P2.A2A3[i]))
          
        } else if (I1.A3[i]==1 & Y1.A3[i]==0) {
          # P2.A3A1, P2.A3A2
          if (is.na(r2.A3A1)){
            P2.A3A1[i] <- 1/2
          } else if (r2.A3A1 >= 0 & r2.A3A1 < e2) {
            P2.A3A1[i] <- e2
          } else if (r2.A3A1 > (1-e2) & r2.A3A1 <= 1) {
            P2.A3A1[i] <- 1-e2
          } else if (e2<= r2.A3A1 & r2.A3A1 <= 1-e2) {
            P2.A3A1[i] <- r2.A3A1
          }
          P2.A3A2[i] <- 1-P2.A3A1[i]
          A_2[i] <- sample(1:2, size = 1, replace = TRUE, prob = c(P2.A3A1[i], P2.A3A2[i]))
        }
        
        I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
        I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
        I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
        Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],
                                        Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],
                                        Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i]), na.rm = T), NA)
        
        Y[i] <- ifelse(Y1[i]==1, Y1[i], Y2[i])
        
        
      }
      
    }
    
    if(n0.burn<=n1.burn & n1.burn<n) {
      # i=n1+1,...,n: update the randomization probability using following different types of method
      for (i in (n1.burn+1):n) {
        if (AR=="AR1") {
          
          # G-computation
          
          # P(Y2.A1A2=1|Y1.A1=0)
          Y2.A1A2hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)
          # P(Y2.A1A3=1|Y1.A1=0)
          Y2.A1A3hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)
          # P(Y2.A2A1=1|Y1.A2=0)
          Y2.A2A1hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)
          # P(Y2.A2A3=1|Y1.A2=0)
          Y2.A2A3hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)
          # P(Y2.A3A1=1|Y1.A3=0)
          Y2.A3A1hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)
          # P(Y2.A3A2=1|Y1.A3=0)
          Y2.A3A2hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)
          
          
          
          if (c=="n.N") {
            c3 <- i/n
          } else if (c=="n.2N") {
            c3 <- i/(2*n)
          } else {
            c3 <- c[3]
          }
          
          r3.A1A2 <- ((Y2.A1A2hat)^c3)/((Y2.A1A2hat)^c3+(Y2.A1A3hat)^c3)
          r3.A1A3 <- 1-r3.A1A2
          r3.A2A1 <- ((Y2.A2A1hat)^c3)/((Y2.A2A1hat)^c3+(Y2.A2A3hat)^c3)
          r3.A2A3 <- 1-r3.A2A1
          r3.A3A1 <- ((Y2.A3A1hat)^c3)/((Y2.A3A1hat)^c3+(Y2.A3A2hat)^c3)
          r3.A3A2 <- 1-r3.A3A1
          
          # AR(c,e)
          if (I1.A1[i]==1 & Y1.A1[i]==0) {
            # P2.A1A2, P2.A1A3
            if (is.na(r3.A1A2)){
              P2.A1A2[i] <- 1/2
            } else if (r3.A1A2 >= 0 & r3.A1A2 < e3) {
              P2.A1A2[i] <- e3
            } else if (r3.A1A2 > (1-e3) & r3.A1A2 <= 1) {
              P2.A1A2[i] <- 1-e3
            } else if (e3<= r3.A1A2 & r3.A1A2 <= 1-e3) {
              P2.A1A2[i] <- r3.A1A2
            }
            P2.A1A3[i] <- 1-P2.A1A2[i]
            A_2[i] <- sample(2:3, size = 1, replace = TRUE, prob = c(P2.A1A2[i], P2.A1A3[i]))
            
          } else if (I1.A2[i]==1 & Y1.A2[i]==0) {
            # P2.A2A1, P2.A2A3
            if (is.na(r3.A2A1)){
              P2.A2A1[i] <- 1/2
            } else if (r3.A2A1 >= 0 & r3.A2A1 < e3) {
              P2.A2A1[i] <- e3
            } else if (r3.A2A1 > (1-e3) & r3.A2A1 <= 1) {
              P2.A2A1[i] <- 1-e3
            } else if (e3<= r3.A2A1 & r3.A2A1 <= 1-e3) {
              P2.A2A1[i] <- r3.A2A1
            }
            P2.A2A3[i] <- 1-P2.A2A1[i]
            A_2[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P2.A2A1[i], P2.A2A3[i]))
            
          } else if (I1.A3[i]==1 & Y1.A3[i]==0) {
            # P2.A3A1, P2.A3A2
            if (is.na(r3.A3A1)){
              P2.A3A1[i] <- 1/2
            } else if (r3.A3A1 >= 0 & r3.A3A1 < e3) {
              P2.A3A1[i] <- e3
            } else if (r3.A3A1 > (1-e3) & r3.A3A1 <= 1) {
              P2.A3A1[i] <- 1-e3
            } else if (e3<= r3.A3A1 & r3.A3A1 <= 1-e3) {
              P2.A3A1[i] <- r3.A3A1
            }
            P2.A3A2[i] <- 1-P2.A3A1[i]
            A_2[i] <- sample(1:2, size = 1, replace = TRUE, prob = c(P2.A3A1[i], P2.A3A2[i]))
          }
          
          I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
          I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
          I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
          Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],
                                          Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],
                                          Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i]), na.rm = T), NA)
          
          Y[i] <- ifelse(Y1[i]==1, Y1[i], Y2[i])
          
        } else if (AR=="AR2") {
          # G-computation
          # P(Y[d(A1, A2)]=1)
          Yd.A1A2hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])+sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)]))/sum(I1.A1[1:(i-1)])*
            sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)
          # P(Y[d(A1, A3)]=1)
          Yd.A1A3hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])+sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)]))/sum(I1.A1[1:(i-1)])*
            sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)
          # P(Y[d(A2, A1)]=1)
          Yd.A2A1hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])+sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)]))/sum(I1.A2[1:(i-1)])*
            sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)
          # P(Y[d(A2, A3)]=1)
          Yd.A2A3hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])+sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)]))/sum(I1.A2[1:(i-1)])*
            sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)
          # P(Y[d(A3, A1)]=1)
          Yd.A3A1hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])+sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)]))/sum(I1.A3[1:(i-1)])*
            sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)
          # P(Y[d(A3, A2)]=1)
          Yd.A3A2hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])+sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)]))/sum(I1.A3[1:(i-1)])*
            sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)
          
          
          if (c=="n.N") {
            c5 <- i/n
          } else if (c=="n.2N") {
            c5 <- i/(2*n)
          } else {
            c5 <- c[5]
          }
          
          
          r5.A1A2 <- (Yd.A1A2hat)^c5/((Yd.A1A2hat)^c5+(Yd.A1A3hat)^c5)
          r5.A1A3 <- 1-r5.A1A2
          r5.A2A1 <- (Yd.A2A1hat)^c5/((Yd.A2A1hat)^c5+(Yd.A2A3hat)^c5)
          r5.A2A3 <- 1-r5.A2A1
          r5.A3A1 <- (Yd.A3A1hat)^c5/((Yd.A3A1hat)^c5+(Yd.A3A2hat)^c5)
          r5.A3A2 <- 1-r5.A3A1
          
          # AR(c,e)
          if (I1.A1[i]==1 & Y1.A1[i]==0) {
            # P2.A1A2, P2.A1A3
            if (is.na(r5.A1A2)){
              P2.A1A2[i] <- 1/2
            } else if (r5.A1A2 >= 0 & r5.A1A2 < e5) {
              P2.A1A2[i] <- e5
            } else if (r5.A1A2 > (1-e5) & r5.A1A2 <= 1) {
              P2.A1A2[i] <- 1-e5
            } else if (e5<= r5.A1A2 & r5.A1A2 <= 1-e5) {
              P2.A1A2[i] <- r5.A1A2
            }
            P2.A1A3[i] <- 1-P2.A1A2[i]
            A_2[i] <- sample(2:3, size = 1, replace = TRUE, prob = c(P2.A1A2[i], P2.A1A3[i]))
            
          } else if (I1.A2[i]==1 & Y1.A2[i]==0) {
            # P2.A2A1, P2.A2A3
            if (is.na(r5.A2A1)){
              P2.A2A1[i] <- 1/2
            } else if (r5.A2A1 >= 0 & r5.A2A1 < e5) {
              P2.A2A1[i] <- e5
            } else if (r5.A2A1 > (1-e5) & r5.A2A1 <= 1) {
              P2.A2A1[i] <- 1-e5
            } else if (e5<= r5.A2A1 & r5.A2A1 <= 1-e5) {
              P2.A2A1[i] <- r5.A2A1
            }
            P2.A2A3[i] <- 1-P2.A2A1[i]
            A_2[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P2.A2A1[i], P2.A2A3[i]))
            
          } else if (I1.A3[i]==1 & Y1.A3[i]==0) {
            # P2.A3A1, P2.A3A2
            if (is.na(r5.A3A1)){
              P2.A3A1[i] <- 1/2
            } else if (r5.A3A1 >= 0 & r5.A3A1 < e5) {
              P2.A3A1[i] <- e5
            } else if (r5.A3A1 > (1-e5) & r5.A3A1 <= 1) {
              P2.A3A1[i] <- 1-e5
            } else if (e5<= r5.A3A1 & r5.A3A1 <= 1-e5) {
              P2.A3A1[i] <- r5.A3A1
            }
            P2.A3A2[i] <- 1-P2.A3A1[i]
            A_2[i] <- sample(1:2, size = 1, replace = TRUE, prob = c(P2.A3A1[i], P2.A3A2[i]))
          }
          
          I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
          I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
          I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
          Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],
                                          Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],
                                          Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i]), na.rm = T), NA)
          
          Y[i] <- ifelse(Y1[i]==1, Y1[i], Y2[i])
          
        }
      }
    }
    
    df <- data.frame(Y1.A1, Y1.A2, Y1.A3, 
                     P1.A1, P1.A2, P1.A3, 
                     Y2.A1A2, Y2.A1A3, Y2.A2A1, Y2.A2A3, Y2.A3A1, Y2.A3A2,
                     A_1, I1.A1, I1.A2, I1.A3, 
                     A_2, I2.A1, I2.A2, I2.A3,
                     P2.A1A2, P2.A1A3, P2.A2A1, P2.A2A3, P2.A3A1, P2.A3A2,
                     Y1, Y2, Y)
    return(df)
  
}


# Input variables for function "data_gene_ra":
## n: sample size; e.g. n=500
## pi1.A1, pi1.A2, pi1.A3: stage 1 true response rate; e.g. pi1.A1=0.5; pi1.A2=0.35; pi1.A3=0.2
## pi2.A1A2, pi2.A1A3, pi2.A2A1, pi2.A2A3, pi2.A3A1, pi2.A3A2: stage 2 true response rate; e.g. pi2.A1A2=0.3; pi2.A1A3=0.4; pi2.A2A1=0.35; pi2.A2A3=0.2; pi2.A3A1=0.25; pi2.A3A2=0.1
## p0.burn: stage 1 burn-in sample proportion; e.g. p0.burn=0.25
## Q: adjusted stage 2 probability for assigning to the inferior treatment; e.g. Q=0.2


data_gene_ra <- function(n,pi1.A1,pi1.A2,pi1.A3,pi2.A1A2,pi2.A1A3,pi2.A2A1,pi2.A2A3,pi2.A3A1,pi2.A3A2,p0.burn,Q) {
  
  ################## Stage I #####################
  
  # stage I potential outcome: Y_1^A1, Y_1^A2, Y_1^A2
  Y1.A1 <- rbinom(n, 1, pi1.A1)
  Y1.A2 <- rbinom(n, 1, pi1.A2)
  Y1.A3 <- rbinom(n, 1, pi1.A3)
  
  # number(n0.burn) and proportion(p0.burn): first stage burn-in sample
  n0.burn <- round(n*p0.burn, 0)
  
  # treatment assignment A in stage I: i=1,...,n.burn equal assign
  A_1 <- sample(1:3, size = n, replace = TRUE, prob = c(1/3,1/3,1/3))
  I1.A1 <- ifelse(A_1==1, 1, 0)
  I1.A2 <- ifelse(A_1==2, 1, 0)
  I1.A3 <- ifelse(A_1==3, 1, 0)
  
  # stage I randomization probabilities
  P1.A1 <- rep(1/3, n)
  P1.A2 <- rep(1/3, n)
  P1.A3 <- rep(1/3, n)
  
  # observed outcome Y1 in stage I
  Y1 <- Y1.A1*I1.A1+Y1.A2*I1.A2+Y1.A3*I1.A3
  
  ################## Stage II #####################
  
  # stage II potential outcome: Y_2^A1A2,Y_2^A1A3, Y_2^A2A1,Y_2^A2A3, Y_2^A3A1,Y_2^A3A2
  Y2.A1A2.raw <- rbinom(n, 1, pi2.A1A2)
  Y2.A1A3.raw <- rbinom(n, 1, pi2.A1A3)
  Y2.A2A1.raw <- rbinom(n, 1, pi2.A2A1)
  Y2.A2A3.raw <- rbinom(n, 1, pi2.A2A3)
  Y2.A3A1.raw <- rbinom(n, 1, pi2.A3A1)
  Y2.A3A2.raw <- rbinom(n, 1, pi2.A3A2)
  
  # if no response in stage I, Y_2^AjAl~Bernouli(pi_2^AjAl)
  Y2.A1A2 <- ifelse(Y1.A1==0, Y2.A1A2.raw, NA)
  Y2.A1A3 <- ifelse(Y1.A1==0, Y2.A1A3.raw, NA)
  Y2.A2A1 <- ifelse(Y1.A2==0, Y2.A2A1.raw, NA)
  Y2.A2A3 <- ifelse(Y1.A2==0, Y2.A2A3.raw, NA)
  Y2.A3A1 <- ifelse(Y1.A3==0, Y2.A3A1.raw, NA)
  Y2.A3A2 <- ifelse(Y1.A3==0, Y2.A3A2.raw, NA)
  
  
  # stage II randomization probabilities
  P2.A1A2 <- rep(NA, n)
  P2.A1A3 <- rep(NA, n)
  P2.A2A1 <- rep(NA, n)
  P2.A2A3 <- rep(NA, n)
  P2.A3A1 <- rep(NA, n)
  P2.A3A2 <- rep(NA, n)
  P2 <- cbind(P2.A1A2, P2.A1A3, P2.A2A1, P2.A2A3, P2.A3A1, P2.A3A2)
  colnames(P2) <- c("P2.A1A2", "P2.A1A3", "P2.A2A1", "P2.A2A3", "P2.A3A1", "P2.A3A2")
  
  # treatment assignment A in stage II
  A_2 <- rep(NA, n)
  I2.A1 <- rep(NA, n)
  I2.A2 <- rep(NA, n)
  I2.A3 <- rep(NA, n)
  
  # observed outcome Y2 in stage II
  Y2 <- rep(NA, n)
  
  # overall observed outcome Y
  Y <- rep(NA, n)
  
  # i=1,...,n0: for those who did not respond in stage I, equally randomized to the rest treatments
  
  for (i in 1:n0.burn) {
    if (I1.A1[i]==1 & Y1.A1[i]==0) {
      P2.A1A2[i] <- 1/2
      P2.A1A3[i] <- 1-P2.A1A2[i]
      A_2[i] <- sample(2:3, size = 1, replace = TRUE, prob = c(P2.A1A2[i], P2.A1A3[i]))
    } else if (I1.A2[i]==1 & Y1.A2[i]==0) {
      P2.A2A1[i] <- 1/2
      P2.A2A3[i] <- 1-P2.A2A1[i]
      A_2[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P2.A2A1[i], P2.A2A3[i]))
    } else if (I1.A3[i]==1 & Y1.A3[i]==0) {
      P2.A3A1[i] <- 1/2
      P2.A3A2[i] <- 1-P2.A3A1[i]
      A_2[i] <- sample(1:2, size = 1, replace = TRUE, prob = c(P2.A3A1[i], P2.A3A2[i]))
    }
    
    I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
    I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
    I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
    
    
    Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],
                                    Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],
                                    Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i]), na.rm = T), NA)
    
    Y[i] <- ifelse(Y1[i]==1, Y1[i], Y2[i])
  }
  
  # first stage response rate from patients 1 to n0
  # P(Y1.A1=1)
  Y1.A1hat <- sum(I1.A1[1:n0.burn]*Y1[1:n0.burn])/sum(I1.A1[1:n0.burn])
  # P(Y1.A2=1)
  Y1.A2hat <- sum(I1.A2[1:n0.burn]*Y1[1:n0.burn])/sum(I1.A2[1:n0.burn])
  # P(Y1.A3=1)
  Y1.A3hat <- sum(I1.A3[1:n0.burn]*Y1[1:n0.burn])/sum(I1.A3[1:n0.burn])
  
  infe <- which.min(c(Y1.A1hat, Y1.A2hat, Y1.A3hat))
  
  if(n0.burn<n) {
    for (i in (n0.burn+1):n) {
      # if receive the inferior treatment in stage I and did not respond
      if (A_1[i]==infe & Y1[i]==0) {
        P2[i,] <- ifelse(names(P2[i,])==paste0("P2.A",infe,"A",c(1,2,3)[-A_1[i]][1]) | names(P2[i,])==paste0("P2.A",infe,"A",c(1,2,3)[-A_1[i]][2]), 1/2, NA)
        A_2[i] <- sample(c(1,2,3)[-infe], size = 1, replace = TRUE, prob = c(1/2, 1/2))
      } # if not receive the inferior treatment in stage I and did not respond
      else if (A_1[i]!=infe & Y1[i]==0) {
        P2[i,which(names(P2[i,])==paste0("P2.A",A_1[i],"A",infe))] = Q
        P2[i,which(names(P2[i,])==paste0("P2.A",A_1[i],"A",c(1,2,3)[-c(A_1[i], infe)]))] = 1-Q
        A_2[i] <- sample(c(infe,c(1,2,3)[-c(A_1[i],infe)]), size=1, replace = TRUE, prob = c(Q, 1-Q))
        
      }
      I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
      I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
      I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
      
      Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],
                                      Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],
                                      Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i]), na.rm = T), NA)
      
      Y[i] <- ifelse(Y1[i]==1, Y1[i], Y2[i])
      
      P2.A1A2[i] <- P2[i,1]
      P2.A1A3[i] <- P2[i,2]
      P2.A2A1[i] <- P2[i,3]
      P2.A2A3[i] <- P2[i,4]
      P2.A3A1[i] <- P2[i,5]
      P2.A3A2[i] <- P2[i,6]
      
    }
    
  }
  
  df <- data.frame(Y1.A1, Y1.A2, Y1.A3, 
                   P1.A1, P1.A2, P1.A3, 
                   Y2.A1A2, Y2.A1A3, Y2.A2A1, Y2.A2A3, Y2.A3A1, Y2.A3A2,
                   A_1, I1.A1, I1.A2, I1.A3, 
                   A_2, I2.A1, I2.A2, I2.A3,
                   P2.A1A2, P2.A1A3, P2.A2A1, P2.A2A3, P2.A3A1, P2.A3A2,
                   Y1, Y2, Y)
  return(df)
  
}


# Input variables for function "monte":
## N: number of Monte Carlo replications; e.g. N=10000
## n: sample size; e.g. n=500
## RA.SMART: whether generate data from RA-SMART design; e.g. RA.SMART=FALSE
## Q: if RA.SMART=TRUE, Q is the adjusted stage 2 probability for assigning to the inferior treatment; e.g. Q=0.2
## pi = c(pi1.A1, pi1.A2, pi1.A3, pi2.A1A2, pi2.A1A3, pi2.A2A1, pi2.A2A3, pi2.A3A1, pi2.A3A2): vector of stage 1 and stage 2 true response rates; e.g. pi=c(0.5, 0.35, 0.2, 0.3, 0.4, 0.35, 0.2, 0.25, 0.1)
## p0.burn: stage 1 burn-in sample proportion; e.g. p0.burn=0.25
## p1.burn: stage 2 burn-in sample proportion; e.g. p1.burn=0.5
## AR: GO-SMART AR-1 or AR-2 design; e.g. AR="AR1"
## c: tuning parameter that controls the degree of dependence of the OAR probabilities on the empirical response rates; e.g. c="n.N" (means c=i/n in the manuscript notation)
## e: randomization constraint parameter, restrict the randomization probability to the interval [e,1-e]; e.g. e=c(0.1, 0.1, 0.1, 0.1, 0.1)
## sce: scenario; e.g. sce="S1"

monte <- function(N,n,RA.SMART,Q,pi,p0.burn,p1.burn,AR,c,e,sce) {
  
  # true response rate
  pi1.A1=pi[1]
  pi1.A2=pi[2]
  pi1.A3=pi[3]
  pi2.A1A2=pi[4]
  pi2.A1A3=pi[5]
  pi2.A2A1=pi[6]
  pi2.A2A3=pi[7]
  pi2.A3A1=pi[8]
  pi2.A3A2=pi[9]
  
  # sample size
  n <- n
  
  # Monte Carlo sample size
  N <- N
  
  # first-stage burn-in sample
  n0 <- n*p0.burn
  
  # treatment assignment
  I1.A1 <- matrix(NA, nrow=N, ncol=n)
  I1.A2 <- matrix(NA, nrow=N, ncol=n)
  I1.A3 <- matrix(NA, nrow=N, ncol=n)
  I2.A1 <- matrix(NA, nrow=N, ncol=n)
  I2.A2 <- matrix(NA, nrow=N, ncol=n)
  I2.A3 <- matrix(NA, nrow=N, ncol=n)
  # observed outcome
  Y1 <- matrix(NA, nrow=N, ncol=n)
  Y2 <- matrix(NA, nrow=N, ncol=n)
  Y <- matrix(NA, nrow=N, ncol=n)
  # counterfactual outcome
  Y1.A1 <- matrix(NA, nrow=N, ncol=n)
  Y1.A2 <- matrix(NA, nrow=N, ncol=n)
  Y1.A3 <- matrix(NA, nrow=N, ncol=n)
  Y2.A1A2 <- matrix(NA, nrow=N, ncol=n)
  Y2.A1A3 <- matrix(NA, nrow=N, ncol=n)
  Y2.A2A1 <- matrix(NA, nrow=N, ncol=n)
  Y2.A2A3 <- matrix(NA, nrow=N, ncol=n)
  Y2.A3A1 <- matrix(NA, nrow=N, ncol=n)
  Y2.A3A2 <- matrix(NA, nrow=N, ncol=n)
  # randomization probability
  P1.A1 <- matrix(NA, nrow=N, ncol=n)
  P1.A2 <- matrix(NA, nrow=N, ncol=n)
  P1.A3 <- matrix(NA, nrow=N, ncol=n)
  P2.A1A2 <- matrix(NA, nrow=N, ncol=n)
  P2.A1A3 <- matrix(NA, nrow=N, ncol=n)
  P2.A2A1 <- matrix(NA, nrow=N, ncol=n)
  P2.A2A3 <- matrix(NA, nrow=N, ncol=n)
  P2.A3A1 <- matrix(NA, nrow=N, ncol=n)
  P2.A3A2 <- matrix(NA, nrow=N, ncol=n)
  
  
  
  
  # IPW estimator
  mu.A1A2hat.ipw <- rep(NA, N)
  mu.A1A3hat.ipw <- rep(NA, N)
  mu.A2A1hat.ipw <- rep(NA, N)
  mu.A2A3hat.ipw <- rep(NA, N)
  mu.A3A1hat.ipw <- rep(NA, N)
  mu.A3A2hat.ipw <- rep(NA, N)
  
  var.A1A2hat.ipw <- rep(NA, N)
  var.A1A3hat.ipw <- rep(NA, N)
  var.A2A1hat.ipw <- rep(NA, N)
  var.A2A3hat.ipw <- rep(NA, N)
  var.A3A1hat.ipw <- rep(NA, N)
  var.A3A2hat.ipw <- rep(NA, N)
  
  
  # IPWN estimator
  mu.A1A2hat.ipwn <- rep(NA, N)
  mu.A1A3hat.ipwn <- rep(NA, N)
  mu.A2A1hat.ipwn <- rep(NA, N)
  mu.A2A3hat.ipwn <- rep(NA, N)
  mu.A3A1hat.ipwn <- rep(NA, N)
  mu.A3A2hat.ipwn <- rep(NA, N)
  
  var.A1A2hat.ipwn <- rep(NA, N)
  var.A1A3hat.ipwn <- rep(NA, N)
  var.A2A1hat.ipwn <- rep(NA, N)
  var.A2A3hat.ipwn <- rep(NA, N)
  var.A3A1hat.ipwn <- rep(NA, N)
  var.A3A2hat.ipwn <- rep(NA, N)
  
  
  
  
  # G-estimator
  mu.A1A2hat.g <- rep(NA, N)
  mu.A1A3hat.g <- rep(NA, N)
  mu.A2A1hat.g <- rep(NA, N)
  mu.A2A3hat.g <- rep(NA, N)
  mu.A3A1hat.g <- rep(NA, N)
  mu.A3A2hat.g <- rep(NA, N)
  
  var.A1A2hat.g <- rep(NA, N)
  var.A1A3hat.g <- rep(NA, N)
  var.A2A1hat.g <- rep(NA, N)
  var.A2A3hat.g <- rep(NA, N)
  var.A3A1hat.g <- rep(NA, N)
  var.A3A2hat.g <- rep(NA, N)
  
  
  # simplest sample mean estimator
  mu.A1A2hat.sm <- rep(NA, N)
  mu.A1A3hat.sm <- rep(NA, N)
  mu.A2A1hat.sm <- rep(NA, N)
  mu.A2A3hat.sm <- rep(NA, N)
  mu.A3A1hat.sm <- rep(NA, N)
  mu.A3A2hat.sm <- rep(NA, N)
  
  var.A1A2hat.sm <- rep(NA, N)
  var.A1A3hat.sm <- rep(NA, N)
  var.A2A1hat.sm <- rep(NA, N)
  var.A2A3hat.sm <- rep(NA, N)
  var.A3A1hat.sm <- rep(NA, N)
  var.A3A2hat.sm <- rep(NA, N)
  
  # estimated optimal/worst DTR
  est.opt.ipw <- rep(NA, N)
  est.opt.ipwn <- rep(NA, N)
  est.opt.g <- rep(NA, N)
  est.opt.sm <- rep(NA, N)
  est.wor.ipw <- rep(NA, N)
  est.wor.ipwn <- rep(NA, N)
  est.wor.g <- rep(NA, N)
  est.wor.sm <- rep(NA, N)
  
  # number of treated in each DTRs
  n.A1A2 <- rep(NA, N)
  n.A1A3 <- rep(NA, N)
  n.A2A1 <- rep(NA, N)
  n.A2A3 <- rep(NA, N)
  n.A3A1 <- rep(NA, N)
  n.A3A2 <- rep(NA, N)
  
  # number of responder in each DTRs
  rn.A1A2 <- rep(NA, N)
  rn.A1A3 <- rep(NA, N)
  rn.A2A1 <- rep(NA, N)
  rn.A2A3 <- rep(NA, N)
  rn.A3A1 <- rep(NA, N)
  rn.A3A2 <- rep(NA, N)
  
  # total number of response in the trial
  NR <- rep(NA, N)
  
  
  # create N Monte Carlo samples
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for (i in 1:N) {
    if (RA.SMART==TRUE) {
      ######### RA-SMART design
      df <- data_gene_ra(n=n,pi1.A1=pi1.A1,pi1.A2=pi1.A2,pi1.A3=pi1.A3,
                         pi2.A1A2=pi2.A1A2,pi2.A1A3=pi2.A1A3,
                         pi2.A2A1=pi2.A2A1,pi2.A2A3=pi2.A2A3,
                         pi2.A3A1=pi2.A3A1,pi2.A3A2=pi2.A3A2,
                         p0.burn=p0.burn, Q=Q)
      
    } else if (RA.SMART==FALSE) {
      ######## SMART and GO-SMART design
      df <- data_gene_c(n=n,pi1.A1=pi1.A1,pi1.A2=pi1.A2,pi1.A3=pi1.A3,
                        pi2.A1A2=pi2.A1A2,pi2.A1A3=pi2.A1A3,
                        pi2.A2A1=pi2.A2A1,pi2.A2A3=pi2.A2A3,
                        pi2.A3A1=pi2.A3A1,pi2.A3A2=pi2.A3A2,
                        p0.burn=p0.burn,p1.burn=p1.burn,AR=AR,
                        e=e, c=c)
      
    } 
    
    # treatment assignment
    I1.A1[i,] <- df$I1.A1
    I1.A2[i,] <- df$I1.A2
    I1.A3[i,] <- df$I1.A3
    I2.A1[i,] <- df$I2.A1
    I2.A2[i,] <- df$I2.A2
    I2.A3[i,] <- df$I2.A3
    # observed outcome
    Y1[i,] <- df$Y1
    Y2[i,] <- df$Y2
    Y[i,] <- df$Y
    # counterfactual outcome
    Y1.A1[i,] <- df$Y1.A1
    Y1.A2[i,] <- df$Y1.A2
    Y1.A3[i,] <- df$Y1.A3
    Y2.A1A2[i,] <- df$Y2.A1A2
    Y2.A1A3[i,] <- df$Y2.A1A3
    Y2.A2A1[i,] <- df$Y2.A2A1
    Y2.A2A3[i,] <- df$Y2.A2A3
    Y2.A3A1[i,] <- df$Y2.A3A1
    Y2.A3A2[i,] <- df$Y2.A3A2
    # randomization probability
    P1.A1[i,] <- df$P1.A1
    P1.A2[i,] <- df$P1.A2
    P1.A3[i,] <- df$P1.A3
    P2.A1A2[i,] <- df$P2.A1A2
    P2.A1A3[i,] <- df$P2.A1A3
    P2.A2A1[i,] <- df$P2.A2A1
    P2.A2A3[i,] <- df$P2.A2A3
    P2.A3A1[i,] <- df$P2.A3A1
    P2.A3A2[i,] <- df$P2.A3A2
    
    # progress indicator
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  # create N Monte Carlo samples
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for (i in 1:N) {
    
    if (RA.SMART==TRUE) {
      ######### RA-SMART design ###########
      
      # most simple sample mean estimator
      mu.A1A2hat.sm[i] <- (sum(Y1[i,]*I1.A1[i,])+sum(Y[i,]*I1.A1[i,]*I2.A2[i,],na.rm = T))/(sum(Y1[i,]*I1.A1[i,])+sum(I1.A1[i,]*I2.A2[i,],na.rm = T))
      mu.A1A3hat.sm[i] <- (sum(Y1[i,]*I1.A1[i,])+sum(Y[i,]*I1.A1[i,]*I2.A3[i,],na.rm = T))/(sum(Y1[i,]*I1.A1[i,])+sum(I1.A1[i,]*I2.A3[i,],na.rm = T))
      mu.A2A1hat.sm[i] <- (sum(Y1[i,]*I1.A2[i,])+sum(Y[i,]*I1.A2[i,]*I2.A1[i,],na.rm = T))/(sum(Y1[i,]*I1.A2[i,])+sum(I1.A2[i,]*I2.A1[i,],na.rm = T))
      mu.A2A3hat.sm[i] <- (sum(Y1[i,]*I1.A2[i,])+sum(Y[i,]*I1.A2[i,]*I2.A3[i,],na.rm = T))/(sum(Y1[i,]*I1.A2[i,])+sum(I1.A2[i,]*I2.A3[i,],na.rm = T))
      mu.A3A1hat.sm[i] <- (sum(Y1[i,]*I1.A3[i,])+sum(Y[i,]*I1.A3[i,]*I2.A1[i,],na.rm = T))/(sum(Y1[i,]*I1.A3[i,])+sum(I1.A3[i,]*I2.A1[i,],na.rm = T))
      mu.A3A2hat.sm[i] <- (sum(Y1[i,]*I1.A3[i,])+sum(Y[i,]*I1.A3[i,]*I2.A2[i,],na.rm = T))/(sum(Y1[i,]*I1.A3[i,])+sum(I1.A3[i,]*I2.A2[i,],na.rm = T))
      
      var.A1A2hat.sm[i] <- mu.A1A2hat.sm[i]*(1-mu.A1A2hat.sm[i])/(sum(Y1[i,]*I1.A1[i,])+sum(I1.A1[i,]*I2.A2[i,],na.rm = T))
      var.A1A3hat.sm[i] <- mu.A1A3hat.sm[i]*(1-mu.A1A3hat.sm[i])/(sum(Y1[i,]*I1.A1[i,])+sum(I1.A1[i,]*I2.A3[i,],na.rm = T))
      var.A2A1hat.sm[i] <- mu.A2A1hat.sm[i]*(1-mu.A2A1hat.sm[i])/(sum(Y1[i,]*I1.A2[i,])+sum(I1.A2[i,]*I2.A1[i,],na.rm = T))
      var.A2A3hat.sm[i] <- mu.A2A3hat.sm[i]*(1-mu.A2A3hat.sm[i])/(sum(Y1[i,]*I1.A2[i,])+sum(I1.A2[i,]*I2.A3[i,],na.rm = T))
      var.A3A1hat.sm[i] <- mu.A3A1hat.sm[i]*(1-mu.A3A1hat.sm[i])/(sum(Y1[i,]*I1.A3[i,])+sum(I1.A3[i,]*I2.A1[i,],na.rm = T))
      var.A3A2hat.sm[i] <- mu.A3A2hat.sm[i]*(1-mu.A3A2hat.sm[i])/(sum(Y1[i,]*I1.A3[i,])+sum(I1.A3[i,]*I2.A2[i,],na.rm = T))
      
      
      # IPW estimator
      mu.A1A2hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,])))/n
      mu.A1A3hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,])))/n
      mu.A2A1hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,])))/n
      mu.A2A3hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,])))/n
      mu.A3A1hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,])))/n
      mu.A3A2hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,])))/n
      
      
      var.A1A2hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,]))-mu.A1A2hat.ipw[i])^2)/(n^2)
      var.A1A3hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,]))-mu.A1A3hat.ipw[i])^2)/(n^2)
      var.A2A1hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,]))-mu.A2A1hat.ipw[i])^2)/(n^2)
      var.A2A3hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,]))-mu.A2A3hat.ipw[i])^2)/(n^2)
      var.A3A1hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,]))-mu.A3A1hat.ipw[i])^2)/(n^2)
      var.A3A2hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,]))-mu.A3A2hat.ipw[i])^2)/(n^2)
      
      
      # IPWN estimator
      w.A1A2 <- Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na((1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,])), 0, (1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,]))
      w.A1A3 <- Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na((1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,])), 0, (1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,]))
      w.A2A1 <- Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na((1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,])), 0, (1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,]))
      w.A2A3 <- Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na((1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,])), 0, (1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,]))
      w.A3A1 <- Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na((1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,])), 0, (1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,]))
      w.A3A2 <- Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na((1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,])), 0, (1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,]))
      
      mu.A1A2hat.ipwn[i] <- mu.A1A2hat.ipw[i]/(sum(w.A1A2)/n)
      mu.A1A3hat.ipwn[i] <- mu.A1A3hat.ipw[i]/(sum(w.A1A3)/n)
      mu.A2A1hat.ipwn[i] <- mu.A2A1hat.ipw[i]/(sum(w.A2A1)/n)
      mu.A2A3hat.ipwn[i] <- mu.A2A3hat.ipw[i]/(sum(w.A2A3)/n)
      mu.A3A1hat.ipwn[i] <- mu.A3A1hat.ipw[i]/(sum(w.A3A1)/n)
      mu.A3A2hat.ipwn[i] <- mu.A3A2hat.ipw[i]/(sum(w.A3A2)/n)
      
      var.A1A2hat.ipwn[i] <- sum((w.A1A2*(Y[i,]-mu.A1A2hat.ipwn[i]))^2)/(n^2)
      var.A1A3hat.ipwn[i] <- sum((w.A1A3*(Y[i,]-mu.A1A3hat.ipwn[i]))^2)/(n^2)
      var.A2A1hat.ipwn[i] <- sum((w.A2A1*(Y[i,]-mu.A2A1hat.ipwn[i]))^2)/(n^2)
      var.A2A3hat.ipwn[i] <- sum((w.A2A3*(Y[i,]-mu.A2A3hat.ipwn[i]))^2)/(n^2)
      var.A3A1hat.ipwn[i] <- sum((w.A3A1*(Y[i,]-mu.A3A1hat.ipwn[i]))^2)/(n^2)
      var.A3A2hat.ipwn[i] <- sum((w.A3A2*(Y[i,]-mu.A3A2hat.ipwn[i]))^2)/(n^2)
      
      
      # G-estimation
      Q_12_tilde <- sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0],na.rm = T)/sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0]))
      Q_13_tilde <- sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0],na.rm = T)/sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0]))
      Q_21_tilde <- sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0],na.rm = T)/sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0]))
      Q_23_tilde <- sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0],na.rm = T)/sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0]))
      Q_31_tilde <- sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0],na.rm = T)/sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0]))
      Q_32_tilde <- sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0],na.rm = T)/sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0]))
      
      pi1.A1.hat.g <- sum(Y1[i,]*I1.A1[i,])/sum(I1.A1[i,])
      pi1.A2.hat.g <- sum(Y1[i,]*I1.A2[i,])/sum(I1.A2[i,])
      pi1.A3.hat.g <- sum(Y1[i,]*I1.A3[i,])/sum(I1.A3[i,])
      
      if (n0<n) {
        
        Q_12_tilde_prime <- sum(I1.A1[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A2[i,(n0+1):n],na.rm = T)/sum(I1.A1[i,(n0+1):n]*(1-Y1[i,(n0+1):n]))
        Q_13_tilde_prime <- sum(I1.A1[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A3[i,(n0+1):n],na.rm = T)/sum(I1.A1[i,(n0+1):n]*(1-Y1[i,(n0+1):n]))
        Q_21_tilde_prime <- sum(I1.A2[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A1[i,(n0+1):n],na.rm = T)/sum(I1.A2[i,(n0+1):n]*(1-Y1[i,(n0+1):n]))
        Q_23_tilde_prime <- sum(I1.A2[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A3[i,(n0+1):n],na.rm = T)/sum(I1.A2[i,(n0+1):n]*(1-Y1[i,(n0+1):n]))
        Q_31_tilde_prime <- sum(I1.A3[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A1[i,(n0+1):n],na.rm = T)/sum(I1.A3[i,(n0+1):n]*(1-Y1[i,(n0+1):n]))
        Q_32_tilde_prime <- sum(I1.A3[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A2[i,(n0+1):n],na.rm = T)/sum(I1.A3[i,(n0+1):n]*(1-Y1[i,(n0+1):n])) 
        
        pi2.A1A2.hat.g <- (sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0]/Q_12_tilde*Y2[i,1:n0],na.rm = T)+
                             sum(I1.A1[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A2[i,(n0+1):n]/Q_12_tilde_prime*Y2[i,(n0+1):n],na.rm = T))/
          (sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0]/Q_12_tilde,na.rm = T)+
             sum(I1.A1[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A2[i,(n0+1):n]/Q_12_tilde_prime,na.rm = T))
        
        pi2.A1A3.hat.g <- (sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0]/Q_13_tilde*Y2[i,1:n0],na.rm = T)+
                             sum(I1.A1[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A3[i,(n0+1):n]/Q_13_tilde_prime*Y2[i,(n0+1):n],na.rm = T))/
          (sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0]/Q_13_tilde,na.rm = T)+
             sum(I1.A1[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A3[i,(n0+1):n]/Q_13_tilde_prime,na.rm = T))
        
        pi2.A2A1.hat.g <- (sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0]/Q_21_tilde*Y2[i,1:n0],na.rm = T)+
                             sum(I1.A2[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A1[i,(n0+1):n]/Q_21_tilde_prime*Y2[i,(n0+1):n],na.rm = T))/
          (sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0]/Q_21_tilde,na.rm = T)+
             sum(I1.A2[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A1[i,(n0+1):n]/Q_21_tilde_prime,na.rm = T))
        
        pi2.A2A3.hat.g <- (sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0]/Q_23_tilde*Y2[i,1:n0],na.rm = T)+
                             sum(I1.A2[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A3[i,(n0+1):n]/Q_23_tilde_prime*Y2[i,(n0+1):n],na.rm = T))/
          (sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0]/Q_23_tilde,na.rm = T)+
             sum(I1.A2[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A3[i,(n0+1):n]/Q_23_tilde_prime,na.rm = T))
        
        pi2.A3A1.hat.g <- (sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0]/Q_31_tilde*Y2[i,1:n0],na.rm = T)+
                             sum(I1.A3[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A1[i,(n0+1):n]/Q_31_tilde_prime*Y2[i,(n0+1):n],na.rm = T))/
          (sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0]/Q_31_tilde,na.rm = T)+
             sum(I1.A3[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A1[i,(n0+1):n]/Q_31_tilde_prime,na.rm = T))
        
        pi2.A3A2.hat.g <- (sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0]/Q_32_tilde*Y2[i,1:n0],na.rm = T)+
                             sum(I1.A3[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A2[i,(n0+1):n]/Q_32_tilde_prime*Y2[i,(n0+1):n],na.rm = T))/
          (sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0]/Q_32_tilde,na.rm = T)+
             sum(I1.A3[i,(n0+1):n]*(1-Y1[i,(n0+1):n])*I2.A2[i,(n0+1):n]/Q_32_tilde_prime,na.rm = T))
        
      } else {
        
        pi2.A1A2.hat.g <- (sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0]/Q_12_tilde*Y2[i,1:n0],na.rm = T))/(sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0]/Q_12_tilde,na.rm = T))
        pi2.A1A3.hat.g <- (sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0]/Q_13_tilde*Y2[i,1:n0],na.rm = T))/(sum(I1.A1[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0]/Q_13_tilde,na.rm = T))
        pi2.A2A1.hat.g <- (sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0]/Q_21_tilde*Y2[i,1:n0],na.rm = T))/(sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0]/Q_21_tilde,na.rm = T))
        pi2.A2A3.hat.g <- (sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0]/Q_23_tilde*Y2[i,1:n0],na.rm = T))/(sum(I1.A2[i,1:n0]*(1-Y1[i,1:n0])*I2.A3[i,1:n0]/Q_23_tilde,na.rm = T))
        pi2.A3A1.hat.g <- (sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0]/Q_31_tilde*Y2[i,1:n0],na.rm = T))/(sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A1[i,1:n0]/Q_31_tilde,na.rm = T))
        pi2.A3A2.hat.g <- (sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0]/Q_32_tilde*Y2[i,1:n0],na.rm = T))/(sum(I1.A3[i,1:n0]*(1-Y1[i,1:n0])*I2.A2[i,1:n0]/Q_32_tilde,na.rm = T))
        
        
      }
      
      mu.A1A2hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A2.hat.g
      mu.A1A3hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A3.hat.g
      mu.A2A1hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A1.hat.g
      mu.A2A3hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A3.hat.g
      mu.A3A1hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A1.hat.g
      mu.A3A2hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A2.hat.g
      
      var.A1A2hat.g[i] <- (1-pi1.A1.hat.g)*(1-pi2.A1A2.hat.g)/(sum(I1.A1[i,])/n)*(pi1.A1.hat.g*(1-pi2.A1A2.hat.g)+(p0.burn/Q_12_tilde+(1-p0.burn)/Q_12_tilde_prime)*pi2.A1A2.hat.g)/n
      var.A1A3hat.g[i] <- (1-pi1.A1.hat.g)*(1-pi2.A1A3.hat.g)/(sum(I1.A1[i,])/n)*(pi1.A1.hat.g*(1-pi2.A1A3.hat.g)+(p0.burn/Q_13_tilde+(1-p0.burn)/Q_13_tilde_prime)*pi2.A1A3.hat.g)/n
      var.A2A1hat.g[i] <- (1-pi1.A2.hat.g)*(1-pi2.A2A1.hat.g)/(sum(I1.A2[i,])/n)*(pi1.A2.hat.g*(1-pi2.A2A1.hat.g)+(p0.burn/Q_21_tilde+(1-p0.burn)/Q_21_tilde_prime)*pi2.A2A1.hat.g)/n
      var.A2A3hat.g[i] <- (1-pi1.A2.hat.g)*(1-pi2.A2A3.hat.g)/(sum(I1.A2[i,])/n)*(pi1.A2.hat.g*(1-pi2.A2A3.hat.g)+(p0.burn/Q_23_tilde+(1-p0.burn)/Q_23_tilde_prime)*pi2.A2A3.hat.g)/n
      var.A3A1hat.g[i] <- (1-pi1.A3.hat.g)*(1-pi2.A3A1.hat.g)/(sum(I1.A3[i,])/n)*(pi1.A3.hat.g*(1-pi2.A3A1.hat.g)+(p0.burn/Q_31_tilde+(1-p0.burn)/Q_31_tilde_prime)*pi2.A3A1.hat.g)/n
      var.A3A2hat.g[i] <- (1-pi1.A3.hat.g)*(1-pi2.A3A2.hat.g)/(sum(I1.A3[i,])/n)*(pi1.A3.hat.g*(1-pi2.A3A2.hat.g)+(p0.burn/Q_32_tilde+(1-p0.burn)/Q_32_tilde_prime)*pi2.A3A2.hat.g)/n
      
      
    } else if (RA.SMART==FALSE) {
      
      ######## SMART or GO-SMART design ########
      
      # most simple sample mean estimator
      mu.A1A2hat.sm[i] <- (sum(Y1[i,]*I1.A1[i,])+sum(Y[i,]*I1.A1[i,]*I2.A2[i,],na.rm = T))/(sum(Y1[i,]*I1.A1[i,])+sum(I1.A1[i,]*I2.A2[i,],na.rm = T))
      mu.A1A3hat.sm[i] <- (sum(Y1[i,]*I1.A1[i,])+sum(Y[i,]*I1.A1[i,]*I2.A3[i,],na.rm = T))/(sum(Y1[i,]*I1.A1[i,])+sum(I1.A1[i,]*I2.A3[i,],na.rm = T))
      mu.A2A1hat.sm[i] <- (sum(Y1[i,]*I1.A2[i,])+sum(Y[i,]*I1.A2[i,]*I2.A1[i,],na.rm = T))/(sum(Y1[i,]*I1.A2[i,])+sum(I1.A2[i,]*I2.A1[i,],na.rm = T))
      mu.A2A3hat.sm[i] <- (sum(Y1[i,]*I1.A2[i,])+sum(Y[i,]*I1.A2[i,]*I2.A3[i,],na.rm = T))/(sum(Y1[i,]*I1.A2[i,])+sum(I1.A2[i,]*I2.A3[i,],na.rm = T))
      mu.A3A1hat.sm[i] <- (sum(Y1[i,]*I1.A3[i,])+sum(Y[i,]*I1.A3[i,]*I2.A1[i,],na.rm = T))/(sum(Y1[i,]*I1.A3[i,])+sum(I1.A3[i,]*I2.A1[i,],na.rm = T))
      mu.A3A2hat.sm[i] <- (sum(Y1[i,]*I1.A3[i,])+sum(Y[i,]*I1.A3[i,]*I2.A2[i,],na.rm = T))/(sum(Y1[i,]*I1.A3[i,])+sum(I1.A3[i,]*I2.A2[i,],na.rm = T))
      
      var.A1A2hat.sm[i] <- mu.A1A2hat.sm[i]*(1-mu.A1A2hat.sm[i])/(sum(Y1[i,]*I1.A1[i,])+sum(I1.A1[i,]*I2.A2[i,],na.rm = T))
      var.A1A3hat.sm[i] <- mu.A1A3hat.sm[i]*(1-mu.A1A3hat.sm[i])/(sum(Y1[i,]*I1.A1[i,])+sum(I1.A1[i,]*I2.A3[i,],na.rm = T))
      var.A2A1hat.sm[i] <- mu.A2A1hat.sm[i]*(1-mu.A2A1hat.sm[i])/(sum(Y1[i,]*I1.A2[i,])+sum(I1.A2[i,]*I2.A1[i,],na.rm = T))
      var.A2A3hat.sm[i] <- mu.A2A3hat.sm[i]*(1-mu.A2A3hat.sm[i])/(sum(Y1[i,]*I1.A2[i,])+sum(I1.A2[i,]*I2.A3[i,],na.rm = T))
      var.A3A1hat.sm[i] <- mu.A3A1hat.sm[i]*(1-mu.A3A1hat.sm[i])/(sum(Y1[i,]*I1.A3[i,])+sum(I1.A3[i,]*I2.A1[i,],na.rm = T))
      var.A3A2hat.sm[i] <- mu.A3A2hat.sm[i]*(1-mu.A3A2hat.sm[i])/(sum(Y1[i,]*I1.A3[i,])+sum(I1.A3[i,]*I2.A2[i,],na.rm = T))
      
      
      # IPW estimator
      mu.A1A2hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,])))/n
      mu.A1A3hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,])))/n
      mu.A2A1hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,])))/n
      mu.A2A3hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,])))/n
      mu.A3A1hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,])))/n
      mu.A3A2hat.ipw[i] <- sum(Y[i,]*Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,])))/n
      
      
      var.A1A2hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,]))-mu.A1A2hat.ipw[i])^2)/(n^2)
      var.A1A3hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,]))-mu.A1A3hat.ipw[i])^2)/(n^2)
      var.A2A1hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,]))-mu.A2A1hat.ipw[i])^2)/(n^2)
      var.A2A3hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,]))-mu.A2A3hat.ipw[i])^2)/(n^2)
      var.A3A1hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,]))-mu.A3A1hat.ipw[i])^2)/(n^2)
      var.A3A2hat.ipw[i] <- sum((Y[i,]*Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na(Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,])), 0, Y[i,]*(1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,]))-mu.A3A2hat.ipw[i])^2)/(n^2)
      
      
      # IPWN estimator
      w.A1A2 <- Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na((1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,])), 0, (1-Y1[i,])*I1.A1[i,]*I2.A2[i,]/(P1.A1[i,]*P2.A1A2[i,]))
      w.A1A3 <- Y1[i,]*I1.A1[i,]/P1.A1[i,]+ifelse(is.na((1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,])), 0, (1-Y1[i,])*I1.A1[i,]*I2.A3[i,]/(P1.A1[i,]*P2.A1A3[i,]))
      w.A2A1 <- Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na((1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,])), 0, (1-Y1[i,])*I1.A2[i,]*I2.A1[i,]/(P1.A2[i,]*P2.A2A1[i,]))
      w.A2A3 <- Y1[i,]*I1.A2[i,]/P1.A2[i,]+ifelse(is.na((1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,])), 0, (1-Y1[i,])*I1.A2[i,]*I2.A3[i,]/(P1.A2[i,]*P2.A2A3[i,]))
      w.A3A1 <- Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na((1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,])), 0, (1-Y1[i,])*I1.A3[i,]*I2.A1[i,]/(P1.A3[i,]*P2.A3A1[i,]))
      w.A3A2 <- Y1[i,]*I1.A3[i,]/P1.A3[i,]+ifelse(is.na((1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,])), 0, (1-Y1[i,])*I1.A3[i,]*I2.A2[i,]/(P1.A3[i,]*P2.A3A2[i,]))
      
      mu.A1A2hat.ipwn[i] <- mu.A1A2hat.ipw[i]/(sum(w.A1A2)/n)
      mu.A1A3hat.ipwn[i] <- mu.A1A3hat.ipw[i]/(sum(w.A1A3)/n)
      mu.A2A1hat.ipwn[i] <- mu.A2A1hat.ipw[i]/(sum(w.A2A1)/n)
      mu.A2A3hat.ipwn[i] <- mu.A2A3hat.ipw[i]/(sum(w.A2A3)/n)
      mu.A3A1hat.ipwn[i] <- mu.A3A1hat.ipw[i]/(sum(w.A3A1)/n)
      mu.A3A2hat.ipwn[i] <- mu.A3A2hat.ipw[i]/(sum(w.A3A2)/n)
      
      var.A1A2hat.ipwn[i] <- sum((w.A1A2*(Y[i,]-mu.A1A2hat.ipwn[i]))^2)/(n^2)
      var.A1A3hat.ipwn[i] <- sum((w.A1A3*(Y[i,]-mu.A1A3hat.ipwn[i]))^2)/(n^2)
      var.A2A1hat.ipwn[i] <- sum((w.A2A1*(Y[i,]-mu.A2A1hat.ipwn[i]))^2)/(n^2)
      var.A2A3hat.ipwn[i] <- sum((w.A2A3*(Y[i,]-mu.A2A3hat.ipwn[i]))^2)/(n^2)
      var.A3A1hat.ipwn[i] <- sum((w.A3A1*(Y[i,]-mu.A3A1hat.ipwn[i]))^2)/(n^2)
      var.A3A2hat.ipwn[i] <- sum((w.A3A2*(Y[i,]-mu.A3A2hat.ipwn[i]))^2)/(n^2)
      
      
      # G-estimator
      pi1.A1.hat.g <- sum(Y1[i,]*I1.A1[i,])/sum(I1.A1[i,])
      pi1.A2.hat.g <- sum(Y1[i,]*I1.A2[i,])/sum(I1.A2[i,])
      pi1.A3.hat.g <- sum(Y1[i,]*I1.A3[i,])/sum(I1.A3[i,])
      pi2.A1A2.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,]*Y2[i,], na.rm = T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,], na.rm = T)
      pi2.A1A3.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,]*Y2[i,], na.rm = T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,], na.rm = T)
      pi2.A2A1.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,]*Y2[i,], na.rm = T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,], na.rm = T)
      pi2.A2A3.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,]*Y2[i,], na.rm = T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,], na.rm = T)
      pi2.A3A1.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,]*Y2[i,], na.rm = T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,], na.rm = T)
      pi2.A3A2.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,]*Y2[i,], na.rm = T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,], na.rm = T)
      
      mu.A1A2hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A2.hat.g
      mu.A1A3hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A3.hat.g
      mu.A2A1hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A1.hat.g
      mu.A2A3hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A3.hat.g
      mu.A3A1hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A1.hat.g
      mu.A3A2hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A2.hat.g
      
      var.A1A2hat.g[i] <- (1-pi2.A1A2.hat.g)^2*pi1.A1.hat.g*(1-pi1.A1.hat.g)/sum(I1.A1[i,])+(1-pi1.A1.hat.g)^2*pi2.A1A2.hat.g*(1-pi2.A1A2.hat.g)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,], na.rm = T)
      var.A1A3hat.g[i] <- (1-pi2.A1A3.hat.g)^2*pi1.A1.hat.g*(1-pi1.A1.hat.g)/sum(I1.A1[i,])+(1-pi1.A1.hat.g)^2*pi2.A1A3.hat.g*(1-pi2.A1A3.hat.g)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,], na.rm = T)
      var.A2A1hat.g[i] <- (1-pi2.A2A1.hat.g)^2*pi1.A2.hat.g*(1-pi1.A2.hat.g)/sum(I1.A2[i,])+(1-pi1.A2.hat.g)^2*pi2.A2A1.hat.g*(1-pi2.A2A1.hat.g)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,], na.rm = T)
      var.A2A3hat.g[i] <- (1-pi2.A2A3.hat.g)^2*pi1.A2.hat.g*(1-pi1.A2.hat.g)/sum(I1.A2[i,])+(1-pi1.A2.hat.g)^2*pi2.A2A3.hat.g*(1-pi2.A2A3.hat.g)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,], na.rm = T)
      var.A3A1hat.g[i] <- (1-pi2.A3A1.hat.g)^2*pi1.A3.hat.g*(1-pi1.A3.hat.g)/sum(I1.A3[i,])+(1-pi1.A3.hat.g)^2*pi2.A3A1.hat.g*(1-pi2.A3A1.hat.g)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,], na.rm = T)
      var.A3A2hat.g[i] <- (1-pi2.A3A2.hat.g)^2*pi1.A3.hat.g*(1-pi1.A3.hat.g)/sum(I1.A3[i,])+(1-pi1.A3.hat.g)^2*pi2.A3A2.hat.g*(1-pi2.A3A2.hat.g)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,], na.rm = T)
      
    }  
    
    # estimated optimal DTR based on IPW estimator
    est.opt.ipw[i] <- which.max(c(mu.A1A2hat.ipw[i], mu.A1A3hat.ipw[i], mu.A2A1hat.ipw[i], mu.A2A3hat.ipw[i], mu.A3A1hat.ipw[i], mu.A3A2hat.ipw[i]))
    # estimated optimal DTR based on IPWN estimator
    est.opt.ipwn[i] <- which.max(c(mu.A1A2hat.ipwn[i], mu.A1A3hat.ipwn[i], mu.A2A1hat.ipwn[i], mu.A2A3hat.ipwn[i], mu.A3A1hat.ipwn[i], mu.A3A2hat.ipwn[i]))
    # estimated optimal DTR based on G-estimator
    est.opt.g[i] <- which.max(c(mu.A1A2hat.g[i], mu.A1A3hat.g[i], mu.A2A1hat.g[i], mu.A2A3hat.g[i], mu.A3A1hat.g[i], mu.A3A2hat.g[i]))
    # estimated optimal DTR based on sample mean estimator
    est.opt.sm[i] <- which.max(c(mu.A1A2hat.sm[i], mu.A1A3hat.sm[i], mu.A2A1hat.sm[i], mu.A2A3hat.sm[i], mu.A3A1hat.sm[i], mu.A3A2hat.sm[i]))
    
    
    # estimated worst DTR based on IPW estimator
    est.wor.ipw[i] <- which.min(c(mu.A1A2hat.ipw[i], mu.A1A3hat.ipw[i], mu.A2A1hat.ipw[i], mu.A2A3hat.ipw[i], mu.A3A1hat.ipw[i], mu.A3A2hat.ipw[i]))
    # estimated worst DTR based on IPWN estimator
    est.wor.ipwn[i] <- which.min(c(mu.A1A2hat.ipwn[i], mu.A1A3hat.ipwn[i], mu.A2A1hat.ipwn[i], mu.A2A3hat.ipwn[i], mu.A3A1hat.ipwn[i], mu.A3A2hat.ipwn[i]))
    # estimated worst DTR based on G-estimator
    est.wor.g[i] <- which.min(c(mu.A1A2hat.g[i], mu.A1A3hat.g[i], mu.A2A1hat.g[i], mu.A2A3hat.g[i], mu.A3A1hat.g[i], mu.A3A2hat.g[i]))
    # estimated worst DTR based on sample mean estimator
    est.wor.sm[i] <- which.min(c(mu.A1A2hat.sm[i], mu.A1A3hat.sm[i], mu.A2A1hat.sm[i], mu.A2A3hat.sm[i], mu.A3A1hat.sm[i], mu.A3A2hat.sm[i]))
    
    # Number of patient treated in each of the 6 DTRs 
    # (response in stage I + not response in stage I and treated in stage II)
    # d(A1,A2)
    n.A1A2[i] <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(I1.A1[i,]*I2.A2[i,], na.rm = T)
    # d(A1,A3)
    n.A1A3[i] <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(I1.A1[i,]*I2.A3[i,], na.rm = T)
    # d(A2,A1)
    n.A2A1[i] <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(I1.A2[i,]*I2.A1[i,], na.rm = T)
    # d(A2,A3)
    n.A2A3[i] <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(I1.A2[i,]*I2.A3[i,], na.rm = T)
    # d(A3,A1)
    n.A3A1[i] <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(I1.A3[i,]*I2.A1[i,], na.rm = T)
    # d(A3,A2)
    n.A3A2[i] <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(I1.A3[i,]*I2.A2[i,], na.rm = T)
    
    # Number of response in each of the 6 DTRs
    # (response in stage I + not response in stage I and response in stage II)
    # d(A1,A2)
    rn.A1A2[i] <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(I1.A1[i,]* I2.A2[i,]*Y[i,], na.rm = T)
    # d(A1,A3)
    rn.A1A3[i] <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(I1.A1[i,]* I2.A3[i,]*Y[i,], na.rm = T)
    # d(A2,A1)
    rn.A2A1[i] <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(I1.A2[i,]* I2.A1[i,]*Y[i,], na.rm = T)
    # d(A2,A3)
    rn.A2A3[i] <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(I1.A2[i,]* I2.A3[i,]*Y[i,], na.rm = T)
    # d(A3,A1)
    rn.A3A1[i] <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(I1.A3[i,]* I2.A1[i,]*Y[i,], na.rm = T)
    # d(A3,A2)
    rn.A3A2[i] <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(I1.A3[i,]* I2.A2[i,]*Y[i,], na.rm = T)
    
    
    # total number of response in the whole trial
    NR[i] <- sum(Y[i,])
    
    # progress indicator
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  out <- data.frame(mu.A1A2hat.sm=mu.A1A2hat.sm, mu.A1A3hat.sm=mu.A1A3hat.sm,
                    mu.A2A1hat.sm=mu.A2A1hat.sm, mu.A2A3hat.sm=mu.A2A3hat.sm,
                    mu.A3A1hat.sm=mu.A3A1hat.sm, mu.A3A2hat.sm=mu.A3A2hat.sm,
                    var.A1A2hat.sm=var.A1A2hat.sm, var.A1A3hat.sm=var.A1A3hat.sm,
                    var.A2A1hat.sm=var.A2A1hat.sm, var.A2A3hat.sm=var.A2A3hat.sm,
                    var.A3A1hat.sm=var.A3A1hat.sm, var.A3A2hat.sm=var.A3A2hat.sm,
                    mu.A1A2hat.ipw=mu.A1A2hat.ipw, mu.A1A3hat.ipw=mu.A1A3hat.ipw,
                    mu.A2A1hat.ipw=mu.A2A1hat.ipw, mu.A2A3hat.ipw=mu.A2A3hat.ipw,
                    mu.A3A1hat.ipw=mu.A3A1hat.ipw, mu.A3A2hat.ipw=mu.A3A2hat.ipw,
                    var.A1A2hat.ipw=var.A1A2hat.ipw, var.A1A3hat.ipw=var.A1A3hat.ipw,
                    var.A2A1hat.ipw=var.A2A1hat.ipw, var.A2A3hat.ipw=var.A2A3hat.ipw,
                    var.A3A1hat.ipw=var.A3A1hat.ipw, var.A3A2hat.ipw=var.A3A2hat.ipw,
                    mu.A1A2hat.ipwn=mu.A1A2hat.ipwn, mu.A1A3hat.ipwn=mu.A1A3hat.ipwn,
                    mu.A2A1hat.ipwn=mu.A2A1hat.ipwn, mu.A2A3hat.ipwn=mu.A2A3hat.ipwn,
                    mu.A3A1hat.ipwn=mu.A3A1hat.ipwn, mu.A3A2hat.ipwn=mu.A3A2hat.ipwn,
                    var.A1A2hat.ipwn=var.A1A2hat.ipwn, var.A1A3hat.ipwn=var.A1A3hat.ipwn,
                    var.A2A1hat.ipwn=var.A2A1hat.ipwn, var.A2A3hat.ipwn=var.A2A3hat.ipwn,
                    var.A3A1hat.ipwn=var.A3A1hat.ipwn, var.A3A2hat.ipwn=var.A3A2hat.ipwn,
                    mu.A1A2hat.g=mu.A1A2hat.g, mu.A1A3hat.g=mu.A1A3hat.g,
                    mu.A2A1hat.g=mu.A2A1hat.g, mu.A2A3hat.g=mu.A2A3hat.g,
                    mu.A3A1hat.g=mu.A3A1hat.g, mu.A3A2hat.g=mu.A3A2hat.g,
                    var.A1A2hat.g=var.A1A2hat.g, var.A1A3hat.g=var.A1A3hat.g,
                    var.A2A1hat.g=var.A2A1hat.g, var.A2A3hat.g=var.A2A3hat.g,
                    var.A3A1hat.g=var.A3A1hat.g, var.A3A2hat.g=var.A3A2hat.g,
                    n.A1A2=n.A1A2, n.A1A3=n.A1A3,
                    n.A2A1=n.A2A1, n.A2A3=n.A2A3,
                    n.A3A1=n.A3A1, n.A3A2=n.A3A2,
                    rn.A1A2=rn.A1A2, rn.A1A3=rn.A1A3,
                    rn.A2A1=rn.A2A1, rn.A2A3=rn.A2A3,
                    rn.A3A1=rn.A3A1, rn.A3A2=rn.A3A2,
                    NR=NR,
                    est.opt.sm,est.opt.ipw,est.opt.ipwn, est.opt.g, est.wor.sm,est.wor.ipw,est.wor.ipwn, est.wor.g)
  
  
  write.csv(out, paste0("./output/dat_", n, "_", sce, "_", AR, "_",p0.burn, "_",p1.burn,"_Q",Q,"_c",c, ".csv"))
  
}


# Input variables for function "eva":
## df: the dataframe that contains all estimators produced by function "monte"
## pi = c(pi1.A1, pi1.A2, pi1.A3, pi2.A1A2, pi2.A1A3, pi2.A2A1, pi2.A2A3, pi2.A3A1, pi2.A3A2): vector of stage 1 and stage 2 true response rate; e.g. pi=c(0.5, 0.35, 0.2, 0.3, 0.4, 0.35, 0.2, 0.25, 0.1)

eva <- function(df,pi) {
  
  pi1.A1=pi[1]
  pi1.A2=pi[2]
  pi1.A3=pi[3]
  pi2.A1A2=pi[4]
  pi2.A1A3=pi[5]
  pi2.A2A1=pi[6]
  pi2.A2A3=pi[7]
  pi2.A3A1=pi[8]
  pi2.A3A2=pi[9]
  
  # True response rate for 6 DTRs
  # d(A1, A2)
  pid.A1A2 <- pi1.A1+(1-pi1.A1)*pi2.A1A2
  # d(A1, A3)
  pid.A1A3 <- pi1.A1+(1-pi1.A1)*pi2.A1A3
  # d(A2, A1)
  pid.A2A1 <- pi1.A2+(1-pi1.A2)*pi2.A2A1
  # d(A2, A3)
  pid.A2A3 <- pi1.A2+(1-pi1.A2)*pi2.A2A3
  # d(A3, A1)
  pid.A3A1 <- pi1.A3+(1-pi1.A3)*pi2.A3A1
  # d(A3, A2)
  pid.A3A2 <- pi1.A3+(1-pi1.A3)*pi2.A3A2
  
  # simple sample mean estimator
  mu.A1A2hat.sm <- df$mu.A1A2hat.sm
  mu.A1A3hat.sm <- df$mu.A1A3hat.sm
  mu.A2A1hat.sm <- df$mu.A2A1hat.sm
  mu.A2A3hat.sm <- df$mu.A2A3hat.sm
  mu.A3A1hat.sm <- df$mu.A3A1hat.sm
  mu.A3A2hat.sm <- df$mu.A3A2hat.sm
  
  var.A1A2hat.sm <- df$var.A1A2hat.sm
  var.A1A3hat.sm <- df$var.A1A3hat.sm
  var.A2A1hat.sm <- df$var.A2A1hat.sm
  var.A2A3hat.sm <- df$var.A2A3hat.sm
  var.A3A1hat.sm <- df$var.A3A1hat.sm
  var.A3A2hat.sm <- df$var.A3A2hat.sm
  
  # IPW estimator
  mu.A1A2hat.ipw <- df$mu.A1A2hat.ipw
  mu.A1A3hat.ipw <- df$mu.A1A3hat.ipw
  mu.A2A1hat.ipw <- df$mu.A2A1hat.ipw
  mu.A2A3hat.ipw <- df$mu.A2A3hat.ipw
  mu.A3A1hat.ipw <- df$mu.A3A1hat.ipw
  mu.A3A2hat.ipw <- df$mu.A3A2hat.ipw
  
  var.A1A2hat.ipw <- df$var.A1A2hat.ipw
  var.A1A3hat.ipw <- df$var.A1A3hat.ipw
  var.A2A1hat.ipw <- df$var.A2A1hat.ipw
  var.A2A3hat.ipw <- df$var.A2A3hat.ipw
  var.A3A1hat.ipw <- df$var.A3A1hat.ipw
  var.A3A2hat.ipw <- df$var.A3A2hat.ipw
  
  # IPWN estimator
  mu.A1A2hat.ipwn <- df$mu.A1A2hat.ipwn
  mu.A1A3hat.ipwn <- df$mu.A1A3hat.ipwn
  mu.A2A1hat.ipwn <- df$mu.A2A1hat.ipwn
  mu.A2A3hat.ipwn <- df$mu.A2A3hat.ipwn
  mu.A3A1hat.ipwn <- df$mu.A3A1hat.ipwn
  mu.A3A2hat.ipwn <- df$mu.A3A2hat.ipwn
  
  var.A1A2hat.ipwn <- df$var.A1A2hat.ipwn
  var.A1A3hat.ipwn <- df$var.A1A3hat.ipwn
  var.A2A1hat.ipwn <- df$var.A2A1hat.ipwn
  var.A2A3hat.ipwn <- df$var.A2A3hat.ipwn
  var.A3A1hat.ipwn <- df$var.A3A1hat.ipwn
  var.A3A2hat.ipwn <- df$var.A3A2hat.ipwn
  
  # G-estimator
  mu.A1A2hat.g <- df$mu.A1A2hat.g
  mu.A1A3hat.g <- df$mu.A1A3hat.g
  mu.A2A1hat.g <- df$mu.A2A1hat.g
  mu.A2A3hat.g <- df$mu.A2A3hat.g
  mu.A3A1hat.g <- df$mu.A3A1hat.g
  mu.A3A2hat.g <- df$mu.A3A2hat.g
  
  var.A1A2hat.g <- df$var.A1A2hat.g
  var.A1A3hat.g <- df$var.A1A3hat.g
  var.A2A1hat.g <- df$var.A2A1hat.g
  var.A2A3hat.g <- df$var.A2A3hat.g
  var.A3A1hat.g <- df$var.A3A1hat.g
  var.A3A2hat.g <- df$var.A3A2hat.g
  
  # estimated optimal/worst DTR
  est.opt.ipw <- df$est.opt.ipw
  est.opt.ipwn <- df$est.opt.ipwn
  est.opt.g <- df$est.opt.g
  est.opt.sm <- df$est.opt.sm
  est.wor.ipw <- df$est.wor.ipw 
  est.wor.ipwn <- df$est.wor.ipwn 
  est.wor.g <- df$est.wor.g
  est.wor.sm <- df$est.wor.sm
  
  # Number of patient treated in each of the 6 DTRs
  n.A1A2 <- df$n.A1A2
  n.A1A3 <- df$n.A1A3
  n.A2A1 <- df$n.A2A1
  n.A2A3 <- df$n.A2A3
  n.A3A1 <- df$n.A3A1
  n.A3A2 <- df$n.A3A2
  
  # Number of patient response in each of the 6 DTRs
  rn.A1A2 <- df$rn.A1A2
  rn.A1A3 <- df$rn.A1A3
  rn.A2A1 <- df$rn.A2A1
  rn.A2A3 <- df$rn.A2A3
  rn.A3A1 <- df$rn.A3A1
  rn.A3A2 <- df$rn.A3A2
  
  # Number if response in the trial
  nr <- df$NR
  
  result <- list()
  
  
  
  ################## 1. bias #########################
  
  bias.mu.A1A2hat.sm <- mean(mu.A1A2hat.sm)-pid.A1A2
  bias.mu.A1A3hat.sm <- mean(mu.A1A3hat.sm)-pid.A1A3
  bias.mu.A2A1hat.sm <- mean(mu.A2A1hat.sm)-pid.A2A1
  bias.mu.A2A3hat.sm <- mean(mu.A2A3hat.sm)-pid.A2A3
  bias.mu.A3A1hat.sm <- mean(mu.A3A1hat.sm)-pid.A3A1
  bias.mu.A3A2hat.sm <- mean(mu.A3A2hat.sm)-pid.A3A2
  
  bias.mu.A1A2hat.ipw <- mean(mu.A1A2hat.ipw)-pid.A1A2
  bias.mu.A1A3hat.ipw <- mean(mu.A1A3hat.ipw)-pid.A1A3
  bias.mu.A2A1hat.ipw <- mean(mu.A2A1hat.ipw)-pid.A2A1
  bias.mu.A2A3hat.ipw <- mean(mu.A2A3hat.ipw)-pid.A2A3
  bias.mu.A3A1hat.ipw <- mean(mu.A3A1hat.ipw)-pid.A3A1
  bias.mu.A3A2hat.ipw <- mean(mu.A3A2hat.ipw)-pid.A3A2
  
  bias.mu.A1A2hat.ipwn <- mean(mu.A1A2hat.ipwn)-pid.A1A2
  bias.mu.A1A3hat.ipwn <- mean(mu.A1A3hat.ipwn)-pid.A1A3
  bias.mu.A2A1hat.ipwn <- mean(mu.A2A1hat.ipwn)-pid.A2A1
  bias.mu.A2A3hat.ipwn <- mean(mu.A2A3hat.ipwn)-pid.A2A3
  bias.mu.A3A1hat.ipwn <- mean(mu.A3A1hat.ipwn)-pid.A3A1
  bias.mu.A3A2hat.ipwn <- mean(mu.A3A2hat.ipwn)-pid.A3A2
  
  bias.mu.A1A2hat.g <- mean(mu.A1A2hat.g)-pid.A1A2
  bias.mu.A1A3hat.g <- mean(mu.A1A3hat.g)-pid.A1A3
  bias.mu.A2A1hat.g <- mean(mu.A2A1hat.g)-pid.A2A1
  bias.mu.A2A3hat.g <- mean(mu.A2A3hat.g)-pid.A2A3
  bias.mu.A3A1hat.g <- mean(mu.A3A1hat.g)-pid.A3A1
  bias.mu.A3A2hat.g <- mean(mu.A3A2hat.g)-pid.A3A2
  
  result$bias.mu.A1A2hat.sm <- round(bias.mu.A1A2hat.sm,3)
  result$bias.mu.A1A3hat.sm <- round(bias.mu.A1A3hat.sm,3)
  result$bias.mu.A2A1hat.sm <- round(bias.mu.A2A1hat.sm,3)
  result$bias.mu.A2A3hat.sm <- round(bias.mu.A2A3hat.sm,3)
  result$bias.mu.A3A1hat.sm <- round(bias.mu.A3A1hat.sm,3)
  result$bias.mu.A3A2hat.sm <- round(bias.mu.A3A2hat.sm,3)
  
  result$bias.mu.A1A2hat.ipw <- round(bias.mu.A1A2hat.ipw,3)
  result$bias.mu.A1A3hat.ipw <- round(bias.mu.A1A3hat.ipw,3)
  result$bias.mu.A2A1hat.ipw <- round(bias.mu.A2A1hat.ipw,3)
  result$bias.mu.A2A3hat.ipw <- round(bias.mu.A2A3hat.ipw,3)
  result$bias.mu.A3A1hat.ipw <- round(bias.mu.A3A1hat.ipw,3)
  result$bias.mu.A3A2hat.ipw <- round(bias.mu.A3A2hat.ipw,3)
  
  result$bias.mu.A1A2hat.ipwn <- round(bias.mu.A1A2hat.ipwn,3)
  result$bias.mu.A1A3hat.ipwn <- round(bias.mu.A1A3hat.ipwn,3)
  result$bias.mu.A2A1hat.ipwn <- round(bias.mu.A2A1hat.ipwn,3)
  result$bias.mu.A2A3hat.ipwn <- round(bias.mu.A2A3hat.ipwn,3)
  result$bias.mu.A3A1hat.ipwn <- round(bias.mu.A3A1hat.ipwn,3)
  result$bias.mu.A3A2hat.ipwn <- round(bias.mu.A3A2hat.ipwn,3)
  
  result$bias.mu.A1A2hat.g <- round(bias.mu.A1A2hat.g,3)
  result$bias.mu.A1A3hat.g <- round(bias.mu.A1A3hat.g,3)
  result$bias.mu.A2A1hat.g <- round(bias.mu.A2A1hat.g,3)
  result$bias.mu.A2A3hat.g <- round(bias.mu.A2A3hat.g,3)
  result$bias.mu.A3A1hat.g <- round(bias.mu.A3A1hat.g,3)
  result$bias.mu.A3A2hat.g <- round(bias.mu.A3A2hat.g,3)
  
  
  ################## 2. Variance ratio #########################
  
  monsd.mu.A1A2hat.sm <- sd(mu.A1A2hat.sm)
  monsd.mu.A1A3hat.sm <- sd(mu.A1A3hat.sm)
  monsd.mu.A2A1hat.sm <- sd(mu.A2A1hat.sm)
  monsd.mu.A2A3hat.sm <- sd(mu.A2A3hat.sm)
  monsd.mu.A3A1hat.sm <- sd(mu.A3A1hat.sm)
  monsd.mu.A3A2hat.sm <- sd(mu.A3A2hat.sm)
  
  avesd.mu.A1A2hat.sm <- mean(sqrt(var.A1A2hat.sm))
  avesd.mu.A1A3hat.sm <- mean(sqrt(var.A1A3hat.sm))
  avesd.mu.A2A1hat.sm <- mean(sqrt(var.A2A1hat.sm))
  avesd.mu.A2A3hat.sm <- mean(sqrt(var.A2A3hat.sm))
  avesd.mu.A3A1hat.sm <- mean(sqrt(var.A3A1hat.sm))
  avesd.mu.A3A2hat.sm <- mean(sqrt(var.A3A2hat.sm))
  
  ratio.mu.A1A2hat.sm <- avesd.mu.A1A2hat.sm/monsd.mu.A1A2hat.sm
  ratio.mu.A1A3hat.sm <- avesd.mu.A1A3hat.sm/monsd.mu.A1A3hat.sm
  ratio.mu.A2A1hat.sm <- avesd.mu.A2A1hat.sm/monsd.mu.A2A1hat.sm
  ratio.mu.A2A3hat.sm <- avesd.mu.A2A3hat.sm/monsd.mu.A2A3hat.sm
  ratio.mu.A3A1hat.sm <- avesd.mu.A3A1hat.sm/monsd.mu.A3A1hat.sm
  ratio.mu.A3A2hat.sm <- avesd.mu.A3A2hat.sm/monsd.mu.A3A2hat.sm
  
  monsd.mu.A1A2hat.ipw <- sd(mu.A1A2hat.ipw)
  monsd.mu.A1A3hat.ipw <- sd(mu.A1A3hat.ipw)
  monsd.mu.A2A1hat.ipw <- sd(mu.A2A1hat.ipw)
  monsd.mu.A2A3hat.ipw <- sd(mu.A2A3hat.ipw)
  monsd.mu.A3A1hat.ipw <- sd(mu.A3A1hat.ipw)
  monsd.mu.A3A2hat.ipw <- sd(mu.A3A2hat.ipw)
  
  avesd.mu.A1A2hat.ipw <- mean(sqrt(var.A1A2hat.ipw))
  avesd.mu.A1A3hat.ipw <- mean(sqrt(var.A1A3hat.ipw))
  avesd.mu.A2A1hat.ipw <- mean(sqrt(var.A2A1hat.ipw))
  avesd.mu.A2A3hat.ipw <- mean(sqrt(var.A2A3hat.ipw))
  avesd.mu.A3A1hat.ipw <- mean(sqrt(var.A3A1hat.ipw))
  avesd.mu.A3A2hat.ipw <- mean(sqrt(var.A3A2hat.ipw))
  
  ratio.mu.A1A2hat.ipw <- avesd.mu.A1A2hat.ipw/monsd.mu.A1A2hat.ipw
  ratio.mu.A1A3hat.ipw <- avesd.mu.A1A3hat.ipw/monsd.mu.A1A3hat.ipw
  ratio.mu.A2A1hat.ipw <- avesd.mu.A2A1hat.ipw/monsd.mu.A2A1hat.ipw
  ratio.mu.A2A3hat.ipw <- avesd.mu.A2A3hat.ipw/monsd.mu.A2A3hat.ipw
  ratio.mu.A3A1hat.ipw <- avesd.mu.A3A1hat.ipw/monsd.mu.A3A1hat.ipw
  ratio.mu.A3A2hat.ipw <- avesd.mu.A3A2hat.ipw/monsd.mu.A3A2hat.ipw
  
  monsd.mu.A1A2hat.ipwn <- sd(mu.A1A2hat.ipwn)
  monsd.mu.A1A3hat.ipwn <- sd(mu.A1A3hat.ipwn)
  monsd.mu.A2A1hat.ipwn <- sd(mu.A2A1hat.ipwn)
  monsd.mu.A2A3hat.ipwn <- sd(mu.A2A3hat.ipwn)
  monsd.mu.A3A1hat.ipwn <- sd(mu.A3A1hat.ipwn)
  monsd.mu.A3A2hat.ipwn <- sd(mu.A3A2hat.ipwn)
  
  avesd.mu.A1A2hat.ipwn <- mean(sqrt(var.A1A2hat.ipwn))
  avesd.mu.A1A3hat.ipwn <- mean(sqrt(var.A1A3hat.ipwn))
  avesd.mu.A2A1hat.ipwn <- mean(sqrt(var.A2A1hat.ipwn))
  avesd.mu.A2A3hat.ipwn <- mean(sqrt(var.A2A3hat.ipwn))
  avesd.mu.A3A1hat.ipwn <- mean(sqrt(var.A3A1hat.ipwn))
  avesd.mu.A3A2hat.ipwn <- mean(sqrt(var.A3A2hat.ipwn))
  
  ratio.mu.A1A2hat.ipwn <- avesd.mu.A1A2hat.ipwn/monsd.mu.A1A2hat.ipwn
  ratio.mu.A1A3hat.ipwn <- avesd.mu.A1A3hat.ipwn/monsd.mu.A1A3hat.ipwn
  ratio.mu.A2A1hat.ipwn <- avesd.mu.A2A1hat.ipwn/monsd.mu.A2A1hat.ipwn
  ratio.mu.A2A3hat.ipwn <- avesd.mu.A2A3hat.ipwn/monsd.mu.A2A3hat.ipwn
  ratio.mu.A3A1hat.ipwn <- avesd.mu.A3A1hat.ipwn/monsd.mu.A3A1hat.ipwn
  ratio.mu.A3A2hat.ipwn <- avesd.mu.A3A2hat.ipwn/monsd.mu.A3A2hat.ipwn
  
  monsd.mu.A1A2hat.g <- sd(mu.A1A2hat.g)
  monsd.mu.A1A3hat.g <- sd(mu.A1A3hat.g)
  monsd.mu.A2A1hat.g <- sd(mu.A2A1hat.g)
  monsd.mu.A2A3hat.g <- sd(mu.A2A3hat.g)
  monsd.mu.A3A1hat.g <- sd(mu.A3A1hat.g)
  monsd.mu.A3A2hat.g <- sd(mu.A3A2hat.g)
  
  avesd.mu.A1A2hat.g <- mean(sqrt(var.A1A2hat.g))
  avesd.mu.A1A3hat.g <- mean(sqrt(var.A1A3hat.g))
  avesd.mu.A2A1hat.g <- mean(sqrt(var.A2A1hat.g))
  avesd.mu.A2A3hat.g <- mean(sqrt(var.A2A3hat.g))
  avesd.mu.A3A1hat.g <- mean(sqrt(var.A3A1hat.g))
  avesd.mu.A3A2hat.g <- mean(sqrt(var.A3A2hat.g))
  
  ratio.mu.A1A2hat.g <- avesd.mu.A1A2hat.g/monsd.mu.A1A2hat.g
  ratio.mu.A1A3hat.g <- avesd.mu.A1A3hat.g/monsd.mu.A1A3hat.g
  ratio.mu.A2A1hat.g <- avesd.mu.A2A1hat.g/monsd.mu.A2A1hat.g
  ratio.mu.A2A3hat.g <- avesd.mu.A2A3hat.g/monsd.mu.A2A3hat.g
  ratio.mu.A3A1hat.g <- avesd.mu.A3A1hat.g/monsd.mu.A3A1hat.g
  ratio.mu.A3A2hat.g <- avesd.mu.A3A2hat.g/monsd.mu.A3A2hat.g
  
  result$monsd.mu.A1A2hat.sm <- round(monsd.mu.A1A2hat.sm*100,3)
  result$monsd.mu.A1A3hat.sm <- round(monsd.mu.A1A3hat.sm*100,3)
  result$monsd.mu.A2A1hat.sm <- round(monsd.mu.A2A1hat.sm*100,3)
  result$monsd.mu.A2A3hat.sm <- round(monsd.mu.A2A3hat.sm*100,3)
  result$monsd.mu.A3A1hat.sm <- round(monsd.mu.A3A1hat.sm*100,3)
  result$monsd.mu.A3A2hat.sm <- round(monsd.mu.A3A2hat.sm*100,3)
  
  result$avesd.mu.A1A2hat.sm <- round(avesd.mu.A1A2hat.sm*100,3)
  result$avesd.mu.A1A3hat.sm <- round(avesd.mu.A1A3hat.sm*100,3)
  result$avesd.mu.A2A1hat.sm <- round(avesd.mu.A2A1hat.sm*100,3)
  result$avesd.mu.A2A3hat.sm <- round(avesd.mu.A2A3hat.sm*100,3)
  result$avesd.mu.A3A1hat.sm <- round(avesd.mu.A3A1hat.sm*100,3)
  result$avesd.mu.A3A2hat.sm <- round(avesd.mu.A3A2hat.sm*100,3)
  
  result$ratio.mu.A1A2hat.sm <- round(ratio.mu.A1A2hat.sm,3)
  result$ratio.mu.A1A3hat.sm <- round(ratio.mu.A1A3hat.sm,3)
  result$ratio.mu.A2A1hat.sm <- round(ratio.mu.A2A1hat.sm,3)
  result$ratio.mu.A2A3hat.sm <- round(ratio.mu.A2A3hat.sm,3)
  result$ratio.mu.A3A1hat.sm <- round(ratio.mu.A3A1hat.sm,3)
  result$ratio.mu.A3A2hat.sm <- round(ratio.mu.A3A2hat.sm,3)
  
  result$monsd.mu.A1A2hat.ipw <- round(monsd.mu.A1A2hat.ipw*100,3)
  result$monsd.mu.A1A3hat.ipw <- round(monsd.mu.A1A3hat.ipw*100,3)
  result$monsd.mu.A2A1hat.ipw <- round(monsd.mu.A2A1hat.ipw*100,3)
  result$monsd.mu.A2A3hat.ipw <- round(monsd.mu.A2A3hat.ipw*100,3)
  result$monsd.mu.A3A1hat.ipw <- round(monsd.mu.A3A1hat.ipw*100,3)
  result$monsd.mu.A3A2hat.ipw <- round(monsd.mu.A3A2hat.ipw*100,3)
  
  result$avesd.mu.A1A2hat.ipw <- round(avesd.mu.A1A2hat.ipw*100,3)
  result$avesd.mu.A1A3hat.ipw <- round(avesd.mu.A1A3hat.ipw*100,3)
  result$avesd.mu.A2A1hat.ipw <- round(avesd.mu.A2A1hat.ipw*100,3)
  result$avesd.mu.A2A3hat.ipw <- round(avesd.mu.A2A3hat.ipw*100,3)
  result$avesd.mu.A3A1hat.ipw <- round(avesd.mu.A3A1hat.ipw*100,3)
  result$avesd.mu.A3A2hat.ipw <- round(avesd.mu.A3A2hat.ipw*100,3)
  
  result$ratio.mu.A1A2hat.ipw <- round(ratio.mu.A1A2hat.ipw,3)
  result$ratio.mu.A1A3hat.ipw <- round(ratio.mu.A1A3hat.ipw,3)
  result$ratio.mu.A2A1hat.ipw <- round(ratio.mu.A2A1hat.ipw,3)
  result$ratio.mu.A2A3hat.ipw <- round(ratio.mu.A2A3hat.ipw,3)
  result$ratio.mu.A3A1hat.ipw <- round(ratio.mu.A3A1hat.ipw,3)
  result$ratio.mu.A3A2hat.ipw <- round(ratio.mu.A3A2hat.ipw,3)
  
  result$monsd.mu.A1A2hat.ipwn <- round(monsd.mu.A1A2hat.ipwn*100,3)
  result$monsd.mu.A1A3hat.ipwn <- round(monsd.mu.A1A3hat.ipwn*100,3)
  result$monsd.mu.A2A1hat.ipwn <- round(monsd.mu.A2A1hat.ipwn*100,3)
  result$monsd.mu.A2A3hat.ipwn <- round(monsd.mu.A2A3hat.ipwn*100,3)
  result$monsd.mu.A3A1hat.ipwn <- round(monsd.mu.A3A1hat.ipwn*100,3)
  result$monsd.mu.A3A2hat.ipwn <- round(monsd.mu.A3A2hat.ipwn*100,3)
  
  result$avesd.mu.A1A2hat.ipwn <- round(avesd.mu.A1A2hat.ipwn*100,3)
  result$avesd.mu.A1A3hat.ipwn <- round(avesd.mu.A1A3hat.ipwn*100,3)
  result$avesd.mu.A2A1hat.ipwn <- round(avesd.mu.A2A1hat.ipwn*100,3)
  result$avesd.mu.A2A3hat.ipwn <- round(avesd.mu.A2A3hat.ipwn*100,3)
  result$avesd.mu.A3A1hat.ipwn <- round(avesd.mu.A3A1hat.ipwn*100,3)
  result$avesd.mu.A3A2hat.ipwn <- round(avesd.mu.A3A2hat.ipwn*100,3)
  
  result$ratio.mu.A1A2hat.ipwn <- round(ratio.mu.A1A2hat.ipwn,3)
  result$ratio.mu.A1A3hat.ipwn <- round(ratio.mu.A1A3hat.ipwn,3)
  result$ratio.mu.A2A1hat.ipwn <- round(ratio.mu.A2A1hat.ipwn,3)
  result$ratio.mu.A2A3hat.ipwn <- round(ratio.mu.A2A3hat.ipwn,3)
  result$ratio.mu.A3A1hat.ipwn <- round(ratio.mu.A3A1hat.ipwn,3)
  result$ratio.mu.A3A2hat.ipwn <- round(ratio.mu.A3A2hat.ipwn,3)
  
  result$monsd.mu.A1A2hat.g <- round(monsd.mu.A1A2hat.g*100,3)
  result$monsd.mu.A1A3hat.g <- round(monsd.mu.A1A3hat.g*100,3)
  result$monsd.mu.A2A1hat.g <- round(monsd.mu.A2A1hat.g*100,3)
  result$monsd.mu.A2A3hat.g <- round(monsd.mu.A2A3hat.g*100,3)
  result$monsd.mu.A3A1hat.g <- round(monsd.mu.A3A1hat.g*100,3)
  result$monsd.mu.A3A2hat.g <- round(monsd.mu.A3A2hat.g*100,3)
  
  result$avesd.mu.A1A2hat.g <- round(avesd.mu.A1A2hat.g*100,3)
  result$avesd.mu.A1A3hat.g <- round(avesd.mu.A1A3hat.g*100,3)
  result$avesd.mu.A2A1hat.g <- round(avesd.mu.A2A1hat.g*100,3)
  result$avesd.mu.A2A3hat.g <- round(avesd.mu.A2A3hat.g*100,3)
  result$avesd.mu.A3A1hat.g <- round(avesd.mu.A3A1hat.g*100,3)
  result$avesd.mu.A3A2hat.g <- round(avesd.mu.A3A2hat.g*100,3)
  
  result$ratio.mu.A1A2hat.g <- round(ratio.mu.A1A2hat.g,3)
  result$ratio.mu.A1A3hat.g <- round(ratio.mu.A1A3hat.g,3)
  result$ratio.mu.A2A1hat.g <- round(ratio.mu.A2A1hat.g,3)
  result$ratio.mu.A2A3hat.g <- round(ratio.mu.A2A3hat.g,3)
  result$ratio.mu.A3A1hat.g <- round(ratio.mu.A3A1hat.g,3)
  result$ratio.mu.A3A2hat.g <- round(ratio.mu.A3A2hat.g,3)
  
  ################## 3. Coverage probability ######################
  
  # 95% Wald type CI
  
  # simple sample mean estimator
  lo.A1A2.sm <- mu.A1A2hat.sm-1.96*sqrt(var.A1A2hat.sm)
  up.A1A2.sm <- mu.A1A2hat.sm+1.96*sqrt(var.A1A2hat.sm)
  lo.A1A3.sm <- mu.A1A3hat.sm-1.96*sqrt(var.A1A3hat.sm)
  up.A1A3.sm <- mu.A1A3hat.sm+1.96*sqrt(var.A1A3hat.sm)
  lo.A2A1.sm <- mu.A2A1hat.sm-1.96*sqrt(var.A2A1hat.sm)
  up.A2A1.sm <- mu.A2A1hat.sm+1.96*sqrt(var.A2A1hat.sm)
  lo.A2A3.sm <- mu.A2A3hat.sm-1.96*sqrt(var.A2A3hat.sm)
  up.A2A3.sm <- mu.A2A3hat.sm+1.96*sqrt(var.A2A3hat.sm)
  lo.A3A1.sm <- mu.A3A1hat.sm-1.96*sqrt(var.A3A1hat.sm)
  up.A3A1.sm <- mu.A3A1hat.sm+1.96*sqrt(var.A3A1hat.sm)
  lo.A3A2.sm <- mu.A3A2hat.sm-1.96*sqrt(var.A3A2hat.sm)
  up.A3A2.sm <- mu.A3A2hat.sm+1.96*sqrt(var.A3A2hat.sm)
  
  cp.A1A2.sm <- mean((lo.A1A2.sm<=pid.A1A2)*(up.A1A2.sm>=pid.A1A2))
  cp.A1A3.sm <- mean((lo.A1A3.sm<=pid.A1A3)*(up.A1A3.sm>=pid.A1A3))
  cp.A2A1.sm <- mean((lo.A2A1.sm<=pid.A2A1)*(up.A2A1.sm>=pid.A2A1))
  cp.A2A3.sm <- mean((lo.A2A3.sm<=pid.A2A3)*(up.A2A3.sm>=pid.A2A3))
  cp.A3A1.sm <- mean((lo.A3A1.sm<=pid.A3A1)*(up.A3A1.sm>=pid.A3A1))
  cp.A3A2.sm <- mean((lo.A3A2.sm<=pid.A3A2)*(up.A3A2.sm>=pid.A3A2))
  
  # IPW estimator
  lo.A1A2.ipw <- mu.A1A2hat.ipw-1.96*sqrt(var.A1A2hat.ipw)
  up.A1A2.ipw <- mu.A1A2hat.ipw+1.96*sqrt(var.A1A2hat.ipw)
  lo.A1A3.ipw <- mu.A1A3hat.ipw-1.96*sqrt(var.A1A3hat.ipw)
  up.A1A3.ipw <- mu.A1A3hat.ipw+1.96*sqrt(var.A1A3hat.ipw)
  lo.A2A1.ipw <- mu.A2A1hat.ipw-1.96*sqrt(var.A2A1hat.ipw)
  up.A2A1.ipw <- mu.A2A1hat.ipw+1.96*sqrt(var.A2A1hat.ipw)
  lo.A2A3.ipw <- mu.A2A3hat.ipw-1.96*sqrt(var.A2A3hat.ipw)
  up.A2A3.ipw <- mu.A2A3hat.ipw+1.96*sqrt(var.A2A3hat.ipw)
  lo.A3A1.ipw <- mu.A3A1hat.ipw-1.96*sqrt(var.A3A1hat.ipw)
  up.A3A1.ipw <- mu.A3A1hat.ipw+1.96*sqrt(var.A3A1hat.ipw)
  lo.A3A2.ipw <- mu.A3A2hat.ipw-1.96*sqrt(var.A3A2hat.ipw)
  up.A3A2.ipw <- mu.A3A2hat.ipw+1.96*sqrt(var.A3A2hat.ipw)
  
  cp.A1A2.ipw <- mean((lo.A1A2.ipw<=pid.A1A2)*(up.A1A2.ipw>=pid.A1A2))
  cp.A1A3.ipw <- mean((lo.A1A3.ipw<=pid.A1A3)*(up.A1A3.ipw>=pid.A1A3))
  cp.A2A1.ipw <- mean((lo.A2A1.ipw<=pid.A2A1)*(up.A2A1.ipw>=pid.A2A1))
  cp.A2A3.ipw <- mean((lo.A2A3.ipw<=pid.A2A3)*(up.A2A3.ipw>=pid.A2A3))
  cp.A3A1.ipw <- mean((lo.A3A1.ipw<=pid.A3A1)*(up.A3A1.ipw>=pid.A3A1))
  cp.A3A2.ipw <- mean((lo.A3A2.ipw<=pid.A3A2)*(up.A3A2.ipw>=pid.A3A2))
  
  # IPWN estimator
  lo.A1A2.ipwn <- mu.A1A2hat.ipwn-1.96*sqrt(var.A1A2hat.ipwn)
  up.A1A2.ipwn <- mu.A1A2hat.ipwn+1.96*sqrt(var.A1A2hat.ipwn)
  lo.A1A3.ipwn <- mu.A1A3hat.ipwn-1.96*sqrt(var.A1A3hat.ipwn)
  up.A1A3.ipwn <- mu.A1A3hat.ipwn+1.96*sqrt(var.A1A3hat.ipwn)
  lo.A2A1.ipwn <- mu.A2A1hat.ipwn-1.96*sqrt(var.A2A1hat.ipwn)
  up.A2A1.ipwn <- mu.A2A1hat.ipwn+1.96*sqrt(var.A2A1hat.ipwn)
  lo.A2A3.ipwn <- mu.A2A3hat.ipwn-1.96*sqrt(var.A2A3hat.ipwn)
  up.A2A3.ipwn <- mu.A2A3hat.ipwn+1.96*sqrt(var.A2A3hat.ipwn)
  lo.A3A1.ipwn <- mu.A3A1hat.ipwn-1.96*sqrt(var.A3A1hat.ipwn)
  up.A3A1.ipwn <- mu.A3A1hat.ipwn+1.96*sqrt(var.A3A1hat.ipwn)
  lo.A3A2.ipwn <- mu.A3A2hat.ipwn-1.96*sqrt(var.A3A2hat.ipwn)
  up.A3A2.ipwn <- mu.A3A2hat.ipwn+1.96*sqrt(var.A3A2hat.ipwn)
  
  cp.A1A2.ipwn <- mean((lo.A1A2.ipwn<=pid.A1A2)*(up.A1A2.ipwn>=pid.A1A2))
  cp.A1A3.ipwn <- mean((lo.A1A3.ipwn<=pid.A1A3)*(up.A1A3.ipwn>=pid.A1A3))
  cp.A2A1.ipwn <- mean((lo.A2A1.ipwn<=pid.A2A1)*(up.A2A1.ipwn>=pid.A2A1))
  cp.A2A3.ipwn <- mean((lo.A2A3.ipwn<=pid.A2A3)*(up.A2A3.ipwn>=pid.A2A3))
  cp.A3A1.ipwn <- mean((lo.A3A1.ipwn<=pid.A3A1)*(up.A3A1.ipwn>=pid.A3A1))
  cp.A3A2.ipwn <- mean((lo.A3A2.ipwn<=pid.A3A2)*(up.A3A2.ipwn>=pid.A3A2))
  
  # G-estimator
  lo.A1A2.g <- mu.A1A2hat.g-1.96*sqrt(var.A1A2hat.g)
  up.A1A2.g <- mu.A1A2hat.g+1.96*sqrt(var.A1A2hat.g)
  lo.A1A3.g <- mu.A1A3hat.g-1.96*sqrt(var.A1A3hat.g)
  up.A1A3.g <- mu.A1A3hat.g+1.96*sqrt(var.A1A3hat.g)
  lo.A2A1.g <- mu.A2A1hat.g-1.96*sqrt(var.A2A1hat.g)
  up.A2A1.g <- mu.A2A1hat.g+1.96*sqrt(var.A2A1hat.g)
  lo.A2A3.g <- mu.A2A3hat.g-1.96*sqrt(var.A2A3hat.g)
  up.A2A3.g <- mu.A2A3hat.g+1.96*sqrt(var.A2A3hat.g)
  lo.A3A1.g <- mu.A3A1hat.g-1.96*sqrt(var.A3A1hat.g)
  up.A3A1.g <- mu.A3A1hat.g+1.96*sqrt(var.A3A1hat.g)
  lo.A3A2.g <- mu.A3A2hat.g-1.96*sqrt(var.A3A2hat.g)
  up.A3A2.g <- mu.A3A2hat.g+1.96*sqrt(var.A3A2hat.g)
  
  cp.A1A2.g <- mean((lo.A1A2.g<=pid.A1A2)*(up.A1A2.g>=pid.A1A2))
  cp.A1A3.g <- mean((lo.A1A3.g<=pid.A1A3)*(up.A1A3.g>=pid.A1A3))
  cp.A2A1.g <- mean((lo.A2A1.g<=pid.A2A1)*(up.A2A1.g>=pid.A2A1))
  cp.A2A3.g <- mean((lo.A2A3.g<=pid.A2A3)*(up.A2A3.g>=pid.A2A3))
  cp.A3A1.g <- mean((lo.A3A1.g<=pid.A3A1)*(up.A3A1.g>=pid.A3A1))
  cp.A3A2.g <- mean((lo.A3A2.g<=pid.A3A2)*(up.A3A2.g>=pid.A3A2))
  
  result$cp.A1A2.sm <- round(cp.A1A2.sm,4)
  result$cp.A1A3.sm <- round(cp.A1A3.sm,4)
  result$cp.A2A1.sm <- round(cp.A2A1.sm,4)
  result$cp.A2A3.sm <- round(cp.A2A3.sm,4)
  result$cp.A3A1.sm <- round(cp.A3A1.sm,4)
  result$cp.A3A2.sm <- round(cp.A3A2.sm,4)
  
  result$cp.A1A2.ipw <- round(cp.A1A2.ipw,4)
  result$cp.A1A3.ipw <- round(cp.A1A3.ipw,4)
  result$cp.A2A1.ipw <- round(cp.A2A1.ipw,4)
  result$cp.A2A3.ipw <- round(cp.A2A3.ipw,4)
  result$cp.A3A1.ipw <- round(cp.A3A1.ipw,4)
  result$cp.A3A2.ipw <- round(cp.A3A2.ipw,4)
  
  result$cp.A1A2.ipwn <- round(cp.A1A2.ipwn,4)
  result$cp.A1A3.ipwn <- round(cp.A1A3.ipwn,4)
  result$cp.A2A1.ipwn <- round(cp.A2A1.ipwn,4)
  result$cp.A2A3.ipwn <- round(cp.A2A3.ipwn,4)
  result$cp.A3A1.ipwn <- round(cp.A3A1.ipwn,4)
  result$cp.A3A2.ipwn <- round(cp.A3A2.ipwn,4)
  
  result$cp.A1A2.g <- round(cp.A1A2.g,4)
  result$cp.A1A3.g <- round(cp.A1A3.g,4)
  result$cp.A2A1.g <- round(cp.A2A1.g,4)
  result$cp.A2A3.g <- round(cp.A2A3.g,4)
  result$cp.A3A1.g <- round(cp.A3A1.g,4)
  result$cp.A3A2.g <- round(cp.A3A2.g,4)
  
  
  ################## 3. # of patients treated in each DTR ######################
  
  N.A1A2 <- mean(n.A1A2)
  N.A1A3 <- mean(n.A1A3)
  N.A2A1 <- mean(n.A2A1)
  N.A2A3 <- mean(n.A2A3)
  N.A3A1 <- mean(n.A3A1)
  N.A3A2 <- mean(n.A3A2)
  N.A1A2.CI <- quantile(n.A1A2, c(0.025, 0.975))
  N.A1A3.CI <- quantile(n.A1A3, c(0.025, 0.975))
  N.A2A1.CI <- quantile(n.A2A1, c(0.025, 0.975))
  N.A2A3.CI <- quantile(n.A2A3, c(0.025, 0.975))
  N.A3A1.CI <- quantile(n.A3A1, c(0.025, 0.975))
  N.A3A2.CI <- quantile(n.A3A2, c(0.025, 0.975))
  
  result$N.A1A2 <- paste0(N.A1A2, "(",N.A1A2.CI[1], ",",N.A1A2.CI[2], ")")
  result$N.A1A3 <- paste0(N.A1A3, "(",N.A1A3.CI[1], ",",N.A1A3.CI[2], ")")
  result$N.A2A1 <- paste0(N.A2A1, "(",N.A2A1.CI[1], ",",N.A2A1.CI[2], ")")
  result$N.A2A3 <- paste0(N.A2A3, "(",N.A2A3.CI[1], ",",N.A2A3.CI[2], ")")
  result$N.A3A1 <- paste0(N.A3A1, "(",N.A3A1.CI[1], ",",N.A3A1.CI[2], ")")
  result$N.A3A2 <- paste0(N.A3A2, "(",N.A3A2.CI[1], ",",N.A3A2.CI[2], ")")
  
  # number of patients treated in the designated optimal DTR
  N.opt <- c(N.A1A2, N.A1A3, N.A2A1, N.A2A3, N.A3A1, N.A3A2)[which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  N.opt.CI.lo <- c(N.A1A2.CI[1], N.A1A3.CI[1], N.A2A1.CI[1], N.A2A3.CI[1], N.A3A1.CI[1], N.A3A2.CI[1])[which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  N.opt.CI.up <- c(N.A1A2.CI[2], N.A1A3.CI[2], N.A2A1.CI[2], N.A2A3.CI[2], N.A3A1.CI[2], N.A3A2.CI[2])[which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  result$N.opt <- paste0(N.opt, "(",N.opt.CI.lo, ",",N.opt.CI.up, ")")
  
  # number of patients treated in the designated worst DTR
  N.wor <- c(N.A1A2, N.A1A3, N.A2A1, N.A2A3, N.A3A1, N.A3A2)[which.min(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  N.wor.CI.lo <- c(N.A1A2.CI[1], N.A1A3.CI[1], N.A2A1.CI[1], N.A2A3.CI[1], N.A3A1.CI[1], N.A3A2.CI[1])[which.min(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  N.wor.CI.up <- c(N.A1A2.CI[2], N.A1A3.CI[2], N.A2A1.CI[2], N.A2A3.CI[2], N.A3A1.CI[2], N.A3A2.CI[2])[which.min(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  result$N.wor <- paste0(N.wor, "(",N.wor.CI.lo, ",",N.wor.CI.up, ")")
  
  ################## 4. # of patients response (success) in each DTR ######################
  
  R.A1A2 <- mean(rn.A1A2)
  R.A1A3 <- mean(rn.A1A3)
  R.A2A1 <- mean(rn.A2A1)
  R.A2A3 <- mean(rn.A2A3)
  R.A3A1 <- mean(rn.A3A1)
  R.A3A2 <- mean(rn.A3A2)
  
  R.A1A2.CI <- quantile(rn.A1A2, c(0.025, 0.975))
  R.A1A3.CI <- quantile(rn.A1A3, c(0.025, 0.975))
  R.A2A1.CI <- quantile(rn.A2A1, c(0.025, 0.975))
  R.A2A3.CI <- quantile(rn.A2A3, c(0.025, 0.975))
  R.A3A1.CI <- quantile(rn.A3A1, c(0.025, 0.975))
  R.A3A2.CI <- quantile(rn.A3A2, c(0.025, 0.975))
  
  # number of patients response in the optimal DTR
  R.opt <- c(R.A1A2, R.A1A3, R.A2A1, R.A2A3, R.A3A1, R.A3A2)[which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  R.opt.CI.lo <- c(R.A1A2.CI[1], R.A1A3.CI[1], R.A2A1.CI[1], R.A2A3.CI[1], R.A3A1.CI[1], R.A3A2.CI[1])[which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  R.opt.CI.up <- c(R.A1A2.CI[2], R.A1A3.CI[2], R.A2A1.CI[2], R.A2A3.CI[2], R.A3A1.CI[2], R.A3A2.CI[2])[which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  result$R.opt <- paste0(R.opt, "(",R.opt.CI.lo, ",",R.opt.CI.up, ")")
  
  
  # number of patients response in the worst DTR
  R.wor <- c(R.A1A2, R.A1A3, R.A2A1, R.A2A3, R.A3A1, R.A3A2)[which.min(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  R.wor.CI.lo <- c(R.A1A2.CI[1], R.A1A3.CI[1], R.A2A1.CI[1], R.A2A3.CI[1], R.A3A1.CI[1], R.A3A2.CI[1])[which.min(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  R.wor.CI.up <- c(R.A1A2.CI[2], R.A1A3.CI[2], R.A2A1.CI[2], R.A2A3.CI[2], R.A3A1.CI[2], R.A3A2.CI[2])[which.min(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2))]
  result$R.wor <- paste0(R.wor, "(",R.wor.CI.lo, ",",R.wor.CI.up, ")")
  
  ################## 5. # of patients response (success) in the trial ######################
  NR <- mean(nr)
  NR.CI <- quantile(nr, c(0.025, 0.975))
  result$NR <- paste0(NR, "(",NR.CI[1], ",",NR.CI[2], ")")
  
  ############################### 6. type I error ######################################
  # the proportion of 10000 MC replicates of which the designated optimal DTR is estimated to be the optimal
  
  if (pid.A1A2==pid.A1A3 & pid.A1A3==pid.A2A1 & pid.A2A1==pid.A2A3 & pid.A2A3==pid.A3A1 & pid.A3A1==pid.A3A2) {
    
    set.seed(1005)
    typeI.sm <- mean(est.opt.sm==sample(1:6, length(est.opt.sm), replace=T))
    result$typeI.sm <- typeI.sm
    
    set.seed(1005)
    typeI.ipw <- mean(est.opt.ipw==sample(1:6, length(est.opt.ipw), replace=T))
    result$typeI.ipw <- typeI.ipw
    
    
    set.seed(1005)
    typeI.ipwn <- mean(est.opt.ipwn==sample(1:6, length(est.opt.ipwn), replace=T))
    result$typeI.ipwn <- typeI.ipwn
    
    set.seed(1005)
    typeI.g <- mean(est.opt.g==sample(1:6, length(est.opt.g), replace=T))
    result$typeI.g <- typeI.g
  }
  
  
  
  ################################# 7. Power #############################################
  
  else {
    ## 1) the proportion of 10000 MC replicates of which the designated optimal DTR is estimated to be the optimal
    pow.sm <- mean(est.opt.sm==which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2)))
    result$pow.sm <- pow.sm
    
    pow.ipw <- mean(est.opt.ipw==which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2)))
    result$pow.ipw <- pow.ipw
    
    pow.ipwn <- mean(est.opt.ipwn==which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2)))
    result$pow.ipwn <- pow.ipwn
    
    pow.g <- mean(est.opt.g==which.max(c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2)))
    result$pow.g <- pow.g
    
    
    
    ## 2) power to identify a set of top two or three within 5% margin, 
    true.para <- c(pid.A1A2, pid.A1A3, pid.A2A1, pid.A2A3, pid.A3A1, pid.A3A2)
    true.max <- which.max(true.para)
    # a function that returns the position of n-th largest
    maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
    thre <- 0.05
    
    ### sample mean estimator
    df.sm <- cbind(mu.A1A2hat.sm, mu.A1A3hat.sm, mu.A2A1hat.sm, mu.A2A3hat.sm, mu.A3A1hat.sm, mu.A3A2hat.sm)
    top.sm <- apply(df.sm, 1, maxn(1))
    second.sm <- apply(df.sm, 1, maxn(2))
    third.sm <- apply(df.sm, 1, maxn(3))
    sm.1 <- sapply(1:length(top.sm),function(i){df.sm[i,top.sm[i]]})
    sm.2 <- sapply(1:length(second.sm),function(i){df.sm[i,second.sm[i]]})
    sm.3 <- sapply(1:length(third.sm),function(i){df.sm[i,third.sm[i]]})
    # within the threshold
    sm2.ind <- ifelse(sm.2>=(sm.1-thre), second.sm, 0)
    sm3.ind <- ifelse(sm.3>=(sm.1-thre), third.sm, 0)
    
    topset.sm <- cbind(top.sm, sm2.ind, sm3.ind)
    
    # power defined as the proportion of c(top.sm, sm2.ind, sm3.ind) cover true.max
    powset.sm <- mean(apply(topset.sm, 1, function(r) any(true.max %in% r)))
    result$powset.sm <- powset.sm
    
    
    ### G-estimator
    df.g <- cbind(mu.A1A2hat.g, mu.A1A3hat.g, mu.A2A1hat.g, mu.A2A3hat.g, mu.A3A1hat.g, mu.A3A2hat.g)
    top.g <- apply(df.g, 1, maxn(1))
    second.g <- apply(df.g, 1, maxn(2))
    third.g <- apply(df.g, 1, maxn(3))
    g.1 <- sapply(1:length(top.g),function(i){df.g[i,top.g[i]]})
    g.2 <- sapply(1:length(second.g),function(i){df.g[i,second.g[i]]})
    g.3 <- sapply(1:length(third.g),function(i){df.g[i,third.g[i]]})
    # within the threshold
    g2.ind <- ifelse(g.2>=(g.1-thre), second.g, 0)
    g3.ind <- ifelse(g.3>=(g.1-thre), third.g, 0)
    
    topset.g <- cbind(top.g, g2.ind, g3.ind)
    
    # power defined as the proportion of c(top.g, g2.ind, g3.ind) cover true.max
    powset.g <- mean(apply(topset.g, 1, function(r) any(true.max %in% r)))
    result$powset.g <- powset.g
    
    
    ### IPW estimator
    df.ipw <- cbind(mu.A1A2hat.ipw, mu.A1A3hat.ipw, mu.A2A1hat.ipw, mu.A2A3hat.ipw, mu.A3A1hat.ipw, mu.A3A2hat.ipw)
    top.ipw <- apply(df.ipw, 1, maxn(1))
    second.ipw <- apply(df.ipw, 1, maxn(2))
    third.ipw <- apply(df.ipw, 1, maxn(3))
    ipw.1 <- sapply(1:length(top.ipw),function(i){df.ipw[i,top.ipw[i]]})
    ipw.2 <- sapply(1:length(second.ipw),function(i){df.ipw[i,second.ipw[i]]})
    ipw.3 <- sapply(1:length(third.ipw),function(i){df.ipw[i,third.ipw[i]]})
    # within the threshold
    ipw2.ind <- ifelse(ipw.2>=(ipw.1-thre), second.ipw, 0)
    ipw3.ind <- ifelse(ipw.3>=(ipw.1-thre), third.ipw, 0)
    
    topset.ipw <- cbind(top.ipw, ipw2.ind, ipw3.ind)
    
    # power defined as the proportion of c(top.ipw, ipw2.ind, ipw3.ind) cover true.max
    powset.ipw <- mean(apply(topset.ipw, 1, function(r) any(true.max %in% r)))
    result$powset.ipw <- powset.ipw
    
    
    ### IPWN estimator
    df.ipwn <- cbind(mu.A1A2hat.ipwn, mu.A1A3hat.ipwn, mu.A2A1hat.ipwn, mu.A2A3hat.ipwn, mu.A3A1hat.ipwn, mu.A3A2hat.ipwn)
    top.ipwn <- apply(df.ipwn, 1, maxn(1))
    second.ipwn <- apply(df.ipwn, 1, maxn(2))
    third.ipwn <- apply(df.ipwn, 1, maxn(3))
    ipwn.1 <- sapply(1:length(top.ipwn),function(i){df.ipwn[i,top.ipwn[i]]})
    ipwn.2 <- sapply(1:length(second.ipwn),function(i){df.ipwn[i,second.ipwn[i]]})
    ipwn.3 <- sapply(1:length(third.ipwn),function(i){df.ipwn[i,third.ipwn[i]]})
    # within the threshold
    ipwn2.ind <- ifelse(ipwn.2>=(ipwn.1-thre), second.ipwn, 0)
    ipwn3.ind <- ifelse(ipwn.3>=(ipwn.1-thre), third.ipwn, 0)
    
    topset.ipwn <- cbind(top.ipwn, ipwn2.ind, ipwn3.ind)
    
    # power defined as the proportion of c(top.ipwn, ipwn2.ind, ipwn3.ind) cover true.max
    powset.ipwn <- mean(apply(topset.ipwn, 1, function(r) any(true.max %in% r)))
    result$powset.ipwn <- powset.ipwn
    
  }
  
  result <- as.data.frame(unlist(result))
  colnames(result) <- NULL
  return(result)
  
}

