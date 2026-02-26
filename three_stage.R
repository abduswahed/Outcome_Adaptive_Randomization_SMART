#######################################################################
# Three-stage version GO-SMART compared with SMART                    #
# Assume 4 potential treatments, T=(A1, A2, A3, A4)                   #
#######################################################################

# This code is for the three-stage version GO-SMART design in comparison with three-stage SMART design (Web Appendix F)
# Three functions included in this code:
## 1. data generation function: "data_gene_three_c" 
## 2. Monte Carlo replication function: "monte_three" 
## 3. Evaluating function: "eva_three"

# Reproduce simulation under the alternative scenario under n=1000 in the sTable 1 in Web Appendix

######## Step1: data generation ##########

data_gene_three_c <- function(n,pi,p0.burn,p1.burn,p2.burn,AR,c,e) {
  
  # true parameters
  pi1.A1 <- pi[1]
  pi1.A2 <- pi[2]
  pi1.A3 <- pi[3]
  pi1.A4 <- pi[4]
  pi2.A1A2 <- pi[5]
  pi2.A1A3 <- pi[6]
  pi2.A1A4 <- pi[7]
  pi2.A2A1 <- pi[8]
  pi2.A2A3 <- pi[9]
  pi2.A2A4 <- pi[10]
  pi2.A3A1 <- pi[11]
  pi2.A3A2 <- pi[12]
  pi2.A3A4 <- pi[13]
  pi2.A4A1 <- pi[14]
  pi2.A4A2 <- pi[15]
  pi2.A4A3 <- pi[16]
  pi3.A1A2A3 <- pi[17]
  pi3.A1A2A4 <- pi[18]
  pi3.A1A3A2 <- pi[19]
  pi3.A1A3A4 <- pi[20]
  pi3.A1A4A2 <- pi[21]
  pi3.A1A4A3 <- pi[22]
  pi3.A2A1A3 <- pi[23]
  pi3.A2A1A4 <- pi[24]
  pi3.A2A3A1 <- pi[25]
  pi3.A2A3A4 <- pi[26]
  pi3.A2A4A1 <- pi[27]
  pi3.A2A4A3 <- pi[28]
  pi3.A3A1A2 <- pi[29]
  pi3.A3A1A4 <- pi[30]
  pi3.A3A2A1 <- pi[31]
  pi3.A3A2A4 <- pi[32]
  pi3.A3A4A1 <- pi[33]
  pi3.A3A4A2 <- pi[34]
  pi3.A4A1A2 <- pi[35]
  pi3.A4A1A3 <- pi[36]
  pi3.A4A2A1 <- pi[37]
  pi3.A4A2A3 <- pi[38]
  pi3.A4A3A1 <- pi[39]
  pi3.A4A3A2 <- pi[40]
  
  ##################################################################################################
  ###################################### Stage I ##################################################
  ##################################################################################################
  
  
  # stage I potential outcomes: Y_1^A1, Y_1^A2, Y_1^A3, Y_1^A4
  Y1.A1 <- rbinom(n, 1, pi1.A1)
  Y1.A2 <- rbinom(n, 1, pi1.A2)
  Y1.A3 <- rbinom(n, 1, pi1.A3)
  Y1.A4 <- rbinom(n, 1, pi1.A4)
  
  # number(n0.burn) and proportion(p0.burn): first stage burn-in sample
  n0.burn <- round(n*p0.burn, 0)
  
  # treatment assignment A in stage I
  A_1 <- rep(NA, n)
  I1.A1 <- rep(NA, n)
  I1.A2 <- rep(NA, n)
  I1.A3 <- rep(NA, n)
  I1.A4 <- rep(NA, n)
  
  ## i=1,...,n0.burn: equally assign
  A_1_n0.burn <- sample(1:4, size = n0.burn, replace = TRUE, prob = c(1/4,1/4,1/4,1/4))
  A_1[1:n0.burn] <- A_1_n0.burn
  I1.A1[1:n0.burn]<- ifelse(A_1_n0.burn==1, 1, 0)
  I1.A2[1:n0.burn]<- ifelse(A_1_n0.burn==2, 1, 0)
  I1.A3[1:n0.burn]<- ifelse(A_1_n0.burn==3, 1, 0)
  I1.A4[1:n0.burn]<- ifelse(A_1_n0.burn==4, 1, 0)
  
  # stage I randomization probabilities
  P1.A1 <- rep(NA, n)
  P1.A2 <- rep(NA, n)
  P1.A3 <- rep(NA, n)
  P1.A4 <- rep(NA, n)
  P1.A1[1:n0.burn] <- rep(1/4, n0.burn)
  P1.A2[1:n0.burn] <- rep(1/4, n0.burn)
  P1.A3[1:n0.burn] <- rep(1/4, n0.burn)
  P1.A4[1:n0.burn] <- rep(1/4, n0.burn)
  
  
  # observed outcome Y1 in stage I
  Y1 <- rep(NA, n)
  Y1[1:n0.burn] <- Y1.A1[1:n0.burn]*I1.A1[1:n0.burn]+Y1.A2[1:n0.burn]*I1.A2[1:n0.burn]+Y1.A3[1:n0.burn]*I1.A3[1:n0.burn]+Y1.A4[1:n0.burn]*I1.A4[1:n0.burn]
  
  ## i=n0+1,...,n: update the randomization probability using all previous data
  if(n0.burn<n) {
    for (i in (n0.burn+1):n) { 
      # P(Y1.A1=1)
      Y1.A1hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])
      # P(Y1.A2=1)
      Y1.A2hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])
      # P(Y1.A3=1)
      Y1.A3hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])
      # P(Y1.A4=1)
      Y1.A4hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])
      
      #AR(c,e)
      if (c=="n.2N"){
        c1 <- i/(2*n)
      } else if (c=="n.N") {
        c1 <- i/n
      } else {
        c1 <- c
      }
        
      r1.A1 <- ((Y1.A1hat)^c1)/((Y1.A1hat)^c1+(Y1.A2hat)^c1+(Y1.A3hat)^c1+(Y1.A4hat)^c1)
      r1.A2 <- ((Y1.A2hat)^c1)/((Y1.A1hat)^c1+(Y1.A2hat)^c1+(Y1.A3hat)^c1+(Y1.A4hat)^c1)
      r1.A3 <- ((Y1.A3hat)^c1)/((Y1.A1hat)^c1+(Y1.A2hat)^c1+(Y1.A3hat)^c1+(Y1.A4hat)^c1)
      r1.A4 <- ((Y1.A4hat)^c1)/((Y1.A1hat)^c1+(Y1.A2hat)^c1+(Y1.A3hat)^c1+(Y1.A4hat)^c1)
      
      
      # AR(c,e)
      e1 <- e[1]
      
      P1.A1[i] <- ifelse(is.na(r1.A1), 1/4,  min(max(c(r1.A1, e1)), 1-e1))
      P1.A2[i] <- ifelse(is.na(r1.A2), 1/4,  min(max(c(r1.A2, e1)), 1-e1))
      P1.A3[i] <- ifelse(is.na(r1.A3), 1/4,  min(max(c(r1.A3, e1)), 1-e1))
      P1.A4[i] <- ifelse(is.na(r1.A4), 1/4,  min(max(c(r1.A4, e1)), 1-e1))
      
      
      # check if the bound used!!
      P1.old <- c(P1.A1[i], P1.A2[i], P1.A3[i], P1.A4[i])
      P1.new <- rep(NA,4)
      # if e1=0.1 used, update the probability
      if(e1 %in% P1.old) {
        P1.new[which(e1==P1.old)] <- P1.old[which(e1==P1.old)]
        P1.new[-which(e1==P1.old)] <- (1-sum(P1.old[which(e1==P1.old)]))*P1.old[-which(e1==P1.old)]/sum(P1.old[-which(e1==P1.old)])    
        
        P1.A1[i] <- P1.new[1]
        P1.A2[i] <- P1.new[2]
        P1.A3[i] <- P1.new[3]
        P1.A4[i] <- P1.new[4]
      } 
      
      A_1[i] <- sample(1:4, size = 1, replace = TRUE, prob = c(P1.A1[i],P1.A2[i],P1.A3[i],P1.A4[i]))
      I1.A1[i]<- ifelse(A_1[i]==1, 1, 0)
      I1.A2[i]<- ifelse(A_1[i]==2, 1, 0)
      I1.A3[i]<- ifelse(A_1[i]==3, 1, 0)
      I1.A4[i]<- ifelse(A_1[i]==4, 1, 0)
      Y1[i] <- Y1.A1[i]*I1.A1[i]+Y1.A2[i]*I1.A2[i]+Y1.A3[i]*I1.A3[i]+Y1.A4[i]*I1.A4[i]
      
    }
  }
  
  
  ##################################################################################################
  ###################################### Stage II ##################################################
  ##################################################################################################
  
  # stage II potential outcome: Y_2^A1A2,Y_2^A1A3,Y_2^A1A4, Y_2^A2A1,Y_2^A2A3, Y_2^A2A4, Y_2^A3A1,Y_2^A3A2, Y_2^A3A4,Y_2^A4A1,Y_2^A4A2, Y_2^A4A3
  Y2.A1A2.raw <- rbinom(n, 1, pi2.A1A2)
  Y2.A1A3.raw <- rbinom(n, 1, pi2.A1A3)
  Y2.A1A4.raw <- rbinom(n, 1, pi2.A1A4)
  Y2.A2A1.raw <- rbinom(n, 1, pi2.A2A1)
  Y2.A2A3.raw <- rbinom(n, 1, pi2.A2A3)
  Y2.A2A4.raw <- rbinom(n, 1, pi2.A2A4)
  Y2.A3A1.raw <- rbinom(n, 1, pi2.A3A1)
  Y2.A3A2.raw <- rbinom(n, 1, pi2.A3A2)
  Y2.A3A4.raw <- rbinom(n, 1, pi2.A3A4)
  Y2.A4A1.raw <- rbinom(n, 1, pi2.A4A1)
  Y2.A4A2.raw <- rbinom(n, 1, pi2.A4A2)
  Y2.A4A3.raw <- rbinom(n, 1, pi2.A4A3)
  
  # if no response in stage I, Y_2^AjAl~Bernouli(pi_2^AjAl)
  Y2.A1A2 <- ifelse(Y1.A1==0, Y2.A1A2.raw, NA)
  Y2.A1A3 <- ifelse(Y1.A1==0, Y2.A1A3.raw, NA)
  Y2.A1A4 <- ifelse(Y1.A1==0, Y2.A1A4.raw, NA)
  Y2.A2A1 <- ifelse(Y1.A2==0, Y2.A2A1.raw, NA)
  Y2.A2A3 <- ifelse(Y1.A2==0, Y2.A2A3.raw, NA)
  Y2.A2A4 <- ifelse(Y1.A2==0, Y2.A2A4.raw, NA)
  Y2.A3A1 <- ifelse(Y1.A3==0, Y2.A3A1.raw, NA)
  Y2.A3A2 <- ifelse(Y1.A3==0, Y2.A3A2.raw, NA)
  Y2.A3A4 <- ifelse(Y1.A3==0, Y2.A3A4.raw, NA)
  Y2.A4A1 <- ifelse(Y1.A4==0, Y2.A4A1.raw, NA)
  Y2.A4A2 <- ifelse(Y1.A4==0, Y2.A4A2.raw, NA)
  Y2.A4A3 <- ifelse(Y1.A4==0, Y2.A4A3.raw, NA)
  
  
  # stage II randomization probabilities
  P2.A1A2 <- rep(NA, n)
  P2.A1A3 <- rep(NA, n)
  P2.A1A4 <- rep(NA, n)
  P2.A2A1 <- rep(NA, n)
  P2.A2A3 <- rep(NA, n)
  P2.A2A4 <- rep(NA, n)
  P2.A3A1 <- rep(NA, n)
  P2.A3A2 <- rep(NA, n)
  P2.A3A4 <- rep(NA, n)
  P2.A4A1 <- rep(NA, n)
  P2.A4A2 <- rep(NA, n)
  P2.A4A3 <- rep(NA, n)
  
  
  # treatment assignment A in stage II
  A_2 <- rep(NA, n)
  I2.A1 <- rep(NA, n)
  I2.A2 <- rep(NA, n)
  I2.A3 <- rep(NA, n)
  I2.A4 <- rep(NA, n)
  
  # observed outcome Y2 in stage II
  Y2 <- rep(NA, n)
  
  # i=1,...,n0: for those who did not respond in stage I, equally randomized to the rest treatments
  if(n0.burn<=n) {
    for (i in 1:n0.burn) {
      if (I1.A1[i]==1 & Y1.A1[i]==0) {
        P2.A1A2[i] <- 1/3
        P2.A1A3[i] <- 1/3
        P2.A1A4[i] <- 1/3
        A_2[i] <- sample(c(2,3,4), size = 1, replace = TRUE, prob = c(P2.A1A2[i], P2.A1A3[i], P2.A1A4[i]))
      } else if (I1.A2[i]==1 & Y1.A2[i]==0) {
        P2.A2A1[i] <- 1/3
        P2.A2A3[i] <- 1/3
        P2.A2A4[i] <- 1/3
        A_2[i] <- sample(c(1,3,4), size = 1, replace = TRUE, prob = c(P2.A2A1[i], P2.A2A3[i], P2.A2A4[i]))
      } else if (I1.A3[i]==1 & Y1.A3[i]==0) {
        P2.A3A1[i] <- 1/3
        P2.A3A2[i] <- 1/3
        P2.A3A4[i] <- 1/3
        A_2[i] <- sample(c(1,2,4), size = 1, replace = TRUE, prob = c(P2.A3A1[i], P2.A3A2[i], P2.A3A4[i]))
      } else if (I1.A4[i]==1 & Y1.A4[i]==0){
        P2.A4A1[i] <- 1/3
        P2.A4A2[i] <- 1/3
        P2.A4A3[i] <- 1/3
        A_2[i] <- sample(c(1,2,3), size = 1, replace = TRUE, prob = c(P2.A4A1[i], P2.A4A2[i], P2.A4A3[i]))
      }
      
      I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
      I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
      I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
      I2.A4[i]<- ifelse(A_2[i]==4, 1, 0)
      
      
      
      Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],Y2.A1A4[i]*I2.A4[i]*I1.A1[i],
                                      Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],Y2.A2A4[i]*I2.A4[i]*I1.A2[i],
                                      Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i],Y2.A3A4[i]*I2.A4[i]*I1.A3[i],
                                      Y2.A4A1[i]*I2.A1[i]*I1.A4[i],Y2.A4A2[i]*I2.A2[i]*I1.A4[i],Y2.A4A3[i]*I2.A3[i]*I1.A4[i]), na.rm = T), NA)
      
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
      # P(Y1.A4=1)
      Y1.A4hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])
      
      if (c=="n.2N") {
        c2 <- i/(2*n)
      } else if (c=="n.N") {
        c2 <- i/n
      } else {
        c2 <- c
      }
      
      
      r2.A1A2 <- ((Y1.A2hat)^c2)/((Y1.A2hat)^c2+(Y1.A3hat)^c2+(Y1.A4hat)^c2)
      r2.A1A3 <- ((Y1.A3hat)^c2)/((Y1.A2hat)^c2+(Y1.A3hat)^c2+(Y1.A4hat)^c2)
      r2.A1A4 <- ((Y1.A4hat)^c2)/((Y1.A2hat)^c2+(Y1.A3hat)^c2+(Y1.A4hat)^c2)
      r2.A2A1 <- ((Y1.A1hat)^c2)/((Y1.A1hat)^c2+(Y1.A3hat)^c2+(Y1.A4hat)^c2)
      r2.A2A3 <- ((Y1.A3hat)^c2)/((Y1.A1hat)^c2+(Y1.A3hat)^c2+(Y1.A4hat)^c2)
      r2.A2A4 <- ((Y1.A4hat)^c2)/((Y1.A1hat)^c2+(Y1.A3hat)^c2+(Y1.A4hat)^c2)
      r2.A3A1 <- ((Y1.A1hat)^c2)/((Y1.A1hat)^c2+(Y1.A2hat)^c2+(Y1.A4hat)^c2)
      r2.A3A2 <- ((Y1.A2hat)^c2)/((Y1.A1hat)^c2+(Y1.A2hat)^c2+(Y1.A4hat)^c2)
      r2.A3A4 <- ((Y1.A4hat)^c2)/((Y1.A1hat)^c2+(Y1.A2hat)^c2+(Y1.A4hat)^c2)
      r2.A4A1 <- ((Y1.A1hat)^c2)/((Y1.A1hat)^c2+(Y1.A2hat)^c2+(Y1.A3hat)^c2)
      r2.A4A2 <- ((Y1.A2hat)^c2)/((Y1.A1hat)^c2+(Y1.A2hat)^c2+(Y1.A3hat)^c2)
      r2.A4A3 <- ((Y1.A3hat)^c2)/((Y1.A1hat)^c2+(Y1.A2hat)^c2+(Y1.A3hat)^c2)
      
      # AR(c,e)
      e2 <- e[2]
      
      if (I1.A1[i]==1 & Y1.A1[i]==0) {
        # P2.A1A2, P2.A1A3, P2.A1A4
        P2.A1A2[i] <- ifelse(is.na(r2.A1A2), 1/3,  min(max(c(r2.A1A2, e2)), 1-e2))
        P2.A1A3[i] <- ifelse(is.na(r2.A1A3), 1/3,  min(max(c(r2.A1A3, e2)), 1-e2))
        P2.A1A4[i] <- ifelse(is.na(r2.A1A4), 1/3,  min(max(c(r2.A1A4, e2)), 1-e2))
        
        # check if the bound used!!
        P2.A1.y.old <- c(P2.A1A2[i], P2.A1A3[i], P2.A1A4[i])
        P2.A1.y.new <- rep(NA,3)
        # if e2=0.1 used, update the probability 
        if(e2 %in% P2.A1.y.old) {
          P2.A1.y.new[which(e2==P2.A1.y.old)] <- P2.A1.y.old[which(e2==P2.A1.y.old)]
          P2.A1.y.new[-which(e2==P2.A1.y.old)] <- (1-sum(P2.A1.y.old[which(e2==P2.A1.y.old)]))*P2.A1.y.old[-which(e2==P2.A1.y.old)]/sum(P2.A1.y.old[-which(e2==P2.A1.y.old)])    
          
          P2.A1A2[i] <- P2.A1.y.new[1]
          P2.A1A3[i] <- P2.A1.y.new[2]
          P2.A1A4[i] <- P2.A1.y.new[3]
        } 
      
        A_2[i] <- sample(c(2,3,4), size = 1, replace = TRUE, prob = c(P2.A1A2[i], P2.A1A3[i], P2.A1A4[i]))
        
      } else if (I1.A2[i]==1 & Y1.A2[i]==0) {
        # P2.A2A1, P2.A2A3, P2.A2A4
        P2.A2A1[i] <- ifelse(is.na(r2.A2A1), 1/3,  min(max(c(r2.A2A1, e2)), 1-e2))
        P2.A2A3[i] <- ifelse(is.na(r2.A2A3), 1/3,  min(max(c(r2.A2A3, e2)), 1-e2))
        P2.A2A4[i] <- ifelse(is.na(r2.A2A4), 1/3,  min(max(c(r2.A2A4, e2)), 1-e2))
        
        # check if the bound was used!!
        P2.A2.y.old <- c(P2.A2A1[i], P2.A2A3[i], P2.A2A4[i])
        P2.A2.y.new <- rep(NA,3)
        # if e2=0.1 used, update the probability 
        if(e2 %in% P2.A2.y.old) {
          P2.A2.y.new[which(e2==P2.A2.y.old)] <- P2.A2.y.old[which(e2==P2.A2.y.old)]
          P2.A2.y.new[-which(e2==P2.A2.y.old)] <- (1-sum(P2.A2.y.old[which(e2==P2.A2.y.old)]))*P2.A2.y.old[-which(e2==P2.A2.y.old)]/sum(P2.A2.y.old[-which(e2==P2.A2.y.old)])    
          
          P2.A2A1[i] <- P2.A2.y.new[1]
          P2.A2A3[i] <- P2.A2.y.new[2]
          P2.A2A4[i] <- P2.A2.y.new[3]
        } 
      
        A_2[i] <- sample(c(1,3,4), size = 1, replace = TRUE, prob = c(P2.A2A1[i], P2.A2A3[i], P2.A2A4[i]))
        
      } else if (I1.A3[i]==1 & Y1.A3[i]==0) {
        # P2.A3A1, P2.A3A2, P2.A3A4
        P2.A3A1[i] <- ifelse(is.na(r2.A3A1), 1/3,  min(max(c(r2.A3A1, e2)), 1-e2))
        P2.A3A2[i] <- ifelse(is.na(r2.A3A2), 1/3,  min(max(c(r2.A3A2, e2)), 1-e2))
        P2.A3A4[i] <- ifelse(is.na(r2.A3A4), 1/3,  min(max(c(r2.A3A4, e2)), 1-e2))
        
        # check if the bound was used!!
        P2.A3.y.old <- c(P2.A3A1[i], P2.A3A2[i], P2.A3A4[i])
        P2.A3.y.new <- rep(NA,3)
        # if e2=0.1 used, update the probability 
        if(e2 %in% P2.A3.y.old) {
          P2.A3.y.new[which(e2==P2.A3.y.old)] <- P2.A3.y.old[which(e2==P2.A3.y.old)]
          P2.A3.y.new[-which(e2==P2.A3.y.old)] <- (1-sum(P2.A3.y.old[which(e2==P2.A3.y.old)]))*P2.A3.y.old[-which(e2==P2.A3.y.old)]/sum(P2.A3.y.old[-which(e2==P2.A3.y.old)])    
          
          P2.A3A1[i] <- P2.A3.y.new[1]
          P2.A3A2[i] <- P2.A3.y.new[2]
          P2.A3A4[i] <- P2.A3.y.new[3]
        } 
        
        A_2[i] <- sample(c(1,2,4), size = 1, replace = TRUE, prob = c(P2.A3A1[i], P2.A3A2[i], P2.A3A4[i]))
        
      } else if (I1.A4[i]==1 & Y1.A4[i]==0) {
        # P2.A4A1, P2.A4A2, P2.A4A3
        P2.A4A1[i] <- ifelse(is.na(r2.A4A1), 1/3,  min(max(c(r2.A4A1, e2)), 1-e2))
        P2.A4A2[i] <- ifelse(is.na(r2.A4A2), 1/3,  min(max(c(r2.A4A2, e2)), 1-e2))
        P2.A4A3[i] <- ifelse(is.na(r2.A4A3), 1/3,  min(max(c(r2.A4A3, e2)), 1-e2))
        
        # check if the bound was used!!
        P2.A4.y.old <- c(P2.A4A1[i], P2.A4A2[i], P2.A4A3[i])
        P2.A4.y.new <- rep(NA,3)
        # if e2=0.1 used, update the probability 
        if(e2 %in% P2.A4.y.old) {
          P2.A4.y.new[which(e2==P2.A4.y.old)] <- P2.A4.y.old[which(e2==P2.A4.y.old)]
          P2.A4.y.new[-which(e2==P2.A4.y.old)] <- (1-sum(P2.A4.y.old[which(e2==P2.A4.y.old)]))*P2.A4.y.old[-which(e2==P2.A4.y.old)]/sum(P2.A4.y.old[-which(e2==P2.A4.y.old)])    
          
          P2.A4A1[i] <- P2.A4.y.new[1]
          P2.A4A2[i] <- P2.A4.y.new[2]
          P2.A4A3[i] <- P2.A4.y.new[3]
        } 
        
        A_2[i] <- sample(c(1,2,3), size = 1, replace = TRUE, prob = c(P2.A4A1[i], P2.A4A2[i], P2.A4A3[i]))
      }
      
      I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
      I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
      I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
      I2.A4[i]<- ifelse(A_2[i]==4, 1, 0)
      
      Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],Y2.A1A4[i]*I2.A4[i]*I1.A1[i],
                                      Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],Y2.A2A4[i]*I2.A4[i]*I1.A2[i],
                                      Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i],Y2.A3A4[i]*I2.A4[i]*I1.A3[i],
                                      Y2.A4A1[i]*I2.A1[i]*I1.A4[i],Y2.A4A2[i]*I2.A2[i]*I1.A4[i],Y2.A4A3[i]*I2.A3[i]*I1.A4[i]), na.rm = T), NA)
      
    
    }
  }
  
  if(n0.burn<=n1.burn & n1.burn<n) {
    # i=n1+1,...,n: update the randomization probability using conditional probability 
    for (i in (n1.burn+1):n) {
      
      # G-computation
      
      # P(Y2.A1A2=1|Y1.A1=0)
      Y2.A1A2hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)
      # P(Y2.A1A3=1|Y1.A1=0)
      Y2.A1A3hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)
      # P(Y2.A1A4=1|Y1.A1=0)
      Y2.A1A4hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)
      
      # P(Y2.A2A1=1|Y1.A2=0)
      Y2.A2A1hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)
      # P(Y2.A2A3=1|Y1.A2=0)
      Y2.A2A3hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)
      # P(Y2.A2A4=1|Y1.A2=0)
      Y2.A2A4hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)
      
      # P(Y2.A3A1=1|Y1.A3=0)
      Y2.A3A1hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)
      # P(Y2.A3A2=1|Y1.A3=0)
      Y2.A3A2hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)
      # P(Y2.A3A4=1|Y1.A3=0)
      Y2.A3A4hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)
      
      # P(Y2.A4A1=1|Y1.A4=0)
      Y2.A4A1hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)
      # P(Y2.A4A2=1|Y1.A4=0)
      Y2.A4A2hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)
      # P(Y2.A4A3=1|Y1.A4=0)
      Y2.A4A3hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)
      
      if (c=="n.2N") {
        c3 <- i/(2*n)
      } else if (c=="n.N") {
        c3 <- i/n
      } else {
        c3 <- c
      }
        
      r3.A1A2 <- ((Y2.A1A2hat)^c3)/((Y2.A1A2hat)^c3+(Y2.A1A3hat)^c3+(Y2.A1A4hat)^c3)
      r3.A1A3 <- ((Y2.A1A3hat)^c3)/((Y2.A1A2hat)^c3+(Y2.A1A3hat)^c3+(Y2.A1A4hat)^c3)
      r3.A1A4 <- ((Y2.A1A4hat)^c3)/((Y2.A1A2hat)^c3+(Y2.A1A3hat)^c3+(Y2.A1A4hat)^c3)
      r3.A2A1 <- ((Y2.A2A1hat)^c3)/((Y2.A2A1hat)^c3+(Y2.A2A3hat)^c3+(Y2.A2A4hat)^c3)
      r3.A2A3 <- ((Y2.A2A3hat)^c3)/((Y2.A2A1hat)^c3+(Y2.A2A3hat)^c3+(Y2.A2A4hat)^c3)
      r3.A2A4 <- ((Y2.A2A4hat)^c3)/((Y2.A2A1hat)^c3+(Y2.A2A3hat)^c3+(Y2.A2A4hat)^c3)
      r3.A3A1 <- ((Y2.A3A1hat)^c3)/((Y2.A3A1hat)^c3+(Y2.A3A2hat)^c3+(Y2.A3A4hat)^c3)
      r3.A3A2 <- ((Y2.A3A2hat)^c3)/((Y2.A3A1hat)^c3+(Y2.A3A2hat)^c3+(Y2.A3A4hat)^c3)
      r3.A3A4 <- ((Y2.A3A4hat)^c3)/((Y2.A3A1hat)^c3+(Y2.A3A2hat)^c3+(Y2.A3A4hat)^c3)
      r3.A4A1 <- ((Y2.A4A1hat)^c3)/((Y2.A4A1hat)^c3+(Y2.A4A2hat)^c3+(Y2.A4A3hat)^c3)
      r3.A4A2 <- ((Y2.A4A2hat)^c3)/((Y2.A4A1hat)^c3+(Y2.A4A2hat)^c3+(Y2.A4A3hat)^c3)
      r3.A4A3 <- ((Y2.A4A3hat)^c3)/((Y2.A4A1hat)^c3+(Y2.A4A2hat)^c3+(Y2.A4A3hat)^c3)
        
        # AR(c,e)
      e3 <- e[3]
      if (I1.A1[i]==1 & Y1.A1[i]==0) {
        # P2.A1A2, P2.A1A3, p2.A1A4
        P2.A1A2[i] <- ifelse(is.na(r3.A1A2), 1/3,  min(max(c(r3.A1A2, e3)), 1-e3))
        P2.A1A3[i] <- ifelse(is.na(r3.A1A3), 1/3,  min(max(c(r3.A1A3, e3)), 1-e3))
        P2.A1A4[i] <- ifelse(is.na(r3.A1A4), 1/3,  min(max(c(r3.A1A4, e3)), 1-e3))
        
        # check if the bound was used!!
        P2.A1.y.old <- c(P2.A1A2[i], P2.A1A3[i], P2.A1A4[i])
        P2.A1.y.new <- rep(NA,3)
        # if e3=0.1 used, update the probability 
        if(e3 %in% P2.A1.y.old) {
          P2.A1.y.new[which(e3==P2.A1.y.old)] <- P2.A1.y.old[which(e3==P2.A1.y.old)]
          P2.A1.y.new[-which(e3==P2.A1.y.old)] <- (1-sum(P2.A1.y.old[which(e3==P2.A1.y.old)]))*P2.A1.y.old[-which(e3==P2.A1.y.old)]/sum(P2.A1.y.old[-which(e3==P2.A1.y.old)])    
          
          P2.A1A2[i] <- P2.A1.y.new[1]
          P2.A1A3[i] <- P2.A1.y.new[2]
          P2.A1A4[i] <- P2.A1.y.new[3]
        } 

        A_2[i] <- sample(c(2,3,4), size = 1, replace = TRUE, prob = c(P2.A1A2[i], P2.A1A3[i], P2.A1A4[i]))
        
      } else if (I1.A2[i]==1 & Y1.A2[i]==0) {
        # P2.A2A1, P2.A2A3, P2.A2A4
        P2.A2A1[i] <- ifelse(is.na(r3.A2A1), 1/3,  min(max(c(r3.A2A1, e3)), 1-e3))
        P2.A2A3[i] <- ifelse(is.na(r3.A2A3), 1/3,  min(max(c(r3.A2A3, e3)), 1-e3))
        P2.A2A4[i] <- ifelse(is.na(r3.A2A4), 1/3,  min(max(c(r3.A2A4, e3)), 1-e3))
        
        # check if the bound was used!!
        P2.A2.y.old <- c(P2.A2A1[i], P2.A2A3[i], P2.A2A4[i])
        P2.A2.y.new <- rep(NA,3)
        # if e3=0.1 used, update the probability 
        if(e3 %in% P2.A2.y.old) {
          P2.A2.y.new[which(e3==P2.A2.y.old)] <- P2.A2.y.old[which(e3==P2.A2.y.old)]
          P2.A2.y.new[-which(e3==P2.A2.y.old)] <- (1-sum(P2.A2.y.old[which(e3==P2.A2.y.old)]))*P2.A2.y.old[-which(e3==P2.A2.y.old)]/sum(P2.A2.y.old[-which(e3==P2.A2.y.old)])    
          
          P2.A2A1[i] <- P2.A2.y.new[1]
          P2.A2A3[i] <- P2.A2.y.new[2]
          P2.A2A4[i] <- P2.A2.y.new[3]
        } 
        
        A_2[i] <- sample(c(1,3,4), size = 1, replace = TRUE, prob = c(P2.A2A1[i], P2.A2A3[i], P2.A2A4[i]))
        
      } else if (I1.A3[i]==1 & Y1.A3[i]==0) {
        # P2.A3A1, P2.A3A2, P2.A3A4
        P2.A3A1[i] <- ifelse(is.na(r3.A3A1), 1/3,  min(max(c(r3.A3A1, e3)), 1-e3))
        P2.A3A2[i] <- ifelse(is.na(r3.A3A2), 1/3,  min(max(c(r3.A3A2, e3)), 1-e3))
        P2.A3A4[i] <- ifelse(is.na(r3.A3A4), 1/3,  min(max(c(r3.A3A4, e3)), 1-e3))
        
        # check if the bound was used!!
        P2.A3.y.old <- c(P2.A3A1[i], P2.A3A2[i], P2.A3A4[i])
        P2.A3.y.new <- rep(NA,3)
        # if e3=0.1 used, update the probability 
        if(e3 %in% P2.A3.y.old) {
          P2.A3.y.new[which(e3==P2.A3.y.old)] <- P2.A3.y.old[which(e3==P2.A3.y.old)]
          P2.A3.y.new[-which(e3==P2.A3.y.old)] <- (1-sum(P2.A3.y.old[which(e3==P2.A3.y.old)]))*P2.A3.y.old[-which(e3==P2.A3.y.old)]/sum(P2.A3.y.old[-which(e3==P2.A3.y.old)])    
          
          P2.A3A1[i] <- P2.A3.y.new[1]
          P2.A3A2[i] <- P2.A3.y.new[2]
          P2.A3A4[i] <- P2.A3.y.new[3]
        } 
        
        A_2[i] <- sample(c(1,2,4), size = 1, replace = TRUE, prob = c(P2.A3A1[i], P2.A3A2[i], P2.A3A4[i]))
        
      } else if (I1.A4[i]==1 & Y1.A4[i]==0) {
        # P2.A4A1, P2.A4A2, P2.A4A3
        P2.A4A1[i] <- ifelse(is.na(r3.A4A1), 1/3,  min(max(c(r3.A4A1, e3)), 1-e3))
        P2.A4A2[i] <- ifelse(is.na(r3.A4A2), 1/3,  min(max(c(r3.A4A2, e3)), 1-e3))
        P2.A4A3[i] <- ifelse(is.na(r3.A4A3), 1/3,  min(max(c(r3.A4A3, e3)), 1-e3))
        
        # check if the bound was used!!
        P2.A4.y.old <- c(P2.A4A1[i], P2.A4A2[i], P2.A4A3[i])
        P2.A4.y.new <- rep(NA,3)
        # if e3=0.1 used, update the probability 
        if(e3 %in% P2.A4.y.old) {
          P2.A4.y.new[which(e3==P2.A4.y.old)] <- P2.A4.y.old[which(e3==P2.A4.y.old)]
          P2.A4.y.new[-which(e3==P2.A4.y.old)] <- (1-sum(P2.A4.y.old[which(e3==P2.A4.y.old)]))*P2.A4.y.old[-which(e3==P2.A4.y.old)]/sum(P2.A4.y.old[-which(e3==P2.A4.y.old)])    
          
          P2.A4A1[i] <- P2.A4.y.new[1]
          P2.A4A2[i] <- P2.A4.y.new[2]
          P2.A4A3[i] <- P2.A4.y.new[3]
        } 
        
        A_2[i] <- sample(c(1,2,3), size = 1, replace = TRUE, prob = c(P2.A4A1[i], P2.A4A2[i], P2.A4A3[i]))
        
      }
      
      I2.A1[i]<- ifelse(A_2[i]==1, 1, 0)
      I2.A2[i]<- ifelse(A_2[i]==2, 1, 0)
      I2.A3[i]<- ifelse(A_2[i]==3, 1, 0)
      I2.A4[i]<- ifelse(A_2[i]==4, 1, 0)
      Y2[i] <- ifelse(Y1[i]==0, sum(c(Y2.A1A2[i]*I2.A2[i]*I1.A1[i],Y2.A1A3[i]*I2.A3[i]*I1.A1[i],Y2.A1A4[i]*I2.A4[i]*I1.A1[i],
                                      Y2.A2A1[i]*I2.A1[i]*I1.A2[i],Y2.A2A3[i]*I2.A3[i]*I1.A2[i],Y2.A2A4[i]*I2.A4[i]*I1.A2[i],
                                      Y2.A3A1[i]*I2.A1[i]*I1.A3[i],Y2.A3A2[i]*I2.A2[i]*I1.A3[i],Y2.A3A4[i]*I2.A4[i]*I1.A3[i],
                                      Y2.A4A1[i]*I2.A1[i]*I1.A4[i],Y2.A4A2[i]*I2.A2[i]*I1.A4[i],Y2.A4A3[i]*I2.A3[i]*I1.A4[i]), na.rm = T), NA)

      
    }
  }
  
  ##################################################################################################
  ###################################### Stage III #################################################
  ##################################################################################################
  
  # stage III potential outcome
  Y3.A1A2A3.raw <- rbinom(n, 1, pi3.A1A2A3)
  Y3.A1A2A4.raw <- rbinom(n, 1, pi3.A1A2A4)
  Y3.A1A3A2.raw <- rbinom(n, 1, pi3.A1A3A2)
  Y3.A1A3A4.raw <- rbinom(n, 1, pi3.A1A3A4)
  Y3.A1A4A2.raw <- rbinom(n, 1, pi3.A1A4A2)
  Y3.A1A4A3.raw <- rbinom(n, 1, pi3.A1A4A3)
  Y3.A2A1A3.raw <- rbinom(n, 1, pi3.A2A1A3)
  Y3.A2A1A4.raw <- rbinom(n, 1, pi3.A2A1A4)
  Y3.A2A3A1.raw <- rbinom(n, 1, pi3.A2A3A1)
  Y3.A2A3A4.raw <- rbinom(n, 1, pi3.A2A3A4)
  Y3.A2A4A1.raw <- rbinom(n, 1, pi3.A2A4A1)
  Y3.A2A4A3.raw <- rbinom(n, 1, pi3.A2A4A3)
  Y3.A3A1A2.raw <- rbinom(n, 1, pi3.A3A1A2)
  Y3.A3A1A4.raw <- rbinom(n, 1, pi3.A3A1A4)
  Y3.A3A2A1.raw <- rbinom(n, 1, pi3.A3A2A1)
  Y3.A3A2A4.raw <- rbinom(n, 1, pi3.A3A2A4)
  Y3.A3A4A1.raw <- rbinom(n, 1, pi3.A3A4A1)
  Y3.A3A4A2.raw <- rbinom(n, 1, pi3.A3A4A2)
  Y3.A4A1A2.raw <- rbinom(n, 1, pi3.A4A1A2)
  Y3.A4A1A3.raw <- rbinom(n, 1, pi3.A4A1A3)
  Y3.A4A2A1.raw <- rbinom(n, 1, pi3.A4A2A1)
  Y3.A4A2A3.raw <- rbinom(n, 1, pi3.A4A2A3)
  Y3.A4A3A1.raw <- rbinom(n, 1, pi3.A4A3A1)
  Y3.A4A3A2.raw <- rbinom(n, 1, pi3.A4A3A2)
  
  # if no response in stage I and stage II, Y_3^AjAlAm~Bernouli(pi_3^AjAlAm)
  Y3.A1A2A3 <- ifelse(Y1.A1==0 & Y2.A1A2==0, Y3.A1A2A3.raw, NA)
  Y3.A1A2A4 <- ifelse(Y1.A1==0 & Y2.A1A2==0, Y3.A1A2A4.raw, NA)
  Y3.A1A3A2 <- ifelse(Y1.A1==0 & Y2.A1A3==0, Y3.A1A3A2.raw, NA)
  Y3.A1A3A4 <- ifelse(Y1.A1==0 & Y2.A1A3==0, Y3.A1A3A4.raw, NA)
  Y3.A1A4A2 <- ifelse(Y1.A1==0 & Y2.A1A4==0, Y3.A1A4A2.raw, NA)
  Y3.A1A4A3 <- ifelse(Y1.A1==0 & Y2.A1A4==0, Y3.A1A4A3.raw, NA)
  Y3.A2A1A3 <- ifelse(Y1.A2==0 & Y2.A2A1==0, Y3.A2A1A3.raw, NA)
  Y3.A2A1A4 <- ifelse(Y1.A2==0 & Y2.A2A1==0, Y3.A2A1A4.raw, NA)
  Y3.A2A3A1 <- ifelse(Y1.A2==0 & Y2.A2A3==0, Y3.A2A3A1.raw, NA)
  Y3.A2A3A4 <- ifelse(Y1.A2==0 & Y2.A2A3==0, Y3.A2A3A4.raw, NA)
  Y3.A2A4A1 <- ifelse(Y1.A2==0 & Y2.A2A4==0, Y3.A2A4A1.raw, NA)
  Y3.A2A4A3 <- ifelse(Y1.A2==0 & Y2.A2A4==0, Y3.A2A4A3.raw, NA)
  Y3.A3A1A2 <- ifelse(Y1.A3==0 & Y2.A3A1==0, Y3.A3A1A2.raw, NA)
  Y3.A3A1A4 <- ifelse(Y1.A3==0 & Y2.A3A1==0, Y3.A3A1A4.raw, NA)
  Y3.A3A2A1 <- ifelse(Y1.A3==0 & Y2.A3A2==0, Y3.A3A2A1.raw, NA)
  Y3.A3A2A4 <- ifelse(Y1.A3==0 & Y2.A3A2==0, Y3.A3A2A4.raw, NA)
  Y3.A3A4A1 <- ifelse(Y1.A3==0 & Y2.A3A4==0, Y3.A3A4A1.raw, NA)
  Y3.A3A4A2 <- ifelse(Y1.A3==0 & Y2.A3A4==0, Y3.A3A4A2.raw, NA)
  Y3.A4A1A2 <- ifelse(Y1.A4==0 & Y2.A4A1==0, Y3.A4A1A2.raw, NA)
  Y3.A4A1A3 <- ifelse(Y1.A4==0 & Y2.A4A1==0, Y3.A4A1A3.raw, NA)
  Y3.A4A2A1 <- ifelse(Y1.A4==0 & Y2.A4A2==0, Y3.A4A2A1.raw, NA)
  Y3.A4A2A3 <- ifelse(Y1.A4==0 & Y2.A4A2==0, Y3.A4A2A3.raw, NA)
  Y3.A4A3A1 <- ifelse(Y1.A4==0 & Y2.A4A3==0, Y3.A4A3A1.raw, NA)
  Y3.A4A3A2 <- ifelse(Y1.A4==0 & Y2.A4A3==0, Y3.A4A3A2.raw, NA)
  
  
  # stage III randomization probabilities
  P3.A1A2A3 <- rep(NA, n)
  P3.A1A2A4 <- rep(NA, n)
  P3.A1A3A2 <- rep(NA, n)
  P3.A1A3A4 <- rep(NA, n)
  P3.A1A4A2 <- rep(NA, n)
  P3.A1A4A3 <- rep(NA, n)
  P3.A2A1A3 <- rep(NA, n)
  P3.A2A1A4 <- rep(NA, n)
  P3.A2A3A1 <- rep(NA, n)
  P3.A2A3A4 <- rep(NA, n)
  P3.A2A4A1 <- rep(NA, n)
  P3.A2A4A3 <- rep(NA, n)
  P3.A3A1A2 <- rep(NA, n)
  P3.A3A1A4 <- rep(NA, n)
  P3.A3A2A1 <- rep(NA, n)
  P3.A3A2A4 <- rep(NA, n)
  P3.A3A4A1 <- rep(NA, n)
  P3.A3A4A2 <- rep(NA, n)
  P3.A4A1A2 <- rep(NA, n)
  P3.A4A1A3 <- rep(NA, n)
  P3.A4A2A1 <- rep(NA, n)
  P3.A4A2A3 <- rep(NA, n)
  P3.A4A3A1 <- rep(NA, n)
  P3.A4A3A2 <- rep(NA, n)
  
  
  # treatment assignment A in stage III
  A_3 <- rep(NA, n)
  I3.A1 <- rep(NA, n)
  I3.A2 <- rep(NA, n)
  I3.A3 <- rep(NA, n)
  I3.A4 <- rep(NA, n)
  
  # observed outcome Y3 in stage III
  Y3 <- rep(NA, n)
  
  # overall observed outcome Y
  Y <- rep(NA, n)
  
  # i=1,...,n0: for those who did not respond in stage I and stage II, equally randomized to the rest treatments
  if(n0.burn<=n) {
    for (i in 1:n0.burn) {
      if (I1.A1[i]==1 & I2.A2[i]==1 & Y1.A1[i]==0 & Y2.A1A2[i]==0) {
        P3.A1A2A3[i] <- 1/2
        P3.A1A2A4[i] <- 1/2
        A_3[i] <- sample(c(3,4), size = 1, replace = TRUE, prob = c(P3.A1A2A3[i],P3.A1A2A4[i]))
      } else if (I1.A1[i]==1 & I2.A3[i]==1 & Y1.A1[i]==0 & Y2.A1A3[i]==0) {
        P3.A1A3A2[i] <- 1/2
        P3.A1A3A4[i] <- 1/2
        A_3[i] <- sample(c(2,4), size = 1, replace = TRUE, prob = c(P3.A1A3A2[i],P3.A1A3A4[i]))
      } else if (I1.A1[i]==1 & I2.A4[i]==1 & Y1.A1[i]==0 & Y2.A1A4[i]==0) {
        P3.A1A4A2[i] <- 1/2
        P3.A1A4A3[i] <- 1/2
        A_3[i] <- sample(c(2,3), size = 1, replace = TRUE, prob = c(P3.A1A4A2[i],P3.A1A4A3[i]))
      } else if (I1.A2[i]==1 & I2.A1[i]==1 & Y1.A2[i]==0 & Y2.A2A1[i]==0) {
        P3.A2A1A3[i] <- 1/2
        P3.A2A1A4[i] <- 1/2
        A_3[i] <- sample(c(3,4), size = 1, replace = TRUE, prob = c(P3.A2A1A3[i],P3.A2A1A4[i]))
      } else if (I1.A2[i]==1 & I2.A3[i]==1 & Y1.A2[i]==0 & Y2.A2A3[i]==0){
        P3.A2A3A1[i] <- 1/2
        P3.A2A3A4[i] <- 1/2
        A_3[i] <- sample(c(1,4), size = 1, replace = TRUE, prob = c(P3.A2A3A1[i],P3.A2A3A4[i]))
      } else if (I1.A2[i]==1 & I2.A4[i]==1 & Y1.A2[i]==0 & Y2.A2A4[i]==0) {
        P3.A2A4A1[i] <- 1/2
        P3.A2A4A3[i] <- 1/2
        A_3[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P3.A2A4A1[i],P3.A2A4A3[i]))
      } else if (I1.A3[i]==1 & I2.A1[i]==1 & Y1.A3[i]==0 & Y2.A3A1[i]==0) {
        P3.A3A1A2[i] <- 1/2
        P3.A3A1A4[i] <- 1/2
        A_3[i] <- sample(c(2,4), size = 1, replace = TRUE, prob = c(P3.A3A1A2[i],P3.A3A1A4[i]))
      } else if (I1.A3[i]==1 & I2.A2[i]==1 & Y1.A3[i]==0 & Y2.A3A2[i]==0) {
        P3.A3A2A1[i] <- 1/2
        P3.A3A2A4[i] <- 1/2
        A_3[i] <- sample(c(1,4), size = 1, replace = TRUE, prob = c(P3.A3A2A1[i],P3.A3A2A4[i]))
      } else if (I1.A3[i]==1 & I2.A4[i]==1 & Y1.A3[i]==0 & Y2.A3A4[i]==0) {
        P3.A3A4A1[i] <- 1/2
        P3.A3A4A2[i] <- 1/2
        A_3[i] <- sample(c(1,2), size = 1, replace = TRUE, prob = c(P3.A3A4A1[i],P3.A3A4A2[i]))
      } else if (I1.A4[i]==1 & I2.A1[i]==1 & Y1.A4[i]==0 & Y2.A4A1[i]==0) {
        P3.A4A1A2[i] <- 1/2
        P3.A4A1A3[i] <- 1/2
        A_3[i] <- sample(c(2,3), size = 1, replace = TRUE, prob = c(P3.A4A1A2[i],P3.A4A1A3[i]))
      } else if (I1.A4[i]==1 & I2.A2[i]==1 & Y1.A4[i]==0 & Y2.A4A2[i]==0) {
        P3.A4A2A1[i] <- 1/2
        P3.A4A2A3[i] <- 1/2
        A_3[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P3.A4A2A1[i],P3.A4A2A3[i]))
      } else if (I1.A4[i]==1 & I2.A3[i]==1 & Y1.A4[i]==0 & Y2.A4A3[i]==0) {
        P3.A4A3A1[i] <- 1/2
        P3.A4A3A2[i] <- 1/2
        A_3[i] <- sample(c(1,2), size = 1, replace = TRUE, prob = c(P3.A4A3A1[i],P3.A4A3A2[i]))
      }
      
      I3.A1[i]<- ifelse(A_3[i]==1, 1, 0)
      I3.A2[i]<- ifelse(A_3[i]==2, 1, 0)
      I3.A3[i]<- ifelse(A_3[i]==3, 1, 0)
      I3.A4[i]<- ifelse(A_3[i]==4, 1, 0)
      
      
      
      Y3[i] <- ifelse(Y1[i]==0 & Y2[i]==0, 
                      sum(c(Y3.A1A2A3[i]*I3.A3[i]*I2.A2[i]*I1.A1[i],Y3.A1A2A4[i]*I3.A4[i]*I2.A2[i]*I1.A1[i],
                            Y3.A1A3A2[i]*I3.A2[i]*I2.A3[i]*I1.A1[i],Y3.A1A3A4[i]*I3.A4[i]*I2.A3[i]*I1.A1[i],
                            Y3.A1A4A2[i]*I3.A2[i]*I2.A4[i]*I1.A1[i],Y3.A1A4A3[i]*I3.A3[i]*I2.A4[i]*I1.A1[i],
                            Y3.A2A1A3[i]*I3.A3[i]*I2.A1[i]*I1.A2[i],Y3.A2A1A4[i]*I3.A4[i]*I2.A1[i]*I1.A2[i],
                            Y3.A2A3A1[i]*I3.A1[i]*I2.A3[i]*I1.A2[i],Y3.A2A3A4[i]*I3.A4[i]*I2.A3[i]*I1.A2[i],
                            Y3.A2A4A1[i]*I3.A1[i]*I2.A4[i]*I1.A2[i],Y3.A2A4A3[i]*I3.A3[i]*I2.A4[i]*I1.A2[i],
                            Y3.A3A1A2[i]*I3.A2[i]*I2.A1[i]*I1.A3[i],Y3.A3A1A4[i]*I3.A4[i]*I2.A1[i]*I1.A3[i],
                            Y3.A3A2A1[i]*I3.A1[i]*I2.A2[i]*I1.A3[i],Y3.A3A2A4[i]*I3.A4[i]*I2.A2[i]*I1.A3[i],
                            Y3.A3A4A1[i]*I3.A1[i]*I2.A4[i]*I1.A3[i],Y3.A3A4A2[i]*I3.A2[i]*I2.A4[i]*I1.A3[i],
                            Y3.A4A1A2[i]*I3.A2[i]*I2.A1[i]*I1.A4[i],Y3.A4A1A3[i]*I3.A3[i]*I2.A1[i]*I1.A4[i],
                            Y3.A4A2A1[i]*I3.A1[i]*I2.A2[i]*I1.A4[i],Y3.A4A2A3[i]*I3.A3[i]*I2.A2[i]*I1.A4[i],
                            Y3.A4A3A1[i]*I3.A1[i]*I2.A3[i]*I1.A4[i],Y3.A4A3A2[i]*I3.A2[i]*I2.A3[i]*I1.A4[i]), na.rm = T), NA)
      
      Y[i] <- ifelse(Y1[i]==1, Y1[i], ifelse(Y2[i]==1, Y2[i], Y3[i]))
      
    }
  }
  
  # number(n2.burn) and proportion(p2.burn): third stage burn-in sample
  n2.burn <- round(n*p2.burn, 0)
  
  if(n0.burn<n2.burn) {
    # i=n0+1,...,n2: update the randomization probability using stage I observed response rate
    for (i in (n0.burn+1):n2.burn) {
      # P(Y1.A1=1)
      Y1.A1hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])
      # P(Y1.A2=1)
      Y1.A2hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])
      # P(Y1.A3=1)
      Y1.A3hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])
      # P(Y1.A4=1)
      Y1.A4hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])
      
      if (c=="n.2N") {
        c4 <- i/(2*n)
      } else if (c=="n.N") {
        c4 <- i/n
      } else {
        c4 <- c
      }
      
      
      r4.A1A2A3 <- ((Y1.A3hat)^c4)/((Y1.A3hat)^c4+(Y1.A4hat)^c4)
      r4.A1A2A4 <- ((Y1.A4hat)^c4)/((Y1.A3hat)^c4+(Y1.A4hat)^c4)
      r4.A1A3A2 <- ((Y1.A2hat)^c4)/((Y1.A2hat)^c4+(Y1.A4hat)^c4)
      r4.A1A3A4 <- ((Y1.A4hat)^c4)/((Y1.A2hat)^c4+(Y1.A4hat)^c4)
      r4.A1A4A2 <- ((Y1.A2hat)^c4)/((Y1.A2hat)^c4+(Y1.A3hat)^c4)
      r4.A1A4A3 <- ((Y1.A3hat)^c4)/((Y1.A2hat)^c4+(Y1.A3hat)^c4)
      r4.A2A1A3 <- ((Y1.A3hat)^c4)/((Y1.A3hat)^c4+(Y1.A4hat)^c4)
      r4.A2A1A4 <- ((Y1.A4hat)^c4)/((Y1.A3hat)^c4+(Y1.A4hat)^c4)
      r4.A2A3A1 <- ((Y1.A1hat)^c4)/((Y1.A1hat)^c4+(Y1.A4hat)^c4)
      r4.A2A3A4 <- ((Y1.A4hat)^c4)/((Y1.A1hat)^c4+(Y1.A4hat)^c4)
      r4.A2A4A1 <- ((Y1.A1hat)^c4)/((Y1.A1hat)^c4+(Y1.A3hat)^c4)
      r4.A2A4A3 <- ((Y1.A3hat)^c4)/((Y1.A1hat)^c4+(Y1.A3hat)^c4)
      r4.A3A1A2 <- ((Y1.A2hat)^c4)/((Y1.A2hat)^c4+(Y1.A4hat)^c4)
      r4.A3A1A4 <- ((Y1.A4hat)^c4)/((Y1.A2hat)^c4+(Y1.A4hat)^c4)
      r4.A3A2A1 <- ((Y1.A1hat)^c4)/((Y1.A1hat)^c4+(Y1.A4hat)^c4)
      r4.A3A2A4 <- ((Y1.A4hat)^c4)/((Y1.A1hat)^c4+(Y1.A4hat)^c4)
      r4.A3A4A1 <- ((Y1.A1hat)^c4)/((Y1.A1hat)^c4+(Y1.A2hat)^c4)
      r4.A3A4A2 <- ((Y1.A2hat)^c4)/((Y1.A1hat)^c4+(Y1.A2hat)^c4)
      r4.A4A1A2 <- ((Y1.A2hat)^c4)/((Y1.A2hat)^c4+(Y1.A3hat)^c4)
      r4.A4A1A3 <- ((Y1.A3hat)^c4)/((Y1.A2hat)^c4+(Y1.A3hat)^c4)
      r4.A4A2A1 <- ((Y1.A1hat)^c4)/((Y1.A1hat)^c4+(Y1.A3hat)^c4)
      r4.A4A2A3 <- ((Y1.A3hat)^c4)/((Y1.A1hat)^c4+(Y1.A3hat)^c4)
      r4.A4A3A1 <- ((Y1.A1hat)^c4)/((Y1.A1hat)^c4+(Y1.A2hat)^c4)
      r4.A4A3A2 <- ((Y1.A2hat)^c4)/((Y1.A1hat)^c4+(Y1.A2hat)^c4)
      
      # AR(c,e)
      e4 <- e[4]
      
      if (I1.A1[i]==1 & I2.A2[i]==1 & Y1.A1[i]==0 & Y2.A1A2[i]==0) {
        P3.A1A2A3[i] <- ifelse(is.na(r4.A1A2A3), 1/2,  min(max(c(r4.A1A2A3, e4)), 1-e4))
        P3.A1A2A4[i] <- 1-P3.A1A2A3[i]
        A_3[i] <- sample(c(3,4), size = 1, replace = TRUE, prob = c(P3.A1A2A3[i],P3.A1A2A4[i]))
      } else if (I1.A1[i]==1 & I2.A3[i]==1 & Y1.A1[i]==0 & Y2.A1A3[i]==0) {
        P3.A1A3A2[i] <- ifelse(is.na(r4.A1A3A2), 1/2,  min(max(c(r4.A1A3A2, e4)), 1-e4))
        P3.A1A3A4[i] <- 1-P3.A1A3A2[i]
        A_3[i] <- sample(c(2,4), size = 1, replace = TRUE, prob = c(P3.A1A3A2[i],P3.A1A3A4[i]))
      } else if (I1.A1[i]==1 & I2.A4[i]==1 & Y1.A1[i]==0 & Y2.A1A4[i]==0) {
        P3.A1A4A2[i] <- ifelse(is.na(r4.A1A4A2), 1/2,  min(max(c(r4.A1A4A2, e4)), 1-e4))
        P3.A1A4A3[i] <- 1-P3.A1A4A2[i]
        A_3[i] <- sample(c(2,3), size = 1, replace = TRUE, prob = c(P3.A1A4A2[i],P3.A1A4A3[i]))
      } else if (I1.A2[i]==1 & I2.A1[i]==1 & Y1.A2[i]==0 & Y2.A2A1[i]==0) {
        P3.A2A1A3[i] <- ifelse(is.na(r4.A2A1A3), 1/2,  min(max(c(r4.A2A1A3, e4)), 1-e4))
        P3.A2A1A4[i] <- 1-P3.A2A1A3[i]
        A_3[i] <- sample(c(3,4), size = 1, replace = TRUE, prob = c(P3.A2A1A3[i],P3.A2A1A4[i]))
      } else if (I1.A2[i]==1 & I2.A3[i]==1 & Y1.A2[i]==0 & Y2.A2A3[i]==0){
        P3.A2A3A1[i] <- ifelse(is.na(r4.A2A3A1), 1/2,  min(max(c(r4.A2A3A1, e4)), 1-e4))
        P3.A2A3A4[i] <- 1-P3.A2A3A1[i]
        A_3[i] <- sample(c(1,4), size = 1, replace = TRUE, prob = c(P3.A2A3A1[i],P3.A2A3A4[i]))
      } else if (I1.A2[i]==1 & I2.A4[i]==1 & Y1.A2[i]==0 & Y2.A2A4[i]==0) {
        P3.A2A4A1[i] <- ifelse(is.na(r4.A2A4A1), 1/2,  min(max(c(r4.A2A4A1, e4)), 1-e4))
        P3.A2A4A3[i] <- 1-P3.A2A4A1[i]
        A_3[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P3.A2A4A1[i],P3.A2A4A3[i]))
      } else if (I1.A3[i]==1 & I2.A1[i]==1 & Y1.A3[i]==0 & Y2.A3A1[i]==0) {
        P3.A3A1A2[i] <- ifelse(is.na(r4.A3A1A2), 1/2,  min(max(c(r4.A3A1A2, e4)), 1-e4))
        P3.A3A1A4[i] <- 1-P3.A3A1A2[i]
        A_3[i] <- sample(c(2,4), size = 1, replace = TRUE, prob = c(P3.A3A1A2[i],P3.A3A1A4[i]))
      } else if (I1.A3[i]==1 & I2.A2[i]==1 & Y1.A3[i]==0 & Y2.A3A2[i]==0) {
        P3.A3A2A1[i] <- ifelse(is.na(r4.A3A2A1), 1/2,  min(max(c(r4.A3A2A1, e4)), 1-e4))
        P3.A3A2A4[i] <- 1-P3.A3A2A1[i] 
        A_3[i] <- sample(c(1,4), size = 1, replace = TRUE, prob = c(P3.A3A2A1[i],P3.A3A2A4[i]))
      } else if (I1.A3[i]==1 & I2.A4[i]==1 & Y1.A3[i]==0 & Y2.A3A4[i]==0) {
        P3.A3A4A1[i] <- ifelse(is.na(r4.A3A4A1), 1/2,  min(max(c(r4.A3A4A1, e4)), 1-e4))
        P3.A3A4A2[i] <- 1-P3.A3A4A1[i]
        A_3[i] <- sample(c(1,2), size = 1, replace = TRUE, prob = c(P3.A3A4A1[i],P3.A3A4A2[i]))
      } else if (I1.A4[i]==1 & I2.A1[i]==1 & Y1.A4[i]==0 & Y2.A4A1[i]==0) {
        P3.A4A1A2[i] <- ifelse(is.na(r4.A4A1A2), 1/2,  min(max(c(r4.A4A1A2, e4)), 1-e4))
        P3.A4A1A3[i] <- 1-P3.A4A1A2[i]
        A_3[i] <- sample(c(2,3), size = 1, replace = TRUE, prob = c(P3.A4A1A2[i],P3.A4A1A3[i]))
      } else if (I1.A4[i]==1 & I2.A2[i]==1 & Y1.A4[i]==0 & Y2.A4A2[i]==0) {
        P3.A4A2A1[i] <- ifelse(is.na(r4.A4A2A1), 1/2,  min(max(c(r4.A4A2A1, e4)), 1-e4))
        P3.A4A2A3[i] <- 1-P3.A4A2A1[i]
        A_3[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P3.A4A2A1[i],P3.A4A2A3[i]))
      } else if (I1.A4[i]==1 & I2.A3[i]==1 & Y1.A4[i]==0 & Y2.A4A3[i]==0) {
        P3.A4A3A1[i] <- ifelse(is.na(r4.A4A3A1), 1/2,  min(max(c(r4.A4A3A1, e4)), 1-e4))
        P3.A4A3A2[i] <- 1-P3.A4A3A1[i]
        A_3[i] <- sample(c(1,2), size = 1, replace = TRUE, prob = c(P3.A4A3A1[i],P3.A4A3A2[i]))
      }
      
      I3.A1[i]<- ifelse(A_3[i]==1, 1, 0)
      I3.A2[i]<- ifelse(A_3[i]==2, 1, 0)
      I3.A3[i]<- ifelse(A_3[i]==3, 1, 0)
      I3.A4[i]<- ifelse(A_3[i]==4, 1, 0)
      
      
      
      Y3[i] <- ifelse(Y1[i]==0 & Y2[i]==0, 
                      sum(c(Y3.A1A2A3[i]*I3.A3[i]*I2.A2[i]*I1.A1[i],Y3.A1A2A4[i]*I3.A4[i]*I2.A2[i]*I1.A1[i],
                            Y3.A1A3A2[i]*I3.A2[i]*I2.A3[i]*I1.A1[i],Y3.A1A3A4[i]*I3.A4[i]*I2.A3[i]*I1.A1[i],
                            Y3.A1A4A2[i]*I3.A2[i]*I2.A4[i]*I1.A1[i],Y3.A1A4A3[i]*I3.A3[i]*I2.A4[i]*I1.A1[i],
                            Y3.A2A1A3[i]*I3.A3[i]*I2.A1[i]*I1.A2[i],Y3.A2A1A4[i]*I3.A4[i]*I2.A1[i]*I1.A2[i],
                            Y3.A2A3A1[i]*I3.A1[i]*I2.A3[i]*I1.A2[i],Y3.A2A3A4[i]*I3.A4[i]*I2.A3[i]*I1.A2[i],
                            Y3.A2A4A1[i]*I3.A1[i]*I2.A4[i]*I1.A2[i],Y3.A2A4A3[i]*I3.A3[i]*I2.A4[i]*I1.A2[i],
                            Y3.A3A1A2[i]*I3.A2[i]*I2.A1[i]*I1.A3[i],Y3.A3A1A4[i]*I3.A4[i]*I2.A1[i]*I1.A3[i],
                            Y3.A3A2A1[i]*I3.A1[i]*I2.A2[i]*I1.A3[i],Y3.A3A2A4[i]*I3.A4[i]*I2.A2[i]*I1.A3[i],
                            Y3.A3A4A1[i]*I3.A1[i]*I2.A4[i]*I1.A3[i],Y3.A3A4A2[i]*I3.A2[i]*I2.A4[i]*I1.A3[i],
                            Y3.A4A1A2[i]*I3.A2[i]*I2.A1[i]*I1.A4[i],Y3.A4A1A3[i]*I3.A3[i]*I2.A1[i]*I1.A4[i],
                            Y3.A4A2A1[i]*I3.A1[i]*I2.A2[i]*I1.A4[i],Y3.A4A2A3[i]*I3.A3[i]*I2.A2[i]*I1.A4[i],
                            Y3.A4A3A1[i]*I3.A1[i]*I2.A3[i]*I1.A4[i],Y3.A4A3A2[i]*I3.A2[i]*I2.A3[i]*I1.A4[i]), na.rm = T), NA)
      
      Y[i] <- ifelse(Y1[i]==1, Y1[i], ifelse(Y2[i]==1, Y2[i], Y3[i]))
    }
  }
  
  if(n0.burn<=n2.burn & n2.burn<n) {
    # i=n2+1,...,n: update the randomization probability using following different types of method
    for (i in (n2.burn+1):n) {
      if (AR=="AR1") {
        # conditional probability
        # G-computation
        
        # P(Y3.A1A2A3=1|Y1.A1=0,Y2.A1A2=0)
        Y3.A1A2A3hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        # P(Y3.A1A2A4=1|Y1.A1=0,Y2.A1A2=0)
        Y3.A1A2A4hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Y3.A1A3A2=1|Y1.A1=0,Y2.A1A3=0)
        Y3.A1A3A2hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        # P(Y3.A1A3A4=1|Y1.A1=0,Y2.A1A3=0)
        Y3.A1A3A4hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Y3.A1A4A2=1|Y1.A1=0,Y2.A1A4=0)
        Y3.A1A4A2hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        # P(Y3.A1A4A3=1|Y1.A1=0,Y2.A1A4=0)
        Y3.A1A4A3hat <- sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        
        # P(Y3.A2A1A3=1|Y1.A2=0,Y2.A2A1=0)
        Y3.A2A1A3hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        # P(Y3.A2A1A4=1|Y1.A2=0,Y2.A2A1=0)
        Y3.A2A1A4hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Y3.A2A3A1=1|Y1.A2=0,Y2.A2A3=0)
        Y3.A2A3A1hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Y3.A2A3A4=1|Y1.A2=0,Y2.A2A3=0)
        Y3.A2A3A4hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Y3.A2A4A1=1|Y1.A2=0,Y2.A2A4=0)
        Y3.A2A4A1hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Y3.A2A4A3=1|Y1.A2=0,Y2.A2A4=0)
        Y3.A2A4A3hat <- sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        
        # P(Y3.A3A1A2=1|Y1.A3=0,Y2.A3A1=0)
        Y3.A3A1A2hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        # P(Y3.A3A1A4=1|Y1.A3=0,Y2.A3A1=0)
        Y3.A3A1A4hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Y3.A3A2A1=1|Y1.A3=0,Y2.A3A2=0)
        Y3.A3A2A1hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Y3.A3A2A1=1|Y1.A3=0,Y2.A3A2=0)
        Y3.A3A2A4hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Y3.A3A4A1=1|Y1.A3=0,Y2.A3A4=0)
        Y3.A3A4A1hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Y3.A3A4A2=1|Y1.A3=0,Y2.A3A4=0)
        Y3.A3A4A2hat <- sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        
        # P(Y3.A4A1A2=1|Y1.A4=0,Y2.A4A1=0)
        Y3.A4A1A2hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        # P(Y3.A4A1A3=1|Y1.A4=0,Y2.A4A1=0)
        Y3.A4A1A3hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        # P(Y3.A4A2A1=1|Y1.A4=0,Y2.A4A2=0)
        Y3.A4A2A1hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Y3.A4A2A3=1|Y1.A4=0,Y2.A4A2=0)
        Y3.A4A2A3hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        # P(Y3.A4A3A1=1|Y1.A4=0,Y2.A4A3=0)
        Y3.A4A3A1hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Y3.A4A3A2=1|Y1.A4=0,Y2.A4A3=0)
        Y3.A4A3A2hat <- sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        
        if (c=="n.2N") {
          c5 <- i/(2*n)
        } else if (c=="n.N") {
          c5 <- i/n
        } else {
          c5 <- c
        }
        
        r5.A1A2A3 <- ((Y3.A1A2A3hat)^c5)/((Y3.A1A2A3hat)^c5+(Y3.A1A2A4hat)^c5)
        r5.A1A2A4 <- ((Y3.A1A2A4hat)^c5)/((Y3.A1A2A3hat)^c5+(Y3.A1A2A4hat)^c5)
        r5.A1A3A2 <- ((Y3.A1A3A2hat)^c5)/((Y3.A1A3A2hat)^c5+(Y3.A1A3A4hat)^c5)
        r5.A1A3A4 <- ((Y3.A1A3A4hat)^c5)/((Y3.A1A3A2hat)^c5+(Y3.A1A3A4hat)^c5)
        r5.A1A4A2 <- ((Y3.A1A4A2hat)^c5)/((Y3.A1A4A2hat)^c5+(Y3.A1A4A3hat)^c5)
        r5.A1A4A3 <- ((Y3.A1A4A3hat)^c5)/((Y3.A1A4A2hat)^c5+(Y3.A1A4A3hat)^c5)
        r5.A2A1A3 <- ((Y3.A2A1A3hat)^c5)/((Y3.A2A1A3hat)^c5+(Y3.A2A1A4hat)^c5)
        r5.A2A1A4 <- ((Y3.A2A1A4hat)^c5)/((Y3.A2A1A3hat)^c5+(Y3.A2A1A4hat)^c5)
        r5.A2A3A1 <- ((Y3.A2A3A1hat)^c5)/((Y3.A2A3A1hat)^c5+(Y3.A2A3A4hat)^c5)
        r5.A2A3A4 <- ((Y3.A2A3A4hat)^c5)/((Y3.A2A3A1hat)^c5+(Y3.A2A3A4hat)^c5)
        r5.A2A4A1 <- ((Y3.A2A4A1hat)^c5)/((Y3.A2A4A1hat)^c5+(Y3.A2A4A3hat)^c5)
        r5.A2A4A3 <- ((Y3.A2A4A3hat)^c5)/((Y3.A2A4A1hat)^c5+(Y3.A2A4A3hat)^c5)
        r5.A3A1A2 <- ((Y3.A3A1A2hat)^c5)/((Y3.A3A1A2hat)^c5+(Y3.A3A1A4hat)^c5)
        r5.A3A1A4 <- ((Y3.A3A1A4hat)^c5)/((Y3.A3A1A2hat)^c5+(Y3.A3A1A4hat)^c5)
        r5.A3A2A1 <- ((Y3.A3A2A1hat)^c5)/((Y3.A3A2A1hat)^c5+(Y3.A3A2A4hat)^c5)
        r5.A3A2A4 <- ((Y3.A3A2A4hat)^c5)/((Y3.A3A2A1hat)^c5+(Y3.A3A2A4hat)^c5)
        r5.A3A4A1 <- ((Y3.A3A4A1hat)^c5)/((Y3.A3A4A1hat)^c5+(Y3.A3A4A2hat)^c5)
        r5.A3A4A2 <- ((Y3.A3A4A2hat)^c5)/((Y3.A3A4A1hat)^c5+(Y3.A3A4A2hat)^c5)
        r5.A4A1A2 <- ((Y3.A4A1A2hat)^c5)/((Y3.A4A1A2hat)^c5+(Y3.A4A1A3hat)^c5)
        r5.A4A1A3 <- ((Y3.A4A1A3hat)^c5)/((Y3.A4A1A2hat)^c5+(Y3.A4A1A3hat)^c5)
        r5.A4A2A1 <- ((Y3.A4A2A1hat)^c5)/((Y3.A4A2A1hat)^c5+(Y3.A4A2A3hat)^c5)
        r5.A4A2A3 <- ((Y3.A4A2A3hat)^c5)/((Y3.A4A2A1hat)^c5+(Y3.A4A2A3hat)^c5)
        r5.A4A3A1 <- ((Y3.A4A3A1hat)^c5)/((Y3.A4A3A1hat)^c5+(Y3.A4A3A2hat)^c5)
        r5.A4A3A2 <- ((Y3.A4A3A2hat)^c5)/((Y3.A4A3A1hat)^c5+(Y3.A4A3A2hat)^c5)
        
        # AR(c,e)
        e5 <- e[5]
        
        if (I1.A1[i]==1 & I2.A2[i]==1 & Y1.A1[i]==0 & Y2.A1A2[i]==0) {
          P3.A1A2A3[i] <- ifelse(is.na(r5.A1A2A3), 1/2,  min(max(c(r5.A1A2A3, e5)), 1-e5))
          P3.A1A2A4[i] <- 1-P3.A1A2A3[i]
          A_3[i] <- sample(c(3,4), size = 1, replace = TRUE, prob = c(P3.A1A2A3[i],P3.A1A2A4[i]))
        } else if (I1.A1[i]==1 & I2.A3[i]==1 & Y1.A1[i]==0 & Y2.A1A3[i]==0) {
          P3.A1A3A2[i] <- ifelse(is.na(r5.A1A3A2), 1/2,  min(max(c(r5.A1A3A2, e5)), 1-e5))
          P3.A1A3A4[i] <- 1-P3.A1A3A2[i]
          A_3[i] <- sample(c(2,4), size = 1, replace = TRUE, prob = c(P3.A1A3A2[i],P3.A1A3A4[i]))
        } else if (I1.A1[i]==1 & I2.A4[i]==1 & Y1.A1[i]==0 & Y2.A1A4[i]==0) {
          P3.A1A4A2[i] <- ifelse(is.na(r5.A1A4A2), 1/2,  min(max(c(r5.A1A4A2, e5)), 1-e5))
          P3.A1A4A3[i] <- 1-P3.A1A4A2[i]
          A_3[i] <- sample(c(2,3), size = 1, replace = TRUE, prob = c(P3.A1A4A2[i],P3.A1A4A3[i]))
        } else if (I1.A2[i]==1 & I2.A1[i]==1 & Y1.A2[i]==0 & Y2.A2A1[i]==0) {
          P3.A2A1A3[i] <- ifelse(is.na(r5.A2A1A3), 1/2,  min(max(c(r5.A2A1A3, e5)), 1-e5))
          P3.A2A1A4[i] <- 1-P3.A2A1A3[i]
          A_3[i] <- sample(c(3,4), size = 1, replace = TRUE, prob = c(P3.A2A1A3[i],P3.A2A1A4[i]))
        } else if (I1.A2[i]==1 & I2.A3[i]==1 & Y1.A2[i]==0 & Y2.A2A3[i]==0){
          P3.A2A3A1[i] <- ifelse(is.na(r5.A2A3A1), 1/2,  min(max(c(r5.A2A3A1, e5)), 1-e5))
          P3.A2A3A4[i] <- 1-P3.A2A3A1[i]
          A_3[i] <- sample(c(1,4), size = 1, replace = TRUE, prob = c(P3.A2A3A1[i],P3.A2A3A4[i]))
        } else if (I1.A2[i]==1 & I2.A4[i]==1 & Y1.A2[i]==0 & Y2.A2A4[i]==0) {
          P3.A2A4A1[i] <- ifelse(is.na(r5.A2A4A1), 1/2,  min(max(c(r5.A2A4A1, e5)), 1-e5))
          P3.A2A4A3[i] <- 1-P3.A2A4A1[i]
          A_3[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P3.A2A4A1[i],P3.A2A4A3[i]))
        } else if (I1.A3[i]==1 & I2.A1[i]==1 & Y1.A3[i]==0 & Y2.A3A1[i]==0) {
          P3.A3A1A2[i] <- ifelse(is.na(r5.A3A1A2), 1/2,  min(max(c(r5.A3A1A2, e5)), 1-e5))
          P3.A3A1A4[i] <- 1-P3.A3A1A2[i]
          A_3[i] <- sample(c(2,4), size = 1, replace = TRUE, prob = c(P3.A3A1A2[i],P3.A3A1A4[i]))
        } else if (I1.A3[i]==1 & I2.A2[i]==1 & Y1.A3[i]==0 & Y2.A3A2[i]==0) {
          P3.A3A2A1[i] <- ifelse(is.na(r5.A3A2A1), 1/2,  min(max(c(r5.A3A2A1, e5)), 1-e5))
          P3.A3A2A4[i] <- 1-P3.A3A2A1[i] 
          A_3[i] <- sample(c(1,4), size = 1, replace = TRUE, prob = c(P3.A3A2A1[i],P3.A3A2A4[i]))
        } else if (I1.A3[i]==1 & I2.A4[i]==1 & Y1.A3[i]==0 & Y2.A3A4[i]==0) {
          P3.A3A4A1[i] <- ifelse(is.na(r5.A3A4A1), 1/2,  min(max(c(r5.A3A4A1, e5)), 1-e5))
          P3.A3A4A2[i] <- 1-P3.A3A4A1[i]
          A_3[i] <- sample(c(1,2), size = 1, replace = TRUE, prob = c(P3.A3A4A1[i],P3.A3A4A2[i]))
        } else if (I1.A4[i]==1 & I2.A1[i]==1 & Y1.A4[i]==0 & Y2.A4A1[i]==0) {
          P3.A4A1A2[i] <- ifelse(is.na(r5.A4A1A2), 1/2,  min(max(c(r5.A4A1A2, e5)), 1-e5))
          P3.A4A1A3[i] <- 1-P3.A4A1A2[i]
          A_3[i] <- sample(c(2,3), size = 1, replace = TRUE, prob = c(P3.A4A1A2[i],P3.A4A1A3[i]))
        } else if (I1.A4[i]==1 & I2.A2[i]==1 & Y1.A4[i]==0 & Y2.A4A2[i]==0) {
          P3.A4A2A1[i] <- ifelse(is.na(r5.A4A2A1), 1/2,  min(max(c(r5.A4A2A1, e5)), 1-e5))
          P3.A4A2A3[i] <- 1-P3.A4A2A1[i]
          A_3[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P3.A4A2A1[i],P3.A4A2A3[i]))
        } else if (I1.A4[i]==1 & I2.A3[i]==1 & Y1.A4[i]==0 & Y2.A4A3[i]==0) {
          P3.A4A3A1[i] <- ifelse(is.na(r5.A4A3A1), 1/2,  min(max(c(r5.A4A3A1, e5)), 1-e5))
          P3.A4A3A2[i] <- 1-P3.A4A3A1[i]
          A_3[i] <- sample(c(1,2), size = 1, replace = TRUE, prob = c(P3.A4A3A1[i],P3.A4A3A2[i]))
        }
        
        I3.A1[i]<- ifelse(A_3[i]==1, 1, 0)
        I3.A2[i]<- ifelse(A_3[i]==2, 1, 0)
        I3.A3[i]<- ifelse(A_3[i]==3, 1, 0)
        I3.A4[i]<- ifelse(A_3[i]==4, 1, 0)
        
        
        
        Y3[i] <- ifelse(Y1[i]==0 & Y2[i]==0, 
                        sum(c(Y3.A1A2A3[i]*I3.A3[i]*I2.A2[i]*I1.A1[i],Y3.A1A2A4[i]*I3.A4[i]*I2.A2[i]*I1.A1[i],
                              Y3.A1A3A2[i]*I3.A2[i]*I2.A3[i]*I1.A1[i],Y3.A1A3A4[i]*I3.A4[i]*I2.A3[i]*I1.A1[i],
                              Y3.A1A4A2[i]*I3.A2[i]*I2.A4[i]*I1.A1[i],Y3.A1A4A3[i]*I3.A3[i]*I2.A4[i]*I1.A1[i],
                              Y3.A2A1A3[i]*I3.A3[i]*I2.A1[i]*I1.A2[i],Y3.A2A1A4[i]*I3.A4[i]*I2.A1[i]*I1.A2[i],
                              Y3.A2A3A1[i]*I3.A1[i]*I2.A3[i]*I1.A2[i],Y3.A2A3A4[i]*I3.A4[i]*I2.A3[i]*I1.A2[i],
                              Y3.A2A4A1[i]*I3.A1[i]*I2.A4[i]*I1.A2[i],Y3.A2A4A3[i]*I3.A3[i]*I2.A4[i]*I1.A2[i],
                              Y3.A3A1A2[i]*I3.A2[i]*I2.A1[i]*I1.A3[i],Y3.A3A1A4[i]*I3.A4[i]*I2.A1[i]*I1.A3[i],
                              Y3.A3A2A1[i]*I3.A1[i]*I2.A2[i]*I1.A3[i],Y3.A3A2A4[i]*I3.A4[i]*I2.A2[i]*I1.A3[i],
                              Y3.A3A4A1[i]*I3.A1[i]*I2.A4[i]*I1.A3[i],Y3.A3A4A2[i]*I3.A2[i]*I2.A4[i]*I1.A3[i],
                              Y3.A4A1A2[i]*I3.A2[i]*I2.A1[i]*I1.A4[i],Y3.A4A1A3[i]*I3.A3[i]*I2.A1[i]*I1.A4[i],
                              Y3.A4A2A1[i]*I3.A1[i]*I2.A2[i]*I1.A4[i],Y3.A4A2A3[i]*I3.A3[i]*I2.A2[i]*I1.A4[i],
                              Y3.A4A3A1[i]*I3.A1[i]*I2.A3[i]*I1.A4[i],Y3.A4A3A2[i]*I3.A2[i]*I2.A3[i]*I1.A4[i]), na.rm = T), NA)
        
        Y[i] <- ifelse(Y1[i]==1, Y1[i], ifelse(Y2[i]==1, Y2[i], Y3[i]))
        
      } else if (AR=="AR2") {
        # DTR response rate
        # G-computation
        # P(Yd.A1A2A3=1)
        Yd.A1A2A3hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*(1-sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T))*
          sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        # P(Yd.A1A2A4=1)
        Yd.A1A2A4hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*(1-sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T))*
          sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Yd.A1A3A2=1)
        Yd.A1A3A2hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*(1-sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T))*
          sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        # P(Yd.A1A3A4=1)
        Yd.A1A3A4hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*(1-sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T))*
          sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Yd.A1A4A2=1)
        Yd.A1A4A2hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*(1-sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T))*
          sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        # P(Yd.A1A4A3=1)
        Yd.A1A4A3hat <- sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)])+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)+
          (1-sum(I1.A1[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A1[1:(i-1)]))*(1-sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T))*
          sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A1[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        
        # P(Yd.A2A1A3=1)
        Yd.A2A1A3hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*(1-sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T))*
          sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        # P(Yd.A2A1A4=1)
        Yd.A2A1A4hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*(1-sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T))*
          sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Yd.A2A3A1=1)
        Yd.A2A3A1hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*(1-sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T))*
          sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Yd.A2A3A4=1)
        Yd.A2A3A4hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*(1-sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T))*
          sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Yd.A2A4A1=1)
        Yd.A2A4A1hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*(1-sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T))*
          sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Yd.A2A4A3=1)
        Yd.A2A4A3hat <- sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)])+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)+
          (1-sum(I1.A2[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A2[1:(i-1)]))*(1-sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T))*
          sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A2[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        
        # P(Yd.A3A1A2=1)
        Yd.A3A1A2hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*(1-sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T))*
          sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        # P(Yd.A3A1A4=1)
        Yd.A3A1A4hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*(1-sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T))*
          sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Yd.A3A2A1=1)
        Yd.A3A2A1hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*(1-sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T))*
          sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Yd.A3A2A1=1)
        Yd.A3A2A4hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*(1-sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T))*
          sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A4[1:(i-1)],na.rm=T)
        # P(Yd.A3A4A1=1)
        Yd.A3A4A1hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*(1-sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T))*
          sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Yd.A3A4A2=1)
        Yd.A3A4A2hat <- sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)])+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T)+
          (1-sum(I1.A3[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A3[1:(i-1)]))*(1-sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)],na.rm=T))*
          sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A3[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A4[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        
        # P(Yd.A4A1A2=1)
        Yd.A4A1A2hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*(1-sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T))*
          sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        # P(Yd.A4A1A3=1)
        Yd.A4A1A3hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T)+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*(1-sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)],na.rm=T))*
          sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A1[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        # P(Yd.A4A2A1=1)
        Yd.A4A2A1hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*(1-sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T))*
          sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Yd.A4A2A3=1)
        Yd.A4A2A3hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T)+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*(1-sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)],na.rm=T))*
          sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A2[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A3[1:(i-1)],na.rm=T)
        # P(Yd.A4A3A1=1)
        Yd.A4A3A1hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*(1-sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T))*
          sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A1[1:(i-1)],na.rm=T)
        # P(Yd.A4A3A2=1)
        Yd.A4A3A2hat <- sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)])+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T)+
          (1-sum(I1.A4[1:(i-1)]*Y1[1:(i-1)])/sum(I1.A4[1:(i-1)]))*(1-sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*Y2[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)],na.rm=T))*
          sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)]*Y3[1:(i-1)],na.rm=T)/sum(I1.A4[1:(i-1)]*(1-Y1[1:(i-1)])*I2.A3[1:(i-1)]*(1-Y2[1:(i-1)])*I3.A2[1:(i-1)],na.rm=T)
        
        if (c=="n.2N") {
          c5 <- i/(2*n)
        } else if (c=="n.N") {
          c5 <- i/n
        } else {
          c5 <- c
        }
        
        r5.A1A2A3 <- ((Yd.A1A2A3hat)^c5)/((Yd.A1A2A3hat)^c5+(Yd.A1A2A4hat)^c5)
        r5.A1A2A4 <- ((Yd.A1A2A4hat)^c5)/((Yd.A1A2A3hat)^c5+(Yd.A1A2A4hat)^c5)
        r5.A1A3A2 <- ((Yd.A1A3A2hat)^c5)/((Yd.A1A3A2hat)^c5+(Yd.A1A3A4hat)^c5)
        r5.A1A3A4 <- ((Yd.A1A3A4hat)^c5)/((Yd.A1A3A2hat)^c5+(Yd.A1A3A4hat)^c5)
        r5.A1A4A2 <- ((Yd.A1A4A2hat)^c5)/((Yd.A1A4A2hat)^c5+(Yd.A1A4A3hat)^c5)
        r5.A1A4A3 <- ((Yd.A1A4A3hat)^c5)/((Yd.A1A4A2hat)^c5+(Yd.A1A4A3hat)^c5)
        r5.A2A1A3 <- ((Yd.A2A1A3hat)^c5)/((Yd.A2A1A3hat)^c5+(Yd.A2A1A4hat)^c5)
        r5.A2A1A4 <- ((Yd.A2A1A4hat)^c5)/((Yd.A2A1A3hat)^c5+(Yd.A2A1A4hat)^c5)
        r5.A2A3A1 <- ((Yd.A2A3A1hat)^c5)/((Yd.A2A3A1hat)^c5+(Yd.A2A3A4hat)^c5)
        r5.A2A3A4 <- ((Yd.A2A3A4hat)^c5)/((Yd.A2A3A1hat)^c5+(Yd.A2A3A4hat)^c5)
        r5.A2A4A1 <- ((Yd.A2A4A1hat)^c5)/((Yd.A2A4A1hat)^c5+(Yd.A2A4A3hat)^c5)
        r5.A2A4A3 <- ((Yd.A2A4A3hat)^c5)/((Yd.A2A4A1hat)^c5+(Yd.A2A4A3hat)^c5)
        r5.A3A1A2 <- ((Yd.A3A1A2hat)^c5)/((Yd.A3A1A2hat)^c5+(Yd.A3A1A4hat)^c5)
        r5.A3A1A4 <- ((Yd.A3A1A4hat)^c5)/((Yd.A3A1A2hat)^c5+(Yd.A3A1A4hat)^c5)
        r5.A3A2A1 <- ((Yd.A3A2A1hat)^c5)/((Yd.A3A2A1hat)^c5+(Yd.A3A2A4hat)^c5)
        r5.A3A2A4 <- ((Yd.A3A2A4hat)^c5)/((Yd.A3A2A1hat)^c5+(Yd.A3A2A4hat)^c5)
        r5.A3A4A1 <- ((Yd.A3A4A1hat)^c5)/((Yd.A3A4A1hat)^c5+(Yd.A3A4A2hat)^c5)
        r5.A3A4A2 <- ((Yd.A3A4A2hat)^c5)/((Yd.A3A4A1hat)^c5+(Yd.A3A4A2hat)^c5)
        r5.A4A1A2 <- ((Yd.A4A1A2hat)^c5)/((Yd.A4A1A2hat)^c5+(Yd.A4A1A3hat)^c5)
        r5.A4A1A3 <- ((Yd.A4A1A3hat)^c5)/((Yd.A4A1A2hat)^c5+(Yd.A4A1A3hat)^c5)
        r5.A4A2A1 <- ((Yd.A4A2A1hat)^c5)/((Yd.A4A2A1hat)^c5+(Yd.A4A2A3hat)^c5)
        r5.A4A2A3 <- ((Yd.A4A2A3hat)^c5)/((Yd.A4A2A1hat)^c5+(Yd.A4A2A3hat)^c5)
        r5.A4A3A1 <- ((Yd.A4A3A1hat)^c5)/((Yd.A4A3A1hat)^c5+(Yd.A4A3A2hat)^c5)
        r5.A4A3A2 <- ((Yd.A4A3A2hat)^c5)/((Yd.A4A3A1hat)^c5+(Yd.A4A3A2hat)^c5)
        
        # AR(c,e)
        e5 <- e[5]
        
        if (I1.A1[i]==1 & I2.A2[i]==1 & Y1.A1[i]==0 & Y2.A1A2[i]==0) {
          P3.A1A2A3[i] <- ifelse(is.na(r5.A1A2A3), 1/2,  min(max(c(r5.A1A2A3, e5)), 1-e5))
          P3.A1A2A4[i] <- 1-P3.A1A2A3[i]
          A_3[i] <- sample(c(3,4), size = 1, replace = TRUE, prob = c(P3.A1A2A3[i],P3.A1A2A4[i]))
        } else if (I1.A1[i]==1 & I2.A3[i]==1 & Y1.A1[i]==0 & Y2.A1A3[i]==0) {
          P3.A1A3A2[i] <- ifelse(is.na(r5.A1A3A2), 1/2,  min(max(c(r5.A1A3A2, e5)), 1-e5))
          P3.A1A3A4[i] <- 1-P3.A1A3A2[i]
          A_3[i] <- sample(c(2,4), size = 1, replace = TRUE, prob = c(P3.A1A3A2[i],P3.A1A3A4[i]))
        } else if (I1.A1[i]==1 & I2.A4[i]==1 & Y1.A1[i]==0 & Y2.A1A4[i]==0) {
          P3.A1A4A2[i] <- ifelse(is.na(r5.A1A4A2), 1/2,  min(max(c(r5.A1A4A2, e5)), 1-e5))
          P3.A1A4A3[i] <- 1-P3.A1A4A2[i]
          A_3[i] <- sample(c(2,3), size = 1, replace = TRUE, prob = c(P3.A1A4A2[i],P3.A1A4A3[i]))
        } else if (I1.A2[i]==1 & I2.A1[i]==1 & Y1.A2[i]==0 & Y2.A2A1[i]==0) {
          P3.A2A1A3[i] <- ifelse(is.na(r5.A2A1A3), 1/2,  min(max(c(r5.A2A1A3, e5)), 1-e5))
          P3.A2A1A4[i] <- 1-P3.A2A1A3[i]
          A_3[i] <- sample(c(3,4), size = 1, replace = TRUE, prob = c(P3.A2A1A3[i],P3.A2A1A4[i]))
        } else if (I1.A2[i]==1 & I2.A3[i]==1 & Y1.A2[i]==0 & Y2.A2A3[i]==0){
          P3.A2A3A1[i] <- ifelse(is.na(r5.A2A3A1), 1/2,  min(max(c(r5.A2A3A1, e5)), 1-e5))
          P3.A2A3A4[i] <- 1-P3.A2A3A1[i]
          A_3[i] <- sample(c(1,4), size = 1, replace = TRUE, prob = c(P3.A2A3A1[i],P3.A2A3A4[i]))
        } else if (I1.A2[i]==1 & I2.A4[i]==1 & Y1.A2[i]==0 & Y2.A2A4[i]==0) {
          P3.A2A4A1[i] <- ifelse(is.na(r5.A2A4A1), 1/2,  min(max(c(r5.A2A4A1, e5)), 1-e5))
          P3.A2A4A3[i] <- 1-P3.A2A4A1[i]
          A_3[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P3.A2A4A1[i],P3.A2A4A3[i]))
        } else if (I1.A3[i]==1 & I2.A1[i]==1 & Y1.A3[i]==0 & Y2.A3A1[i]==0) {
          P3.A3A1A2[i] <- ifelse(is.na(r5.A3A1A2), 1/2,  min(max(c(r5.A3A1A2, e5)), 1-e5))
          P3.A3A1A4[i] <- 1-P3.A3A1A2[i]
          A_3[i] <- sample(c(2,4), size = 1, replace = TRUE, prob = c(P3.A3A1A2[i],P3.A3A1A4[i]))
        } else if (I1.A3[i]==1 & I2.A2[i]==1 & Y1.A3[i]==0 & Y2.A3A2[i]==0) {
          P3.A3A2A1[i] <- ifelse(is.na(r5.A3A2A1), 1/2,  min(max(c(r5.A3A2A1, e5)), 1-e5))
          P3.A3A2A4[i] <- 1-P3.A3A2A1[i] 
          A_3[i] <- sample(c(1,4), size = 1, replace = TRUE, prob = c(P3.A3A2A1[i],P3.A3A2A4[i]))
        } else if (I1.A3[i]==1 & I2.A4[i]==1 & Y1.A3[i]==0 & Y2.A3A4[i]==0) {
          P3.A3A4A1[i] <- ifelse(is.na(r5.A3A4A1), 1/2,  min(max(c(r5.A3A4A1, e5)), 1-e5))
          P3.A3A4A2[i] <- 1-P3.A3A4A1[i]
          A_3[i] <- sample(c(1,2), size = 1, replace = TRUE, prob = c(P3.A3A4A1[i],P3.A3A4A2[i]))
        } else if (I1.A4[i]==1 & I2.A1[i]==1 & Y1.A4[i]==0 & Y2.A4A1[i]==0) {
          P3.A4A1A2[i] <- ifelse(is.na(r5.A4A1A2), 1/2,  min(max(c(r5.A4A1A2, e5)), 1-e5))
          P3.A4A1A3[i] <- 1-P3.A4A1A2[i]
          A_3[i] <- sample(c(2,3), size = 1, replace = TRUE, prob = c(P3.A4A1A2[i],P3.A4A1A3[i]))
        } else if (I1.A4[i]==1 & I2.A2[i]==1 & Y1.A4[i]==0 & Y2.A4A2[i]==0) {
          P3.A4A2A1[i] <- ifelse(is.na(r5.A4A2A1), 1/2,  min(max(c(r5.A4A2A1, e5)), 1-e5))
          P3.A4A2A3[i] <- 1-P3.A4A2A1[i]
          A_3[i] <- sample(c(1,3), size = 1, replace = TRUE, prob = c(P3.A4A2A1[i],P3.A4A2A3[i]))
        } else if (I1.A4[i]==1 & I2.A3[i]==1 & Y1.A4[i]==0 & Y2.A4A3[i]==0) {
          P3.A4A3A1[i] <- ifelse(is.na(r5.A4A3A1), 1/2,  min(max(c(r5.A4A3A1, e5)), 1-e5))
          P3.A4A3A2[i] <- 1-P3.A4A3A1[i]
          A_3[i] <- sample(c(1,2), size = 1, replace = TRUE, prob = c(P3.A4A3A1[i],P3.A4A3A2[i]))
        }
        
        I3.A1[i]<- ifelse(A_3[i]==1, 1, 0)
        I3.A2[i]<- ifelse(A_3[i]==2, 1, 0)
        I3.A3[i]<- ifelse(A_3[i]==3, 1, 0)
        I3.A4[i]<- ifelse(A_3[i]==4, 1, 0)
        
        
        
        Y3[i] <- ifelse(Y1[i]==0 & Y2[i]==0, 
                        sum(c(Y3.A1A2A3[i]*I3.A3[i]*I2.A2[i]*I1.A1[i],Y3.A1A2A4[i]*I3.A4[i]*I2.A2[i]*I1.A1[i],
                              Y3.A1A3A2[i]*I3.A2[i]*I2.A3[i]*I1.A1[i],Y3.A1A3A4[i]*I3.A4[i]*I2.A3[i]*I1.A1[i],
                              Y3.A1A4A2[i]*I3.A2[i]*I2.A4[i]*I1.A1[i],Y3.A1A4A3[i]*I3.A3[i]*I2.A4[i]*I1.A1[i],
                              Y3.A2A1A3[i]*I3.A3[i]*I2.A1[i]*I1.A2[i],Y3.A2A1A4[i]*I3.A4[i]*I2.A1[i]*I1.A2[i],
                              Y3.A2A3A1[i]*I3.A1[i]*I2.A3[i]*I1.A2[i],Y3.A2A3A4[i]*I3.A4[i]*I2.A3[i]*I1.A2[i],
                              Y3.A2A4A1[i]*I3.A1[i]*I2.A4[i]*I1.A2[i],Y3.A2A4A3[i]*I3.A3[i]*I2.A4[i]*I1.A2[i],
                              Y3.A3A1A2[i]*I3.A2[i]*I2.A1[i]*I1.A3[i],Y3.A3A1A4[i]*I3.A4[i]*I2.A1[i]*I1.A3[i],
                              Y3.A3A2A1[i]*I3.A1[i]*I2.A2[i]*I1.A3[i],Y3.A3A2A4[i]*I3.A4[i]*I2.A2[i]*I1.A3[i],
                              Y3.A3A4A1[i]*I3.A1[i]*I2.A4[i]*I1.A3[i],Y3.A3A4A2[i]*I3.A2[i]*I2.A4[i]*I1.A3[i],
                              Y3.A4A1A2[i]*I3.A2[i]*I2.A1[i]*I1.A4[i],Y3.A4A1A3[i]*I3.A3[i]*I2.A1[i]*I1.A4[i],
                              Y3.A4A2A1[i]*I3.A1[i]*I2.A2[i]*I1.A4[i],Y3.A4A2A3[i]*I3.A3[i]*I2.A2[i]*I1.A4[i],
                              Y3.A4A3A1[i]*I3.A1[i]*I2.A3[i]*I1.A4[i],Y3.A4A3A2[i]*I3.A2[i]*I2.A3[i]*I1.A4[i]), na.rm = T), NA)
        
        Y[i] <- ifelse(Y1[i]==1, Y1[i], ifelse(Y2[i]==1, Y2[i], Y3[i]))
      }
    }
  }
  
  df <- data.frame(Y1.A1, Y1.A2, Y1.A3, Y1.A4, 
                   P1.A1, P1.A2, P1.A3, P1.A4, 
                   Y2.A1A2, Y2.A1A3, Y2.A1A4, Y2.A2A1, Y2.A2A3, Y2.A2A4, Y2.A3A1, Y2.A3A2, Y2.A3A4, Y2.A4A1, Y2.A4A2, Y2.A4A3,
                   A_1, I1.A1, I1.A2, I1.A3, I1.A4, 
                   A_2, I2.A1, I2.A2, I2.A3, I2.A4,
                   A_3, I3.A1, I3.A2, I3.A3, I3.A4,
                   P2.A1A2, P2.A1A3, P2.A1A4, P2.A2A1, P2.A2A3, P2.A2A4, P2.A3A1, P2.A3A2, P2.A3A4, P2.A4A1, P2.A4A2, P2.A4A3,
                   Y3.A1A2A3, Y3.A1A2A4, Y3.A1A3A2, Y3.A1A3A4, Y3.A1A4A2, Y3.A1A4A3, 
                   Y3.A2A1A3, Y3.A2A1A4, Y3.A2A3A1, Y3.A2A3A4, Y3.A2A4A1, Y3.A2A4A3, 
                   Y3.A3A1A2, Y3.A3A1A4, Y3.A3A2A1, Y3.A3A2A4, Y3.A3A4A1, Y3.A3A4A2, 
                   Y3.A4A1A2, Y3.A4A1A3, Y3.A4A2A1, Y3.A4A2A3, Y3.A4A3A1, Y3.A4A3A2,
                   P3.A1A2A3, P3.A1A2A4, P3.A1A3A2, P3.A1A3A4, P3.A1A4A2, P3.A1A4A3, 
                   P3.A2A1A3, P3.A2A1A4, P3.A2A3A1, P3.A2A3A4, P3.A2A4A1, P3.A2A4A3, 
                   P3.A3A1A2, P3.A3A1A4, P3.A3A2A1, P3.A3A2A4, P3.A3A4A1, P3.A3A4A2, 
                   P3.A4A1A2, P3.A4A1A3, P3.A4A2A1, P3.A4A2A3, P3.A4A3A1, P3.A4A3A2,
                   Y1, Y2, Y3, Y)
  return(df)
} 

  
######## Step2:  Monte Carlo replication function #############

monte_three <- function(N,n,pi,p0.burn,p1.burn,p2.burn,AR,c,e,sce) {
  # true response rate
  pi1.A1 <- pi[1]
  pi1.A2 <- pi[2]
  pi1.A3 <- pi[3]
  pi1.A4 <- pi[4]
  pi2.A1A2 <- pi[5]
  pi2.A1A3 <- pi[6]
  pi2.A1A4 <- pi[7]
  pi2.A2A1 <- pi[8]
  pi2.A2A3 <- pi[9]
  pi2.A2A4 <- pi[10]
  pi2.A3A1 <- pi[11]
  pi2.A3A2 <- pi[12]
  pi2.A3A4 <- pi[13]
  pi2.A4A1 <- pi[14]
  pi2.A4A2 <- pi[15]
  pi2.A4A3 <- pi[16]
  pi3.A1A2A3 <- pi[17]
  pi3.A1A2A4 <- pi[18]
  pi3.A1A3A2 <- pi[19]
  pi3.A1A3A4 <- pi[20]
  pi3.A1A4A2 <- pi[21]
  pi3.A1A4A3 <- pi[22]
  pi3.A2A1A3 <- pi[23]
  pi3.A2A1A4 <- pi[24]
  pi3.A2A3A1 <- pi[25]
  pi3.A2A3A4 <- pi[26]
  pi3.A2A4A1 <- pi[27]
  pi3.A2A4A3 <- pi[28]
  pi3.A3A1A2 <- pi[29]
  pi3.A3A1A4 <- pi[30]
  pi3.A3A2A1 <- pi[31]
  pi3.A3A2A4 <- pi[32]
  pi3.A3A4A1 <- pi[33]
  pi3.A3A4A2 <- pi[34]
  pi3.A4A1A2 <- pi[35]
  pi3.A4A1A3 <- pi[36]
  pi3.A4A2A1 <- pi[37]
  pi3.A4A2A3 <- pi[38]
  pi3.A4A3A1 <- pi[39]
  pi3.A4A3A2 <- pi[40]
  
  # sample size
  n <- n
  
  # Monte Carlo sample size
  N <- N
  
  # stage 1 burn-in sample
  n0 <- n*p0.burn
  # stage 2 burn-in sample
  n1 <- n*p1.burn
  # stage 3 burn-in sample
  n2 <- n*p2.burn
  
  # treatment assignment
  I1.A1 <- matrix(NA, nrow=N, ncol=n)
  I1.A2 <- matrix(NA, nrow=N, ncol=n)
  I1.A3 <- matrix(NA, nrow=N, ncol=n)
  I1.A4 <- matrix(NA, nrow=N, ncol=n)
  I2.A1 <- matrix(NA, nrow=N, ncol=n)
  I2.A2 <- matrix(NA, nrow=N, ncol=n)
  I2.A3 <- matrix(NA, nrow=N, ncol=n)
  I2.A4 <- matrix(NA, nrow=N, ncol=n)
  I3.A1 <- matrix(NA, nrow=N, ncol=n)
  I3.A2 <- matrix(NA, nrow=N, ncol=n)
  I3.A3 <- matrix(NA, nrow=N, ncol=n)
  I3.A4 <- matrix(NA, nrow=N, ncol=n)
  
  # observed outcome
  Y1 <- matrix(NA, nrow=N, ncol=n)
  Y2 <- matrix(NA, nrow=N, ncol=n)
  Y3 <- matrix(NA, nrow=N, ncol=n)
  Y <- matrix(NA, nrow=N, ncol=n)
  
  # counterfactual outcome
  Y1.A1 <- matrix(NA, nrow=N, ncol=n)
  Y1.A2 <- matrix(NA, nrow=N, ncol=n)
  Y1.A3 <- matrix(NA, nrow=N, ncol=n)
  Y1.A4 <- matrix(NA, nrow=N, ncol=n)
  Y2.A1A2 <- matrix(NA, nrow=N, ncol=n)
  Y2.A1A3 <- matrix(NA, nrow=N, ncol=n)
  Y2.A1A4 <- matrix(NA, nrow=N, ncol=n)
  Y2.A2A1 <- matrix(NA, nrow=N, ncol=n)
  Y2.A2A3 <- matrix(NA, nrow=N, ncol=n)
  Y2.A2A4 <- matrix(NA, nrow=N, ncol=n)
  Y2.A3A1 <- matrix(NA, nrow=N, ncol=n)
  Y2.A3A2 <- matrix(NA, nrow=N, ncol=n)
  Y2.A3A4 <- matrix(NA, nrow=N, ncol=n)
  Y2.A4A1 <- matrix(NA, nrow=N, ncol=n)
  Y2.A4A2 <- matrix(NA, nrow=N, ncol=n)
  Y2.A4A3 <- matrix(NA, nrow=N, ncol=n)
  Y3.A1A2A3 <- matrix(NA, nrow=N, ncol=n)
  Y3.A1A2A4 <- matrix(NA, nrow=N, ncol=n)
  Y3.A1A3A2 <- matrix(NA, nrow=N, ncol=n)
  Y3.A1A3A4 <- matrix(NA, nrow=N, ncol=n)
  Y3.A1A4A2 <- matrix(NA, nrow=N, ncol=n)
  Y3.A1A4A3 <- matrix(NA, nrow=N, ncol=n)
  Y3.A2A1A3 <- matrix(NA, nrow=N, ncol=n)
  Y3.A2A1A4 <- matrix(NA, nrow=N, ncol=n)
  Y3.A2A3A1 <- matrix(NA, nrow=N, ncol=n)
  Y3.A2A3A4 <- matrix(NA, nrow=N, ncol=n)
  Y3.A2A4A1 <- matrix(NA, nrow=N, ncol=n)
  Y3.A2A4A3 <- matrix(NA, nrow=N, ncol=n)
  Y3.A3A1A2 <- matrix(NA, nrow=N, ncol=n)
  Y3.A3A1A4 <- matrix(NA, nrow=N, ncol=n)
  Y3.A3A2A1 <- matrix(NA, nrow=N, ncol=n)
  Y3.A3A2A4 <- matrix(NA, nrow=N, ncol=n)
  Y3.A3A4A1 <- matrix(NA, nrow=N, ncol=n)
  Y3.A3A4A2 <- matrix(NA, nrow=N, ncol=n)
  Y3.A4A1A2 <- matrix(NA, nrow=N, ncol=n)
  Y3.A4A1A3 <- matrix(NA, nrow=N, ncol=n)
  Y3.A4A2A1 <- matrix(NA, nrow=N, ncol=n)
  Y3.A4A2A3 <- matrix(NA, nrow=N, ncol=n)
  Y3.A4A3A1 <- matrix(NA, nrow=N, ncol=n)
  Y3.A4A3A2 <- matrix(NA, nrow=N, ncol=n)
  
  
  # G-estimator
  mu.A1A2A3hat.g <- rep(NA, N)
  mu.A1A2A4hat.g <- rep(NA, N)
  mu.A1A3A2hat.g <- rep(NA, N)
  mu.A1A3A4hat.g <- rep(NA, N)
  mu.A1A4A2hat.g <- rep(NA, N)
  mu.A1A4A3hat.g <- rep(NA, N)
  mu.A2A1A3hat.g <- rep(NA, N)
  mu.A2A1A4hat.g <- rep(NA, N)
  mu.A2A3A1hat.g <- rep(NA, N)
  mu.A2A3A4hat.g <- rep(NA, N)
  mu.A2A4A1hat.g <- rep(NA, N)
  mu.A2A4A3hat.g <- rep(NA, N)
  mu.A3A1A2hat.g <- rep(NA, N)
  mu.A3A1A4hat.g <- rep(NA, N)
  mu.A3A2A1hat.g <- rep(NA, N)
  mu.A3A2A4hat.g <- rep(NA, N)
  mu.A3A4A1hat.g <- rep(NA, N)
  mu.A3A4A2hat.g <- rep(NA, N)
  mu.A4A1A2hat.g <- rep(NA, N)
  mu.A4A1A3hat.g <- rep(NA, N)
  mu.A4A2A1hat.g <- rep(NA, N)
  mu.A4A2A3hat.g <- rep(NA, N)
  mu.A4A3A1hat.g <- rep(NA, N)
  mu.A4A3A2hat.g <- rep(NA, N)
  
  
  # estimated optimal/worst DTR
  est.opt.g <- rep(NA, N)
  est.wor.g <- rep(NA, N)
  
  
  # number of treated in each DTRs
  n.A1A2A3 <- rep(NA, N)
  n.A1A2A4 <- rep(NA, N)
  n.A1A3A2 <- rep(NA, N)
  n.A1A3A4 <- rep(NA, N)
  n.A1A4A2 <- rep(NA, N)
  n.A1A4A3 <- rep(NA, N)
  n.A2A1A3 <- rep(NA, N)
  n.A2A1A4 <- rep(NA, N)
  n.A2A3A1 <- rep(NA, N)
  n.A2A3A4 <- rep(NA, N)
  n.A2A4A1 <- rep(NA, N)
  n.A2A4A3 <- rep(NA, N)
  n.A3A1A2 <- rep(NA, N)
  n.A3A1A4 <- rep(NA, N)
  n.A3A2A1 <- rep(NA, N)
  n.A3A2A4 <- rep(NA, N)
  n.A3A4A1 <- rep(NA, N)
  n.A3A4A2 <- rep(NA, N)
  n.A4A1A2 <- rep(NA, N)
  n.A4A1A3 <- rep(NA, N)
  n.A4A2A1 <- rep(NA, N)
  n.A4A2A3 <- rep(NA, N)
  n.A4A3A1 <- rep(NA, N)
  n.A4A3A2 <- rep(NA, N)
  
  # number of responder in each DTRs
  rn.A1A2A3 <- rep(NA, N)
  rn.A1A2A4 <- rep(NA, N)
  rn.A1A3A2 <- rep(NA, N)
  rn.A1A3A4 <- rep(NA, N)
  rn.A1A4A2 <- rep(NA, N)
  rn.A1A4A3 <- rep(NA, N)
  rn.A2A1A3 <- rep(NA, N)
  rn.A2A1A4 <- rep(NA, N)
  rn.A2A3A1 <- rep(NA, N)
  rn.A2A3A4 <- rep(NA, N)
  rn.A2A4A1 <- rep(NA, N)
  rn.A2A4A3 <- rep(NA, N)
  rn.A3A1A2 <- rep(NA, N)
  rn.A3A1A4 <- rep(NA, N)
  rn.A3A2A1 <- rep(NA, N)
  rn.A3A2A4 <- rep(NA, N)
  rn.A3A4A1 <- rep(NA, N)
  rn.A3A4A2 <- rep(NA, N)
  rn.A4A1A2 <- rep(NA, N)
  rn.A4A1A3 <- rep(NA, N)
  rn.A4A2A1 <- rep(NA, N)
  rn.A4A2A3 <- rep(NA, N)
  rn.A4A3A1 <- rep(NA, N)
  rn.A4A3A2 <- rep(NA, N)
  
  # total number of response in the trial
  NR <- rep(NA, N)
  
  # create N Monte Carlo samples
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for (i in 1:N) {
    # SMART & GO-SMART  
    df <- data_gene_three_c(n=n,pi=pi,p0.burn=p0.burn,p1.burn=p1.burn,p2.burn=p2.burn,AR=AR,c=c,e=e)
    
    # treatment assignment
    I1.A1[i,] <- df$I1.A1
    I1.A2[i,] <- df$I1.A2
    I1.A3[i,] <- df$I1.A3
    I1.A4[i,] <- df$I1.A4
    I2.A1[i,] <- df$I2.A1
    I2.A2[i,] <- df$I2.A2
    I2.A3[i,] <- df$I2.A3
    I2.A4[i,] <- df$I2.A4
    I3.A1[i,] <- df$I3.A1
    I3.A2[i,] <- df$I3.A2
    I3.A3[i,] <- df$I3.A3
    I3.A4[i,] <- df$I3.A4
    
    # observed outcome
    Y1[i,] <- df$Y1
    Y2[i,] <- df$Y2
    Y3[i,] <- df$Y3
    Y[i,] <- df$Y
    
    # G-estimator
    pi1.A1.hat.g <- sum(Y1[i,]*I1.A1[i,])/sum(I1.A1[i,])
    pi1.A2.hat.g <- sum(Y1[i,]*I1.A2[i,])/sum(I1.A2[i,])
    pi1.A3.hat.g <- sum(Y1[i,]*I1.A3[i,])/sum(I1.A3[i,])
    pi1.A4.hat.g <- sum(Y1[i,]*I1.A4[i,])/sum(I1.A4[i,])
    
    pi2.A1A2.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,]*Y2[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,],na.rm=T)
    pi2.A1A3.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,]*Y2[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,],na.rm=T)
    pi2.A1A4.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A4[i,]*Y2[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A4[i,],na.rm=T)
    pi2.A2A1.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,]*Y2[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,],na.rm=T)
    pi2.A2A3.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,]*Y2[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,],na.rm=T)
    pi2.A2A4.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A4[i,]*Y2[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A4[i,],na.rm=T)
    pi2.A3A1.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,]*Y2[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,],na.rm=T)
    pi2.A3A2.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,]*Y2[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,],na.rm=T)
    pi2.A3A4.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A4[i,]*Y2[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A4[i,],na.rm=T)
    pi2.A4A1.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A1[i,]*Y2[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A1[i,],na.rm=T)
    pi2.A4A2.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A2[i,]*Y2[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A2[i,],na.rm=T)
    pi2.A4A3.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A3[i,]*Y2[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A3[i,],na.rm=T)
    
    pi3.A1A2A3.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A3[i,]*Y3[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A3[i,],na.rm=T)
    pi3.A1A2A4.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A4[i,]*Y3[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A4[i,],na.rm=T)
    pi3.A1A3A2.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A2[i,]*Y3[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A2[i,],na.rm=T)
    pi3.A1A3A4.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A4[i,]*Y3[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A4[i,],na.rm=T)
    pi3.A1A4A2.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A2[i,]*Y3[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A2[i,],na.rm=T)
    pi3.A1A4A3.hat.g <- sum(I1.A1[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A3[i,]*Y3[i,],na.rm=T)/sum(I1.A1[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A3[i,],na.rm=T)
    pi3.A2A1A3.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A3[i,]*Y3[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A3[i,],na.rm=T)
    pi3.A2A1A4.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A4[i,]*Y3[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A4[i,],na.rm=T)
    pi3.A2A3A1.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A1[i,]*Y3[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A1[i,],na.rm=T)
    pi3.A2A3A4.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A4[i,]*Y3[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A4[i,],na.rm=T)
    pi3.A2A4A1.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A1[i,]*Y3[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A1[i,],na.rm=T)
    pi3.A2A4A3.hat.g <- sum(I1.A2[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A3[i,]*Y3[i,],na.rm=T)/sum(I1.A2[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A3[i,],na.rm=T)
    pi3.A3A1A2.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A2[i,]*Y3[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A2[i,],na.rm=T)
    pi3.A3A1A4.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A4[i,]*Y3[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A4[i,],na.rm=T)
    pi3.A3A2A1.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A1[i,]*Y3[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A1[i,],na.rm=T)
    pi3.A3A2A4.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A4[i,]*Y3[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A4[i,],na.rm=T)
    pi3.A3A4A1.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A1[i,]*Y3[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A1[i,],na.rm=T)
    pi3.A3A4A2.hat.g <- sum(I1.A3[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A2[i,]*Y3[i,],na.rm=T)/sum(I1.A3[i,]*(1-Y1[i,])*I2.A4[i,]*(1-Y2[i,])*I3.A2[i,],na.rm=T)
    pi3.A4A1A2.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A2[i,]*Y3[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A2[i,],na.rm=T)
    pi3.A4A1A3.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A3[i,]*Y3[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A1[i,]*(1-Y2[i,])*I3.A3[i,],na.rm=T)
    pi3.A4A2A1.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A1[i,]*Y3[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A1[i,],na.rm=T)
    pi3.A4A2A3.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A3[i,]*Y3[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A2[i,]*(1-Y2[i,])*I3.A3[i,],na.rm=T)
    pi3.A4A3A1.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A1[i,]*Y3[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A1[i,],na.rm=T)
    pi3.A4A3A2.hat.g <- sum(I1.A4[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A2[i,]*Y3[i,],na.rm=T)/sum(I1.A4[i,]*(1-Y1[i,])*I2.A3[i,]*(1-Y2[i,])*I3.A2[i,],na.rm=T)
    
    mu.A1A2A3hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A2.hat.g+(1-pi1.A1.hat.g)*(1-pi2.A1A2.hat.g)*pi3.A1A2A3.hat.g
    mu.A1A2A4hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A2.hat.g+(1-pi1.A1.hat.g)*(1-pi2.A1A2.hat.g)*pi3.A1A2A4.hat.g
    mu.A1A3A2hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A3.hat.g+(1-pi1.A1.hat.g)*(1-pi2.A1A3.hat.g)*pi3.A1A3A2.hat.g
    mu.A1A3A4hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A3.hat.g+(1-pi1.A1.hat.g)*(1-pi2.A1A3.hat.g)*pi3.A1A3A4.hat.g
    mu.A1A4A2hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A4.hat.g+(1-pi1.A1.hat.g)*(1-pi2.A1A4.hat.g)*pi3.A1A4A2.hat.g
    mu.A1A4A3hat.g[i] <- pi1.A1.hat.g+(1-pi1.A1.hat.g)*pi2.A1A4.hat.g+(1-pi1.A1.hat.g)*(1-pi2.A1A4.hat.g)*pi3.A1A4A3.hat.g
    
    mu.A2A1A3hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A1.hat.g+(1-pi1.A2.hat.g)*(1-pi2.A2A1.hat.g)*pi3.A2A1A3.hat.g
    mu.A2A1A4hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A1.hat.g+(1-pi1.A2.hat.g)*(1-pi2.A2A1.hat.g)*pi3.A2A1A4.hat.g
    mu.A2A3A1hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A3.hat.g+(1-pi1.A2.hat.g)*(1-pi2.A2A3.hat.g)*pi3.A2A3A1.hat.g
    mu.A2A3A4hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A3.hat.g+(1-pi1.A2.hat.g)*(1-pi2.A2A3.hat.g)*pi3.A2A3A4.hat.g
    mu.A2A4A1hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A4.hat.g+(1-pi1.A2.hat.g)*(1-pi2.A2A4.hat.g)*pi3.A2A4A1.hat.g
    mu.A2A4A3hat.g[i] <- pi1.A2.hat.g+(1-pi1.A2.hat.g)*pi2.A2A4.hat.g+(1-pi1.A2.hat.g)*(1-pi2.A2A4.hat.g)*pi3.A2A4A3.hat.g
    
    mu.A3A1A2hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A1.hat.g+(1-pi1.A3.hat.g)*(1-pi2.A3A1.hat.g)*pi3.A3A1A2.hat.g
    mu.A3A1A4hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A1.hat.g+(1-pi1.A3.hat.g)*(1-pi2.A3A1.hat.g)*pi3.A3A1A4.hat.g
    mu.A3A2A1hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A2.hat.g+(1-pi1.A3.hat.g)*(1-pi2.A3A2.hat.g)*pi3.A3A2A1.hat.g
    mu.A3A2A4hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A2.hat.g+(1-pi1.A3.hat.g)*(1-pi2.A3A2.hat.g)*pi3.A3A2A4.hat.g
    mu.A3A4A1hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A4.hat.g+(1-pi1.A3.hat.g)*(1-pi2.A3A4.hat.g)*pi3.A3A4A1.hat.g
    mu.A3A4A2hat.g[i] <- pi1.A3.hat.g+(1-pi1.A3.hat.g)*pi2.A3A4.hat.g+(1-pi1.A3.hat.g)*(1-pi2.A3A4.hat.g)*pi3.A3A4A2.hat.g
    
    mu.A4A1A2hat.g[i] <- pi1.A4.hat.g+(1-pi1.A4.hat.g)*pi2.A4A1.hat.g+(1-pi1.A4.hat.g)*(1-pi2.A4A1.hat.g)*pi3.A4A1A2.hat.g
    mu.A4A1A3hat.g[i] <- pi1.A4.hat.g+(1-pi1.A4.hat.g)*pi2.A4A1.hat.g+(1-pi1.A4.hat.g)*(1-pi2.A4A1.hat.g)*pi3.A4A1A3.hat.g
    mu.A4A2A1hat.g[i] <- pi1.A4.hat.g+(1-pi1.A4.hat.g)*pi2.A4A2.hat.g+(1-pi1.A4.hat.g)*(1-pi2.A4A2.hat.g)*pi3.A4A2A1.hat.g
    mu.A4A2A3hat.g[i] <- pi1.A4.hat.g+(1-pi1.A4.hat.g)*pi2.A4A2.hat.g+(1-pi1.A4.hat.g)*(1-pi2.A4A2.hat.g)*pi3.A4A2A3.hat.g
    mu.A4A3A1hat.g[i] <- pi1.A4.hat.g+(1-pi1.A4.hat.g)*pi2.A4A3.hat.g+(1-pi1.A4.hat.g)*(1-pi2.A4A3.hat.g)*pi3.A4A3A1.hat.g
    mu.A4A3A2hat.g[i] <- pi1.A4.hat.g+(1-pi1.A4.hat.g)*pi2.A4A3.hat.g+(1-pi1.A4.hat.g)*(1-pi2.A4A3.hat.g)*pi3.A4A3A2.hat.g
    
    # estimated optimal DTR based on G-estimator
    est.opt.g[i] <- which.max(c(mu.A1A2A3hat.g[i], mu.A1A2A4hat.g[i], 
                                mu.A1A3A2hat.g[i], mu.A1A3A4hat.g[i],
                                mu.A1A4A2hat.g[i], mu.A1A4A3hat.g[i],
                                mu.A2A1A3hat.g[i], mu.A2A1A4hat.g[i],
                                mu.A2A3A1hat.g[i], mu.A2A3A4hat.g[i],
                                mu.A2A4A1hat.g[i], mu.A2A4A3hat.g[i],
                                mu.A3A1A2hat.g[i], mu.A3A1A4hat.g[i],
                                mu.A3A2A1hat.g[i], mu.A3A2A4hat.g[i],
                                mu.A3A4A1hat.g[i], mu.A3A4A2hat.g[i],
                                mu.A4A1A2hat.g[i], mu.A4A1A3hat.g[i],
                                mu.A4A2A1hat.g[i], mu.A4A2A3hat.g[i],
                                mu.A4A3A1hat.g[i], mu.A4A3A2hat.g[i]))
    
    est.wor.g[i] <- which.min(c(mu.A1A2A3hat.g[i], mu.A1A2A4hat.g[i], 
                                mu.A1A3A2hat.g[i], mu.A1A3A4hat.g[i],
                                mu.A1A4A2hat.g[i], mu.A1A4A3hat.g[i],
                                mu.A2A1A3hat.g[i], mu.A2A1A4hat.g[i],
                                mu.A2A3A1hat.g[i], mu.A2A3A4hat.g[i],
                                mu.A2A4A1hat.g[i], mu.A2A4A3hat.g[i],
                                mu.A3A1A2hat.g[i], mu.A3A1A4hat.g[i],
                                mu.A3A2A1hat.g[i], mu.A3A2A4hat.g[i],
                                mu.A3A4A1hat.g[i], mu.A3A4A2hat.g[i],
                                mu.A4A1A2hat.g[i], mu.A4A1A3hat.g[i],
                                mu.A4A2A1hat.g[i], mu.A4A2A3hat.g[i],
                                mu.A4A3A1hat.g[i], mu.A4A3A2hat.g[i]))
    
    # Number of patient treated in each of the 24 DTRs 
    # (response in stage I + no response in stage I and response in stage II + no response in stage I & II and treated in stage III)
    n.A1A2A3 <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(Y2[i,]*I1.A1[i,]*I2.A2[i,], na.rm = T)+sum(I1.A1[i,]*I2.A2[i,]*I3.A3[i,], na.rm = T)
    n.A1A2A4 <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(Y2[i,]*I1.A1[i,]*I2.A2[i,], na.rm = T)+sum(I1.A1[i,]*I2.A2[i,]*I3.A4[i,], na.rm = T)
    n.A1A3A2 <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(Y2[i,]*I1.A1[i,]*I2.A3[i,], na.rm = T)+sum(I1.A1[i,]*I2.A3[i,]*I3.A2[i,], na.rm = T)
    n.A1A3A4 <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(Y2[i,]*I1.A1[i,]*I2.A3[i,], na.rm = T)+sum(I1.A1[i,]*I2.A3[i,]*I3.A4[i,], na.rm = T)
    n.A1A4A2 <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(Y2[i,]*I1.A1[i,]*I2.A4[i,], na.rm = T)+sum(I1.A1[i,]*I2.A4[i,]*I3.A2[i,], na.rm = T)
    n.A1A4A3 <- sum(Y1[i,]*I1.A1[i,], na.rm = T)+sum(Y2[i,]*I1.A1[i,]*I2.A4[i,], na.rm = T)+sum(I1.A1[i,]*I2.A4[i,]*I3.A3[i,], na.rm = T)
    n.A2A1A3 <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(Y2[i,]*I1.A2[i,]*I2.A1[i,], na.rm = T)+sum(I1.A2[i,]*I2.A1[i,]*I3.A3[i,], na.rm = T)
    n.A2A1A4 <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(Y2[i,]*I1.A2[i,]*I2.A1[i,], na.rm = T)+sum(I1.A2[i,]*I2.A1[i,]*I3.A4[i,], na.rm = T)
    n.A2A3A1 <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(Y2[i,]*I1.A2[i,]*I2.A3[i,], na.rm = T)+sum(I1.A2[i,]*I2.A3[i,]*I3.A1[i,], na.rm = T)
    n.A2A3A4 <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(Y2[i,]*I1.A2[i,]*I2.A3[i,], na.rm = T)+sum(I1.A2[i,]*I2.A3[i,]*I3.A4[i,], na.rm = T)
    n.A2A4A1 <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(Y2[i,]*I1.A2[i,]*I2.A4[i,], na.rm = T)+sum(I1.A2[i,]*I2.A4[i,]*I3.A1[i,], na.rm = T)
    n.A2A4A3 <- sum(Y1[i,]*I1.A2[i,], na.rm = T)+sum(Y2[i,]*I1.A2[i,]*I2.A4[i,], na.rm = T)+sum(I1.A2[i,]*I2.A4[i,]*I3.A3[i,], na.rm = T)
    n.A3A1A2 <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(Y2[i,]*I1.A3[i,]*I2.A1[i,], na.rm = T)+sum(I1.A3[i,]*I2.A1[i,]*I3.A2[i,], na.rm = T)
    n.A3A1A4 <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(Y2[i,]*I1.A3[i,]*I2.A1[i,], na.rm = T)+sum(I1.A3[i,]*I2.A1[i,]*I3.A4[i,], na.rm = T)
    n.A3A2A1 <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(Y2[i,]*I1.A3[i,]*I2.A2[i,], na.rm = T)+sum(I1.A3[i,]*I2.A2[i,]*I3.A1[i,], na.rm = T)
    n.A3A2A4 <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(Y2[i,]*I1.A3[i,]*I2.A2[i,], na.rm = T)+sum(I1.A3[i,]*I2.A2[i,]*I3.A4[i,], na.rm = T)
    n.A3A4A1 <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(Y2[i,]*I1.A3[i,]*I2.A4[i,], na.rm = T)+sum(I1.A3[i,]*I2.A4[i,]*I3.A1[i,], na.rm = T)
    n.A3A4A2 <- sum(Y1[i,]*I1.A3[i,], na.rm = T)+sum(Y2[i,]*I1.A3[i,]*I2.A4[i,], na.rm = T)+sum(I1.A3[i,]*I2.A4[i,]*I3.A2[i,], na.rm = T)
    n.A4A1A2 <- sum(Y1[i,]*I1.A4[i,], na.rm = T)+sum(Y2[i,]*I1.A4[i,]*I2.A1[i,], na.rm = T)+sum(I1.A4[i,]*I2.A1[i,]*I3.A2[i,], na.rm = T)
    n.A4A1A3 <- sum(Y1[i,]*I1.A4[i,], na.rm = T)+sum(Y2[i,]*I1.A4[i,]*I2.A1[i,], na.rm = T)+sum(I1.A4[i,]*I2.A1[i,]*I3.A3[i,], na.rm = T)
    n.A4A2A1 <- sum(Y1[i,]*I1.A4[i,], na.rm = T)+sum(Y2[i,]*I1.A4[i,]*I2.A2[i,], na.rm = T)+sum(I1.A4[i,]*I2.A2[i,]*I3.A1[i,], na.rm = T)
    n.A4A2A3 <- sum(Y1[i,]*I1.A4[i,], na.rm = T)+sum(Y2[i,]*I1.A4[i,]*I2.A2[i,], na.rm = T)+sum(I1.A4[i,]*I2.A2[i,]*I3.A3[i,], na.rm = T)
    n.A4A3A1 <- sum(Y1[i,]*I1.A4[i,], na.rm = T)+sum(Y2[i,]*I1.A4[i,]*I2.A3[i,], na.rm = T)+sum(I1.A4[i,]*I2.A3[i,]*I3.A1[i,], na.rm = T)
    n.A4A3A2 <- sum(Y1[i,]*I1.A4[i,], na.rm = T)+sum(Y2[i,]*I1.A4[i,]*I2.A3[i,], na.rm = T)+sum(I1.A4[i,]*I2.A3[i,]*I3.A2[i,], na.rm = T)
    
    
    # total number of response in the trial
    NR[i] <- sum(Y[i,])
    
    
    # progress indicator
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  out <- data.frame(mu.A1A2A3hat.g=mu.A1A2A3hat.g, mu.A1A2A4hat.g=mu.A1A2A4hat.g, 
                    mu.A1A3A2hat.g=mu.A1A3A2hat.g, mu.A1A3A4hat.g=mu.A1A3A4hat.g,
                    mu.A1A4A2hat.g=mu.A1A4A2hat.g, mu.A1A4A3hat.g=mu.A1A4A3hat.g,
                    mu.A2A1A3hat.g=mu.A2A1A3hat.g, mu.A2A1A4hat.g=mu.A2A1A4hat.g,
                    mu.A2A3A1hat.g=mu.A2A3A1hat.g, mu.A2A3A4hat.g=mu.A2A3A4hat.g,
                    mu.A2A4A1hat.g=mu.A2A4A1hat.g, mu.A2A4A3hat.g=mu.A2A4A3hat.g,
                    mu.A3A1A2hat.g=mu.A3A1A2hat.g, mu.A3A1A4hat.g=mu.A3A1A4hat.g,
                    mu.A3A2A1hat.g=mu.A3A2A1hat.g, mu.A3A2A4hat.g=mu.A3A2A4hat.g,
                    mu.A3A4A1hat.g=mu.A3A4A1hat.g, mu.A3A4A2hat.g=mu.A3A4A2hat.g,
                    mu.A4A1A2hat.g=mu.A4A1A2hat.g, mu.A4A1A3hat.g=mu.A4A1A3hat.g,
                    mu.A4A2A1hat.g=mu.A4A2A1hat.g, mu.A4A2A3hat.g=mu.A4A2A3hat.g,
                    mu.A4A3A1hat.g=mu.A4A3A1hat.g, mu.A4A3A2hat.g=mu.A4A3A2hat.g,
                    n.A1A2A3=n.A1A2A3, n.A1A2A4=n.A1A2A4, 
                    n.A1A3A2=n.A1A3A2, n.A1A3A4=n.A1A3A4,
                    n.A1A4A2=n.A1A4A2, n.A1A4A3=n.A1A4A3,
                    n.A2A1A3=n.A2A1A3, n.A2A1A4=n.A2A1A4,
                    n.A2A3A1=n.A2A3A1, n.A2A3A4=n.A2A3A4,
                    n.A2A4A1=n.A2A4A1, n.A2A4A3=n.A2A4A3,
                    n.A3A1A2=n.A3A1A2, n.A3A1A4=n.A3A1A4,
                    n.A3A2A1=n.A3A2A1, n.A3A2A4=n.A3A2A4,
                    n.A3A4A1=n.A3A4A1, n.A3A4A2=n.A3A4A2,
                    n.A4A1A2=n.A4A1A2, n.A4A1A3=n.A4A1A3,
                    n.A4A2A1=n.A4A2A1, n.A4A2A3=n.A4A2A3,
                    n.A4A3A1=n.A4A3A1, n.A4A3A2=n.A4A3A2,
                    NR=NR,
                    est.opt.g,est.wor.g)
  
  write.csv(out, paste0("./output/dat_three_", n, "_",p0.burn, "_",p1.burn, "_",p2.burn, "_", AR, "_", sce, "_c",c,".csv"))
  

}

######## Step3: Evaluating the operating characteristics ##############

eva_three <- function(df,pi,H0) {
  # true response rate
  pi1.A1 <- pi[1]
  pi1.A2 <- pi[2]
  pi1.A3 <- pi[3]
  pi1.A4 <- pi[4]
  pi2.A1A2 <- pi[5]
  pi2.A1A3 <- pi[6]
  pi2.A1A4 <- pi[7]
  pi2.A2A1 <- pi[8]
  pi2.A2A3 <- pi[9]
  pi2.A2A4 <- pi[10]
  pi2.A3A1 <- pi[11]
  pi2.A3A2 <- pi[12]
  pi2.A3A4 <- pi[13]
  pi2.A4A1 <- pi[14]
  pi2.A4A2 <- pi[15]
  pi2.A4A3 <- pi[16]
  pi3.A1A2A3 <- pi[17]
  pi3.A1A2A4 <- pi[18]
  pi3.A1A3A2 <- pi[19]
  pi3.A1A3A4 <- pi[20]
  pi3.A1A4A2 <- pi[21]
  pi3.A1A4A3 <- pi[22]
  pi3.A2A1A3 <- pi[23]
  pi3.A2A1A4 <- pi[24]
  pi3.A2A3A1 <- pi[25]
  pi3.A2A3A4 <- pi[26]
  pi3.A2A4A1 <- pi[27]
  pi3.A2A4A3 <- pi[28]
  pi3.A3A1A2 <- pi[29]
  pi3.A3A1A4 <- pi[30]
  pi3.A3A2A1 <- pi[31]
  pi3.A3A2A4 <- pi[32]
  pi3.A3A4A1 <- pi[33]
  pi3.A3A4A2 <- pi[34]
  pi3.A4A1A2 <- pi[35]
  pi3.A4A1A3 <- pi[36]
  pi3.A4A2A1 <- pi[37]
  pi3.A4A2A3 <- pi[38]
  pi3.A4A3A1 <- pi[39]
  pi3.A4A3A2 <- pi[40]
  
  # True response rate for 24 DTRs
  pid.A1A2A3 <- pi1.A1+(1-pi1.A1)*pi2.A1A2+(1-pi1.A1)*(1-pi2.A1A2)*pi3.A1A2A3
  pid.A1A2A4 <- pi1.A1+(1-pi1.A1)*pi2.A1A2+(1-pi1.A1)*(1-pi2.A1A2)*pi3.A1A2A4
  pid.A1A3A2 <- pi1.A1+(1-pi1.A1)*pi2.A1A3+(1-pi1.A1)*(1-pi2.A1A3)*pi3.A1A3A2
  pid.A1A3A4 <- pi1.A1+(1-pi1.A1)*pi2.A1A3+(1-pi1.A1)*(1-pi2.A1A3)*pi3.A1A3A4
  pid.A1A4A2 <- pi1.A1+(1-pi1.A1)*pi2.A1A4+(1-pi1.A1)*(1-pi2.A1A4)*pi3.A1A4A2
  pid.A1A4A3 <- pi1.A1+(1-pi1.A1)*pi2.A1A4+(1-pi1.A1)*(1-pi2.A1A4)*pi3.A1A4A3
  pid.A2A1A3 <- pi1.A2+(1-pi1.A2)*pi2.A2A1+(1-pi1.A2)*(1-pi2.A2A1)*pi3.A2A1A3
  pid.A2A1A4 <- pi1.A2+(1-pi1.A2)*pi2.A2A1+(1-pi1.A2)*(1-pi2.A2A1)*pi3.A2A1A4
  pid.A2A3A1 <- pi1.A2+(1-pi1.A2)*pi2.A2A3+(1-pi1.A2)*(1-pi2.A2A3)*pi3.A2A3A1
  pid.A2A3A4 <- pi1.A2+(1-pi1.A2)*pi2.A2A3+(1-pi1.A2)*(1-pi2.A2A3)*pi3.A2A3A4
  pid.A2A4A1 <- pi1.A2+(1-pi1.A2)*pi2.A2A4+(1-pi1.A2)*(1-pi2.A2A4)*pi3.A2A4A1
  pid.A2A4A3 <- pi1.A2+(1-pi1.A2)*pi2.A2A4+(1-pi1.A2)*(1-pi2.A2A4)*pi3.A2A4A3
  pid.A3A1A2 <- pi1.A3+(1-pi1.A3)*pi2.A3A1+(1-pi1.A3)*(1-pi2.A3A1)*pi3.A3A1A2
  pid.A3A1A4 <- pi1.A3+(1-pi1.A3)*pi2.A3A1+(1-pi1.A3)*(1-pi2.A3A1)*pi3.A3A1A4
  pid.A3A2A1 <- pi1.A3+(1-pi1.A3)*pi2.A3A2+(1-pi1.A3)*(1-pi2.A3A2)*pi3.A3A2A1
  pid.A3A2A4 <- pi1.A3+(1-pi1.A3)*pi2.A3A2+(1-pi1.A3)*(1-pi2.A3A2)*pi3.A3A2A4
  pid.A3A4A1 <- pi1.A3+(1-pi1.A3)*pi2.A3A4+(1-pi1.A3)*(1-pi2.A3A4)*pi3.A3A4A1
  pid.A3A4A2 <- pi1.A3+(1-pi1.A3)*pi2.A3A4+(1-pi1.A3)*(1-pi2.A3A4)*pi3.A3A4A2
  pid.A4A1A2 <- pi1.A4+(1-pi1.A4)*pi2.A4A1+(1-pi1.A4)*(1-pi2.A4A1)*pi3.A4A1A2
  pid.A4A1A3 <- pi1.A4+(1-pi1.A4)*pi2.A4A1+(1-pi1.A4)*(1-pi2.A4A1)*pi3.A4A1A3
  pid.A4A2A1 <- pi1.A4+(1-pi1.A4)*pi2.A4A2+(1-pi1.A4)*(1-pi2.A4A2)*pi3.A4A2A1
  pid.A4A2A3 <- pi1.A4+(1-pi1.A4)*pi2.A4A2+(1-pi1.A4)*(1-pi2.A4A2)*pi3.A4A2A3
  pid.A4A3A1 <- pi1.A4+(1-pi1.A4)*pi2.A4A3+(1-pi1.A4)*(1-pi2.A4A3)*pi3.A4A3A1
  pid.A4A3A2 <- pi1.A4+(1-pi1.A4)*pi2.A4A3+(1-pi1.A4)*(1-pi2.A4A3)*pi3.A4A3A2
  
  # G-estimator
  mu.A1A2A3hat.g <- df$mu.A1A2A3hat.g[!is.na(df$mu.A1A2A3hat.g)]
  mu.A1A2A4hat.g <- df$mu.A1A2A4hat.g[!is.na(df$mu.A1A2A4hat.g)] 
  mu.A1A3A2hat.g <- df$mu.A1A3A2hat.g[!is.na(df$mu.A1A3A2hat.g)] 
  mu.A1A3A4hat.g <- df$mu.A1A3A4hat.g[!is.na(df$mu.A1A3A4hat.g)] 
  mu.A1A4A2hat.g <- df$mu.A1A4A2hat.g[!is.na(df$mu.A1A4A2hat.g)] 
  mu.A1A4A3hat.g <- df$mu.A1A4A3hat.g[!is.na(df$mu.A1A4A3hat.g)] 
  mu.A2A1A3hat.g <- df$mu.A2A1A3hat.g[!is.na(df$mu.A2A1A3hat.g)]
  mu.A2A1A4hat.g <- df$mu.A2A1A4hat.g[!is.na(df$mu.A2A1A4hat.g)] 
  mu.A2A3A1hat.g <- df$mu.A2A3A1hat.g[!is.na(df$mu.A2A3A1hat.g)] 
  mu.A2A3A4hat.g <- df$mu.A2A3A4hat.g[!is.na(df$mu.A2A3A4hat.g)] 
  mu.A2A4A1hat.g <- df$mu.A2A4A1hat.g[!is.na(df$mu.A2A4A1hat.g)] 
  mu.A2A4A3hat.g <- df$mu.A2A4A3hat.g[!is.na(df$mu.A2A4A3hat.g)] 
  mu.A3A1A2hat.g <- df$mu.A3A1A2hat.g[!is.na(df$mu.A3A1A2hat.g)] 
  mu.A3A1A4hat.g <- df$mu.A3A1A4hat.g[!is.na(df$mu.A3A1A4hat.g)] 
  mu.A3A2A1hat.g <- df$mu.A3A2A1hat.g[!is.na(df$mu.A3A2A1hat.g)] 
  mu.A3A2A4hat.g <- df$mu.A3A2A4hat.g[!is.na(df$mu.A3A2A4hat.g)] 
  mu.A3A4A1hat.g <- df$mu.A3A4A1hat.g[!is.na(df$mu.A3A4A1hat.g)] 
  mu.A3A4A2hat.g <- df$mu.A3A4A2hat.g[!is.na(df$mu.A3A4A2hat.g)] 
  mu.A4A1A2hat.g <- df$mu.A4A1A2hat.g[!is.na(df$mu.A4A1A2hat.g)] 
  mu.A4A1A3hat.g <- df$mu.A4A1A3hat.g[!is.na(df$mu.A4A1A3hat.g)] 
  mu.A4A2A1hat.g <- df$mu.A4A2A1hat.g[!is.na(df$mu.A4A2A1hat.g)] 
  mu.A4A2A3hat.g <- df$mu.A4A2A3hat.g[!is.na(df$mu.A4A2A3hat.g)] 
  mu.A4A3A1hat.g <- df$mu.A4A3A1hat.g[!is.na(df$mu.A4A3A1hat.g)] 
  mu.A4A3A2hat.g <- df$mu.A4A3A2hat.g[!is.na(df$mu.A4A3A2hat.g)] 
  
  
  # estimated optimal/worst DTR
  est.opt.g <- df$est.opt.g
  est.wor.g <- df$est.wor.g
  
  
  # Number of patient treated in each of the 24 DTRs
  n.A1A2A3 <- df$n.A1A2A3
  n.A1A2A4 <- df$n.A1A2A4
  n.A1A3A2 <- df$n.A1A3A2
  n.A1A3A4 <- df$n.A1A3A4
  n.A1A4A2 <- df$n.A1A4A2
  n.A1A4A3 <- df$n.A1A4A3
  n.A2A1A3 <- df$n.A2A1A3
  n.A2A1A4 <- df$n.A2A1A4
  n.A2A3A1 <- df$n.A2A3A1
  n.A2A3A4 <- df$n.A2A3A4
  n.A2A4A1 <- df$n.A2A4A1
  n.A2A4A3 <- df$n.A2A4A3
  n.A3A1A2 <- df$n.A3A1A2
  n.A3A1A4 <- df$n.A3A1A4
  n.A3A2A1 <- df$n.A3A2A1
  n.A3A2A4 <- df$n.A3A2A4
  n.A3A4A1 <- df$n.A3A4A1
  n.A3A4A2 <- df$n.A3A4A2
  n.A4A1A2 <- df$n.A4A1A2
  n.A4A1A3 <- df$n.A4A1A3
  n.A4A2A1 <- df$n.A4A2A1
  n.A4A2A3 <- df$n.A4A2A3
  n.A4A3A1 <- df$n.A4A3A1
  n.A4A3A2 <- df$n.A4A3A2
  
  
  # number of response in the trial
  nr <- df$NR
  
  
  
  
  result <- list()
  
  
  
  ################## 1. bias #########################
  
  bias.mu.A1A2A3hat.g <- mean(mu.A1A2A3hat.g)-pid.A1A2A3
  bias.mu.A1A2A4hat.g <- mean(mu.A1A2A4hat.g)-pid.A1A2A4
  bias.mu.A1A3A2hat.g <- mean(mu.A1A3A2hat.g)-pid.A1A3A2
  bias.mu.A1A3A4hat.g <- mean(mu.A1A3A4hat.g)-pid.A1A3A4
  bias.mu.A1A4A2hat.g <- mean(mu.A1A4A2hat.g)-pid.A1A4A2
  bias.mu.A1A4A3hat.g <- mean(mu.A1A4A3hat.g)-pid.A1A4A3
  
  bias.mu.A2A1A3hat.g <- mean(mu.A2A1A3hat.g)-pid.A2A1A3
  bias.mu.A2A1A4hat.g <- mean(mu.A2A1A4hat.g)-pid.A2A1A4
  bias.mu.A2A3A1hat.g <- mean(mu.A2A3A1hat.g)-pid.A2A3A1
  bias.mu.A2A3A4hat.g <- mean(mu.A2A3A4hat.g)-pid.A2A3A4
  bias.mu.A2A4A1hat.g <- mean(mu.A2A4A1hat.g)-pid.A2A4A1
  bias.mu.A2A4A3hat.g <- mean(mu.A2A4A3hat.g)-pid.A2A4A3
  
  bias.mu.A3A1A2hat.g <- mean(mu.A3A1A2hat.g)-pid.A3A1A2
  bias.mu.A3A1A4hat.g <- mean(mu.A3A1A4hat.g)-pid.A3A1A4
  bias.mu.A3A2A1hat.g <- mean(mu.A3A2A1hat.g)-pid.A3A2A1
  bias.mu.A3A2A4hat.g <- mean(mu.A3A2A4hat.g)-pid.A3A2A4
  bias.mu.A3A4A1hat.g <- mean(mu.A3A4A1hat.g)-pid.A3A4A1
  bias.mu.A3A4A2hat.g <- mean(mu.A3A4A2hat.g)-pid.A3A4A2
  
  bias.mu.A4A1A2hat.g <- mean(mu.A4A1A2hat.g)-pid.A4A1A2
  bias.mu.A4A1A3hat.g <- mean(mu.A4A1A3hat.g)-pid.A4A1A3
  bias.mu.A4A2A1hat.g <- mean(mu.A4A2A1hat.g)-pid.A4A2A1
  bias.mu.A4A2A3hat.g <- mean(mu.A4A2A3hat.g)-pid.A4A2A3
  bias.mu.A4A3A1hat.g <- mean(mu.A4A3A1hat.g)-pid.A4A3A1
  bias.mu.A4A3A2hat.g <- mean(mu.A4A3A2hat.g)-pid.A4A3A2
  
  
  result$bias.mu.A1A2A3hat.g <- round(bias.mu.A1A2A3hat.g,3)
  result$bias.mu.A1A2A4hat.g <- round(bias.mu.A1A2A4hat.g,3)
  result$bias.mu.A1A3A2hat.g <- round(bias.mu.A1A3A2hat.g,3)
  result$bias.mu.A1A3A4hat.g <- round(bias.mu.A1A3A4hat.g,3)
  result$bias.mu.A1A4A2hat.g <- round(bias.mu.A1A4A2hat.g,3)
  result$bias.mu.A1A4A3hat.g <- round(bias.mu.A1A4A3hat.g,3)
  
  result$bias.mu.A2A1A3hat.g <- round(bias.mu.A2A1A3hat.g,3)
  result$bias.mu.A2A1A4hat.g <- round(bias.mu.A2A1A4hat.g,3)
  result$bias.mu.A2A3A1hat.g <- round(bias.mu.A2A3A1hat.g,3)
  result$bias.mu.A2A3A4hat.g <- round(bias.mu.A2A3A4hat.g,3)
  result$bias.mu.A2A4A1hat.g <- round(bias.mu.A2A4A1hat.g,3)
  result$bias.mu.A2A4A3hat.g <- round(bias.mu.A2A4A3hat.g,3)
  
  result$bias.mu.A3A1A2hat.g <- round(bias.mu.A3A1A2hat.g,3)
  result$bias.mu.A3A1A4hat.g <- round(bias.mu.A3A1A4hat.g,3)
  result$bias.mu.A3A2A1hat.g <- round(bias.mu.A3A2A1hat.g,3)
  result$bias.mu.A3A2A4hat.g <- round(bias.mu.A3A2A4hat.g,3)
  result$bias.mu.A3A4A1hat.g <- round(bias.mu.A3A4A1hat.g,3)
  result$bias.mu.A3A4A2hat.g <- round(bias.mu.A3A4A2hat.g,3)
  
  result$bias.mu.A4A1A2hat.g <- round(bias.mu.A4A1A2hat.g,3)
  result$bias.mu.A4A1A3hat.g <- round(bias.mu.A4A1A3hat.g,3)
  result$bias.mu.A4A2A1hat.g <- round(bias.mu.A4A2A1hat.g,3)
  result$bias.mu.A4A2A3hat.g <- round(bias.mu.A4A2A3hat.g,3)
  result$bias.mu.A4A3A1hat.g <- round(bias.mu.A4A3A1hat.g,3)
  result$bias.mu.A4A3A2hat.g <- round(bias.mu.A4A3A2hat.g,3)
  
  
  
  ################## 2. # of patients treated in each DTR ######################
  
  N.A1A2A3 <- mean(n.A1A2A3)
  N.A1A2A4 <- mean(n.A1A2A4)
  N.A1A3A2 <- mean(n.A1A3A2)
  N.A1A3A4 <- mean(n.A1A3A4)
  N.A1A4A2 <- mean(n.A1A4A2)
  N.A1A4A3 <- mean(n.A1A4A3)
  N.A2A1A3 <- mean(n.A2A1A3)
  N.A2A1A4 <- mean(n.A2A1A4)
  N.A2A3A1 <- mean(n.A2A3A1)
  N.A2A3A4 <- mean(n.A2A3A4)
  N.A2A4A1 <- mean(n.A2A4A1)
  N.A2A4A3 <- mean(n.A2A4A3)
  N.A3A1A2 <- mean(n.A3A1A2)
  N.A3A1A4 <- mean(n.A3A1A4)
  N.A3A2A1 <- mean(n.A3A2A1)
  N.A3A2A4 <- mean(n.A3A2A4)
  N.A3A4A1 <- mean(n.A3A4A1)
  N.A3A4A2 <- mean(n.A3A4A2)
  N.A4A1A2 <- mean(n.A4A1A2)
  N.A4A1A3 <- mean(n.A4A1A3)
  N.A4A2A1 <- mean(n.A4A2A1)
  N.A4A2A3 <- mean(n.A4A2A3)
  N.A4A3A1 <- mean(n.A4A3A1)
  N.A4A3A2 <- mean(n.A4A3A2)
  
  
  N.A1A2A3.CI <- quantile(n.A1A2A3, c(0.025, 0.975))
  N.A1A2A4.CI <- quantile(n.A1A2A4, c(0.025, 0.975))
  N.A1A3A2.CI <- quantile(n.A1A3A2, c(0.025, 0.975))
  N.A1A3A4.CI <- quantile(n.A1A3A4, c(0.025, 0.975))
  N.A1A4A2.CI <- quantile(n.A1A4A2, c(0.025, 0.975))
  N.A1A4A3.CI <- quantile(n.A1A4A3, c(0.025, 0.975))
  N.A2A1A3.CI <- quantile(n.A2A1A3, c(0.025, 0.975))
  N.A2A1A4.CI <- quantile(n.A2A1A4, c(0.025, 0.975))
  N.A2A3A1.CI <- quantile(n.A2A3A1, c(0.025, 0.975))
  N.A2A3A4.CI <- quantile(n.A2A3A4, c(0.025, 0.975))
  N.A2A4A1.CI <- quantile(n.A2A4A1, c(0.025, 0.975))
  N.A2A4A3.CI <- quantile(n.A2A4A3, c(0.025, 0.975))
  N.A3A1A2.CI <- quantile(n.A3A1A2, c(0.025, 0.975))
  N.A3A1A4.CI <- quantile(n.A3A1A4, c(0.025, 0.975))
  N.A3A2A1.CI <- quantile(n.A3A2A1, c(0.025, 0.975))
  N.A3A2A4.CI <- quantile(n.A3A2A4, c(0.025, 0.975))
  N.A3A4A1.CI <- quantile(n.A3A4A1, c(0.025, 0.975))
  N.A3A4A2.CI <- quantile(n.A3A4A2, c(0.025, 0.975))
  N.A4A1A2.CI <- quantile(n.A4A1A2, c(0.025, 0.975))
  N.A4A1A3.CI <- quantile(n.A4A1A3, c(0.025, 0.975))
  N.A4A2A1.CI <- quantile(n.A4A2A1, c(0.025, 0.975))
  N.A4A2A3.CI <- quantile(n.A4A2A3, c(0.025, 0.975))
  N.A4A3A1.CI <- quantile(n.A4A3A1, c(0.025, 0.975))
  N.A4A3A2.CI <- quantile(n.A4A3A2, c(0.025, 0.975))
  
  
  result$N.A1A2A3 <- paste0(N.A1A2A3, "(",N.A1A2A3.CI[1], ",",N.A1A2A3.CI[2], ")")
  result$N.A1A2A4 <- paste0(N.A1A2A4, "(",N.A1A2A4.CI[1], ",",N.A1A2A4.CI[2], ")")
  result$N.A1A3A2 <- paste0(N.A1A3A2, "(",N.A1A3A2.CI[1], ",",N.A1A3A2.CI[2], ")")
  result$N.A1A3A4 <- paste0(N.A1A3A4, "(",N.A1A3A4.CI[1], ",",N.A1A3A4.CI[2], ")")
  result$N.A1A4A2 <- paste0(N.A1A4A2, "(",N.A1A4A2.CI[1], ",",N.A1A4A2.CI[2], ")")
  result$N.A1A4A3 <- paste0(N.A1A4A3, "(",N.A1A4A3.CI[1], ",",N.A1A4A3.CI[2], ")")
  
  result$N.A2A1A3 <- paste0(N.A2A1A3, "(",N.A2A1A3.CI[1], ",",N.A2A1A3.CI[2], ")")
  result$N.A2A1A4 <- paste0(N.A2A1A4, "(",N.A2A1A4.CI[1], ",",N.A2A1A4.CI[2], ")")
  result$N.A2A3A1 <- paste0(N.A2A3A1, "(",N.A2A3A1.CI[1], ",",N.A2A3A1.CI[2], ")")
  result$N.A2A3A4 <- paste0(N.A2A3A4, "(",N.A2A3A4.CI[1], ",",N.A2A3A4.CI[2], ")")
  result$N.A2A4A1 <- paste0(N.A2A4A1, "(",N.A2A4A1.CI[1], ",",N.A2A4A1.CI[2], ")")
  result$N.A2A4A3 <- paste0(N.A2A4A3, "(",N.A2A4A3.CI[1], ",",N.A2A4A3.CI[2], ")")
  
  result$N.A3A1A2 <- paste0(N.A3A1A2, "(",N.A3A1A2.CI[1], ",",N.A3A1A2.CI[2], ")")
  result$N.A3A1A4 <- paste0(N.A3A1A4, "(",N.A3A1A4.CI[1], ",",N.A3A1A4.CI[2], ")")
  result$N.A3A2A1 <- paste0(N.A3A2A1, "(",N.A3A2A1.CI[1], ",",N.A3A2A1.CI[2], ")")
  result$N.A3A2A4 <- paste0(N.A3A2A4, "(",N.A3A2A4.CI[1], ",",N.A3A2A4.CI[2], ")")
  result$N.A3A4A1 <- paste0(N.A3A4A1, "(",N.A3A4A1.CI[1], ",",N.A3A4A1.CI[2], ")")
  result$N.A3A4A2 <- paste0(N.A3A4A2, "(",N.A3A4A2.CI[1], ",",N.A3A4A2.CI[2], ")")
  
  result$N.A4A1A2 <- paste0(N.A4A1A2, "(",N.A4A1A2.CI[1], ",",N.A4A1A2.CI[2], ")")
  result$N.A4A1A3 <- paste0(N.A4A1A3, "(",N.A4A1A3.CI[1], ",",N.A4A1A3.CI[2], ")")
  result$N.A4A2A1 <- paste0(N.A4A2A1, "(",N.A4A2A1.CI[1], ",",N.A4A2A1.CI[2], ")")
  result$N.A4A2A3 <- paste0(N.A4A2A3, "(",N.A4A2A3.CI[1], ",",N.A4A2A3.CI[2], ")")
  result$N.A4A3A1 <- paste0(N.A4A3A1, "(",N.A4A3A1.CI[1], ",",N.A4A3A1.CI[2], ")")
  result$N.A4A3A2 <- paste0(N.A4A3A2, "(",N.A4A3A2.CI[1], ",",N.A4A3A2.CI[2], ")")
  
  
  
  # number of patients treated in the designated optimal DTR
  true.para <- c(pid.A1A2A3,pid.A1A2A4,pid.A1A3A2,pid.A1A3A4,pid.A1A4A2,pid.A1A4A3,
                 pid.A2A1A3,pid.A2A1A4,pid.A2A3A1,pid.A2A3A4,pid.A2A4A1,pid.A2A4A3,
                 pid.A3A1A2,pid.A3A1A4,pid.A3A2A1,pid.A3A2A4,pid.A3A4A1,pid.A3A4A2,
                 pid.A4A1A2,pid.A4A1A3,pid.A4A2A1,pid.A4A2A3,pid.A4A3A1,pid.A4A3A2)
  
  N.opt <- max(c(N.A1A2A3,N.A1A2A4,N.A1A3A2,N.A1A3A4,N.A1A4A2,N.A1A4A3,
                 N.A2A1A3,N.A2A1A4,N.A2A3A1,N.A2A3A4,N.A2A4A1,N.A2A4A3,
                 N.A3A1A2,N.A3A1A4,N.A3A2A1,N.A3A2A4,N.A3A4A1,N.A3A4A2,
                 N.A4A1A2,N.A4A1A3,N.A4A2A1,N.A4A2A3,N.A4A3A1,N.A4A3A2)[which(true.para==max(true.para))])
  N.opt.CI.lo <- max(c(N.A1A2A3.CI[1],N.A1A2A4.CI[1],N.A1A3A2.CI[1],N.A1A3A4.CI[1],N.A1A4A2.CI[1],N.A1A4A3.CI[1],
                       N.A2A1A3.CI[1],N.A2A1A4.CI[1],N.A2A3A1.CI[1],N.A2A3A4.CI[1],N.A2A4A1.CI[1],N.A2A4A3.CI[1],
                       N.A3A1A2.CI[1],N.A3A1A4.CI[1],N.A3A2A1.CI[1],N.A3A2A4.CI[1],N.A3A4A1.CI[1],N.A3A4A2.CI[1],
                       N.A4A1A2.CI[1],N.A4A1A3.CI[1],N.A4A2A1.CI[1],N.A4A2A3.CI[1],N.A4A3A1.CI[1],N.A4A3A2.CI[1])[which(true.para==max(true.para))])
  
  N.opt.CI.up <- max(c(N.A1A2A3.CI[2],N.A1A2A4.CI[2],N.A1A3A2.CI[2],N.A1A3A4.CI[2],N.A1A4A2.CI[2],N.A1A4A3.CI[2],
                       N.A2A1A3.CI[2],N.A2A1A4.CI[2],N.A2A3A1.CI[2],N.A2A3A4.CI[2],N.A2A4A1.CI[2],N.A2A4A3.CI[2],
                       N.A3A1A2.CI[2],N.A3A1A4.CI[2],N.A3A2A1.CI[2],N.A3A2A4.CI[2],N.A3A4A1.CI[2],N.A3A4A2.CI[2],
                       N.A4A1A2.CI[2],N.A4A1A3.CI[2],N.A4A2A1.CI[2],N.A4A2A3.CI[2],N.A4A3A1.CI[2],N.A4A3A2.CI[2])[which(true.para==max(true.para))])
  result$N.opt <- paste0(N.opt, "(",N.opt.CI.lo, ",",N.opt.CI.up, ")")
  
  # number of patients treated in the designated worst DTR
  N.wor <- min(c(N.A1A2A3,N.A1A2A4,N.A1A3A2,N.A1A3A4,N.A1A4A2,N.A1A4A3,
                 N.A2A1A3,N.A2A1A4,N.A2A3A1,N.A2A3A4,N.A2A4A1,N.A2A4A3,
                 N.A3A1A2,N.A3A1A4,N.A3A2A1,N.A3A2A4,N.A3A4A1,N.A3A4A2,
                 N.A4A1A2,N.A4A1A3,N.A4A2A1,N.A4A2A3,N.A4A3A1,N.A4A3A2)[which(true.para==min(true.para))])
  N.wor.CI.lo <- min(c(N.A1A2A3.CI[1],N.A1A2A4.CI[1],N.A1A3A2.CI[1],N.A1A3A4.CI[1],N.A1A4A2.CI[1],N.A1A4A3.CI[1],
                       N.A2A1A3.CI[1],N.A2A1A4.CI[1],N.A2A3A1.CI[1],N.A2A3A4.CI[1],N.A2A4A1.CI[1],N.A2A4A3.CI[1],
                       N.A3A1A2.CI[1],N.A3A1A4.CI[1],N.A3A2A1.CI[1],N.A3A2A4.CI[1],N.A3A4A1.CI[1],N.A3A4A2.CI[1],
                       N.A4A1A2.CI[1],N.A4A1A3.CI[1],N.A4A2A1.CI[1],N.A4A2A3.CI[1],N.A4A3A1.CI[1],N.A4A3A2.CI[1])[which(true.para==min(true.para))])
  
  N.wor.CI.up <- min(c(N.A1A2A3.CI[2],N.A1A2A4.CI[2],N.A1A3A2.CI[2],N.A1A3A4.CI[2],N.A1A4A2.CI[2],N.A1A4A3.CI[2],
                       N.A2A1A3.CI[2],N.A2A1A4.CI[2],N.A2A3A1.CI[2],N.A2A3A4.CI[2],N.A2A4A1.CI[2],N.A2A4A3.CI[2],
                       N.A3A1A2.CI[2],N.A3A1A4.CI[2],N.A3A2A1.CI[2],N.A3A2A4.CI[2],N.A3A4A1.CI[2],N.A3A4A2.CI[2],
                       N.A4A1A2.CI[2],N.A4A1A3.CI[2],N.A4A2A1.CI[2],N.A4A2A3.CI[2],N.A4A3A1.CI[2],N.A4A3A2.CI[2])[which(true.para==min(true.para))])
  
  result$N.wor <- paste0(N.wor, "(",N.wor.CI.lo, ",",N.wor.CI.up, ")")
  
  ################## 3. # of response in the trial ######################
  
  NR <- mean(nr)
  NR.CI <- quantile(nr, c(0.025, 0.975))
  result$NR <- paste0(NR, "(",NR.CI[1], ",",NR.CI[2], ")")
  
  
  ############################### 4. Type I error ######################################
  # the proportion of 10000 MC replicates of which the designated optimal DTR is estimated to be the optimal
  
  if (H0==TRUE) {
    
    
    set.seed(1006)
    
    typeI.g <- mean(est.opt.g==sample(1:24, length(est.opt.g), replace=T))
    result$typeI.g <- typeI.g
    

  }
  
  
  
  ################################# 5. Power ############################################
  # the proportion of 10000 MC replicates of which the designated optimal DTR is estimated to be the optimal
  else {
    
    true.para <- c(pid.A1A2A3,pid.A1A2A4,pid.A1A3A2,pid.A1A3A4,pid.A1A4A2,pid.A1A4A3,
                   pid.A2A1A3,pid.A2A1A4,pid.A2A3A1,pid.A2A3A4,pid.A2A4A1,pid.A2A4A3,
                   pid.A3A1A2,pid.A3A1A4,pid.A3A2A1,pid.A3A2A4,pid.A3A4A1,pid.A3A4A2,
                   pid.A4A1A2,pid.A4A1A3,pid.A4A2A1,pid.A4A2A3,pid.A4A3A1,pid.A4A3A2)
    
    pow.g <- mean(est.opt.g %in% which(true.para==max(true.para)))
    
    result$pow.g <- pow.g
    
    
    
    ## 2) power to identify a set of top two or three within 5% margin
    true.max <- which(true.para==max(true.para))
    # a function that returns the position of n-th largest
    maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
    thre <- 0.05
    
    ### G-estimator
    df.g <- cbind(mu.A1A2A3hat.g, mu.A1A2A4hat.g, mu.A1A3A2hat.g, mu.A1A3A4hat.g, mu.A1A4A2hat.g, mu.A1A4A3hat.g, 
                  mu.A2A1A3hat.g, mu.A2A1A4hat.g, mu.A2A3A1hat.g, mu.A2A3A4hat.g, mu.A2A4A1hat.g, mu.A2A4A3hat.g, 
                  mu.A3A1A2hat.g, mu.A3A1A4hat.g, mu.A3A2A1hat.g, mu.A3A2A4hat.g, mu.A3A4A1hat.g, mu.A3A4A2hat.g, 
                  mu.A4A1A2hat.g, mu.A4A1A3hat.g, mu.A4A2A1hat.g, mu.A4A2A3hat.g, mu.A4A3A1hat.g, mu.A4A3A2hat.g)
    
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
    
  }
  
  result <- as.data.frame(unlist(result))
  colnames(result) <- NULL
  return(result)
  
}



############################### Alternative scenario ############################################

pi <- c(0.5,0.3,0.2,0.4, 0.5,0.4,rep(0.35,10),0.4,0.4,0.3,0.3,rep(0.2,20))

#################### n=1000 ########################

###### SMART: p0.burn=1, p1.burn=0, p2.burn=0
set.seed(1027)
monte_three(N=10000,n=1000,pi=pi,p0.burn=1,p1.burn=0,p2.burn=0,AR=NA,c=NA,e=NA,sce="S1")

##### GO-SMART: AR1, p0.burn=0.25, c="i/n"
set.seed(1027)
monte_three(N=10000,n=1000,pi=pi,p0.burn=0.25,p1.burn=0.5,p2.burn=0.75,AR="AR1",c="n.N",e=c(0.1,0.1,0.1,0.1,0.1),sce="S1")

##### GO-SMART: AR2, p0.burn=0.25, c="i/n"
set.seed(1027)
monte_three(N=10000,n=1000,pi=pi,p0.burn=0.25,p1.burn=0.5,p2.burn=0.75,AR="AR2",c="n.N",e=c(0.1,0.1,0.1,0.1,0.1),sce="S1")


# import the datasets produced by monte_three functions above
dat_1000_pi1_SMART <- read.csv("./output/dat_three_1000_1_0_0_NA_S1_cNA.csv")
dat_1000_pi1_AR1_0.25_cn.N <- read.csv("./output/dat_three_1000_0.25_0.5_0.75_AR1_S1_cn.N.csv")
dat_1000_pi1_AR2_0.25_cn.N <- read.csv("./output/dat_three_1000_0.25_0.5_0.75_AR2_S1_cn.N.csv")

eva_1000_pi1_SMART <- eva_three(df=dat_1000_pi1_SMART,pi=pi, H0=F)
eva_1000_pi1_AR1_0.25_cn.N <- eva_three(df=dat_1000_pi1_AR1_0.25_cn.N,pi=pi, H0=F)
eva_1000_pi1_AR2_0.25_cn.N <- eva_three(df=dat_1000_pi1_AR2_0.25_cn.N,pi=pi, H0=F)
eva_1000_S1 <- cbind(eva_1000_pi1_SMART, eva_1000_pi1_AR1_0.25_cn.N, eva_1000_pi1_AR2_0.25_cn.N)

#View(eva_1000_S1)
