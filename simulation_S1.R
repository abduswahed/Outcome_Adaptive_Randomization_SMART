#####################################################################################
##################################  Scenario S1  ####################################
#####################################################################################

# This code reproduces simulation scenario S1 under n=600 in the manuscript
# Be sure to load "fun.R" file in your own location before running this code
# Three steps included in this code:
## Step 1. simulate 10,000 datasets under SMART design, RA-SMART design, and GO-SMART (AR-1 and AR-2) design, and save the datasets in the output folder;
## Step 2. import the datasets produced in step 1, and do analysis;
## Step 3. visualize the evaluation results. Figures in the manuscript and Web Appendix are reproduced in this step.

source("./fun.R")


######################################################################################
################### Step 1: simulate 10,000 datasets #################################
####### For SMART design, RA-SMART design, and GO-SMART (AR-1 and AR-2) designs ######
######################################################################################

# pi1.A1=0.5; pi1.A2=0.35; pi1.A3=0.2; 
# pi2.A1A2=0.3; pi2.A1A3=0.4; pi2.A2A1=0.35; pi2.A2A3=0.2; pi2.A3A1=0.25; pi2.A3A2=0.1
# n=600 

pi1 <- c(0.5, 0.35, 0.2, 0.3, 0.4, 0.35, 0.2, 0.25, 0.1)

###### SMART: p0.burn=1, p1.burn=0 (Q, c and e input NA)
set.seed(1027)
monte(N=10000,n=600,RA.SMART=F,Q=NA,pi=pi1,p0.burn=1,p1.burn=0,AR=NA,c=NA,e=c(NA, NA, NA, NA, NA), sce="S1")


###### RA-SMART: p0.burn=0.25, Q=0.2 (p1.burn, c and e input NA)
set.seed(1027)
monte(N=10000,n=600,RA.SMART=T,Q=0.2,pi=pi1,p0.burn=0.25,p1.burn=NA,AR=NA,c=NA,e=c(NA, NA, NA, NA, NA),sce="S1")


##### GO-SMART: p0.burn=0.25, p1.burn=0.5, c="n.N" (meaning c=i/n in the manuscript)
AR=c("AR1", "AR2")
for (i in 1:length(AR)) {
  set.seed(1027)
  monte(N=10000,n=600,RA.SMART=F,Q=NA,pi=pi1,p0.burn=0.25,p1.burn=0.5,AR=AR[i],c="n.N",e=rep(0.1,5),sce="S1")
}


###### RA-SMART: p0.burn=0.5, Q=0.2
set.seed(1027)
monte(N=10000,n=600,RA.SMART=T,Q=0.2,pi=pi1,p0.burn=0.5,p1.burn=NA,AR=NA,c=NA,e=c(NA, NA, NA, NA, NA), sce="S1")
  

##### GO-SMART: p0.burn=0.5, p1.burn=0.75, c="n.N" (meaning c=i/n in the manuscript)
AR=c("AR1", "AR2")
for (i in 1:length(AR)) {
  set.seed(1027)
  monte(N=10000,n=600,RA.SMART=F,Q=NA,pi=pi1,p0.burn=0.5,p1.burn=0.75,AR=AR[i],c="n.N",e=rep(0.1,5),sce="S1")
}


######################################################################################
######################### Step 2: Analyze the data  ##################################
################ Bias, variance ratio, coverage probability  #########################
################## number of patients treated in each DTR  ###########################
################## number of patients response in each DTR  ##########################
################## number of patients response in the trial  #########################
########################### Type I error and power  ##################################
######################################################################################


# import the data produced by the "monte" function above
dat_600_pi1_SMART <- read.csv("./output/dat_600_S1_NA_1_0_QNA_cNA.csv")
dat_600_pi1_AR1_0.25_0.5_cn.N <- read.csv("./output/dat_600_S1_AR1_0.25_0.5_QNA_cn.N.csv")
dat_600_pi1_AR2_0.25_0.5_cn.N <- read.csv("./output/dat_600_S1_AR2_0.25_0.5_QNA_cn.N.csv")
dat_600_pi1_AR1_0.5_0.75_cn.N <- read.csv("./output/dat_600_S1_AR1_0.5_0.75_QNA_cn.N.csv")
dat_600_pi1_AR2_0.5_0.75_cn.N <- read.csv("./output/dat_600_S1_AR2_0.5_0.75_QNA_cn.N.csv")
dat_600_pi1_RA_SMART_0.25_Q0.2 <-read.csv("./output/dat_600_S1_NA_0.25_NA_Q0.2_cNA.csv")
dat_600_pi1_RA_SMART_0.5_Q0.2 <-read.csv("./output/dat_600_S1_NA_0.5_NA_Q0.2_cNA.csv")


# evaluate the operating characteristics

pi1 <- c(0.5, 0.35, 0.2, 0.3, 0.4, 0.35, 0.2, 0.25, 0.1)

# SMART
eva_600_pi1_SMART <- eva(df=dat_600_pi1_SMART,pi=pi1)

# p0.burn=0.25,p1.burn=0.5
# RA-SMART, Q=0.2
eva_600_pi1_RA_SMART_0.25_Q0.2 <- eva(df=dat_600_pi1_RA_SMART_0.25_Q0.2,pi=pi1)

# GO-SMART, c=n/N
eva_600_pi1_AR1_0.25_0.5_cn.N <- eva(df=dat_600_pi1_AR1_0.25_0.5_cn.N,pi=pi1)
eva_600_pi1_AR2_0.25_0.5_cn.N <- eva(df=dat_600_pi1_AR2_0.25_0.5_cn.N,pi=pi1)

# p0.burn=0.5,p1.burn=0.75
# RA-SMART, Q=0.2
eva_600_pi1_RA_SMART_0.5_Q0.2 <- eva(df=dat_600_pi1_RA_SMART_0.5_Q0.2,pi=pi1)

# GO-SMART, c=n/N
eva_600_pi1_AR1_0.5_0.75_cn.N <- eva(df=dat_600_pi1_AR1_0.5_0.75_cn.N,pi=pi1)
eva_600_pi1_AR2_0.5_0.75_cn.N <- eva(df=dat_600_pi1_AR2_0.5_0.75_cn.N,pi=pi1)

eva_600_S1_cn.N <- cbind(eva_600_pi1_SMART,
                         eva_600_pi1_RA_SMART_0.25_Q0.2,
                         eva_600_pi1_AR1_0.25_0.5_cn.N,
                         eva_600_pi1_AR2_0.25_0.5_cn.N,
                         eva_600_pi1_RA_SMART_0.5_Q0.2,
                         eva_600_pi1_AR1_0.5_0.75_cn.N,
                         eva_600_pi1_AR2_0.5_0.75_cn.N)
colnames(eva_600_S1_cn.N) <- c("SMART", "RA_SMART_Q0.2_0.25", 
                                "GO_SMART_AR1_0.25_0.5_cn.N","GO_SMART_AR2_0.25_0.5_cn.N",
                                "RA_SMART_Q0.2_0.5", 
                                "GO_SMART_AR1_0.5_0.75_cn.N","GO_SMART_AR2_0.5_0.75_cn.N")


######################################################################################
######################### Step 3: Visualize the results  #############################
######################################################################################

library(stringr)
library(tidyverse)
library(ggplot2)
library(latex2exp)


cbp1 <- c("#56B4E9", "#0072B2", "#E69F00", "#D55E00", "#999999")

## 1. Coverage probability (reproduce S1 panel in Figure 2)
df <- eva_600_S1_cn.N
cp <- df[str_detect(row.names(df), "cp"),]
cp$est <- rep(c("SM", "IPRW", "NIPRW", "G"),each=6)
cp$dtr <- rep(c("d(A1,A2)", "d(A1,A3)", "d(A2,A1)", "d(A2,A3)", "d(A3,A1)", "d(A3,A2)"),4)
cp <- cp %>% 
  gather(key = design, value=cp, SMART:GO_SMART_AR2_0.5_0.75_cn.N)
cp$p0 <- c(rep(1,24),rep(0.25,72), rep(0.5,72))
cp <- cp %>% 
  mutate(design=ifelse(substr(design,1,2)=="RA", substr(design,1,13), 
                       ifelse(substr(design,1,2)=="GO", substr(design,1,12), design)),
         design_big=ifelse(substr(design,1,2)=="RA", "RA-SMART", 
                           ifelse(substr(design,1,2)=="GO", "GO-SMART", "SMART")),
         ar = ifelse(design=="SMART", "SMART", substr(design,10,nchar(design))),
         design_name=ifelse(design=="SMART", "SMART", 
                            paste0(design_big, ": ", ar)))

cp <- cp %>% 
  mutate(design_name=ifelse(design_name=="RA-SMART: Q0.2", "RA-SMART", design_name),
         design_name=ifelse(design_name=="GO-SMART: AR1", "GO-SMART AR-1",
                            ifelse(design_name=="GO-SMART: AR2", "GO-SMART AR-2", design_name)),
         design_name=factor(design_name, levels=c("GO-SMART AR-1", "GO-SMART AR-2", "RA-SMART", "SMART")))

theme_set(theme_bw())
cp %>% 
  ggplot(aes(x=design_name, y=as.numeric(cp), color=design_name, shape=as.character(p0))) +
  geom_point(alpha=1, size=2.5)+
  facet_grid(dtr~est)+
  geom_hline(aes(yintercept=0.95), color="black",linetype = "dashed")+
  labs(
    x="Design",
    y="Coverage probability",
    shape=TeX("$\\p_0$"),
    color="Design",
    title = "S1"
  )+
  #scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values = cbp1)+
  #scale_color_brewer(palette = "Paired")+
  theme(#legend.position = "bottom",
    #legend.box = "vertical",
    axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_y_continuous(limits=c(as.numeric(min(cp$cp))-0.01,max(c(0.96, as.numeric(max(cp$cp))+0.05))))



## 2. Bias (reproduce S1 panel in sFigure 3 in Web Appendix)
df <- eva_600_S1_cn.N
bias <- df[str_detect(row.names(df), "bias"),]
bias$est <- rep(c("SM", "IPRW", "NIPRW", "G"),each=6)
bias$dtr <- rep(c("d(A1,A2)", "d(A1,A3)", "d(A2,A1)", "d(A2,A3)", "d(A3,A1)", "d(A3,A2)"),4)
bias <- bias %>% 
  gather(key = design, value=bias, SMART:GO_SMART_AR2_0.5_0.75_cn.N)
bias$p0 <- c(rep(1,24),rep(0.25,72), rep(0.5,72))
bias <- bias %>% 
  mutate(design=ifelse(substr(design,1,2)=="RA", substr(design,1,13), 
                       ifelse(substr(design,1,2)=="GO", substr(design,1,12), design)),
         design_big=ifelse(substr(design,1,2)=="RA", "RA-SMART", 
                           ifelse(substr(design,1,2)=="GO", "GO-SMART", "SMART")),
         ar = ifelse(design=="SMART", "SMART", substr(design,10,nchar(design))),
         design_name=ifelse(design=="SMART", "SMART", 
                            paste0(design_big, ": ", ar)))

bias <- bias %>% 
  mutate(design_name=ifelse(design_name=="RA-SMART: Q0.2", "RA-SMART", design_name),
         design_name=ifelse(design_name=="GO-SMART: AR1", "GO-SMART AR-1",
                            ifelse(design_name=="GO-SMART: AR2", "GO-SMART AR-2", design_name)),
         design_name=factor(design_name, levels=c("GO-SMART AR-1", "GO-SMART AR-2", "RA-SMART", "SMART")))


theme_set(theme_bw())
bias %>% 
  ggplot(aes(x=design_name, y=as.numeric(bias), color=design_name, shape=as.character(p0))) +
  geom_point(alpha=1, size=3)+
  facet_grid(dtr~est)+
  labs(
    x="Design",
    y="Bias",
    color="Design",
    shape=TeX("$\\p_0$"),
    title = "S1"
  )+
  scale_color_manual(values = cbp1)+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))



## 3. Type I error & Power (reproduce S1 panel in Figure 3)
df <- eva_600_S1_cn.N
pow <- df[str_detect(row.names(df), "pow"),]
pow$est <- rep(c("SM", "IPRW", "NIPRW", "G"),2)
pow$set <- rep(c("No", "Yes"),each=4)
pow <- pow %>% 
  gather(key = design, value=power, SMART:GO_SMART_AR2_0.5_0.75_cn.N)
pow$p0 <- c(rep(1,8),rep(0.25,24), rep(0.5,24))
pow <- pow %>% 
  filter(est != "SM") %>% 
  filter(set == "No") %>% 
  mutate(design=ifelse(substr(design,1,2)=="RA", substr(design,1,13), 
                       ifelse(substr(design,1,2)=="GO", substr(design,1,12), design)),
         design_big=ifelse(substr(design,1,2)=="RA", "RA-SMART", 
                           ifelse(substr(design,1,2)=="GO", "GO-SMART", "SMART")),
         ar = ifelse(design=="SMART", "SMART", substr(design,10,nchar(design))),
         design_name=ifelse(design=="SMART", "SMART", 
                            paste0(design_big, ": ", ar)),
         power=as.numeric(power),
         power=power*100) 

pow <- pow %>% 
  mutate(design_name=ifelse(design_name=="RA-SMART: Q0.2", "RA-SMART", design_name),
         design_name=ifelse(design_name=="GO-SMART: AR1", "GO-SMART AR-1",
                            ifelse(design_name=="GO-SMART: AR2", "GO-SMART AR-2", design_name)),
         design_name=factor(design_name, levels=c("GO-SMART AR-1", "GO-SMART AR-2", "RA-SMART", "SMART")))

theme_set(theme_bw())
pow %>% 
  ggplot(aes(x=design_name, y=as.numeric(power), color=design_name, shape=as.character(p0))) +
  geom_point(size=3)+
  facet_grid(~est)+
  labs(
    shape=TeX("$\\p_0$"),
    x="Design",
    y="Power (%)",
    color="Design",
    title = "S1"
  )+
  scale_color_manual(values = cbp1)+
  theme(#legend.position = "bottom",
    #     legend.box = "vertical",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits=c(60,100),breaks = c(seq(60,100,5)))


## 4. Power set plot (reproduce similar plot as in the right panel of sFigure 4 in Web Appendix (for sample size 600, not 300))
df <- eva_600_S1_cn.N
pow <- df[str_detect(row.names(df), "pow"),]
pow$est <- rep(c("SM", "IPRW", "NIPRW", "G"),2)
pow$set <- rep(c("No", "Yes"),each=4)
pow <- pow %>% 
  gather(key = design, value=power, SMART:GO_SMART_AR2_0.5_0.75_cn.N)
pow$p0 <- c(rep(1,8),rep(0.25,24), rep(0.5,24))
pow <- pow %>% 
  filter(est != "SM") %>% 
  filter(set == "Yes") %>% 
  mutate(design=ifelse(substr(design,1,2)=="RA", substr(design,1,13), 
                       ifelse(substr(design,1,2)=="GO", substr(design,1,12), design)),
         design_big=ifelse(substr(design,1,2)=="RA", "RA-SMART", 
                           ifelse(substr(design,1,2)=="GO", "GO-SMART", "SMART")),
         ar = ifelse(design=="SMART", "SMART", substr(design,10,nchar(design))),
         design_name=ifelse(design=="SMART", "SMART", 
                            paste0(design_big, ": ", ar)),
         power=as.numeric(power),
         power=power*100) 

pow <- pow %>% 
  mutate(design_name=ifelse(design_name=="RA-SMART: Q0.2", "RA-SMART", design_name),
         design_name=ifelse(design_name=="GO-SMART: AR1", "GO-SMART AR-1",
                            ifelse(design_name=="GO-SMART: AR2", "GO-SMART AR-2", design_name)),
         design_name=factor(design_name, levels=c("GO-SMART AR-1", "GO-SMART AR-2", "RA-SMART", "SMART")))


theme_set(theme_bw())
pow %>% 
  ggplot(aes(x=design_name, y=as.numeric(power), color=design_name, shape=as.character(p0))) +
  geom_point(size=3)+
  facet_grid(~est)+
  labs(
    shape=TeX("$\\p_0$"),
    x="Design",
    y="Power (%)",
    color="Design",
    title = "Identifying the best set of DTRs"
  )+
  scale_color_manual(values = cbp1)+
  theme(#legend.position = "bottom",
    #     legend.box = "vertical",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits=c(60,100),breaks = c(seq(55,100,5)))



##  5. Average number of patients treated with the optimal and worst DTRs (reproduce similar plot as in Figure 4)
df <- eva_600_S1_cn.N
num.opt <- df[str_detect(row.names(df), "N.opt"),]
num.wor <- df[str_detect(row.names(df), "N.wor"),]
num <- rbind(num.opt,num.wor)
num$cho <- c("Optimal DTR", "Worst DTR")
num <- num %>% 
  gather(key = design, value=number, SMART:GO_SMART_AR2_0.5_0.75_cn.N)
num$p0 <- c(rep(1,2),rep(0.25,6), rep(0.5,6))
num <- num %>% 
  mutate(design=ifelse(substr(design,1,2)=="RA", substr(design,1,13), 
                       ifelse(substr(design,1,2)=="GO", substr(design,1,12), design)),
         design_big=ifelse(substr(design,1,2)=="RA", "RA-SMART", 
                           ifelse(substr(design,1,2)=="GO", "GO-SMART", "SMART")),
         ar = ifelse(design=="SMART", "SMART", substr(design,10,nchar(design))),
         design_name=ifelse(design=="SMART", "SMART", 
                            paste0(design_big, ": ", ar))) %>% 
  mutate(number=substr(number, 1, str_locate(number, "\\(")[,1]-1),
         number=as.numeric(number))

num <- num %>% 
  mutate(design_name=ifelse(design_name=="RA-SMART: Q0.2", "RA-SMART", design_name),
         design_name=ifelse(design_name=="GO-SMART: AR1", "GO-SMART AR-1",
                            ifelse(design_name=="GO-SMART: AR2", "GO-SMART AR-2", design_name)),
         design_name=factor(design_name, levels=c("GO-SMART AR-1", "GO-SMART AR-2", "RA-SMART", "SMART")))

theme_set(theme_bw())
num %>% 
  ggplot(aes(x=as.character(p0), y=number, fill=design_name, group=as.factor(design_name))) +
  geom_bar(stat="identity", position=position_dodge2(preserve = "single"))+
  facet_grid(~cho)+
  labs(
    x=TeX("$\\p_0$"),
    y="Mean number of patients treated",
    fill="Design",
    title = "S1"
  )+
  scale_fill_manual(values = cbp1)+
  scale_color_brewer(palette = "Blues")+
  scale_y_continuous(limits=c(0,220), breaks = c(seq(0,220,10)))


## 6. Total number of patients responding in the trial (reproduce S1 panel in sFigure 5 in Web Appendix)
df <- eva_600_S1_cn.N
rep <- df[str_detect(row.names(df), "NR"),]
rep <- rep %>% 
  gather(key = design, value=number, SMART:GO_SMART_AR2_0.5_0.75_cn.N)
rep$p0 <- c(rep(1,1),rep(0.25,3), rep(0.5,3))
rep <- rep %>% 
  mutate(design=ifelse(substr(design,1,2)=="RA", substr(design,1,13), 
                       ifelse(substr(design,1,2)=="GO", substr(design,1,12), design)),
         design_big=ifelse(substr(design,1,2)=="RA", "RA-SMART", 
                           ifelse(substr(design,1,2)=="GO", "GO-SMART", "SMART")),
         ar = ifelse(design=="SMART", "SMART", substr(design,10,nchar(design))),
         design_name=ifelse(design=="SMART", "SMART", 
                            paste0(design_big, ": ", ar))) %>% 
  mutate(number=substr(number, 1, str_locate(number, "\\(")[,1]-1),
         number=as.numeric(number))

rep <- rep %>% 
  mutate(design_name=ifelse(design_name=="RA-SMART: Q0.2", "RA-SMART", design_name),
         design_name=ifelse(design_name=="GO-SMART: AR1", "GO-SMART AR-1",
                            ifelse(design_name=="GO-SMART: AR2", "GO-SMART AR-2", design_name)),
         design_name=factor(design_name, levels=c("GO-SMART AR-1", "GO-SMART AR-2", "RA-SMART", "SMART")))


theme_set(theme_bw())
rep %>% 
  ggplot(aes(x=as.character(p0), y=number, fill=design_name, group=as.factor(design_name))) +
  geom_bar(stat="identity", position=position_dodge2(preserve = "single"))+
  labs(
    x=TeX("$\\p_0$"),
    y="Total number of patients response",
    fill="Design",
    title = "S1"
  )+
  scale_fill_manual(values = cbp1)+
  scale_y_continuous(limits=c(0,400), breaks = c(seq(0,400,50)))


