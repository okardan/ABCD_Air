# "Neighborhood air pollution is negatively associated with neurocognitive maturation in early adolescence"
# Omid Kardan, Chacriya Sereeyothin, Kathryn E. Schertz, Mike Angstadt, Alexander S. Weigard, Marc G. Berman, Monica D. Rosenberg

# script by Omid Kardan 
# contact omidk@med.umich.edu
# Produces Regression Tables 1, 2, S1, and S2; Produces the Violin plots in Fig 1 and 2; Produces the combined margin+scatterplots in Fig 3; Performs mediation analysis

# Requires dataforMR.csv and nihcpm_scores.csv which are respectively produced by Kardan_et_al_ABCD_pull_variables.R and Kardan_et_al_Apply_ccCPM.m.

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ppcor)
library(fastDummies)
library(sjPlot)
library(corrplot)
library(lme4)
library(plyr)
library(mediation)

setwd("E:/Omid/Air Pollution") 

dataforMR <- read.csv('dataforMR.csv')
dat0 <- dataforMR
dat0 <- dat0 %>% mutate(NBK_y2 = (tfmri_nb_all_beh_c0b_rate_y2 + tfmri_nb_all_beh_c2b_rate_y2)/2)
dat0<- dat0 %>% mutate(NBK_bl = (tfmri_nb_all_beh_c0b_rate_b + tfmri_nb_all_beh_c2b_rate_b)/2)
dat0_nbk_complete <- filter(dat0,NBK_bl+NBK_y2 !="NaN") # n = 7278 with n-back performance in both years

#white <- filter(dat0_nbk_complete, White ==1) # 4024
#Black <- filter(dat0_nbk_complete, Black ==1) # 949
#Hisp <- filter(dat0_nbk_complete, Hispanic ==1) # 1423
#Male <- filter(dat0_nbk_complete, Male ==1) # 3850

t.test(dat0_nbk_complete$NBK_y2,dat0_nbk_complete$NBK_bl, paired = TRUE) # t(7277) = 45.80, M = 7.3%
(sd(dat0_nbk_complete$NBK_y2-dat0_nbk_complete$NBK_bl)) # sd = 13.5%

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PLOTS fig 1 violin plot @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

twoRels_both <- dat0_nbk_complete[,c('NBK_bl','NBK_y2')]


twoRels_both$obs <- 1:nrow(twoRels_both)
d2 <- tidyr::gather(twoRels_both, time, value, -obs)
grad2 <- as.data.frame(cbind(c(mean(twoRels_both$NBK_bl) , mean(twoRels_both$NBK_y2) ),
                             c(1,2)  )   )

ggplot(d2, aes(time, value*100)) + 
  geom_violin(draw_quantiles = c(.5)) +
  geom_point() +
  geom_line(aes(group = obs),alpha = .05,color = "#00AFBB") +
  scale_x_discrete(limits = c('NBK_bl', 'NBK_y2'))+
  theme_minimal(base_size = 21)+
  geom_line(data = grad2,aes(x = V2,
                             y = V1*100, size=1),size=2, color = "red")+
  
  labs(y = "n-back Accuracy (%)", x = 'Age (ABCD Release)', color="", title = "") 
ggsave('violinS_nbk_bl_y2.png', dpi = 600)

####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## making Area Deprivation Index Principal Component
dat <- dat0_nbk_complete %>% dplyr::select(-c("subjectkey", "site_id_l_y2","site_id_l_b","site_id_l_y1"))
dat2<-as.data.frame(scale(dat))
tmp <- dat2[,c("reshist_addr1_adi_edu_l_b" ,
               "reshist_addr1_adi_edu_h_b", 
               "reshist_addr1_adi_work_c_b",
               "reshist_addr1_adi_income_b",
               "reshist_addr1_adi_in_dis_b",
               "reshist_addr1_adi_home_v_b",
               "reshist_addr1_adi_rent_b",
               "reshist_addr1_adi_mortg_b",
               "reshist_addr1_adi_home_o_b",
               "reshist_addr1_adi_unemp_b",
               "reshist_addr1_adi_pov_b",
               "reshist_addr1_adi_b138_b",
               "reshist_addr1_adi_sp_b" ,
               "reshist_addr1_adi_ncar_b_log", 
               "reshist_addr1_adi_ntel_b_log", 
               "reshist_addr1_adi_nplumb_b_log",
               "reshist_addr1_adi_crowd_b_log")]
tmp <- tmp[complete.cases(tmp),]
tt<- prcomp(tmp, center = TRUE, scale = TRUE)
summary(tt)

datCompleteADI <- filter(dat0_nbk_complete, !is.na(reshist_addr1_adi_edu_l_b) & !is.na(reshist_addr1_adi_edu_h_b) & !is.na(reshist_addr1_adi_work_c_b) & !is.na(reshist_addr1_adi_income_b)
                         & !is.na(reshist_addr1_adi_in_dis_b) & !is.na(reshist_addr1_adi_home_v_b) & !is.na(reshist_addr1_adi_rent_b) & !is.na(reshist_addr1_adi_mortg_b) & !is.na(reshist_addr1_adi_home_o_b) & !is.na(reshist_addr1_adi_unemp_b) & !is.na(reshist_addr1_adi_pov_b) & !is.na(reshist_addr1_adi_b138_b)
                         & !is.na(reshist_addr1_adi_sp_b) & !is.na(reshist_addr1_adi_ncar_b_log) & !is.na(reshist_addr1_adi_ntel_b_log) & !is.na(reshist_addr1_adi_nplumb_b_log) & !is.na(reshist_addr1_adi_crowd_b_log))

datCompleteADI$AreaDeprivationIndex <- tt$x[,1] # First PC will beincluded as covariate in the regressions
datCompleteADI$AreaDeprivationIndex2 <- tt$x[,2] # second PC was not related to cog performance so not included in regressions as covariate
datComplete_env <- filter(datCompleteADI,NBK_y2+ NBK_bl +  Income_y1 + HighestEd_y1 + FamilyConflict+ Male +
                            White + Black + Hispanic + Asian + Neurocog_BPPC1_GenAbility + Neurocog_BPPC2_ExecFunc+ Neurocog_BPPC3_LearnMem +
                            NeighSafety + NeighCrime + AreaDeprivationIndex + PM25+reshist_addr1_no2_b+ reshist_addr1_popdensity_b_log !="NaN") # 5256 with full covariates data

#white <- filter(datComplete_env, White ==1) # 3095
#Black <- filter(datComplete_env, Black ==1) # 570
#Hisp <- filter(datComplete_env, Hispanic ==1) # 974
#Male <- filter(datComplete_env, Male ==1) # 2774

datComplete_env$PM25_z <- scale(datComplete_env$PM25)
datComplete_env$NO2_z <- scale(datComplete_env$reshist_addr1_no2_b)
datComplete_env <- datComplete_env %>% mutate(AirPollution = (PM25_z + NO2_z)/2)

##
dat3a <- datComplete_env %>% dplyr::select(-c("subjectkey", "site_id_l_b", "site_id_l_y2","site_id_l_y1"))
dat3 <-as.data.frame(scale(dat3a)) # z-scoring variables so regression coefficients are standardized betas
dat3$subjectkey<- datComplete_env$subjectkey
dat3$site_id_l_y2<- datComplete_env$site_id_l_y2
dat3$site_id_l_y1<- datComplete_env$site_id_l_y1
dat3$site_id_l_b<- datComplete_env$site_id_l_b
dat3$delta_nbk_nonz <- dat3a$NBK_y2 - dat3a$NBK_bl # non-zscored change in n-back for use in scatterplots

### linear regressions

dat3$delta_nbk <- dat3$NBK_y2 - dat3$NBK_bl
cor.test(dat3$delta_nbk,dat3$AirPollution)  # r = -.055, t = -4.00
dat3$delta_nbk_z <- scale(dat3$delta_nbk)   # z-scoring variables so regression coefficients are standardized betas
dat3$AirPollution_z <- scale(dat3$AirPollution)

lmdeltanbk_adj <- lm(delta_nbk_z ~ NBK_bl+Neurocog_BPPC1_GenAbility + Neurocog_BPPC2_ExecFunc+ Neurocog_BPPC3_LearnMem +Male +
                       White + Black + Hispanic + Asian +  Income_y1 + HighestEd_y1 + FamilyConflict+
                       reshist_addr1_popdensity_b_log+ NeighSafety + NeighCrime + AreaDeprivationIndex +AirPollution_z  
                     , data=dat3)  # use delta_nbk_nonz to get estimate for intercept / sd() to standardize beta = 0.56462
summary(lmdeltanbk_adj)
tab_model(lmdeltanbk_adj, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, 
          string.se = "SE",dv.labels = c("Adjusted change in n-back Accuracy"),file = "~/Table1") #### Table 1  ########


lmdeltanbk <- lm(delta_nbk ~  Neurocog_BPPC1_GenAbility + Neurocog_BPPC2_ExecFunc+ Neurocog_BPPC3_LearnMem +Male +
                   White + Black + Hispanic + Asian +  Income_y1 + HighestEd_y1 + FamilyConflict+
                   reshist_addr1_popdensity_b_log+ NeighSafety + NeighCrime + AreaDeprivationIndex +AirPollution_z
           , data=dat3)
summary(lmdeltanbk)
tab_model(lmdeltanbk, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, 
          string.se = "SE",dv.labels = c("Change in n-back Accuracy"),file = "~/TableS1")  # Supplementary Table S1

#@@@@@@@@@@@@@ ccCPM Network Strength @@@@@@@@@@@##########
brainscores <- read.csv('nihcpm_scores.csv')
brainscores <- filter(brainscores, !is.na(nihcpm_score_bl) & !is.na(nihcpm_score_yr2))  # n = 1159 with ccCPM score at both timepoints
brainscores$nihcpm_score_bl_z <- brainscores$nihcpm_score_bl/sd(brainscores$nihcpm_score_bl) # unit variance for interpretability 
brainscores$nihcpm_score_yr2_z <- brainscores$nihcpm_score_yr2/sd(brainscores$nihcpm_score_yr2)


brainscores_demog <- merge(dataforMR,brainscores,by.x = "subjectkey"
                           ,by.y = "subs",all.y = TRUE)  
#white <- filter(brainscores_demog, White ==1) # 754
#Black <- filter(brainscores_demog, Black ==1) # 78
#Hisp <- filter(brainscores_demog, Hispanic ==1) # 200
#Male <- filter(brainscores_demog, Male ==1) # 554


t.test(brainscores$nihcpm_score_yr2_z,brainscores$nihcpm_score_bl_z, paired = TRUE) # t(1158) = 5.32, M = .17
(sd(brainscores$nihcpm_score_yr2_z - brainscores$nihcpm_score_bl_z)) # sd = 1.1

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PLOTS fig 2 violins @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

twoRels_both <- brainscores[,c('nihcpm_score_bl_z','nihcpm_score_yr2_z')]


twoRels_both$obs <- 1:nrow(twoRels_both)
d2 <- tidyr::gather(twoRels_both, time, value, -obs)
grad2 <- as.data.frame(cbind(c(mean(twoRels_both$nihcpm_score_bl_z) , mean(twoRels_both$nihcpm_score_yr2_z) ),
                             c(1,2)  )   )

ggplot(d2, aes(time, value*1)) + 
  geom_violin(draw_quantiles = c(.5)) +
  geom_point() +
  geom_line(aes(group = obs),alpha = .05,color = "#8F00FF") +
  scale_x_discrete(limits = c('nihcpm_score_bl_z', 'nihcpm_score_yr2_z'))+
  theme_minimal(base_size = 21)+
  geom_line(data = grad2,aes(x = V2,
                             y = V1*1, size=1),size=2, color = "red")+
  
  labs(y = "ccCPM strength (Z)", x = 'Age (ABCD Release)', color="", title = "") 
ggsave('violinS_ccCPM_bl_y2.png', dpi = 600)

####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#regressions forTable 2

df.dat3.withbrain <- merge(dat3,brainscores,by.x = "subjectkey"
                           ,by.y = "subs",all.y = TRUE)
FDs <- read.csv('~/fMRI Data/ABCD_nback_both.csv') # Frame Displacement (FD) info for fMRI runs
df.dat3.withbrain0 <- merge(df.dat3.withbrain,FDs,by.x = "subjectkey"
                            ,by.y = "subjectkey",all.x = TRUE)

df.dat3.withbrain1 <- filter(df.dat3.withbrain0,meanFD.2yr + meanFD.bl+ nihcpm_score_yr2+ nihcpm_score_bl+ NBK_bl +  Income_y1 + HighestEd_y1 + FamilyConflict+ Male +
                               White + Black + Hispanic + Asian + Neurocog_BPPC1_GenAbility + Neurocog_BPPC2_ExecFunc+ Neurocog_BPPC3_LearnMem +
                               NeighSafety + NeighCrime + AreaDeprivationIndex + AirPollution + reshist_addr1_popdensity_b_log !="NaN") # n = 878 with complete covariates data

#white <- filter(df.dat3.withbrain1, White >0) # 595
#Black <- filter(df.dat3.withbrain1, Black >0) # 54
#Hisp <- filter(df.dat3.withbrain1, Hispanic >0) # 155
#Male <- filter(df.dat3.withbrain1, Male >0) # 418


df.dat3.withbrain1 <- df.dat3.withbrain1 %>% mutate(meanFD = (meanFD.2yr + meanFD.bl)/2)
df.dat3.withbrain1$delta_cpm <- df.dat3.withbrain1$nihcpm_score_yr2 - df.dat3.withbrain1$nihcpm_score_bl # change in ccCPM strength
df.dat3.withbrain1$delta_cpm <- scale(df.dat3.withbrain1$delta_cpm) # z-score change in ccCPM
df.dat3.withbrain1$meanFD <- scale(df.dat3.withbrain1$meanFD) # z-scoring variables so regression coefficients are standardized betas
cor.test(df.dat3.withbrain1$delta_cpm,df.dat3.withbrain1$AirPollution) # r = -.042, t = -1.25, p = .213
pcor(cbind(df.dat3.withbrain1$delta_cpm,df.dat3.withbrain1$AirPollution,df.dat3.withbrain1$meanFD)) # r = -.041, t = -1.22, p = .223 adjusted for FD


lmdeltacpm_adj <- lm(delta_cpm ~  nihcpm_score_bl_z+meanFD+ Neurocog_BPPC1_GenAbility + Neurocog_BPPC2_ExecFunc+ Neurocog_BPPC3_LearnMem +Male +
                       White + Black + Hispanic + Asian +  Income_y1 + HighestEd_y1 + FamilyConflict+
                       reshist_addr1_popdensity_b_log+ NeighSafety + NeighCrime + AreaDeprivationIndex +AirPollution_z
            , data=df.dat3.withbrain1)
summary(lmdeltacpm_adj)
tab_model(lmdeltacpm_adj, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, 
          string.se = "SE",dv.labels = c("Change in cognitive_cpm_score"),file = "~/Table2")  ##### Table 2 ######


lmdeltacpm <- lm(delta_cpm ~  meanFD+ Neurocog_BPPC1_GenAbility + Neurocog_BPPC2_ExecFunc+ Neurocog_BPPC3_LearnMem +Male +
                   White + Black + Hispanic + Asian +  Income_y1 + HighestEd_y1 + FamilyConflict+
                   reshist_addr1_popdensity_b_log+ NeighSafety + NeighCrime + AreaDeprivationIndex +AirPollution_z
                     , data=df.dat3.withbrain1)
summary(lmdeltacpm)
tab_model(lmdeltacpm, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, 
          string.se = "SE",dv.labels = c("Change in cognitive_cpm_score"),file = "~/TableS2") # Supp Table S2


# corr values in the Fig 3 scatterplots
cor.test(df.dat3.withbrain1$delta_nbk,df.dat3.withbrain1$delta_cpm)
pcor(cbind(df.dat3.withbrain1$delta_nbk,df.dat3.withbrain1$delta_cpm,df.dat3.withbrain1$meanFD)) # r = .157 adjusted for FD
pcor(cbind(df.dat3.withbrain1$nihcpm_score_bl_z,df.dat3.withbrain1$NBK_bl,df.dat3.withbrain1$meanFD)) # r = .234 adjusted for FD
pcor(cbind(df.dat3.withbrain1$nihcpm_score_yr2_z,df.dat3.withbrain1$NBK_y2,df.dat3.withbrain1$meanFD)) # r = .226 adjusted for FD
##@@@@@@@@@@@@@@@@@@@@@ plot Fig 3 inset scatterplot and scatterplots with margin distributions @@@@@@@@@@@@@@@@@
#install.packages("ggExtra")
dat0_nbk_complete_s <- dat0_nbk_complete[,c('subjectkey','NBK_bl','NBK_y2')] 
df.dat3.withbrain1_s <- df.dat3.withbrain1[,c('subjectkey','nihcpm_score_bl_z','nihcpm_score_yr2_z')] 
dff <- merge(dat0_nbk_complete_s,df.dat3.withbrain1_s,by = "subjectkey"
             ,all.y = TRUE)
dff$delta_nbk <- dff$NBK_y2 - dff$NBK_bl
dff$delta_cpm <- dff$nihcpm_score_yr2_z - dff$nihcpm_score_bl_z 
m1<- mean(dff$delta_nbk)
m2<- mean(dff$delta_cpm)
mm <- as.data.frame(rbind(c(m1,m2),c(0,0)))  
p1 <- ggplot(data=dff) + 
  aes(x = delta_cpm, y = delta_nbk*100) + 
  geom_point(alpha=.09, size = 2,color = "#00AFBB") +
  geom_point(data = mm,aes(x = V2,
                            y = V1*100, size=2),shape = c(4,3), size=c(3,1), color = c("red","black"), inherit.aes = FALSE)+
  geom_smooth(method = "lm", aes(x = delta_cpm, y = delta_nbk*100), size=1, color = "#00AFBB") +
  #  annotate("text", label = "r_adj = .157, p < .001", x = 0, y = 45, size=5) +
  theme_classic(base_size = 12) + coord_cartesian(xlim = c(-3.4,3.4), ylim = c(-30,30))+
    labs(y = "\u0394 Accuracy (%)", x = "\u0394 ccCPM Strength (Z)",  title = " ") 
p1
ggsave('~/fig3_inset_scatter.png', dpi = 600)

library(ggExtra)
dat0_nbk_complete_s <- dat0_nbk_complete[,c('subjectkey','NBK_bl','NBK_y2')] 
df.dat3.withbrain1_s <- df.dat3.withbrain1[,c('subjectkey','nihcpm_score_bl_z','nihcpm_score_yr2_z')] 
dff <- merge(dat0_nbk_complete_s,df.dat3.withbrain1_s,by = "subjectkey"
                           ,all.y = TRUE)

nbl<- mean(dff$NBK_bl)
ny2<- mean(dff$NBK_y2)
cbl<- mean(dff$nihcpm_score_bl_z)
cy2<- mean(dff$nihcpm_score_yr2_z)
grad <- as.data.frame(cbind(c(nbl,ny2),c(cbl,cy2)))

dff_bl <- dff %>% dplyr::select(subjectkey, NBK_bl, nihcpm_score_bl_z) %>% dplyr::rename(NBK = NBK_bl, cpm = nihcpm_score_bl_z)
dff_bl$age = 1
dff_y2 <- dff %>% dplyr::select(subjectkey, NBK_y2, nihcpm_score_yr2_z) %>% dplyr::rename(NBK = NBK_y2, cpm = nihcpm_score_yr2_z)
dff_y2$age = 2
dff2 <- rbind(dff_bl, dff_y2)
dff2$age <- as.factor(dff2$age)

p12 <- ggplot(dff2,aes(x = cpm, y = (NBK)*100, color= age, group=age)) + 
  geom_point(alpha=.20, size = 4) +
  geom_smooth(method = "lm")+
  geom_line(data = grad,aes(x = V2,
                            y = V1*100, size=1),size=1.5, color = "black", arrow = arrow(length = unit(0.02, "npc")), inherit.aes = FALSE)+
  #geom_point(data = grad,aes(x = V2,
   #                          y = V1*100, size=1),shape = 3, size=3, color = "red", inherit.aes = FALSE)+
 # scale_color_manual(values = c("#E7B800","#A94064"), aesthetics = "color") +
  scale_color_manual(values = c("#E7B800","#FD6A02"), aesthetics = "color") +
  theme_classic(base_size = 16) + coord_cartesian(xlim = c(-2,5), ylim = c(0,103))+
  theme(legend.position = "none") +
  labs(y = "n-back Accuracy (%)", x = 'Composite Cognitive Network Strength (z)', color="", title = "") 

p12
ggsave(file = "~/fig3_alt.png", ggMarginal(p12, groupColour = TRUE, groupFill = TRUE), width = 7, height = 7, dpi = 600)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# mediation analysis
m1<-lm(delta_cpm ~  nihcpm_score_bl + meanFD+ Income_y1 + HighestEd_y1 + FamilyConflict+ Male +
         White + Black + Hispanic + Asian + Neurocog_BPPC1_GenAbility + Neurocog_BPPC2_ExecFunc+ Neurocog_BPPC3_LearnMem +
         NeighSafety + NeighCrime + AreaDeprivationIndex  + reshist_addr1_popdensity_b_log
       ,data = df.dat3.withbrain1)
delta_cpm_res <- (m1$residuals - mean(m1$residuals))/sd(m1$residuals)
df.dat3.withbrain1$delta_cpm_res <- delta_cpm_res  # delta ccCPM residulized for the covariates

m2<-lm(delta_nbk ~  NBK_bl  + Income_y1 + HighestEd_y1 + FamilyConflict+ Male +
         White + Black + Hispanic + Asian + Neurocog_BPPC1_GenAbility + Neurocog_BPPC2_ExecFunc+ Neurocog_BPPC3_LearnMem +
         NeighSafety + NeighCrime + AreaDeprivationIndex  + reshist_addr1_popdensity_b_log,data = df.dat3.withbrain1)
delta_nbk_res <- (m2$residuals - mean(m2$residuals))/sd(m2$residuals)
df.dat3.withbrain1$delta_nbk_res <- delta_nbk_res # delta Acc residulized for the covariates
# perform mediation analysis
tot.fit <- lm(delta_nbk_res ~ AirPollution ,   data = df.dat3.withbrain1)
summary(tot.fit)  # total effect
 
med.fit <- lm(delta_cpm_res ~  AirPollution , data = df.dat3.withbrain1)
summary(med.fit)
out.fit <- lm(delta_nbk_res ~ delta_cpm_res+ AirPollution ,   data = df.dat3.withbrain1)
summary(out.fit)
med23.out <- mediate(med.fit, out.fit, treat = "AirPollution", mediator = "delta_cpm_res",
                     control.value = 0, treat.value = 1,sims = 5000 )
summary(med23.out)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
