# Neighborhood air pollution is negatively associated with neurocognitive maturation in early adolescence
# contact: Omid Kardan, email: omidk@med.umich.edu
# Script to compile the variables from ABCD study data tables (run before the Matlab PLS scripts)
# Also produces the Table 1 in the manuscript (demographic breakdown of the groups and comparison)
# Requires abcd_p_demo.csv, abcd_y_lt.csv, ph_p_pds.csv, ph_y_pds.csv, mri_y_adm_info.csv,
# nc_y_nihtb.csv, led_l_no2.csv, led_l_pm25.csv, led_l_o3.csv, led_l_adi.csv
# and the pscyhopathology p-factor values from Brislin et al., 2021


library(dplyr)
library(ggplot2)
library(tidyverse)
library(ppcor)
library(fastDummies)
library(sjPlot)
library(corrplot)
library(lme4)
library(lmerTest)
library(plyr)
library(mediation)
library(effectsize)

setwd('~/New Manuscript')

dat_list <- read.csv('~/abcd_sub_event_list.csv') # a csv file containing all ABCD subids repeated in 4 rows per sub (one row for each year with years labeled as 'eventname')

##################### demog and race/ethnicity #############################
datafordemo0 <- read.csv('~/ABCD_Release5/abcd-general/abcd_p_demo.csv')
datafordemo <- merge(dat_list,datafordemo0[,c('src_subject_id','eventname','demo_comb_income_v2_l',
                                              'demo_prnt_ed_v2_l','demo_prtnr_ed_v2_l',
                                              'race_ethnicity','demo_sex_v2')],
                     by = c('src_subject_id','eventname'), all.x = TRUE)
datafordemo$Income <- datafordemo$demo_comb_income_v2_l
datafordemo$Income[datafordemo$Income==999] <- NA
datafordemo$Income[datafordemo$Income==777] <- NA
datafordemo$Income_cat = factor( datafordemo$Income, levels= 1:10, 
                                 labels = c("5000", "8500", "14000", "20500", "30000",
                                            "42500", "62500", "87500", "150000", "200000") )
datafordemo$HighestEdParent <- datafordemo$demo_prnt_ed_v2_l 
datafordemo$HighestEdParent <- as.numeric(datafordemo$HighestEdParent)
datafordemo$HighestEdParent[datafordemo$HighestEdParent == 999] <- NA
datafordemo$HighestEdParent[datafordemo$HighestEdParent == 777] <- NA
datafordemo$HighestEdPartner <- datafordemo$demo_prtnr_ed_v2_l
datafordemo$HighestEdPartner <- as.numeric(datafordemo$HighestEdPartner)
datafordemo$HighestEdPartner[datafordemo$HighestEdPartner == 999] <- NA
datafordemo$HighestEdPartner[datafordemo$HighestEdPartner == 777] <- NA
#retain the highest education out of the two parents/partners
datafordemo$HighestEd <- pmax(datafordemo$HighestEdParent, datafordemo$HighestEdPartner, na.rm =TRUE)
datafordemo$Male_bin = ifelse(datafordemo$demo_sex_v2 == 1 | datafordemo$demo_sex_v2 == 3,1,0) # 3 intsex_male and no intsex_fem
#dummy code race
datafordemo <- datafordemo %>% 
  dummy_cols(ignore_na = TRUE, select_columns = c("race_ethnicity"))
colnames(datafordemo)[15] = c("White") 
colnames(datafordemo)[16] = c("Black")
colnames(datafordemo)[17] = c("Hispanic")
colnames(datafordemo)[18] = c("Asian")
colnames(datafordemo)[19] = c("Other")

datafordem <- datafordemo[,c('src_subject_id','eventname','Income','HighestEd','Male_bin',
                             'White','Black','Hispanic','Asian','Other')]


##################### age and site and puberty #
dataforage0 <- read.csv('~/ABCD_Release5/abcd-general/abcd_y_lt.csv')
dataforage <- merge(dat_list,dataforage0[,c('src_subject_id','eventname','site_id_l',
                                            'interview_age','rel_family_id')],
                    by = c('src_subject_id','eventname'), all.x = TRUE)

dataforPubDev_p <- read.csv('~/ABCD_Release5/physical-health/ph_p_pds.csv') # pds_p_ss_male_category_2 and pds_p_ss_female_category_2
dataforPubDev_p <- dataforPubDev_p %>% mutate(pds_p_ss = coalesce(pds_p_ss_female_category_2, pds_p_ss_male_category_2))

dataforPubDev_y <- read.csv('~/ABCD_Release5/physical-health/ph_y_pds.csv')
dataforPubDev_y <- dataforPubDev_y %>% mutate(pds_y_ss = coalesce(pds_y_ss_female_category_2, pds_y_ss_male_cat_2))

dataforPubDev0 <- merge(dataforPubDev_p[,c("subid","eventname","pds_p_ss")], dataforPubDev_y[,c("subid","eventname","pds_y_ss")], by = c("subid","eventname"))
dataforPubDev0$pds_ss <- rowMeans(dataforPubDev0[,c("pds_p_ss","pds_y_ss")], na.rm = TRUE)

dataforPubDev <- merge(dat_list,dataforPubDev0[,c('subid','eventname','pds_ss')],
                       by = c('subid','eventname'), all.x = TRUE)

dataforscanman0 <- read.csv('~/ABCD_Release5/mri_y_adm_info.csv')
dataforscanman <- merge(dat_list,dataforscanman0[,c('src_subject_id','eventname','mri_info_manufacturersmn')],
                        by = c('src_subject_id','eventname'), all.x = TRUE)
dataforscanman <- dataforscanman %>% 
  dummy_cols(ignore_na = TRUE, select_columns = c("mri_info_manufacturersmn"))
colnames(dataforscanman)[5] = c("Achieva") 
colnames(dataforscanman)[6] = c("Discovery")
colnames(dataforscanman)[7] = c("Ingenia")
colnames(dataforscanman)[8] = c("Orchestra")
colnames(dataforscanman)[9] = c("Prisma")
colnames(dataforscanman)[10] = c("Pfit")
colnames(dataforscanman)[11] = c("Premier")
colnames(dataforscanman)[12] = c("UHP")

dataforscan <- dataforscanman[,c('src_subject_id','eventname','Achieva','Discovery','Ingenia',
                                 'Orchestra','Prisma','Pfit','Premier','UHP')]

##################### cog #
cog0 <- read.csv('~/abcd_mrinback02_R4.csv')
cog0$nback0_acc <- cog0$The.rate.of.correct.responses.to.0.back.stimuli.during.run.1.and.run.2
cog0$nback2_acc <- cog0$The.rate.of.correct.responses.to.2.back.stimuli.during.run.1.and.run.2
cog_nbk <- merge(dat_list,cog0[,c('src_subject_id','eventname','nback0_acc','nback2_acc')],
                 by = c('src_subject_id','eventname'), all.x = TRUE)


cog_nihtbx0 <- read.csv('~/ABCD_Release5/neurocognition/nc_y_nihtb.csv')

cog_nihtbx <- merge(dat_list,cog_nihtbx0[,c('src_subject_id','eventname','nihtbx_picvocab_uncorrected',	'nihtbx_reading_uncorrected',	
                                            'nihtbx_picture_uncorrected',	'nihtbx_flanker_uncorrected',	'nihtbx_pattern_uncorrected')],
                    by = c('src_subject_id','eventname'), all.x = TRUE)

##################### environment air pol #############################

dataforNO20 <- read.csv('~/R5_linked-external-data/led_l_no2.csv')
dataforNO2 <- merge(dat_list,dataforNO20[,c('src_subject_id','eventname','reshist_addr1_no2_2016_aavg','reshist_addr1_no2_2016_max','reshist_addr1_no2_2016_min')],
                    by = c('src_subject_id','eventname'), all.x = TRUE)


dataforPM250 <- read.csv('~/R5_linked-external-data/led_l_pm25.csv')
dataforPM25 <- merge(dat_list,dataforPM250[,c('src_subject_id','eventname','reshist_addr1_pm252016aa')],
                     by = c('src_subject_id','eventname'), all.x = TRUE)

dataforO30 <- read.csv('~/R5_linked-external-data/led_l_o3.csv')
dataforO3 <- merge(dat_list,dataforO30[,c('src_subject_id','eventname','reshist_addr1_o3_2016_max')],
                     by = c('src_subject_id','eventname'), all.x = TRUE)

##################### environment ADI #############################
adi0 <- read.csv('~/R5_linked-external-data/led_l_adi.csv')
adi <- merge(dat_list,adi0[,c('src_subject_id','eventname','reshist_addr1_adi_wsum')],
             by = c('src_subject_id','eventname'), all.x = TRUE)

##################### psychopathology #############################
# pscyhopathology p-factor values from Brislin et al., 2021

pfac0 <- read.csv('~/ABCD_lavaan_pfactor_time2.csv')
pfac <- merge(dat_list,pfac0[,c('src_subject_id','eventname','PF10_lavaan','PF10_INT_lavaan',	'PF10_EXT_lavaan')],
              by = c('src_subject_id','eventname'), all.x = TRUE)

########## combine all data 
df_total_withcog <- datafordem %>% 
  full_join(dataforage, by = c('src_subject_id','eventname'))%>%
  full_join(dataforPubDev, by = c('src_subject_id','eventname'))%>%
  full_join(dataforscan, by = c('src_subject_id','eventname'))%>%
  full_join(dataforNO2, by = c('src_subject_id','eventname'))%>%
  full_join(dataforPM25, by = c('src_subject_id','eventname'))%>%
  full_join(dataforO3, by = c('src_subject_id','eventname'))%>%  
  full_join(adi, by = c('src_subject_id','eventname'))%>%
  full_join(pfac, by = c('src_subject_id','eventname'))%>%
  full_join(cog_nbk, by = c('src_subject_id','eventname'))%>%
  full_join(cog_nihtbx, by = c('src_subject_id','eventname'))
write.csv(df_total_withcog,'~/df_total_withcog2025.csv')  

##### spread the ses and save as seperate years
df_total_withcog <- read.csv('~/New Manuscript/df_total_withcog2025.csv') 

dy <- df_total_withcog %>% mutate(Income=as.numeric(Income),HighestEd=as.numeric(HighestEd), Male_bin=as.numeric(Male_bin), White=as.numeric(White),
                                  Black=as.numeric(Black), Hispanic=as.numeric(Hispanic), Asian=as.numeric(Asian), Other=as.numeric(Other))

dy <- dy %>% group_by(src_subject_id) %>% mutate_at(vars(Income, HighestEd, Male_bin, White, Black, Hispanic, Asian, Other),
                                                    ~replace_na(.,mean(.,na.rm=TRUE))) # average the SES measures if they are available for multiple years


write.csv(dy[dy$eventname == 'baseline_year_1_arm_1',],'~/New Manuscript/df_filled_demog_withcog_Y0_more.csv') 
write.csv(dy[dy$eventname == '1_year_follow_up_y_arm_1',],'~/New Manuscript/df_filled_demog_withcog_Y1_more.csv') 
write.csv(dy[dy$eventname == '2_year_follow_up_y_arm_1',],'~/New Manuscript/df_filled_demog_withcog_Y2_more.csv') 

####################                                            ##########################
####################            ###########################
###################                                 ##########################
df_filled_demog_withcog_Y0 <- read.csv('~/New Manuscript/df_filled_demog_withcog_Y0_more.csv')
df_filled_demog_withcog_Y0$subid <- df_filled_demog_withcog_Y0$subid.x
df_filled_demog_withcog_Y0$picvocab <- df_filled_demog_withcog_Y0$nihtbx_picvocab_uncorrected
df_filled_demog_withcog_Y0$read <- df_filled_demog_withcog_Y0$nihtbx_reading_uncorrected
df_filled_demog_withcog_Y0$picture <- df_filled_demog_withcog_Y0$nihtbx_picture_uncorrected
df_filled_demog_withcog_Y0$flanker <- df_filled_demog_withcog_Y0$nihtbx_flanker_uncorrected
df_filled_demog_withcog_Y0$pattern <- df_filled_demog_withcog_Y0$nihtbx_pattern_uncorrected

df_filled_demog_withcog_Y2 <- read.csv('~/New Manuscript/df_filled_demog_withcog_Y2_more.csv')
df_filled_demog_withcog_Y2$subid <- df_filled_demog_withcog_Y2$subid.x
df_filled_demog_withcog_Y2$rel_family_id <- df_filled_demog_withcog_Y0$rel_family_id
df_filled_demog_withcog_Y2$Achieva <- df_filled_demog_withcog_Y0$Achieva
df_filled_demog_withcog_Y2$Discovery <- df_filled_demog_withcog_Y0$Discovery
df_filled_demog_withcog_Y2$Ingenia <- df_filled_demog_withcog_Y0$Ingenia
df_filled_demog_withcog_Y2$Orchestra <- df_filled_demog_withcog_Y0$Orchestra
df_filled_demog_withcog_Y2$Prisma <- df_filled_demog_withcog_Y0$Prisma
df_filled_demog_withcog_Y2$Pfit <- df_filled_demog_withcog_Y0$Pfit
df_filled_demog_withcog_Y2$Premier <- df_filled_demog_withcog_Y0$Premier
df_filled_demog_withcog_Y2$UHP <- df_filled_demog_withcog_Y0$UHP
df_filled_demog_withcog_Y2$picvocab <- df_filled_demog_withcog_Y2$nihtbx_picvocab_uncorrected
df_filled_demog_withcog_Y2$read <- df_filled_demog_withcog_Y2$nihtbx_reading_uncorrected
df_filled_demog_withcog_Y2$picture <- df_filled_demog_withcog_Y2$nihtbx_picture_uncorrected
df_filled_demog_withcog_Y2$flanker <- df_filled_demog_withcog_Y2$nihtbx_flanker_uncorrected
df_filled_demog_withcog_Y2$pattern <- df_filled_demog_withcog_Y2$nihtbx_pattern_uncorrected

# This part requires AFC values calculated from rsFC using method in Kardan et al., 2025, Dev. Cog. Neuro. as well as 
# head motion in the rest-state runs and mri_y_smr_thk_dsk for cort thickness 

wbdiff_b <- read.csv('Z:/abcd/UMich_connectomeDev/manuscript/mod_sumrs_baseline_rest_mc2024_allruns.csv')
varlist  <- c('subid','interview_age','Income','HighestEd',
              'Male_bin','White','Black','Hispanic','Asian','Other',
              'Achieva','Discovery','Ingenia','Orchestra','Prisma','Pfit','Premier','UHP',
              'site_id_l','rel_family_id','pds_ss',
              'picvocab','read','picture','flanker','pattern',
              'nback0_acc','nback2_acc','PF10_lavaan','reshist_addr1_popdensity',
              'reshist_addr1_adi_wsum','reshist_addr1_no2_2016_aavg','reshist_addr1_no2_2016_min','reshist_addr1_no2_2016_max','reshist_addr1_pm252016aa','reshist_addr1_o3_2016_max')
wbdiff_all_b <- merge(wbdiff_b,df_filled_demog_withcog_Y0[,varlist],
                      by.x = 'subids', by.y = 'subid', all.x = TRUE)  # 

wbdiff_all_b <- wbdiff_all_b[!is.nan(wbdiff_all_b$wbdiff_adult),] 

wbdiff_y2 <- read.csv('Z:/abcd/UMich_connectomeDev/manuscript/mod_sumrs_2Y_rest_mc2024_allruns.csv')
wbdiff_all_y2 <- merge(wbdiff_y2,df_filled_demog_withcog_Y2[,varlist],
                       by.x = 'subids', by.y = 'subid', all.x = TRUE)  # 

wbdiff_all_y2 <- wbdiff_all_y2[!is.nan(wbdiff_all_y2$wbdiff_adult),] 



mot_us <- read.csv('~/abcd_rsFC_FDs.csv') # mean FD in the rest fMRI runs used when calcualting AFC
wbdiff_all_bt <- merge(wbdiff_all_b,SU_us[,c('subkey','FD_y0','FD_y2','pFD_y0','pFD_y2')],
                       by.x = 'subids', by.y = 'subkey', all.x = TRUE)  # 
wbdiff_all_y2t <- merge(wbdiff_all_y2,SU_us[,c('subkey','FD_y0','FD_y2','pFD_y0','pFD_y2')],
                        by.x = 'subids', by.y = 'subkey', all.x = TRUE)  # 

dat1 <- wbdiff_all_bt %>% mutate(mat_score = wbdiff_adult - wbdiff_baby, age = interview_age/12, age_q = (scale(interview_age))^2)
dat2 <- wbdiff_all_y2t %>% mutate(mat_score = wbdiff_adult - wbdiff_baby, age = interview_age/12, age_q = (scale(interview_age))^2)

dat1$NIH5 <- rowMeans(dat1[,c('picvocab','read','picture','flanker','pattern')], na.rm=T)
dat1$NBK <- rowMeans(dat1[,c('nback0_acc','nback2_acc')], na.rm=T)

dat2$NIH5 <- rowMeans(dat2[,c('picvocab','read','picture','flanker','pattern')], na.rm=T)
dat2$NBK <- rowMeans(dat2[,c('nback0_acc','nback2_acc')], na.rm=T)
dat1<- dat1 %>% mutate(agegroup="Y0")
dat2 <- dat2 %>% mutate(agegroup="Y2")

thkness <- read.csv('~/New Manuscript/mri_y_smr_thk_dsk.csv')  # cortical thickness added
thkness0 <- thkness[thkness$eventname == 'baseline_year_1_arm_1',]
thkness2<-thkness[thkness$eventname == '2_year_follow_up_y_arm_1',] 
dat10<- merge(dat1,thkness0[,c('subids','smri_thick_cdk_mean')],
             by.x = 'subids', by.y = 'subids', all.x = TRUE)  # 
dat20<- merge(dat2,thkness2[,c('subids','smri_thick_cdk_mean')],
             by.x = 'subids', by.y = 'subids', all.x = TRUE)  # 


##########################################
write.csv(dat10,'~/df_filled_beh_Y0_more.csv')
write.csv(dat20,'~/df_filled_beh_Y2_more.csv')


#### Table 1 in the ms (run after the matlab files as it requires the high pollution and matched control groups)
all_dat <- read.csv('~/df_filled_demog_withcog_Y0_more.csv')
comp_dat<- read.csv('~/Sub_samples.csv')  # this is generated by the matlab files and contain the group labels
dats<- merge(all_dat,comp_dat,
              by.x = 'subid.x', by.y = 'usable_subs', all.x = TRUE)
length(dats$subid.x[dats$reshist_addr1_pm252016aa >9 & !is.na(dats$reshist_addr1_pm252016aa )])/length(dats$subid.x[!is.na(dats$reshist_addr1_pm252016aa )])
length(dats$subid.x[dats$reshist_addr1_pm252016aa >9 & !is.na(dats$high_pm )])/length(dats$subid.x[!is.na(dats$high_pm )])


length(dats$subid.x[dats$reshist_addr1_o3_2016_max >70 & !is.na(dats$reshist_addr1_o3_2016_max )])/length(dats$subid.x[!is.na(dats$reshist_addr1_o3_2016_max)])
length(dats$subid.x[dats$reshist_addr1_o3_2016_max >70 & !is.na(dats$high_o3 )])/length(dats$subid.x[!is.na(dats$high_o3 )])

length(dats$subid.x[dats$reshist_addr1_no2_2016_aavg >40 & !is.na(dats$reshist_addr1_no2_2016_aavg ) ])
  

length(dats$subid.x[dats$reshist_addr1_o3_2016_max >70 & dats$reshist_addr1_pm252016aa >9 & !is.na(dats$high_o3 )])

subsamp0 <- dats<- merge(all_dat,comp_dat,
                        by.x = 'subid.x', by.y = 'usable_subs', all.y = TRUE)

all_daty2 <- read.csv('~/df_filled_demog_withcog_Y2_more.csv')
subsamp<- merge(subsamp0,all_daty2[,c('subid.x','interview_age')],
              by.x = 'subid.x', by.y = 'subid.x', all.x = TRUE)

## pm groups
mean(subsamp$HighestEd[ subsamp$high_pm ==1], na.rm =TRUE) # 
sd(subsamp$HighestEd[ subsamp$high_pm ==1], na.rm = TRUE) # 
mean(subsamp$HighestEd[ subsamp$low_pm ==1], na.rm =TRUE) # 
sd(subsamp$HighestEd[ subsamp$low_pm ==1], na.rm = TRUE) # 
cohens_d(subsamp$HighestEd[subsamp$high_pm ==1],subsamp$HighestEd[subsamp$low_pm ==1])

mean(subsamp$HighestEd[ subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7], na.rm =TRUE) # 
sd(subsamp$HighestEd[ subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7], na.rm = TRUE) # 
cohens_d(subsamp$HighestEd[subsamp$high_pm ==1],subsamp$HighestEd[subsamp$high_pm ==0  & subsamp$reshist_addr1_pm252016aa<7])

# levels= 1:10, labels = c("5000", "8500", "14000", "20500", "30000","42500", 7:"62500", 8:"87500", "150000", "200000") )
mean(subsamp$Income[ subsamp$high_pm ==1], na.rm = TRUE) # 
62500 + (7.085049-7)*(87500-62500)
mean(subsamp$Income[ subsamp$low_pm ==1], na.rm = TRUE) # 
62500 + (7.482314-7)*(87500-62500)
cohens_d(subsamp$Income[subsamp$high_pm ==1],subsamp$Income[subsamp$low_pm ==1])
mean(subsamp$Income[ subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7], na.rm = TRUE) # 
87500 + (8.300464-8)*(150000-87500)
min(subsamp$Income[ subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7], na.rm = TRUE)
cohens_d(subsamp$Income[subsamp$high_pm ==1],subsamp$Income[subsamp$high_pm ==0  & subsamp$reshist_addr1_pm252016aa<7])

length(subsamp$subid.x[subsamp$White == 1 & subsamp$high_pm ==1])/348
length(subsamp$subid.x[subsamp$White == 1 & subsamp$low_pm ==1])/279 # 
cohens_d(subsamp$White[subsamp$high_pm ==1],subsamp$White[subsamp$low_pm ==1])
length(subsamp$subid.x[subsamp$White == 1 & subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])/1381
cohens_d(subsamp$White[subsamp$high_pm ==1],subsamp$White[subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])

length(subsamp$subid.x[subsamp$Hispanic == 1 & subsamp$high_pm ==1])/348
length(subsamp$subid.x[subsamp$Hispanic == 1 & subsamp$low_pm ==1])/279 # 
cohens_d(subsamp$Hispanic[subsamp$high_pm ==1],subsamp$Hispanic[subsamp$low_pm ==1])
length(subsamp$subid.x[subsamp$Hispanic == 1 & subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])/1381
cohens_d(subsamp$Hispanic[subsamp$high_pm ==1],subsamp$Hispanic[subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])

length(subsamp$subid.x[subsamp$Black == 1 & subsamp$high_pm ==1])/348
length(subsamp$subid.x[subsamp$Black == 1 & subsamp$low_pm ==1])/279 #
cohens_d(subsamp$Black[subsamp$high_pm ==1],subsamp$Black[subsamp$low_pm ==1])
length(subsamp$subid.x[subsamp$Black == 1 & subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])/1381
cohens_d(subsamp$Black[subsamp$high_pm ==1],subsamp$Black[subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])

length(subsamp$subid.x[subsamp$Asian == 1 & subsamp$high_pm ==1])/348
length(subsamp$subid.x[subsamp$Asian == 1 & subsamp$low_pm ==1])/279 #
cohens_d(subsamp$Asian[subsamp$high_pm ==1],subsamp$Asian[subsamp$low_pm ==1])
length(subsamp$subid.x[subsamp$Asian == 1 & subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])/1381
cohens_d(subsamp$Asian[subsamp$high_pm ==1],subsamp$Asian[subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])

length(subsamp$subid.x[subsamp$Other == 1 & subsamp$high_pm ==1])/348
length(subsamp$subid.x[subsamp$Other == 1 & subsamp$low_pm ==1])/279 #
cohens_d(subsamp$Other[subsamp$high_pm ==1],subsamp$Other[subsamp$low_pm ==1])
length(subsamp$subid.x[subsamp$Other == 1 & subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])/1381
cohens_d(subsamp$Other[subsamp$high_pm ==1],subsamp$Other[subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])

length(subsamp$subid.x[subsamp$Male_bin == 1 & subsamp$high_pm ==1])/348
length(subsamp$subid.x[subsamp$Male_bin == 1 & subsamp$low_pm ==1])/279 # 
cohens_d(subsamp$Male_bin[subsamp$high_pm ==1],subsamp$Male_bin[subsamp$low_pm ==1])
length(subsamp$subid.x[subsamp$Male_bin == 1 & subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])/1381
cohens_d(subsamp$Male_bin[subsamp$high_pm ==1],subsamp$Male_bin[subsamp$high_pm ==0 & subsamp$low_pm ==0 & subsamp$reshist_addr1_pm252016aa<7])

mean(subsamp$interview_age.x[ subsamp$high_pm ==1], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.x[ subsamp$high_pm ==1], na.rm = TRUE)/12 # 
mean(subsamp$interview_age.x[ subsamp$low_pm ==1], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.x[ subsamp$low_pm ==1], na.rm = TRUE)/12 # 
cohens_d(subsamp$interview_age.x[subsamp$high_pm ==1],subsamp$interview_age.x[subsamp$low_pm ==1])
mean(subsamp$interview_age.x[ subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.x[ subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7], na.rm = TRUE)/12 # 
cohens_d(subsamp$interview_age.x[subsamp$high_pm ==1],subsamp$interview_age.x[subsamp$high_pm ==0  & subsamp$reshist_addr1_pm252016aa<7])

mean(subsamp$interview_age.y[ subsamp$high_pm ==1], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.y[ subsamp$high_pm ==1], na.rm = TRUE)/12 # 
mean(subsamp$interview_age.y[ subsamp$low_pm ==1], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.y[ subsamp$low_pm ==1], na.rm = TRUE)/12 # 
cohens_d(subsamp$interview_age.y[subsamp$high_pm ==1],subsamp$interview_age.y[subsamp$low_pm ==1])
mean(subsamp$interview_age.y[ subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.y[ subsamp$high_pm ==0 & subsamp$reshist_addr1_pm252016aa<7], na.rm = TRUE)/12 # 
cohens_d(subsamp$interview_age.y[subsamp$high_pm ==1],subsamp$interview_age.y[subsamp$high_pm ==0  & subsamp$reshist_addr1_pm252016aa<7])

## O3 groups
mean(subsamp$HighestEd[ subsamp$high_o3 ==1], na.rm =TRUE) # 
sd(subsamp$HighestEd[ subsamp$high_o3 ==1], na.rm = TRUE) # 
mean(subsamp$HighestEd[ subsamp$low_o3 ==1], na.rm =TRUE) # 
sd(subsamp$HighestEd[ subsamp$low_o3 ==1], na.rm = TRUE) # 
cohens_d(subsamp$HighestEd[subsamp$high_o3 ==1],subsamp$HighestEd[subsamp$low_o3 ==1])

mean(subsamp$HighestEd[ subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60], na.rm =TRUE) # 
sd(subsamp$HighestEd[ subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60], na.rm = TRUE) # 
cohens_d(subsamp$HighestEd[subsamp$high_o3 ==1],subsamp$HighestEd[subsamp$high_o3 ==0  & subsamp$reshist_addr1_o3_2016_max<60])

# levels= 1:10, labels = c("5000", "8500", "14000", "20500", "30000","42500", 7:"62500", 8:"87500", "150000", "200000") )
mean(subsamp$Income[ subsamp$high_o3 ==1], na.rm = TRUE) # 
87500 + (8.240182-8)*(150000-87500)
mean(subsamp$Income[ subsamp$low_o3 ==1], na.rm = TRUE) # 
87500 + (8.277521-8)*(150000-87500)
cohens_d(subsamp$Income[subsamp$high_o3 ==1],subsamp$Income[subsamp$low_o3 ==1])
mean(subsamp$Income[ subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60], na.rm = TRUE) # 
87500 + (8.04652-8)*(150000-87500)
min(subsamp$Income[ subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60], na.rm = TRUE)
cohens_d(subsamp$Income[subsamp$high_o3 ==1],subsamp$Income[subsamp$high_o3 ==0  & subsamp$reshist_addr1_o3_2016_max<60])

length(subsamp$subid.x[subsamp$White == 1 & subsamp$high_o3 ==1])/355
length(subsamp$subid.x[subsamp$White == 1 & subsamp$low_o3 ==1])/324 # 
cohens_d(subsamp$White[subsamp$high_o3 ==1],subsamp$White[subsamp$low_o3 ==1])
length(subsamp$subid.x[subsamp$White == 1 & subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])/978
cohens_d(subsamp$White[subsamp$high_o3 ==1],subsamp$White[subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])

length(subsamp$subid.x[subsamp$Hispanic == 1 & subsamp$high_o3 ==1])/355
length(subsamp$subid.x[subsamp$Hispanic == 1 & subsamp$low_o3 ==1])/324 # 
cohens_d(subsamp$Hispanic[subsamp$high_o3 ==1],subsamp$Hispanic[subsamp$low_o3 ==1])
length(subsamp$subid.x[subsamp$Hispanic == 1 & subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])/978
cohens_d(subsamp$Hispanic[subsamp$high_o3 ==1],subsamp$Hispanic[subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])

length(subsamp$subid.x[subsamp$Black == 1 & subsamp$high_o3 ==1])/355
length(subsamp$subid.x[subsamp$Black == 1 & subsamp$low_o3 ==1])/324 #
cohens_d(subsamp$Black[subsamp$high_o3 ==1],subsamp$Black[subsamp$low_o3 ==1])
length(subsamp$subid.x[subsamp$Black == 1 & subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])/978
cohens_d(subsamp$Black[subsamp$high_o3 ==1],subsamp$Black[subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])

length(subsamp$subid.x[subsamp$Asian == 1 & subsamp$high_o3 ==1])/355
length(subsamp$subid.x[subsamp$Asian == 1 & subsamp$low_o3 ==1])/324 #
cohens_d(subsamp$Asian[subsamp$high_o3 ==1],subsamp$Asian[subsamp$low_o3 ==1])
length(subsamp$subid.x[subsamp$Asian == 1 & subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])/978
cohens_d(subsamp$Asian[subsamp$high_o3 ==1],subsamp$Asian[subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])

length(subsamp$subid.x[subsamp$Other == 1 & subsamp$high_o3 ==1])/355
length(subsamp$subid.x[subsamp$Other == 1 & subsamp$low_o3 ==1])/324 #
cohens_d(subsamp$Other[subsamp$high_o3 ==1],subsamp$Other[subsamp$low_o3 ==1])
length(subsamp$subid.x[subsamp$Other == 1 & subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])/978
cohens_d(subsamp$Other[subsamp$high_o3 ==1],subsamp$Other[subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])

length(subsamp$subid.x[subsamp$Male_bin == 1 & subsamp$high_o3 ==1])/355
length(subsamp$subid.x[subsamp$Male_bin == 1 & subsamp$low_o3 ==1])/324 # 
cohens_d(subsamp$Male_bin[subsamp$high_o3 ==1],subsamp$Male_bin[subsamp$low_o3 ==1])
length(subsamp$subid.x[subsamp$Male_bin == 1 & subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])/978
cohens_d(subsamp$Male_bin[subsamp$high_o3 ==1],subsamp$Male_bin[subsamp$high_o3 ==0 & subsamp$low_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])

length(subsamp$subid.x[subsamp$Male_bin == 0 & subsamp$high_o3 ==1])/355
length(subsamp$subid.x[subsamp$Male_bin == 0 & subsamp$low_o3 ==1])/324 # 
cohens_d(-1*subsamp$Male_bin[subsamp$high_o3 ==1],-1*subsamp$Male_bin[subsamp$low_o3 ==1])
length(subsamp$subid.x[subsamp$Male_bin == 0 & subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])/978
cohens_d(-1*subsamp$Male_bin[subsamp$high_o3 ==1],-1*subsamp$Male_bin[subsamp$high_o3 ==0 & subsamp$low_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60])

mean(subsamp$interview_age.x[ subsamp$high_o3 ==1], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.x[ subsamp$high_o3 ==1], na.rm = TRUE)/12 # 
mean(subsamp$interview_age.x[ subsamp$low_o3 ==1], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.x[ subsamp$low_o3 ==1], na.rm = TRUE)/12 # 
cohens_d(subsamp$interview_age.x[subsamp$high_o3 ==1],subsamp$interview_age.x[subsamp$low_o3 ==1])
mean(subsamp$interview_age.x[ subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.x[ subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60], na.rm = TRUE)/12 # 
cohens_d(subsamp$interview_age.x[subsamp$high_o3 ==1],subsamp$interview_age.x[subsamp$high_o3 ==0  & subsamp$reshist_addr1_o3_2016_max<60])

mean(subsamp$interview_age.y[ subsamp$high_o3 ==1], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.y[ subsamp$high_o3 ==1], na.rm = TRUE)/12 # 
mean(subsamp$interview_age.y[ subsamp$low_o3 ==1], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.y[ subsamp$low_o3 ==1], na.rm = TRUE)/12 # 
cohens_d(subsamp$interview_age.y[subsamp$high_o3 ==1],subsamp$interview_age.y[subsamp$low_o3 ==1])
mean(subsamp$interview_age.y[ subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60], na.rm =TRUE)/12 # 
sd(subsamp$interview_age.y[ subsamp$high_o3 ==0 & subsamp$reshist_addr1_o3_2016_max<60], na.rm = TRUE)/12 # 
cohens_d(subsamp$interview_age.y[subsamp$high_o3 ==1],subsamp$interview_age.y[subsamp$high_o3 ==0  & subsamp$reshist_addr1_o3_2016_max<60])


