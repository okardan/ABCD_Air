# Script for pulling different variables from the ABCD curated text files across Releases 2, 3, and 4
# including variables used in "Neighborhood air pollution is negatively associated with neurocognitive maturation in early adolescence"
# Omid Kardan, Chacriya Sereeyothin, Kathryn E. Schertz, Mike Angstadt, Alexander S. Weigard, Marc G. Berman, Monica D. Rosenberg

# script by Omid Kardan and Chacriya Sereeyothin
# contact omidk@med.umich.edu
# produces dataforMR.csv which is used in the regressions R script Kardan_et_al_Regressions_and_violins.R

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ppcor)
library(fastDummies)
library(sjPlot)
library(corrplot)
library(lme4)
library(plyr)

######################################################################################
############################## LOAD AND CLEAN VARIABLES ############################## 
######################################################################################

setwd("E:/Omid/Air Pollution") 

dataDir2.0 <- "ABCD_textfiles/ABCDCuratedRelease2.0.1_textfiles/" 
dataDir3.0 <- "ABCD_textfiles/ABCDCuratedRelease3.0_textfiles/"
dataDir4.0 <- "ABCD_textfiles/ABCDCuratedRelease4.0_textfiles/"

mysum  <- function(x)sum(x,na.rm = any(!is.na(x)))
mymean <- function(x)mean(x,na.rm = any(!is.na(x))) 

######### Read files #########  text files can be downloaded from https://nda.nih.gov/data_dictionary.html 
# Alternatively, you can copy and paste the variable names into the NIMH data analysis and exploration portal (DEAP) at
# https://deap.nimhda.org/  and download them as .csv

# each text file is read from the latest release with the file
Demographics <- read.delim(paste(dataDir2.0,"abcddemo01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)
Screener     <- read.delim(paste(dataDir3.0,"abcd_screen01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)
RAChecklist  <- read.delim(paste(dataDir3.0,"abcd_ra01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)
ScannerID    <- read.delim(paste(dataDir4.0,"abcd_mri01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1) #1 fewer column than Release 3

SiteID       <- read.delim(paste(dataDir4.0,"abcd_lt01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1) # +2-1 = 1 more column than Release 3

Family       <- read.delim(paste(dataDir3.0,"acspsw03.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)

NIH_toolbox  <- read.delim(paste(dataDir4.0,"abcd_tbss01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1) #1 fewer column than Release 3

Pearson      <- read.delim(paste(dataDir4.0,"abcd_ps01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1) #1 fewer column than Release 3

CashChoice   <- read.delim(paste(dataDir3.0,"cct01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)
LittleMan    <- read.delim(paste(dataDir2.0,"lmtp201.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)

Nback        <- read.delim(paste(dataDir4.0,"abcd_mrinback02.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1) #1 fewer column than Release 3

RecMem       <- read.delim(paste(dataDir3.0,"mribrec02.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)
SST          <- read.delim(paste(dataDir2.0,"abcd_sst02.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)
MID          <- read.delim(paste(dataDir3.0,"abcd_mid02.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)

## add environmental variables
Environment  <- read.delim(paste(dataDir4.0,"abcd_rhds01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)  # went from 280 columns in R3 to 702 columns in R4

PReport      <- read.delim(paste(dataDir4.0,"abcd_sscep01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1) #1 fewer column than Release 3

YReport      <- read.delim(paste(dataDir4.0,"abcd_sscey01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1) # went from 62 columns in R3 to 66 columns in R4

NeighSafeY   <- read.delim(paste(dataDir3.0,"abcd_nsc01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1)

## add income, race, education level (of parent)
Parent       <- read.delim(paste(dataDir4.0,"abcd_lpds01.txt",sep=""), na.strings=c(""," ","NA"), stringsAsFactors=FALSE) %>% slice(-1) # went from 130 columns in R3 to 160 columns in R4

## add general baseline neurocognitive BPPCA values from Thompson et al. 2018 https://github.com/ABCD-STUDY/BPPCA_Neurocognition_Thompson_etal_2018
NeuroCogPC <- read_csv("ABCD_textfiles/abcd_data_final.csv")  # 
NeuroCogPC <- data.frame(NeuroCogPC)

colnames(NeuroCogPC)[1]    <-c("subjectkey")

#change format of subjectkey in NeuroCogPC data frame to match the subjectkey format in the text files (NDAR_INV)
split_at_index <- function(value, index) {
  first <- substr(value, 1, index) 
  second <- substr(value, index+1, nchar(value))
  newformat <- paste(first,second, sep = "_")
  return(newformat)
}

for(i in 1:nrow(NeuroCogPC)) {
  NeuroCogPC[i,1] <- split_at_index(NeuroCogPC[i,1], 4)
}


######## Transform Income & Education Variables #########
#income (household)
Parent$Income <- Parent$demo_comb_income_v2_l
Parent$Income[Parent$Income==999] <- NA
Parent$Income[Parent$Income==777] <- NA


Parent$Income = factor( Parent$Income, levels= 1:10, 
                        labels = c("5000", "8500", "14000", "20500", "30000",
                                   "42500", "62500", "87500", "150000", "200000") )

#parent education
Parent$HighestEdParent <- Parent$demo_prnt_ed_v2_l # in R4 there's also Parent$demo_prnt_ed_v2_2yr_l
Parent$HighestEdParent <- as.numeric(Parent$HighestEdParent)
Parent$HighestEdParent[Parent$HighestEdParent == 999] <- NA
Parent$HighestEdParent[Parent$HighestEdParent == 777] <- NA


Parent$HighestEdPartner <- Parent$demo_prtnr_ed_v2_l
Parent$HighestEdPartner <- as.numeric(Parent$HighestEdPartner)
Parent$HighestEdPartner[Parent$HighestEdPartner == 999] <- NA
Parent$HighestEdPartner[Parent$HighestEdPartner == 777] <- NA

#retain the highest education out of the two parents/partners
Parent$HighestEd <- pmax(Parent$HighestEdParent, Parent$HighestEdPartner, na.rm =TRUE)

######### Segment Event-Specific Data (baseline/year 1/year 2) #########

######### Baseline ############
Demographics_b <- Demographics[which(Demographics$eventname=="baseline_year_1_arm_1"),]
#Screener     <- Screener[which(Screener$eventname=="baseline_year_1_arm_1"),]
RAChecklist_b  <- RAChecklist[which(RAChecklist$eventname=="baseline_year_1_arm_1"),]
colnames(RAChecklist_b)[10:29] <- paste(colnames(RAChecklist_b)[10:29], "b", sep="_")

ScannerID_b    <- ScannerID[which(ScannerID$eventname=="baseline_year_1_arm_1"),]
colnames(ScannerID_b)[10:16] <- paste(colnames(ScannerID_b)[10:16], "b", sep="_")

SiteID_b       <- SiteID[ which(SiteID$eventname=="baseline_year_1_arm_1"),]
colnames(SiteID_b)[10] <- paste(colnames(SiteID_b)[10], "b", sep="_")

Family_b       <- Family[ which(Family$eventname=="baseline_year_1_arm_1"),]
colnames(Family_b)[10:32] <- paste(colnames(Family_b)[10:32], "b", sep="_")

NIH_toolbox_b  <- NIH_toolbox[ which(NIH_toolbox$eventname=="baseline_year_1_arm_1"),]
colnames(NIH_toolbox_b)[10:106] <- paste(colnames(NIH_toolbox_b)[10:106], "b", sep="_")

Pearson_b      <- Pearson[ which(Pearson$eventname=="baseline_year_1_arm_1"),]
colnames(Pearson_b)[10:72] <- paste(colnames(Pearson_b)[10:72], "b", sep="_")

CashChoice_b   <- CashChoice[ which(CashChoice$eventname=="baseline_year_1_arm_1"),]
colnames(CashChoice_b)[9] <- paste(colnames(CashChoice_b)[9], "b", sep="_")

LittleMan_b    <- LittleMan[ which(LittleMan$eventname=="baseline_year_1_arm_1"),]
colnames(LittleMan_b)[9:58] <- paste(colnames(LittleMan_b)[9:58], "b", sep="_")

Nback_b        <- Nback[ which(Nback$eventname=="baseline_year_1_arm_1"),]
colnames(Nback_b)[10:590] <- paste(colnames(Nback_b)[10:590], "b", sep="_")

RecMem_b       <- RecMem[ which(RecMem$eventname=="baseline_year_1_arm_1"),]
colnames(RecMem_b)[10:55] <- paste(colnames(RecMem_b)[10:55], "b", sep="_")

SST_b          <- SST[ which(SST$eventname=="baseline_year_1_arm_1"),]
colnames(SST_b)[10:100] <- paste(colnames(SST_b)[10:100], "b", sep="_")

MID_b          <- MID[ which(MID$eventname=="baseline_year_1_arm_1"),]
colnames(MID_b)[10:250] <- paste(colnames(MID_b)[10:250], "b", sep="_")

Environment_b  <- Environment[ which(Environment$eventname=="baseline_year_1_arm_1"),]
colnames(Environment_b)[10:702] <- paste(colnames(Environment_b)[10:702], "b", sep="_")

PReport_b      <- PReport[ which(PReport$eventname=="baseline_year_1_arm_1"),]
colnames(PReport_b)[10:84] <- paste(colnames(PReport_b)[10:84], "b", sep="_")

YReport_b      <- YReport[ which(YReport$eventname=="baseline_year_1_arm_1"),]
colnames(YReport_b)[10:66] <- paste(colnames(YReport_b)[10:66], "b", sep="_")

NeighSafeY_b   <- NeighSafeY[ which(NeighSafeY$eventname=="baseline_year_1_arm_1"),]
colnames(NeighSafeY_b)[10] <- paste(colnames(NeighSafeY_b)[10], "b", sep="_")


######### Year 1 ############
SiteID_y1       <- SiteID[ which(SiteID$eventname=="1_year_follow_up_y_arm_1"),]
colnames(SiteID_y1)[10] <- paste(colnames(SiteID_y1)[10], "y1", sep="_")

Family_y1       <- Family[ which(Family$eventname=="1_year_follow_up_y_arm_1"),]
colnames(Family_y1)[10:32] <- paste(colnames(Family_y1)[10:32], "y1", sep="_")

Environment_y1  <- Environment[ which(Environment$eventname=="1_year_follow_up_y_arm_1"),]
colnames(Environment_y1)[10:702] <- paste(colnames(Environment_y1)[10:702], "y1", sep="_")

PReport_y1      <- PReport[ which(PReport$eventname=="1_year_follow_up_y_arm_1"),]
colnames(PReport_y1)[10:84] <- paste(colnames(PReport_y1)[10:84], "y1", sep="_")

YReport_y1      <- YReport[ which(YReport$eventname=="1_year_follow_up_y_arm_1"),]
colnames(YReport_y1)[10:66] <- paste(colnames(YReport_y1)[10:66], "y1", sep="_")

NeighSafeY_y1   <- NeighSafeY[ which(NeighSafeY$eventname=="1_year_follow_up_y_arm_1"),]
colnames(NeighSafeY_y1)[10] <- paste(colnames(NeighSafeY_y1)[10], "y1", sep="_")

Parent_y1       <-Parent[ which(Parent$eventname=="1_year_follow_up_y_arm_1"),]
colnames(Parent_y1)[10:164] <- paste(colnames(Parent_y1)[10:164], "y1", sep="_")


#files withOUT 1 year follow up data: Demographics, Screener, Cash Choice, Little Man, MID, Nback, NIH_toolbox, Pearson, 
#RAChecklist, RecMem, ScannerID, SST


######### Year 2 ############
RAChecklist_y2  <- RAChecklist[which(RAChecklist$eventname=="2_year_follow_up_y_arm_1"),]
colnames(RAChecklist_y2)[10:29] <- paste(colnames(RAChecklist_y2)[10:29], "y2", sep="_")

ScannerID_y2    <- ScannerID[which(ScannerID$eventname=="2_year_follow_up_y_arm_1"),]
colnames(ScannerID_y2)[10:16] <- paste(colnames(ScannerID_y2)[10:16], "y2", sep="_")

SiteID_y2       <- SiteID[ which(SiteID$eventname=="2_year_follow_up_y_arm_1"),]
colnames(SiteID_y2)[10] <- paste(colnames(SiteID_y2)[10], "y2", sep="_")

NIH_toolbox_y2  <- NIH_toolbox[ which(NIH_toolbox$eventname=="2_year_follow_up_y_arm_1"),]
colnames(NIH_toolbox_y2)[10:106] <- paste(colnames(NIH_toolbox_y2)[10:106], "y2", sep="_")

Pearson_y2      <- Pearson[ which(Pearson$eventname=="2_year_follow_up_y_arm_1"),]
colnames(Pearson_y2)[10:72] <- paste(colnames(Pearson_y2)[10:72], "y2", sep="_")

Nback_y2        <- Nback[ which(Nback$eventname=="2_year_follow_up_y_arm_1"),]
colnames(Nback_y2)[10:590] <- paste(colnames(Nback_y2)[10:590], "y2", sep="_")

RecMem_y2       <- RecMem[ which(RecMem$eventname=="2_year_follow_up_y_arm_1"),]
colnames(RecMem_y2)[10:55] <- paste(colnames(RecMem_y2)[10:55], "y2", sep="_")

MID_y2          <- MID[ which(MID$eventname=="2_year_follow_up_y_arm_1"),]
colnames(MID_y2)[10:250] <- paste(colnames(MID_y2)[10:250], "y2", sep="_")

Environment_y2  <- Environment[ which(Environment$eventname=="2_year_follow_up_y_arm_1"),]
colnames(Environment_y2)[10:702] <- paste(colnames(Environment_y2)[10:702], "y2", sep="_")

PReport_y2      <- PReport[ which(PReport$eventname=="2_year_follow_up_y_arm_1"),]
colnames(PReport_y2)[10:84] <- paste(colnames(PReport_y2)[10:84], "y2", sep="_")

YReport_y2      <- YReport[ which(YReport$eventname=="2_year_follow_up_y_arm_1"),]
colnames(YReport_y2)[10:66] <- paste(colnames(YReport_y2)[10:66], "y2", sep="_")

NeighSafeY_y2   <- NeighSafeY[ which(NeighSafeY$eventname=="2_year_follow_up_y_arm_1"),]
colnames(NeighSafeY_y2)[10] <- paste(colnames(NeighSafeY_y2)[10], "y2", sep="_")

Parent_y2       <-Parent[ which(Parent$eventname=="2_year_follow_up_y_arm_1"),]
colnames(Parent_y2)[10:164] <- paste(colnames(Parent_y2)[10:164], "y2", sep="_")


#files withOUT 2 year follow up data: Demographics, Screener, Cash Choice, Family, Little Man, SST


######### Define relevant variables #########
cols.NeuroCogPC     <-c("Neurocog_BPPC1_GenAbility", "Neurocog_BPPC2_ExecFunc","Neurocog_BPPC3_LearnMem")

##### Baseline #####
cols.Demographics_b <- c("interview_age","gender")  # changed to sex later
cols.Screener       <- c("scrn_asd","scrn_medcond_other","scrn_epls","scrn_seizure","scrn_commondx")
cols.RAChecklist_b  <- c("ra_scan_check_list_rcom_b","ra_scan_cl_mid_scan_lap_b","ra_scan_check_list_vemorc_b","ra_scan_cl_nbac_scan_lap_b","ra_scan_check_list_sstrc_b","ra_scan_cl_sst_scan_lap_b")
cols.ScannerID_b    <- c("mri_info_deviceserialnumber_b")
cols.Family_b       <- c("rel_relationship_b","rel_family_id_b","race_ethnicity_b")
cols.SiteID_b       <- c("site_id_l_b")
cols.NIH_toolbox_b  <- c("nihtbx_picvocab_uncorrected_b","nihtbx_flanker_uncorrected_b","nihtbx_list_uncorrected_b","nihtbx_cardsort_uncorrected_b","nihtbx_pattern_uncorrected_b","nihtbx_picture_uncorrected_b","nihtbx_reading_uncorrected_b","nihtbx_fluidcomp_uncorrected_b","nihtbx_cryst_uncorrected_b","nihtbx_totalcomp_uncorrected_b")
cols.Pearson_b      <- c("pea_wiscv_tss_b","pea_ravlt_sd_trial_i_tc_b","pea_ravlt_sd_trial_ii_tc_b","pea_ravlt_sd_trial_iii_tc_b","pea_ravlt_sd_trial_iv_tc_b","pea_ravlt_sd_trial_v_tc_b","pea_ravlt_sd_trial_i_tr_b","pea_ravlt_sd_trial_ii_tr_b","pea_ravlt_sd_trial_iii_tr_b","pea_ravlt_sd_trial_iv_tr_b","pea_ravlt_sd_trial_v_tr_b","pea_ravlt_sd_trial_i_ti_b","pea_ravlt_sd_trial_ii_ti_b","pea_ravlt_sd_trial_iii_ti_b","pea_ravlt_sd_trial_iv_ti_b","pea_ravlt_sd_trial_v_ti_b","pea_ravlt_sd_listb_tc_b","pea_ravlt_sd_listb_tr_b","pea_ravlt_sd_listb_ti_b","pea_ravlt_sd_trial_vi_tc_b","pea_ravlt_sd_trial_vi_tr_b","pea_ravlt_sd_trial_vi_ti_b","pea_ravlt_ld_trial_vii_tc_b","pea_ravlt_ld_trial_vii_tr_b","pea_ravlt_ld_trial_vii_ti_b")
cols.CashChoice_b   <- c("cash_choice_task_b")
cols.LittleMan_b    <- c("lmt_scr_efficiency_b","lmt_scr_perc_correct_b","lmt_scr_rt_correct_b")
cols.Nback_b        <- c("tfmri_nback_beh_switchflag_b","tfmri_nback_beh_performflag_b","tfmri_nb_all_beh_ctotal_mrt_b","tfmri_nb_all_beh_ctotal_stdrt_b","tfmri_nb_all_beh_c0b_rate_b","tfmri_nb_all_beh_c0bnf_rate_b","tfmri_nb_all_beh_c0bngf_rate_b","tfmri_nb_all_beh_c0bp_rate_b","tfmri_nb_all_beh_c0bpf_rate_b","tfmri_nb_all_beh_c2b_rate_b","tfmri_nb_all_beh_c2bnf_rate_b","tfmri_nb_all_beh_c2bngf_rate_b","tfmri_nb_all_beh_c2bp_rate_b","tfmri_nb_all_beh_c2bpf_rate_b","tfmri_nb_all_beh_cnf_rate_b","tfmri_nb_all_beh_cngf_rate_b","tfmri_nb_all_beh_cpf_rate_b","tfmri_nb_all_beh_cplace_rate_b","tfmri_nb_all_beh_ctotal_rate_b","tfmri_nb_r1_beh_c0b_rate_b","tfmri_nb_r2_beh_c0b_rate_b","tfmri_nb_r1_beh_c2b_rate_b","tfmri_nb_r2_beh_c2b_rate_b")
cols.RecMem_b       <- c("tfmri_rec_beh_switchflag_b","tfmri_rec_all_beh_posface_br_b","tfmri_rec_all_beh_posf_dpr_b","tfmri_rec_all_beh_neutface_br_b","tfmri_rec_all_beh_neutf_dp_b","tfmri_rec_all_beh_negface_br_b","tfmri_rec_all_beh_negf_dp_b","tfmri_rec_all_beh_place_br_b","tfmri_rec_all_beh_place_dp_b")

cols.MID_b          <- c("tfmri_mid_beh_switchflag_b","tfmri_mid_beh_performflag_b","tfmri_mid_all_beh_srwpfb_rate_b","tfmri_mid_all_beh_lrwpfb_rate_b","tfmri_mid_all_beh_slpfb_rate_b","tfmri_mid_all_beh_llpfb_rate_b","tfmri_mid_all_beh_ntpfb_rate_b","tfmri_mid_r1_beh_t_earnings_b","tfmri_mid_r2_beh_t_earnings_b","tfmri_mid_all_beh_t_earnings_b","tfmri_mid_all_beh_t_nt_b","tfmri_mid_all_beh_srwpfb_nt_b","tfmri_mid_all_beh_lrwpfb_nt_b","tfmri_mid_all_beh_slpfb_nt_b","tfmri_mid_all_beh_llpfb_nt_b","tfmri_mid_all_beh_ntpfb_nt_b","tfmri_mid_all_beh_srwpfb_mrt_b","tfmri_mid_all_beh_lrwpfb_mrt_b","tfmri_mid_all_beh_slpfb_mrt_b","tfmri_mid_all_beh_llpfb_mrt_b","tfmri_mid_all_beh_ntpfb_mrt_b")

#define PM2.5, NO2 exposure, walkability, proximity to roads, crime rate, drugs, ADI, family conflict
cols.Environment_b  <- c("reshist_addr1_pm25_b","reshist_addr1_no2_b","reshist_addr1_o3_2016_annavg_b","reshist_addr1_walkindex_b","reshist_addr1_proxrd_b","reshist_addr1_popdensity_b", "reshist_addr1_d1a_b", "reshist_addr1_p1vlnt_b","reshist_addr1_drugtot_b","reshist_addr1_drgsale_b","reshist_addr1_drgposs_b","reshist_addr1_dui_b","reshist_addr1_mjsale_b","reshist_addr1_adi_edu_l_b", "reshist_addr1_adi_edu_h_b","reshist_addr1_adi_work_c_b","reshist_addr1_adi_income_b", "reshist_addr1_adi_in_dis_b","reshist_addr1_adi_home_v_b", "reshist_addr1_adi_rent_b", "reshist_addr1_adi_mortg_b", "reshist_addr1_adi_home_o_b", "reshist_addr1_adi_unemp_b", "reshist_addr1_adi_pov_b","reshist_addr1_adi_b138_b", "reshist_addr1_adi_sp_b","reshist_addr1_adi_ncar_b", "reshist_addr1_adi_ntel_b","reshist_addr1_adi_nplumb_b","reshist_addr1_adi_crowd_b")

#define parent reports, youth reports, youth reports of neighborhood safety
cols.PReport_b      <- c("nsc_p_ss_mean_3_items_b","fes_p_ss_fc_pr_b")
cols.YReport_b      <- c("pmq_y_ss_mean_b","fes_y_ss_fc_pr_b","crpbi_y_ss_parent_b","srpf_y_ss_ses_b","srpf_y_ss_iiss_b","srpf_y_ss_dfs_b")
cols.NeighSafeY_b   <- c("neighborhood_crime_y_b")

##### Year 1 #####

cols.SiteID_y1       <- c("site_id_l_y1")
cols.Family_y1       <- c("rel_relationship_y1","rel_family_id_y1","race_ethnicity_y1")
cols.Environment_y1  <- c("reshist_addr1_walkindex_y1","reshist_addr1_d1a_y1", "reshist_addr1_popdensity_y1", "reshist_addr1_p1vlnt_y1","reshist_addr1_drugtot_y1","reshist_addr1_drgsale_y1","reshist_addr1_drgposs_y1","reshist_addr1_dui_y1","reshist_addr1_mjsale_y1","reshist_addr1_adi_edu_l_y1", "reshist_addr1_adi_edu_h_y1","reshist_addr1_adi_work_c_y1","reshist_addr1_adi_income_y1", "reshist_addr1_adi_in_dis_y1","reshist_addr1_adi_home_v_y1", "reshist_addr1_adi_rent_y1", "reshist_addr1_adi_mortg_y1", "reshist_addr1_adi_home_o_y1", "reshist_addr1_adi_unemp_y1", "reshist_addr1_adi_pov_y1","reshist_addr1_adi_b138_y1", "reshist_addr1_adi_sp_y1","reshist_addr1_adi_ncar_y1", "reshist_addr1_adi_ntel_y1","reshist_addr1_adi_nplumb_y1","reshist_addr1_adi_crowd_y1")
cols.PReport_y1      <- c("nsc_p_ss_mean_3_items_y1","fes_p_ss_fc_pr_y1")
cols.YReport_y1      <- c("pmq_y_ss_mean_y1","fes_y_ss_fc_pr_y1","crpbi_y_ss_parent_y1","srpf_y_ss_ses_y1","srpf_y_ss_iiss_y1","srpf_y_ss_dfs_y1")
cols.NeighSafeY_y1   <- c("neighborhood_crime_y_y1")
cols.Parent_y1       <- c("Income_y1", "HighestEd_y1")

##### Year 2 #####

cols.RAChecklist_y2  <- c("ra_scan_check_list_rcom_y2","ra_scan_cl_mid_scan_lap_y2","ra_scan_check_list_vemorc_y2","ra_scan_cl_nbac_scan_lap_y2","ra_scan_check_list_sstrc_y2","ra_scan_cl_sst_scan_lap_y2")
cols.ScannerID_y2    <- c("mri_info_deviceserialnumber_y2")
cols.SiteID_y2       <- c("site_id_l_y2")
cols.NIH_toolbox_y2  <- c("nihtbx_picvocab_uncorrected_y2","nihtbx_flanker_uncorrected_y2","nihtbx_list_uncorrected_y2","nihtbx_cardsort_uncorrected_y2","nihtbx_pattern_uncorrected_y2","nihtbx_picture_uncorrected_y2","nihtbx_reading_uncorrected_y2","nihtbx_fluidcomp_uncorrected_y2","nihtbx_cryst_uncorrected_y2","nihtbx_totalcomp_uncorrected_y2")
cols.Pearson_y2      <- c("pea_wiscv_tss_y2","pea_ravlt_sd_trial_i_tc_y2","pea_ravlt_sd_trial_ii_tc_y2","pea_ravlt_sd_trial_iii_tc_y2","pea_ravlt_sd_trial_iv_tc_y2","pea_ravlt_sd_trial_v_tc_y2","pea_ravlt_sd_trial_i_tr_y2","pea_ravlt_sd_trial_ii_tr_y2","pea_ravlt_sd_trial_iii_tr_y2","pea_ravlt_sd_trial_iv_tr_y2","pea_ravlt_sd_trial_v_tr_y2","pea_ravlt_sd_trial_i_ti_y2","pea_ravlt_sd_trial_ii_ti_y2","pea_ravlt_sd_trial_iii_ti_y2","pea_ravlt_sd_trial_iv_ti_y2","pea_ravlt_sd_trial_v_ti_y2","pea_ravlt_sd_listb_tc_y2","pea_ravlt_sd_listb_tr_y2","pea_ravlt_sd_listb_ti_y2","pea_ravlt_sd_trial_vi_tc_y2","pea_ravlt_sd_trial_vi_tr_y2","pea_ravlt_sd_trial_vi_ti_y2","pea_ravlt_ld_trial_vii_tc_y2","pea_ravlt_ld_trial_vii_tr_y2","pea_ravlt_ld_trial_vii_ti_y2")
cols.Nback_y2        <- c("tfmri_nback_beh_switchflag_y2","tfmri_nback_beh_performflag_y2","tfmri_nb_all_beh_ctotal_mrt_y2","tfmri_nb_all_beh_ctotal_stdrt_y2","tfmri_nb_all_beh_c0b_rate_y2","tfmri_nb_all_beh_c0bnf_rate_y2","tfmri_nb_all_beh_c0bngf_rate_y2","tfmri_nb_all_beh_c0bp_rate_y2","tfmri_nb_all_beh_c0bpf_rate_y2","tfmri_nb_all_beh_c2b_rate_y2","tfmri_nb_all_beh_c2bnf_rate_y2","tfmri_nb_all_beh_c2bngf_rate_y2","tfmri_nb_all_beh_c2bp_rate_y2","tfmri_nb_all_beh_c2bpf_rate_y2","tfmri_nb_all_beh_cnf_rate_y2","tfmri_nb_all_beh_cngf_rate_y2","tfmri_nb_all_beh_cpf_rate_y2","tfmri_nb_all_beh_cplace_rate_y2","tfmri_nb_all_beh_ctotal_rate_y2","tfmri_nb_r1_beh_c0b_rate_y2","tfmri_nb_r2_beh_c0b_rate_y2","tfmri_nb_r1_beh_c2b_rate_y2","tfmri_nb_r2_beh_c2b_rate_y2")
cols.RecMem_y2       <- c("tfmri_rec_beh_switchflag_y2","tfmri_rec_all_beh_posface_br_y2","tfmri_rec_all_beh_posf_dpr_y2","tfmri_rec_all_beh_neutface_br_y2","tfmri_rec_all_beh_neutf_dp_y2","tfmri_rec_all_beh_negface_br_y2","tfmri_rec_all_beh_negf_dp_y2","tfmri_rec_all_beh_place_br_y2","tfmri_rec_all_beh_place_dp_y2")
cols.MID_y2          <- c("tfmri_mid_beh_switchflag_y2","tfmri_mid_beh_performflag_y2","tfmri_mid_all_beh_srwpfb_rate_y2","tfmri_mid_all_beh_lrwpfb_rate_y2","tfmri_mid_all_beh_slpfb_rate_y2","tfmri_mid_all_beh_llpfb_rate_y2","tfmri_mid_all_beh_ntpfb_rate_y2","tfmri_mid_r1_beh_t_earnings_y2","tfmri_mid_r2_beh_t_earnings_y2","tfmri_mid_all_beh_t_earnings_y2","tfmri_mid_all_beh_t_nt_y2","tfmri_mid_all_beh_srwpfb_nt_y2","tfmri_mid_all_beh_lrwpfb_nt_y2","tfmri_mid_all_beh_slpfb_nt_y2","tfmri_mid_all_beh_llpfb_nt_y2","tfmri_mid_all_beh_ntpfb_nt_y2","tfmri_mid_all_beh_srwpfb_mrt_y2","tfmri_mid_all_beh_lrwpfb_mrt_y2","tfmri_mid_all_beh_slpfb_mrt_y2","tfmri_mid_all_beh_llpfb_mrt_y2","tfmri_mid_all_beh_ntpfb_mrt_y2")
cols.Environment_y2  <- c("reshist_addr1_walkindex_y2","reshist_addr1_d1a_y2", "reshist_addr1_popdensity_y2", "reshist_addr1_p1vlnt_y2","reshist_addr1_drugtot_y2","reshist_addr1_drgsale_y2","reshist_addr1_drgposs_y2","reshist_addr1_dui_y2","reshist_addr1_mjsale_y2","reshist_addr1_adi_edu_l_y2", "reshist_addr1_adi_edu_h_y2","reshist_addr1_adi_work_c_y2","reshist_addr1_adi_income_y2", "reshist_addr1_adi_in_dis_y2","reshist_addr1_adi_home_v_y2", "reshist_addr1_adi_rent_y2", "reshist_addr1_adi_mortg_y2", "reshist_addr1_adi_home_o_y2", "reshist_addr1_adi_unemp_y2", "reshist_addr1_adi_pov_y2","reshist_addr1_adi_b138_y2", "reshist_addr1_adi_sp_y2","reshist_addr1_adi_ncar_y2", "reshist_addr1_adi_ntel_y2","reshist_addr1_adi_nplumb_y2","reshist_addr1_adi_crowd_y2")
cols.PReport_y2      <- c("nsc_p_ss_mean_3_items_y2","fes_p_ss_fc_pr_y2")
cols.YReport_y2      <- c("pmq_y_ss_mean_y2","fes_y_ss_fc_pr_y2","srpf_y_ss_ses_y2","srpf_y_ss_iiss_y2","srpf_y_ss_dfs_y2")
cols.NeighSafeY_y2   <- c("neighborhood_crime_y_y2")
cols.Parent_y2       <- c("Income_y2")


######### Retain relevant variables for baseline #########
NeuroCogPC          <-unique(subset(NeuroCogPC,      select = c("subjectkey", cols.NeuroCogPC)))

Demographics_b      <- unique(subset(Demographics_b, select = c("subjectkey", cols.Demographics_b)))
Screener            <- unique(subset(Screener,       select = c("subjectkey", cols.Screener)))
RAChecklist_b       <- unique(subset(RAChecklist_b,  select = c("subjectkey", cols.RAChecklist_b)))
ScannerID_b         <- unique(subset(ScannerID_b,    select = c("subjectkey", cols.ScannerID_b)))
Family_b            <- unique(subset(Family_b,       select = c("subjectkey", cols.Family_b)))
SiteID_b            <- unique(subset(SiteID_b,       select = c("subjectkey", cols.SiteID_b)))
NIH_toolbox_b       <- unique(subset(NIH_toolbox_b,  select = c("subjectkey", cols.NIH_toolbox_b)))
Pearson_b           <- unique(subset(Pearson_b,      select = c("subjectkey", cols.Pearson_b)))
CashChoice_b        <- unique(subset(CashChoice_b,   select = c("subjectkey", cols.CashChoice_b)))
LittleMan_b         <- unique(subset(LittleMan_b,    select = c("subjectkey", cols.LittleMan_b)))
Nback_b             <- unique(subset(Nback_b,        select = c("subjectkey", cols.Nback_b)))
RecMem_b            <- unique(subset(RecMem_b,       select = c("subjectkey", cols.RecMem_b)))
#SST_b               <- unique(subset(SST_b,          select = c("subjectkey", cols.SST_b)))
MID_b               <- unique(subset(MID_b,          select = c("subjectkey", cols.MID_b)))
Environment_b       <- unique(subset(Environment_b,  select = c("subjectkey", cols.Environment_b)))
PReport_b           <- unique(subset(PReport_b,      select = c("subjectkey", cols.PReport_b)))
YReport_b           <- unique(subset(YReport_b,      select = c("subjectkey", cols.YReport_b)))
NeighSafeY_b        <- unique(subset(NeighSafeY_b,   select = c("subjectkey", cols.NeighSafeY_b)))

######### Retain relevant variables for Y1 #########
SiteID_y1            <- unique(subset(SiteID_y1,       select = c("subjectkey", cols.SiteID_y1)))
Family_y1            <- unique(subset(Family_y1,       select = c("subjectkey", cols.Family_y1)))
Environment_y1       <- unique(subset(Environment_y1,  select = c("subjectkey", cols.Environment_y1)))
PReport_y1           <- unique(subset(PReport_y1,      select = c("subjectkey", cols.PReport_y1)))
YReport_y1           <- unique(subset(YReport_y1,      select = c("subjectkey", cols.YReport_y1)))
NeighSafeY_y1        <- unique(subset(NeighSafeY_y1,   select = c("subjectkey", cols.NeighSafeY_y1)))
Parent_y1            <- unique(subset(Parent_y1,       select = c("subjectkey", cols.Parent_y1)))


######### Retain relevant variables for Y2 #########
RAChecklist_y2       <- unique(subset(RAChecklist_y2,  select = c("subjectkey", cols.RAChecklist_y2)))
ScannerID_y2         <- unique(subset(ScannerID_y2,    select = c("subjectkey", cols.ScannerID_y2)))
SiteID_y2            <- unique(subset(SiteID_y2,       select = c("subjectkey", cols.SiteID_y2)))
NIH_toolbox_y2       <- unique(subset(NIH_toolbox_y2,  select = c("subjectkey", cols.NIH_toolbox_y2)))
Pearson_y2           <- unique(subset(Pearson_y2,      select = c("subjectkey", cols.Pearson_y2)))
Nback_y2             <- unique(subset(Nback_y2,        select = c("subjectkey", cols.Nback_y2)))
RecMem_y2            <- unique(subset(RecMem_y2,       select = c("subjectkey", cols.RecMem_y2)))
MID_y2               <- unique(subset(MID_y2,          select = c("subjectkey", cols.MID_y2)))
Environment_y2       <- unique(subset(Environment_y2,  select = c("subjectkey", cols.Environment_y2)))
PReport_y2           <- unique(subset(PReport_y2,      select = c("subjectkey", cols.PReport_y2)))
YReport_y2           <- unique(subset(YReport_y2,      select = c("subjectkey", cols.YReport_y2)))
NeighSafeY_y2        <- unique(subset(NeighSafeY_y2,   select = c("subjectkey", cols.NeighSafeY_y2)))
Parent_y2            <- unique(subset(Parent_y2,       select = c("subjectkey", cols.Parent_y2)))


######### Convert variables to numeric #########

Demographics_b[, "interview_age"]     <- lapply("interview_age",    function(x) as.numeric(Demographics_b[[x]]))
Screener[, cols.Screener]             <- lapply(cols.Screener,      function(x) as.numeric(Screener[[x]]))
RAChecklist_b[, cols.RAChecklist_b]   <- lapply(cols.RAChecklist_b, function(x) as.numeric(RAChecklist_b[[x]]))
Family_b[, cols.Family_b]             <- lapply(cols.Family_b,      function(x) as.numeric(Family_b[[x]]))
NIH_toolbox_b[, cols.NIH_toolbox_b]   <- lapply(cols.NIH_toolbox_b, function(x) as.numeric(NIH_toolbox_b[[x]]))
Pearson_b[, cols.Pearson_b]           <- lapply(cols.Pearson_b,     function(x) as.numeric(Pearson_b[[x]]))
CashChoice_b[, cols.CashChoice_b]     <- lapply(cols.CashChoice_b,  function(x) as.numeric(CashChoice_b[[x]]))
LittleMan_b[, cols.LittleMan_b]       <- lapply(cols.LittleMan_b,   function(x) as.numeric(LittleMan_b[[x]]))
Nback_b[, cols.Nback_b]               <- lapply(cols.Nback_b,       function(x) as.numeric(Nback_b[[x]]))
RecMem_b[, cols.RecMem_b]             <- lapply(cols.RecMem_b,      function(x) as.numeric(RecMem_b[[x]]))
#SST_b[, cols.SST_b]                   <- lapply(cols.SST_b,         function(x) as.numeric(SST_b[[x]]))
MID_b[, cols.MID_b]                   <- lapply(cols.MID_b,         function(x) as.numeric(MID_b[[x]]))
Environment_b[, cols.Environment_b]   <- lapply(cols.Environment_b, function(x) as.numeric(Environment_b[[x]]))
PReport_b[, cols.PReport_b]           <- lapply(cols.PReport_b,     function(x) as.numeric(PReport_b[[x]]))
YReport_b[, cols.YReport_b]           <- lapply(cols.YReport_b,     function(x) as.numeric(YReport_b[[x]]))
NeighSafeY_b[, cols.NeighSafeY_b]     <- lapply(cols.NeighSafeY_b,  function(x) as.numeric(NeighSafeY_b[[x]]))


Family_y1[, cols.Family_y1]           <- lapply(cols.Family_y1,      function(x) as.numeric(Family_y1[[x]]))
Environment_y1[, cols.Environment_y1] <- lapply(cols.Environment_y1, function(x) as.numeric(Environment_y1[[x]]))
PReport_y1[, cols.PReport_y1]         <- lapply(cols.PReport_y1,     function(x) as.numeric(PReport_y1[[x]]))
YReport_y1[, cols.YReport_y1]         <- lapply(cols.YReport_y1,     function(x) as.numeric(YReport_y1[[x]]))
NeighSafeY_y1[, cols.NeighSafeY_y1]   <- lapply(cols.NeighSafeY_y1,  function(x) as.numeric(NeighSafeY_y1[[x]]))
Parent_y1[, cols.Parent_y1]           <- lapply(cols.Parent_y1,      function(x) as.numeric(Parent_y1[[x]]))


NIH_toolbox_y2[, cols.NIH_toolbox_y2] <- lapply(cols.NIH_toolbox_y2, function(x) as.numeric(NIH_toolbox_y2[[x]]))
RAChecklist_y2[, cols.RAChecklist_y2] <- lapply(cols.RAChecklist_y2, function(x) as.numeric(RAChecklist_y2[[x]]))
Pearson_y2[, cols.Pearson_y2]         <- lapply(cols.Pearson_y2,     function(x) as.numeric(Pearson_y2[[x]]))
Nback_y2[, cols.Nback_y2]             <- lapply(cols.Nback_y2,       function(x) as.numeric(Nback_y2[[x]]))
RecMem_y2[, cols.RecMem_y2]           <- lapply(cols.RecMem_y2,      function(x) as.numeric(RecMem_y2[[x]]))
MID_y2[, cols.MID_y2]                 <- lapply(cols.MID_y2,         function(x) as.numeric(MID_y2[[x]]))
Environment_y2[, cols.Environment_y2] <- lapply(cols.Environment_y2, function(x) as.numeric(Environment_y2[[x]]))
PReport_y2[, cols.PReport_y2]         <- lapply(cols.PReport_y2,     function(x) as.numeric(PReport_y2[[x]]))
YReport_y2[, cols.YReport_y2]         <- lapply(cols.YReport_y2,     function(x) as.numeric(YReport_y2[[x]]))
NeighSafeY_y2[, cols.NeighSafeY_y2]   <- lapply(cols.NeighSafeY_y2,  function(x) as.numeric(NeighSafeY_y2[[x]]))
Parent_y2[, cols.Parent_y2]           <- lapply(cols.Parent_y2,      function(x) as.numeric(Parent_y2[[x]]))


######### Add performance measure columns #########
RecMem_b$overall_dprime_b             <- apply(RecMem_b[c('tfmri_rec_all_beh_posf_dpr_b', 'tfmri_rec_all_beh_neutf_dp_b', 'tfmri_rec_all_beh_negf_dp_b', 'tfmri_rec_all_beh_place_dp_b')], 1, mymean)
RecMem_y2$overall_dprime_y2           <- apply(RecMem_y2[c('tfmri_rec_all_beh_posf_dpr_y2', 'tfmri_rec_all_beh_neutf_dp_y2', 'tfmri_rec_all_beh_negf_dp_y2', 'tfmri_rec_all_beh_place_dp_y2')], 1, mymean)
MID_b$mean_earnings_b                 <- apply(MID_b[c('tfmri_mid_r1_beh_t_earnings_b', 'tfmri_mid_r2_beh_t_earnings_b')], 1, mymean)
MID_y2$mean_earnings_y2               <- apply(MID_y2[c('tfmri_mid_r1_beh_t_earnings_y2', 'tfmri_mid_r2_beh_t_earnings_y2')], 1, mymean)

######### Remove cash choice option 3 ("don't know") #########
CashChoice_b$cash_choice_task_no3_b <- CashChoice_b$cash_choice_task_b
CashChoice_b$cash_choice_task_no3_b[CashChoice_b$cash_choice_task_no3_b == 3] <- NA


######### Merge, clean, crop data #########
data.merge_b  <- Reduce(function(x,y) merge(x = x, y = y, by = "subjectkey", all.x = TRUE, all.y = TRUE), list(Demographics_b, Screener, RAChecklist_b, ScannerID_b, SiteID_b, Family_b, NIH_toolbox_b, Pearson_b, CashChoice_b, LittleMan_b, Nback_b, RecMem_b, MID_b, Environment_b, PReport_b, YReport_b, NeighSafeY_b, NeuroCogPC)) # OK removed SST_b
data.merge_y1 <- Reduce(function(x,y) merge(x = x, y = y, by = "subjectkey", all.x = TRUE, all.y = TRUE), list(SiteID_y1, Family_y1, Environment_y1, PReport_y1, YReport_y1, NeighSafeY_y1, Parent_y1))
data.merge_y2 <- Reduce(function(x,y) merge(x = x, y = y, by = "subjectkey", all.x = TRUE, all.y = TRUE), list(RAChecklist_y2, ScannerID_y2, SiteID_y2, NIH_toolbox_y2, Pearson_y2, Nback_y2, RecMem_y2, MID_y2, Environment_y2, PReport_y2, YReport_y2, NeighSafeY_y2, Parent_y2))

data.merge_mega <- Reduce(function(x,y) merge(x = x, y = y, by = "subjectkey", all.x = TRUE, all.y = TRUE), list(Demographics_b, Screener, RAChecklist_b, ScannerID_b, SiteID_b, Family_b, NIH_toolbox_b, Pearson_b, CashChoice_b, LittleMan_b, Nback_b, RecMem_b, MID_b, Environment_b, PReport_b, YReport_b, NeighSafeY_b,SiteID_y1, Family_y1, Environment_y1, PReport_y1, YReport_y1, NeighSafeY_y1,RAChecklist_y2, ScannerID_y2, SiteID_y2, NIH_toolbox_y2, Pearson_y2, Nback_y2, RecMem_y2, MID_y2, Environment_y2, PReport_y2, YReport_y2, NeighSafeY_y2))  # OK removed SST_b

#add variable names for correlation matrix
data.crop_b   <- data.merge_b[ which(data.merge_b$scrn_asd==0 & (data.merge_b$scrn_epls!=1 | is.na(data.merge_b$scrn_epls))), ] #remove participants with autism and epilepsy

data.crop_b  <- subset(data.crop_b, select = c(subjectkey, rel_family_id_b, site_id_l_b, ra_scan_cl_mid_scan_lap_b, ra_scan_cl_nbac_scan_lap_b, ra_scan_cl_sst_scan_lap_b, interview_age, gender, race_ethnicity_b, tfmri_nback_beh_performflag_b, tfmri_mid_beh_performflag_b, nihtbx_list_uncorrected_b, nihtbx_picvocab_uncorrected_b, nihtbx_reading_uncorrected_b, tfmri_nb_all_beh_c2b_rate_b, pea_wiscv_tss_b, nihtbx_picture_uncorrected_b, pea_ravlt_sd_trial_vi_tc_b, pea_ravlt_ld_trial_vii_tc_b, tfmri_nb_all_beh_c0b_rate_b, nihtbx_cardsort_uncorrected_b, nihtbx_flanker_uncorrected_b, lmt_scr_efficiency_b, overall_dprime_b, nihtbx_pattern_uncorrected_b, mean_earnings_b, cash_choice_task_no3_b,        reshist_addr1_pm25_b, reshist_addr1_no2_b, reshist_addr1_o3_2016_annavg_b, reshist_addr1_walkindex_b, reshist_addr1_proxrd_b, reshist_addr1_popdensity_b, reshist_addr1_d1a_b, reshist_addr1_p1vlnt_b, reshist_addr1_drugtot_b, reshist_addr1_drgsale_b, reshist_addr1_drgposs_b, reshist_addr1_dui_b, reshist_addr1_mjsale_b, reshist_addr1_adi_edu_l_b, reshist_addr1_adi_edu_h_b, reshist_addr1_adi_work_c_b, reshist_addr1_adi_income_b, reshist_addr1_adi_in_dis_b, reshist_addr1_adi_home_v_b, reshist_addr1_adi_rent_b,reshist_addr1_adi_mortg_b, reshist_addr1_adi_home_o_b, reshist_addr1_adi_unemp_b, reshist_addr1_adi_pov_b, reshist_addr1_adi_b138_b, reshist_addr1_adi_sp_b, reshist_addr1_adi_ncar_b, reshist_addr1_adi_ntel_b, reshist_addr1_adi_nplumb_b, reshist_addr1_adi_crowd_b,
                                               neighborhood_crime_y_b, nsc_p_ss_mean_3_items_b, fes_p_ss_fc_pr_b, pmq_y_ss_mean_b, fes_y_ss_fc_pr_b, crpbi_y_ss_parent_b, srpf_y_ss_ses_b, srpf_y_ss_iiss_b, srpf_y_ss_dfs_b, Neurocog_BPPC1_GenAbility, Neurocog_BPPC2_ExecFunc,Neurocog_BPPC3_LearnMem)) # OK removed tfmri_sst_beh_performflag_b and tfmri_sst_all_beh_total_meanrt_inv_b
##
data.crop_y1 <- subset(data.merge_y1, select = c(subjectkey, race_ethnicity_y1, rel_family_id_y1, site_id_l_y1, reshist_addr1_walkindex_y1,reshist_addr1_popdensity_y1, reshist_addr1_d1a_y1, reshist_addr1_p1vlnt_y1, reshist_addr1_drugtot_y1, reshist_addr1_drgsale_y1, reshist_addr1_drgposs_y1, reshist_addr1_dui_y1, reshist_addr1_mjsale_y1, reshist_addr1_adi_edu_l_y1, reshist_addr1_adi_edu_h_y1, reshist_addr1_adi_work_c_y1, reshist_addr1_adi_income_y1, reshist_addr1_adi_in_dis_y1, reshist_addr1_adi_home_v_y1, reshist_addr1_adi_rent_y1, reshist_addr1_adi_mortg_y1, reshist_addr1_adi_home_o_y1, reshist_addr1_adi_unemp_y1, reshist_addr1_adi_pov_y1, reshist_addr1_adi_b138_y1, reshist_addr1_adi_sp_y1, reshist_addr1_adi_ncar_y1, reshist_addr1_adi_ntel_y1, reshist_addr1_adi_nplumb_y1, reshist_addr1_adi_crowd_y1,
                                                 neighborhood_crime_y_y1, nsc_p_ss_mean_3_items_y1, fes_p_ss_fc_pr_y1, pmq_y_ss_mean_y1, fes_y_ss_fc_pr_y1, crpbi_y_ss_parent_y1, srpf_y_ss_ses_y1, srpf_y_ss_iiss_y1, srpf_y_ss_dfs_y1, Income_y1, HighestEd_y1))
##
data.crop_y2 <- subset(data.merge_y2, select = c(subjectkey, site_id_l_y2, ra_scan_cl_mid_scan_lap_y2, ra_scan_cl_nbac_scan_lap_y2, ra_scan_cl_sst_scan_lap_y2, tfmri_nback_beh_performflag_y2, tfmri_mid_beh_performflag_y2, nihtbx_list_uncorrected_y2, nihtbx_picvocab_uncorrected_y2, nihtbx_reading_uncorrected_y2, tfmri_nb_all_beh_c2b_rate_y2, nihtbx_picture_uncorrected_y2, pea_ravlt_sd_trial_vi_tc_y2, pea_ravlt_ld_trial_vii_tc_y2, tfmri_nb_all_beh_c0b_rate_y2, nihtbx_cardsort_uncorrected_y2, nihtbx_flanker_uncorrected_y2, overall_dprime_y2, nihtbx_pattern_uncorrected_y2, mean_earnings_y2,
                                                 reshist_addr1_walkindex_y2, reshist_addr1_popdensity_y2, reshist_addr1_d1a_y2, reshist_addr1_p1vlnt_y2, reshist_addr1_drugtot_y2, reshist_addr1_drgsale_y2, reshist_addr1_drgposs_y2, reshist_addr1_dui_y2, reshist_addr1_mjsale_y2, reshist_addr1_adi_edu_l_y2, reshist_addr1_adi_edu_h_y2, reshist_addr1_adi_work_c_y2, reshist_addr1_adi_income_y2, reshist_addr1_adi_in_dis_y2, reshist_addr1_adi_home_v_y2, reshist_addr1_adi_rent_y2, reshist_addr1_adi_mortg_y2,reshist_addr1_adi_home_o_y2, reshist_addr1_adi_unemp_y2, reshist_addr1_adi_pov_y2, reshist_addr1_adi_b138_y2, reshist_addr1_adi_sp_y2, reshist_addr1_adi_ncar_y2, reshist_addr1_adi_ntel_y2, reshist_addr1_adi_nplumb_y2, reshist_addr1_adi_crowd_y2,
                                                 neighborhood_crime_y_y2, nsc_p_ss_mean_3_items_y2, fes_p_ss_fc_pr_y2, pmq_y_ss_mean_y2, fes_y_ss_fc_pr_y2, srpf_y_ss_ses_y2, srpf_y_ss_iiss_y2, srpf_y_ss_dfs_y2, Income_y2))

data.crop_mega <- data.merge_mega[ which(data.merge_mega$scrn_asd==0 & (data.merge_mega$scrn_epls!=1 | is.na(data.merge_mega$scrn_epls))), ]

data.crop_mega <- subset(data.crop_mega, select = c(subjectkey, rel_family_id_b, site_id_l_b, ra_scan_cl_mid_scan_lap_b, ra_scan_cl_nbac_scan_lap_b, ra_scan_cl_sst_scan_lap_b, interview_age, gender, tfmri_nback_beh_performflag_b, tfmri_mid_beh_performflag_b, nihtbx_list_uncorrected_b, nihtbx_picvocab_uncorrected_b, nihtbx_reading_uncorrected_b, tfmri_nb_all_beh_c2b_rate_b, pea_wiscv_tss_b, nihtbx_picture_uncorrected_b, pea_ravlt_sd_trial_vi_tc_b, pea_ravlt_ld_trial_vii_tc_b, tfmri_nb_all_beh_c0b_rate_b, nihtbx_cardsort_uncorrected_b, nihtbx_flanker_uncorrected_b, lmt_scr_efficiency_b, overall_dprime_b, nihtbx_pattern_uncorrected_b, mean_earnings_b, cash_choice_task_no3_b, 
                                                    reshist_addr1_pm25_b, reshist_addr1_no2_b, reshist_addr1_o3_2016_annavg_b, reshist_addr1_walkindex_b, reshist_addr1_proxrd_b, reshist_addr1_d1a_b, reshist_addr1_p1vlnt_b, reshist_addr1_drugtot_b, reshist_addr1_drgsale_b, reshist_addr1_drgposs_b, reshist_addr1_dui_b, reshist_addr1_mjsale_b, reshist_addr1_adi_edu_l_b, reshist_addr1_adi_edu_h_b, reshist_addr1_adi_work_c_b, reshist_addr1_adi_income_b, reshist_addr1_adi_in_dis_b, reshist_addr1_adi_home_v_b, reshist_addr1_adi_rent_b,reshist_addr1_adi_mortg_b, reshist_addr1_adi_home_o_b, reshist_addr1_adi_unemp_b, reshist_addr1_adi_pov_b, reshist_addr1_adi_b138_b, reshist_addr1_adi_sp_b, reshist_addr1_adi_ncar_b, reshist_addr1_adi_ntel_b, reshist_addr1_adi_nplumb_b, reshist_addr1_adi_crowd_b,
                                                    neighborhood_crime_y_b, nsc_p_ss_mean_3_items_b, fes_p_ss_fc_pr_b, pmq_y_ss_mean_b, fes_y_ss_fc_pr_b, crpbi_y_ss_parent_b, srpf_y_ss_ses_b, srpf_y_ss_iiss_b, srpf_y_ss_dfs_b,
                                                    reshist_addr1_walkindex_y1, reshist_addr1_d1a_y1, reshist_addr1_p1vlnt_y1, reshist_addr1_drugtot_y1, reshist_addr1_drgsale_y1, reshist_addr1_drgposs_y1, reshist_addr1_dui_y1, reshist_addr1_mjsale_y1, reshist_addr1_adi_edu_l_y1, reshist_addr1_adi_edu_h_y1, reshist_addr1_adi_work_c_y1, reshist_addr1_adi_income_y1, reshist_addr1_adi_in_dis_y1, reshist_addr1_adi_home_v_y1, reshist_addr1_adi_rent_y1, reshist_addr1_adi_mortg_y1, reshist_addr1_adi_home_o_y1, reshist_addr1_adi_unemp_y1, reshist_addr1_adi_pov_y1, reshist_addr1_adi_b138_y1, reshist_addr1_adi_sp_y1, reshist_addr1_adi_ncar_y1, reshist_addr1_adi_ntel_y1, reshist_addr1_adi_nplumb_y1, reshist_addr1_adi_crowd_y1,
                                                    neighborhood_crime_y_y1, nsc_p_ss_mean_3_items_y1, fes_p_ss_fc_pr_y1, pmq_y_ss_mean_y1, fes_y_ss_fc_pr_y1, crpbi_y_ss_parent_y1, srpf_y_ss_ses_y1, srpf_y_ss_iiss_y1, srpf_y_ss_dfs_y1, 
                                                    nihtbx_list_uncorrected_y2, nihtbx_picvocab_uncorrected_y2, nihtbx_reading_uncorrected_y2, tfmri_nb_all_beh_c2b_rate_y2, nihtbx_picture_uncorrected_y2, pea_ravlt_sd_trial_vi_tc_y2, pea_ravlt_ld_trial_vii_tc_y2, tfmri_nb_all_beh_c0b_rate_y2, nihtbx_cardsort_uncorrected_y2, nihtbx_flanker_uncorrected_y2, overall_dprime_y2, nihtbx_pattern_uncorrected_y2, mean_earnings_y2,
                                                    reshist_addr1_walkindex_y2, reshist_addr1_d1a_y2, reshist_addr1_p1vlnt_y2, reshist_addr1_drugtot_y2, reshist_addr1_drgsale_y2, reshist_addr1_drgposs_y2, reshist_addr1_dui_y2, reshist_addr1_mjsale_y2, reshist_addr1_adi_edu_l_y2, reshist_addr1_adi_edu_h_y2, reshist_addr1_adi_work_c_y2, reshist_addr1_adi_income_y2, reshist_addr1_adi_in_dis_y2, reshist_addr1_adi_home_v_y2, reshist_addr1_adi_rent_y2, reshist_addr1_adi_mortg_y2,reshist_addr1_adi_home_o_y2, reshist_addr1_adi_unemp_y2, reshist_addr1_adi_pov_y2, reshist_addr1_adi_b138_y2, reshist_addr1_adi_sp_y2, reshist_addr1_adi_ncar_y2, reshist_addr1_adi_ntel_y2, reshist_addr1_adi_nplumb_y2, reshist_addr1_adi_crowd_y2,
                                                    neighborhood_crime_y_y2, nsc_p_ss_mean_3_items_y2, fes_p_ss_fc_pr_y2, pmq_y_ss_mean_y2, fes_y_ss_fc_pr_y2, srpf_y_ss_ses_y2, srpf_y_ss_iiss_y2, srpf_y_ss_dfs_y2))

############################### Dataset for MULTIPLE REGRESSIONS ################################ 
dataforMR  <- Reduce(function(x,y) merge(x = x, y = y, by = "subjectkey", all.x = TRUE, all.y = TRUE), list(data.crop_b,data.crop_y1,data.crop_y2))

#### transform and clean data ####

#dummy code sex; Male = 1, Female = 0
dataforMR$gender <- ifelse(dataforMR$gender == "M", 1, 0)
colnames(dataforMR)[8] <- c("Male")

#dummy code race
dataforMR <- dataforMR %>% 
  dummy_cols(ignore_na = TRUE, select_columns = c("race_ethnicity_b"))
colnames(dataforMR)[163] = c("White")  
colnames(dataforMR)[164] = c("Black")
colnames(dataforMR)[165] = c("Hispanic")
colnames(dataforMR)[166] = c("Asian")
colnames(dataforMR)[167] = c("Other")

#air quality: (NO2 + PM2.5 indexes) -- more positive values indicate worse air quality
dataforMR$NeighAirQuality <- -(dataforMR$reshist_addr1_no2_b + dataforMR$reshist_addr1_pm25_b + reshist_addr1_o3_2016_annavg_b)
dataforMR$PM25            <- dataforMR$reshist_addr1_pm25_b
dataforMR$NO2            <- dataforMR$reshist_addr1_no2_b
dataforMR$O3            <- dataforMR$reshist_addr1_o3_2016_annavg_b
#youth reports of neigh safety + adult reports of neigh safety
dataforMR$NeighSafety     <- -(dataforMR$neighborhood_crime_y_b + dataforMR$nsc_p_ss_mean_3_items_b)
#log + 1 of adultviolentcrime + drugsale 
dataforMR$NeighCrime      <- (log1p(dataforMR$reshist_addr1_p1vlnt_b) + log1p(dataforMR$reshist_addr1_drgsale_b))
#youth reports of fam conflict + adult reports of fam conflict
dataforMR$FamilyConflict  <- (dataforMR$fes_y_ss_fc_pr_b + dataforMR$fes_p_ss_fc_pr_b)

#log +1 transform NCar, NTel, NPlumb, crowd, proximity to roads, population density
dataforMR$reshist_addr1_adi_ncar_b_log   <- log1p(dataforMR$reshist_addr1_adi_ncar_b)
dataforMR$reshist_addr1_adi_ntel_b_log   <- log1p(dataforMR$reshist_addr1_adi_ntel_b)
dataforMR$reshist_addr1_adi_nplumb_b_log <- log1p(dataforMR$reshist_addr1_adi_nplumb_b)
dataforMR$reshist_addr1_adi_crowd_b_log  <- log1p(dataforMR$reshist_addr1_adi_crowd_b)
dataforMR$reshist_addr1_proxrd_b_log     <- log1p(dataforMR$reshist_addr1_proxrd_b)
dataforMR$reshist_addr1_popdensity_b_log <- log1p(dataforMR$reshist_addr1_popdensity_b)

#change any cols with NaN to NA
dataforMR$reshist_addr1_adi_ncar_b_log[is.nan(dataforMR$reshist_addr1_adi_ncar_b_log)] <- NA
dataforMR$reshist_addr1_adi_ntel_b_log[is.nan(dataforMR$reshist_addr1_adi_ntel_b_log)] <- NA
dataforMR$reshist_addr1_adi_nplumb_b_log[is.nan(dataforMR$reshist_addr1_adi_nplumb_b_log)] <- NA
dataforMR$reshist_addr1_adi_crowd_b_log[is.nan(dataforMR$reshist_addr1_adi_crowd_b_log)] <- NA
dataforMR$reshist_addr1_proxrd_b_log[is.nan(dataforMR$reshist_addr1_proxrd_b_log)] <- NA
dataforMR$reshist_addr1_popdensity_b_log[is.nan(dataforMR$reshist_addr1_popdensity_b_log)] <- NA
dataforMR$NeighCrime[is.nan(dataforMR$NeighCrime)] <- NA


#set infinite values to NA
dataforMR$reshist_addr1_adi_ncar_b_log[is.infinite(dataforMR$reshist_addr1_adi_ncar_b_log)] <- NA
dataforMR$reshist_addr1_adi_ntel_b_log[is.infinite(dataforMR$reshist_addr1_adi_ntel_b_log)] <- NA
dataforMR$reshist_addr1_adi_nplumb_b_log[is.infinite(dataforMR$reshist_addr1_adi_nplumb_b_log)] <- NA
dataforMR$reshist_addr1_adi_crowd_b_log[is.infinite(dataforMR$reshist_addr1_adi_crowd_log)] <- NA
dataforMR$reshist_addr1_proxrd_b_log[is.infinite(dataforMR$reshist_addr1_proxrd_b_log)] <- NA
dataforMR$reshist_addr1_popdensity_b_log[is.infinite(dataforMR$reshist_addr1_popdensity_b_log)] <- NA
dataforMR$NeighCrime[is.infinite(dataforMR$NeighCrime)] <- NA

write.csv(dataforMR,file = "dataforMR.csv")  # this csv is used in the regressions script