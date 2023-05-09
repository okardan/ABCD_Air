% "Neighborhood air pollution is negatively associated with neurocognitive maturation in early adolescence"
% Omid Kardan, Chacriya Sereeyothin, Kathryn E. Schertz, Mike Angstadt, Alexander S. Weigard, Marc G. Berman, Monica D. Rosenberg

% script by Omid Kardan 
% contact omidk@med.umich.edu
% Produces nihcpm_scores.csv which is used in the regressions R script Kardan_et_al_Regressions_and_violins.R

% Requires pre-processed Shen-268 Parcellated fMRI timeseries and ccCPM by ABCD sites (download 
% trained ccCPMs for each left-out site from https://nda.nih.gov/study.html?id=1849)


%% Read baselineY1 and y2followup clean shen timeseries
clear all
folder_base = '~\fMRI Data\Timeseries\baselineYear1Arm1\'; % folder containing pre-processed parcel timeseries for n-back task fMRI
folder_2y = '~\fMRI Data\Timeseries\2YearFollowUpYArm1\'; % folder containing pre-processed parcel timeseries for n-back task fMRI

abcd_shen = cell(1166,1);
T = readtable('~\fMRI Data\ABCD_nback_both.csv'); % file containing Frame Displacement (FD) for n-back task fMRI

for k =1:1166
subid = T.subjectid{k};
ptseries_base_nbk1 = importdata([folder_base,subid,'_baselineYear1Arm1_run-01_clean_ts.txt']);
ptseries_base_nbk2 = importdata([folder_base,subid,'_baselineYear1Arm1_run-02_clean_ts.txt']);

ptseries_2y_nbk1 = importdata([folder_2y,subid,'_2YearFollowUpYArm1_nback_run-01_clean_ts.txt']);
ptseries_2y_nbk2 = importdata([folder_2y,subid,'_2YearFollowUpYArm1_nback_run-02_clean_ts.txt']);

abcd_shen{k}.subid = subid;
abcd_shen{k}.pt.bl_nbk1 = [ptseries_base_nbk1];
abcd_shen{k}.pt.bl_nbk2 = [ptseries_base_nbk2];
abcd_shen{k}.pt.yr2_nbk1 = [ ptseries_2y_nbk1];
abcd_shen{k}.pt.yr2_nbk2 = [ptseries_2y_nbk2];
abcd_shen{k}.meanFD = [T.meanFD_bl(k)  T.meanFD_2yr(k)];
abcd_shen{k}.site = [T.site_id_l_bl{k}  T.site_id_l_2yr{k}];
abcd_shen{k}.nbkacc = [ [T.acc0b_bl(k); T.acc2b_bl(k)]  [T.acc0b_2yr(k);  T.acc2b_2yr(k)] ];
k
end
%%  make ccCPM strength (nih-toolbox cpm from Kardan et al, 2022, PLOS Biology) scores

nih_cpm =cell(22,1);
for ss=1:22
    sm = ss;
    if ss ==17 | ss==19 % Sites 17 and 19 were not included in Kardan et al 2022 so using any of the site ccCPMs is ok for them
        sm = 2;
    end
    clear mdl
    mdl = NaN(268,268);
    try
   load(['~\BrainWeights\CPM_cogComp_masks\trndMdlcpm01_',num2str(sm),'.mat']); %(download 
   % trained ccCPMs for each left-out site from https://nda.nih.gov/study.html?id=1849)
    catch
    end
    nih_cpm{ss} = mdl;
end

sites =[]; subs =[]; nihcpm_score_bl=[]; nihcpm_score_yr2=[];
for k = 1:1166
    subs = [subs; string(abcd_shen{k}.subid)];
    sites = [sites; str2num(abcd_shen{k}.site(5:6))];
    
    smod = nih_cpm{ str2num(abcd_shen{k}.site(5:6)) };
    corr_mat1 = atanh(corr(abcd_shen{k, 1}.pt.bl_nbk1));
    corr_mat2 = atanh(corr(abcd_shen{k, 1}.pt.bl_nbk2));
    
    corr_mat3 = atanh(corr(abcd_shen{k, 1}.pt.yr2_nbk1));
    corr_mat4 = atanh(corr(abcd_shen{k, 1}.pt.yr2_nbk2));
    
    nihcpm_score_bl = [nihcpm_score_bl; nanmean(nanmean(corr_mat1.*smod)) + nanmean(nanmean(corr_mat2.*smod))];
    nihcpm_score_yr2 = [nihcpm_score_yr2; nanmean(nanmean(corr_mat3.*smod)) + nanmean(nanmean(corr_mat4.*smod))];
    k
end
T2 = table(subs,sites,nihcpm_score_bl,nihcpm_score_yr2);
writetable(T2,'nihcpm_scores.csv');  % used in the regressions R script Kardan_et_al_Regressions_and_violins.R
    
    
