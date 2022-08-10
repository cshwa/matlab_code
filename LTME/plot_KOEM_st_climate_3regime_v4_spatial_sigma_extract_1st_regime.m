close all; clear; clc;
sigsig=[2;3;];
%%
for ixx = 1:2
    clearvars -except ixx sigsig
    
sig = sigsig(ixx);

cd F:\ROMS\roms_tools\Run\
start

cd D:\장기생태\Dynamic\KOEM
% koem data
% load koem_timeseires_monthly_3sig.mat % estimate 3sigma from spatially all pt.
% load koem_timeseires_monthly_gy_only_2sig.mat % estimate 2sigma from spatially gy pt. from 'plot_KOEM_st_for_data_process_step5_fix_yearly_monthly_gy_3sig.m'
% load koem_timeseires_monthly_gy_only_3sig.mat % estimate 3sigma from spatially gy pt. from 'plot_KOEM_st_for_data_process_step5_fix_yearly_monthly_gy_3sig.m'
load(['koem_timeseires_monthly_gy_only_',num2str(sig),'sig.mat']);


% spatial region index
sp_gy = [4,22,28,29,30,32,33,34,35];
sp_s_gy = [38, 39, 42, 41, 44];
sp_gm = [6,9,10];% gamak bay
sp_e_gy = [48, 49];
sp_jj = [36, 37, 45, 46]; %jinju bay

% make eom_d
k=0
for i = 2001:2001
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end
t_tick_pre=sum(eom_d,2);

for i = 1:size(eom_d,1)
    for j = 1:size(eom_d,2)
        eom_d_each(i,j) = sum(eom_d(i,1:j));
    end
end


% obs_full_tdx = 49:60;
% full time climatological mean
% 1997~2018 : 22yr
% 1st regime : 1997 ~ 2003.12 (from koem nh4 shift)
% 2nd regime : 2004.01 ~ 2009.12 (from koem no3 shift)
% 3rd regime : 2010.01 ~ 2018.12 (from koem no3 shift)

% yr10_filter(1:168) = 1;
% yr10_filter(168+1:size(regime_no3,2)) = 1;
find_t_ind = {'2003-12'};
find_date=find(strcmp(ref_date,find_t_ind) == 1);
yr10_filter(1:find_date) = 1;
yr10_filter(find_date+1:length(ref_date)) = NaN;

for i = 1:65
regime_no3(i,:) = regime_no3(i,:) .* yr10_filter; 
regime_nh4(i,:) = regime_nh4(i,:) .* yr10_filter; 
regime_do(i,:) = regime_do(i,:) .* yr10_filter; 
regime_chl(i,:) = regime_chl(i,:) .* yr10_filter; 
regime_temp(i,:) = regime_temp(i,:) .* yr10_filter; 
regime_salt(i,:) = regime_salt(i,:) .* yr10_filter; 
regime_no3_b(i,:) = regime_no3_b(i,:) .* yr10_filter; 
regime_nh4_b(i,:) = regime_nh4_b(i,:) .* yr10_filter; 
regime_do_b(i,:) = regime_do_b(i,:) .* yr10_filter; 
regime_chl_b(i,:) = regime_chl_b(i,:) .* yr10_filter; 
regime_temp_b(i,:) = regime_temp_b(i,:) .* yr10_filter; 
regime_salt_b(i,:) = regime_salt_b(i,:) .* yr10_filter; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  spatial mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % spatial region index
sp_gy = [4,22,28,29,30,32,33,34,35];
sp_s_gy = [38, 39, 42, 41, 44];
sp_gm = [6,9,10];% gamak bay
sp_e_gy = [48, 49];
sp_jj = [36, 37, 45, 46]; %jinju bay


%% obs 
%% full time obs spatial mean

%sp_gy
    obs_gy_no3_ft=squeeze(regime_no3(sp_gy,:));
    obs_gy_no3_b_ft=squeeze(regime_no3_b(sp_gy,:));
    obs_gy_nh4_ft=squeeze(regime_nh4(sp_gy,:));
    obs_gy_nh4_b_ft=squeeze(regime_nh4_b(sp_gy,:));
    obs_gy_chl_ft=squeeze(regime_chl(sp_gy,:));
    obs_gy_chl_b_ft=squeeze(regime_chl_b(sp_gy,:));
    obs_gy_temp_ft=squeeze(regime_temp(sp_gy,:));
    obs_gy_temp_b_ft=squeeze(regime_temp_b(sp_gy,:));
    obs_gy_salt_ft=squeeze(regime_salt(sp_gy,:));
    obs_gy_salt_b_ft=squeeze(regime_salt_b(sp_gy,:));
    obs_gy_do_ft=squeeze(regime_do(sp_gy,:));
    obs_gy_do_b_ft=squeeze(regime_do_b(sp_gy,:));  


%sp_s_gy
    obs_sgy_no3_ft=squeeze(regime_no3(sp_s_gy,:));
    obs_sgy_no3_b_ft=squeeze(regime_no3_b(sp_s_gy,:));
    obs_sgy_nh4_ft=squeeze(regime_nh4(sp_s_gy,:));
    obs_sgy_nh4_b_ft=squeeze(regime_nh4_b(sp_s_gy,:));
    obs_sgy_chl_ft=squeeze(regime_chl(sp_s_gy,:));
    obs_sgy_chl_b_ft=squeeze(regime_chl_b(sp_s_gy,:));
    obs_sgy_temp_ft=squeeze(regime_temp(sp_s_gy,:));
    obs_sgy_temp_b_ft=squeeze(regime_temp_b(sp_s_gy,:));
    obs_sgy_salt_ft=squeeze(regime_salt(sp_s_gy,:));
    obs_sgy_salt_b_ft=squeeze(regime_salt_b(sp_s_gy,:));
    obs_sgy_do_ft=squeeze(regime_do(sp_s_gy,:));
    obs_sgy_do_b_ft=squeeze(regime_do_b(sp_s_gy,:));  


% sp_e_gy
    obs_egy_no3_ft=squeeze(regime_no3(sp_e_gy,:));
    obs_egy_no3_b_ft=squeeze(regime_no3_b(sp_e_gy,:));
    obs_egy_nh4_ft=squeeze(regime_nh4(sp_e_gy,:));
    obs_egy_nh4_b_ft=squeeze(regime_nh4_b(sp_e_gy,:));
    obs_egy_chl_ft=squeeze(regime_chl(sp_e_gy,:));
    obs_egy_chl_b_ft=squeeze(regime_chl_b(sp_e_gy,:));
    obs_egy_temp_ft=squeeze(regime_temp(sp_e_gy,:));
    obs_egy_temp_b_ft=squeeze(regime_temp_b(sp_e_gy,:));
    obs_egy_salt_ft=squeeze(regime_salt(sp_e_gy,:));
    obs_egy_salt_b_ft=squeeze(regime_salt_b(sp_e_gy,:));
    obs_egy_do_ft=squeeze(regime_do(sp_e_gy,:));
    obs_egy_do_b_ft=squeeze(regime_do_b(sp_e_gy,:));  

%     sp_jj
    obs_jj_no3_ft=squeeze(regime_no3(sp_jj,:));
    obs_jj_no3_b_ft=squeeze(regime_no3_b(sp_jj,:));
    obs_jj_nh4_ft=squeeze(regime_nh4(sp_jj,:));
    obs_jj_nh4_b_ft=squeeze(regime_nh4_b(sp_jj,:));
    obs_jj_chl_ft=squeeze(regime_chl(sp_jj,:));
    obs_jj_chl_b_ft=squeeze(regime_chl_b(sp_jj,:));
    obs_jj_temp_ft=squeeze(regime_temp(sp_jj,:));
    obs_jj_temp_b_ft=squeeze(regime_temp_b(sp_jj,:));
    obs_jj_salt_ft=squeeze(regime_salt(sp_jj,:));
    obs_jj_salt_b_ft=squeeze(regime_salt_b(sp_jj,:));
    obs_jj_do_ft=squeeze(regime_do(sp_jj,:));
    obs_jj_do_b_ft=squeeze(regime_do_b(sp_jj,:));  

for i = 1:12
     if i==1
        %% time & spatial std 
        clearvars len_temp
        len_temp=size(obs_gy_no3_ft(:,i:12:end),1)*size(obs_gy_no3_ft(:,i:12:end),2);
        obs_std_gy_no3_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_no3_ft(:,i:12:end),1,len_temp));
        obs_std_gy_no3_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_no3_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_nh4_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_nh4_ft(:,i:12:end),1,len_temp));
        obs_std_gy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_nh4_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_chl_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_chl_ft(:,i:12:end),1,len_temp));
        obs_std_gy_chl_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_chl_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_temp_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_temp_ft(:,i:12:end),1,len_temp));
        obs_std_gy_temp_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_temp_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_salt_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_salt_ft(:,i:12:end),1,len_temp));
        obs_std_gy_salt_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_salt_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_do_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_do_ft(:,i:12:end),1,len_temp));
        obs_std_gy_do_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_gy_do_b_ft(:,i:12:end),1,len_temp));
        clearvars len_temp
        len_temp=size(obs_sgy_no3_ft(:,i:12:end),1)*size(obs_sgy_no3_ft(:,i:12:end),2);
        obs_std_sgy_no3_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_no3_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_no3_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_no3_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_nh4_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_nh4_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_nh4_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_chl_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_chl_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_chl_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_chl_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_temp_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_temp_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_temp_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_temp_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_salt_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_salt_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_salt_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_salt_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_do_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_do_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_do_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_do_b_ft(:,i:12:end),1,len_temp));
        clearvars len_temp
        len_temp=size(obs_egy_no3_ft(:,i:12:end),1)*size(obs_egy_no3_ft(:,i:12:end),2);
        obs_std_egy_no3_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_no3_ft(:,i:12:end),1,len_temp));
        obs_std_egy_no3_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_no3_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_nh4_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_nh4_ft(:,i:12:end),1,len_temp));
        obs_std_egy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_nh4_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_chl_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_chl_ft(:,i:12:end),1,len_temp));
        obs_std_egy_chl_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_chl_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_temp_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_temp_ft(:,i:12:end),1,len_temp));
        obs_std_egy_temp_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_temp_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_salt_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_salt_ft(:,i:12:end),1,len_temp));
        obs_std_egy_salt_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_salt_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_do_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_do_ft(:,i:12:end),1,len_temp));
        obs_std_egy_do_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_egy_do_b_ft(:,i:12:end),1,len_temp));
        clearvars len_temp
        len_temp=size(obs_jj_no3_ft(:,i:12:end),1)*size(obs_jj_no3_ft(:,i:12:end),2);
        obs_std_jj_no3_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_no3_ft(:,i:12:end),1,len_temp));
        obs_std_jj_no3_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_no3_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_nh4_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_nh4_ft(:,i:12:end),1,len_temp));
        obs_std_jj_nh4_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_nh4_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_chl_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_chl_ft(:,i:12:end),1,len_temp));
        obs_std_jj_chl_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_chl_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_temp_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_temp_ft(:,i:12:end),1,len_temp));
        obs_std_jj_temp_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_temp_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_salt_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_salt_ft(:,i:12:end),1,len_temp));
        obs_std_jj_salt_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_salt_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_do_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_do_ft(:,i:12:end),1,len_temp));
        obs_std_jj_do_b_ft(1:eom_d_each(1,i))=nanstd(reshape(obs_jj_do_b_ft(:,i:12:end),1,len_temp));

       %% time std
        spm_obs_std_gy_no3_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_no3_ft(:,i:12:end),1));
        spm_obs_std_gy_no3_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_no3_b_ft(:,i:12:end),1));
        spm_obs_std_gy_nh4_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_nh4_ft(:,i:12:end),1));
        spm_obs_std_gy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_nh4_b_ft(:,i:12:end),1));
        spm_obs_std_gy_chl_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_chl_ft(:,i:12:end),1));
        spm_obs_std_gy_chl_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_chl_b_ft(:,i:12:end),1));
        spm_obs_std_gy_temp_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_temp_ft(:,i:12:end),1));
        spm_obs_std_gy_temp_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_temp_b_ft(:,i:12:end),1));
        spm_obs_std_gy_salt_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_salt_ft(:,i:12:end),1));
        spm_obs_std_gy_salt_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_salt_b_ft(:,i:12:end),1));
        spm_obs_std_gy_do_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_do_ft(:,i:12:end),1));
        spm_obs_std_gy_do_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_do_b_ft(:,i:12:end),1));

        spm_obs_std_sgy_no3_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_no3_ft(:,i:12:end),1));
        spm_obs_std_sgy_no3_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_no3_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_nh4_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_nh4_ft(:,i:12:end),1));
        spm_obs_std_sgy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_nh4_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_chl_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_chl_ft(:,i:12:end),1));
        spm_obs_std_sgy_chl_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_chl_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_temp_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_temp_ft(:,i:12:end),1));
        spm_obs_std_sgy_temp_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_temp_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_salt_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_salt_ft(:,i:12:end),1));
        spm_obs_std_sgy_salt_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_salt_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_do_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_do_ft(:,i:12:end),1));
        spm_obs_std_sgy_do_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_do_b_ft(:,i:12:end),1));

        spm_obs_std_egy_no3_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_no3_ft(:,i:12:end),1));
        spm_obs_std_egy_no3_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_no3_b_ft(:,i:12:end),1));
        spm_obs_std_egy_nh4_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_nh4_ft(:,i:12:end),1));
        spm_obs_std_egy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_nh4_b_ft(:,i:12:end),1));
        spm_obs_std_egy_chl_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_chl_ft(:,i:12:end),1));
        spm_obs_std_egy_chl_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_chl_b_ft(:,i:12:end),1));
        spm_obs_std_egy_temp_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_temp_ft(:,i:12:end),1));
        spm_obs_std_egy_temp_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_temp_b_ft(:,i:12:end),1));
        spm_obs_std_egy_salt_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_salt_ft(:,i:12:end),1));
        spm_obs_std_egy_salt_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_salt_b_ft(:,i:12:end),1));
        spm_obs_std_egy_do_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_do_ft(:,i:12:end),1));
        spm_obs_std_egy_do_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_do_b_ft(:,i:12:end),1));

        spm_obs_std_jj_no3_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_no3_ft(:,i:12:end),1));
        spm_obs_std_jj_no3_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_no3_b_ft(:,i:12:end),1));
        spm_obs_std_jj_nh4_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_nh4_ft(:,i:12:end),1));
        spm_obs_std_jj_nh4_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_nh4_b_ft(:,i:12:end),1));
        spm_obs_std_jj_chl_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_chl_ft(:,i:12:end),1));
        spm_obs_std_jj_chl_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_chl_b_ft(:,i:12:end),1));
        spm_obs_std_jj_temp_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_temp_ft(:,i:12:end),1));
        spm_obs_std_jj_temp_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_temp_b_ft(:,i:12:end),1));
        spm_obs_std_jj_salt_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_salt_ft(:,i:12:end),1));
        spm_obs_std_jj_salt_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_salt_b_ft(:,i:12:end),1));
        spm_obs_std_jj_do_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_do_ft(:,i:12:end),1));
        spm_obs_std_jj_do_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_do_b_ft(:,i:12:end),1));

       %% spatial std
        tm_obs_std_gy_no3_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_no3_ft(:,i:12:end),2));
        tm_obs_std_gy_no3_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_no3_b_ft(:,i:12:end),2));
        tm_obs_std_gy_nh4_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_nh4_ft(:,i:12:end),2));
        tm_obs_std_gy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_nh4_b_ft(:,i:12:end),2));
        tm_obs_std_gy_chl_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_chl_ft(:,i:12:end),2));
        tm_obs_std_gy_chl_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_chl_b_ft(:,i:12:end),2));
        tm_obs_std_gy_temp_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_temp_ft(:,i:12:end),2));
        tm_obs_std_gy_temp_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_temp_b_ft(:,i:12:end),2));
        tm_obs_std_gy_salt_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_salt_ft(:,i:12:end),2));
        tm_obs_std_gy_salt_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_salt_b_ft(:,i:12:end),2));
        tm_obs_std_gy_do_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_do_ft(:,i:12:end),2));
        tm_obs_std_gy_do_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_do_b_ft(:,i:12:end),2));

        tm_obs_std_sgy_no3_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_no3_ft(:,i:12:end),2));
        tm_obs_std_sgy_no3_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_no3_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_nh4_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_nh4_ft(:,i:12:end),2));
        tm_obs_std_sgy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_nh4_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_chl_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_chl_ft(:,i:12:end),2));
        tm_obs_std_sgy_chl_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_chl_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_temp_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_temp_ft(:,i:12:end),2));
        tm_obs_std_sgy_temp_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_temp_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_salt_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_salt_ft(:,i:12:end),2));
        tm_obs_std_sgy_salt_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_salt_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_do_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_do_ft(:,i:12:end),2));
        tm_obs_std_sgy_do_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_do_b_ft(:,i:12:end),2));

        tm_obs_std_egy_no3_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_no3_ft(:,i:12:end),2));
        tm_obs_std_egy_no3_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_no3_b_ft(:,i:12:end),2));
        tm_obs_std_egy_nh4_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_nh4_ft(:,i:12:end),2));
        tm_obs_std_egy_nh4_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_nh4_b_ft(:,i:12:end),2));
        tm_obs_std_egy_chl_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_chl_ft(:,i:12:end),2));
        tm_obs_std_egy_chl_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_chl_b_ft(:,i:12:end),2));
        tm_obs_std_egy_temp_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_temp_ft(:,i:12:end),2));
        tm_obs_std_egy_temp_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_temp_b_ft(:,i:12:end),2));
        tm_obs_std_egy_salt_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_salt_ft(:,i:12:end),2));
        tm_obs_std_egy_salt_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_salt_b_ft(:,i:12:end),2));
        tm_obs_std_egy_do_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_do_ft(:,i:12:end),2));
        tm_obs_std_egy_do_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_do_b_ft(:,i:12:end),2));

        tm_obs_std_jj_no3_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_no3_ft(:,i:12:end),2));
        tm_obs_std_jj_no3_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_no3_b_ft(:,i:12:end),2));
        tm_obs_std_jj_nh4_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_nh4_ft(:,i:12:end),2));
        tm_obs_std_jj_nh4_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_nh4_b_ft(:,i:12:end),2));
        tm_obs_std_jj_chl_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_chl_ft(:,i:12:end),2));
        tm_obs_std_jj_chl_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_chl_b_ft(:,i:12:end),2));
        tm_obs_std_jj_temp_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_temp_ft(:,i:12:end),2));
        tm_obs_std_jj_temp_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_temp_b_ft(:,i:12:end),2));
        tm_obs_std_jj_salt_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_salt_ft(:,i:12:end),2));
        tm_obs_std_jj_salt_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_salt_b_ft(:,i:12:end),2));
        tm_obs_std_jj_do_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_do_ft(:,i:12:end),2));
        tm_obs_std_jj_do_b_ft(1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_do_b_ft(:,i:12:end),2));

        
        %% mean
        obm_gy_no3_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_no3_ft(:,i:12:end)));
        obm_gy_no3_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_no3_b_ft(:,i:12:end)));
        obm_gy_nh4_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_nh4_ft(:,i:12:end)));
        obm_gy_nh4_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_nh4_b_ft(:,i:12:end)));
        obm_gy_chl_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_chl_ft(:,i:12:end)));
        obm_gy_chl_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_chl_b_ft(:,i:12:end)));
        obm_gy_temp_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_temp_ft(:,i:12:end)));
        obm_gy_temp_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_temp_b_ft(:,i:12:end)));
        obm_gy_salt_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_salt_ft(:,i:12:end)));
        obm_gy_salt_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_salt_b_ft(:,i:12:end)));
        obm_gy_do_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_do_ft(:,i:12:end)));
        obm_gy_do_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_do_b_ft(:,i:12:end)));

        obm_sgy_no3_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_no3_ft(:,i:12:end)));
        obm_sgy_no3_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_no3_b_ft(:,i:12:end)));
        obm_sgy_nh4_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_nh4_ft(:,i:12:end)));
        obm_sgy_nh4_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_nh4_b_ft(:,i:12:end)));
        obm_sgy_chl_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_chl_ft(:,i:12:end)));
        obm_sgy_chl_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_chl_b_ft(:,i:12:end)));
        obm_sgy_temp_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_temp_ft(:,i:12:end)));
        obm_sgy_temp_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_temp_b_ft(:,i:12:end)));
        obm_sgy_salt_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_salt_ft(:,i:12:end)));
        obm_sgy_salt_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_salt_b_ft(:,i:12:end)));
        obm_sgy_do_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_do_ft(:,i:12:end)));
        obm_sgy_do_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_do_b_ft(:,i:12:end)));

        obm_egy_no3_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_no3_ft(:,i:12:end)));
        obm_egy_no3_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_no3_b_ft(:,i:12:end)));
        obm_egy_nh4_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_nh4_ft(:,i:12:end)));
        obm_egy_nh4_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_nh4_b_ft(:,i:12:end)));
        obm_egy_chl_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_chl_ft(:,i:12:end)));
        obm_egy_chl_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_chl_b_ft(:,i:12:end)));
        obm_egy_temp_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_temp_ft(:,i:12:end)));
        obm_egy_temp_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_temp_b_ft(:,i:12:end)));
        obm_egy_salt_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_salt_ft(:,i:12:end)));
        obm_egy_salt_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_salt_b_ft(:,i:12:end)));
        obm_egy_do_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_do_ft(:,i:12:end)));
        obm_egy_do_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_do_b_ft(:,i:12:end)));

        obm_jj_no3_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_no3_ft(:,i:12:end)));
        obm_jj_no3_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_no3_b_ft(:,i:12:end)));
        obm_jj_nh4_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_nh4_ft(:,i:12:end)));
        obm_jj_nh4_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_nh4_b_ft(:,i:12:end)));
        obm_jj_chl_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_chl_ft(:,i:12:end)));
        obm_jj_chl_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_chl_b_ft(:,i:12:end)));
        obm_jj_temp_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_temp_ft(:,i:12:end)));
        obm_jj_temp_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_temp_b_ft(:,i:12:end)));
        obm_jj_salt_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_salt_ft(:,i:12:end)));
        obm_jj_salt_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_salt_b_ft(:,i:12:end)));
        obm_jj_do_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_do_ft(:,i:12:end)));
        obm_jj_do_b_ft(1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_do_b_ft(:,i:12:end)));

     else
        %% time & spatial std
        clearvars len_temp
        len_temp=size(obs_gy_no3_ft(:,i:12:end),1)*size(obs_gy_no3_ft(:,i:12:end),2);
        obs_std_gy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_no3_ft(:,i:12:end),1,len_temp));
        obs_std_gy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_no3_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_nh4_ft(:,i:12:end),1,len_temp));
        obs_std_gy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_nh4_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_chl_ft(:,i:12:end),1,len_temp));
        obs_std_gy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_chl_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_temp_ft(:,i:12:end),1,len_temp));
        obs_std_gy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_temp_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_salt_ft(:,i:12:end),1,len_temp));
        obs_std_gy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_salt_b_ft(:,i:12:end),1,len_temp));
        obs_std_gy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_do_ft(:,i:12:end),1,len_temp));
        obs_std_gy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_gy_do_b_ft(:,i:12:end),1,len_temp));
        clearvars len_temp
        len_temp=size(obs_sgy_no3_ft(:,i:12:end),1)*size(obs_sgy_no3_ft(:,i:12:end),2);
        obs_std_sgy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_no3_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_no3_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_nh4_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_nh4_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_chl_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_chl_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_temp_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_temp_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_salt_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_salt_b_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_do_ft(:,i:12:end),1,len_temp));
        obs_std_sgy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_sgy_do_b_ft(:,i:12:end),1,len_temp));
        clearvars len_temp
        len_temp=size(obs_egy_no3_ft(:,i:12:end),1)*size(obs_egy_no3_ft(:,i:12:end),2);
        obs_std_egy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_no3_ft(:,i:12:end),1,len_temp));
        obs_std_egy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_no3_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_nh4_ft(:,i:12:end),1,len_temp));
        obs_std_egy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_nh4_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_chl_ft(:,i:12:end),1,len_temp));
        obs_std_egy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_chl_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_temp_ft(:,i:12:end),1,len_temp));
        obs_std_egy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_temp_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_salt_ft(:,i:12:end),1,len_temp));
        obs_std_egy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_salt_b_ft(:,i:12:end),1,len_temp));
        obs_std_egy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_do_ft(:,i:12:end),1,len_temp));
        obs_std_egy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_egy_do_b_ft(:,i:12:end),1,len_temp));
        clearvars len_temp
        len_temp=size(obs_jj_no3_ft(:,i:12:end),1)*size(obs_jj_no3_ft(:,i:12:end),2);
        obs_std_jj_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_no3_ft(:,i:12:end),1,len_temp));
        obs_std_jj_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_no3_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_nh4_ft(:,i:12:end),1,len_temp));
        obs_std_jj_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_nh4_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_chl_ft(:,i:12:end),1,len_temp));
        obs_std_jj_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_chl_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_temp_ft(:,i:12:end),1,len_temp));
        obs_std_jj_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_temp_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_salt_ft(:,i:12:end),1,len_temp));
        obs_std_jj_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_salt_b_ft(:,i:12:end),1,len_temp));
        obs_std_jj_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_do_ft(:,i:12:end),1,len_temp));
        obs_std_jj_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(reshape(obs_jj_do_b_ft(:,i:12:end),1,len_temp));
        
         %% time std
        spm_obs_std_gy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_no3_ft(:,i:12:end),1));
        spm_obs_std_gy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_no3_b_ft(:,i:12:end),1));
        spm_obs_std_gy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_nh4_ft(:,i:12:end),1));
        spm_obs_std_gy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_nh4_b_ft(:,i:12:end),1));
        spm_obs_std_gy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_chl_ft(:,i:12:end),1));
        spm_obs_std_gy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_chl_b_ft(:,i:12:end),1));
        spm_obs_std_gy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_temp_ft(:,i:12:end),1));
        spm_obs_std_gy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_temp_b_ft(:,i:12:end),1));
        spm_obs_std_gy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_salt_ft(:,i:12:end),1));
        spm_obs_std_gy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_salt_b_ft(:,i:12:end),1));
        spm_obs_std_gy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_do_ft(:,i:12:end),1));
        spm_obs_std_gy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_do_b_ft(:,i:12:end),1));

        spm_obs_std_sgy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_no3_ft(:,i:12:end),1));
        spm_obs_std_sgy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_no3_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_nh4_ft(:,i:12:end),1));
        spm_obs_std_sgy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_nh4_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_chl_ft(:,i:12:end),1));
        spm_obs_std_sgy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_chl_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_temp_ft(:,i:12:end),1));
        spm_obs_std_sgy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_temp_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_salt_ft(:,i:12:end),1));
        spm_obs_std_sgy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_salt_b_ft(:,i:12:end),1));
        spm_obs_std_sgy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_do_ft(:,i:12:end),1));
        spm_obs_std_sgy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_do_b_ft(:,i:12:end),1));

        spm_obs_std_egy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_no3_ft(:,i:12:end),1));
        spm_obs_std_egy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_no3_b_ft(:,i:12:end),1));
        spm_obs_std_egy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_nh4_ft(:,i:12:end),1));
        spm_obs_std_egy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_nh4_b_ft(:,i:12:end),1));
        spm_obs_std_egy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_chl_ft(:,i:12:end),1));
        spm_obs_std_egy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_chl_b_ft(:,i:12:end),1));
        spm_obs_std_egy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_temp_ft(:,i:12:end),1));
        spm_obs_std_egy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_temp_b_ft(:,i:12:end),1));
        spm_obs_std_egy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_salt_ft(:,i:12:end),1));
        spm_obs_std_egy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_salt_b_ft(:,i:12:end),1));
        spm_obs_std_egy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_do_ft(:,i:12:end),1));
        spm_obs_std_egy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_do_b_ft(:,i:12:end),1));

        spm_obs_std_jj_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_no3_ft(:,i:12:end),1));
        spm_obs_std_jj_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_no3_b_ft(:,i:12:end),1));
        spm_obs_std_jj_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_nh4_ft(:,i:12:end),1));
        spm_obs_std_jj_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_nh4_b_ft(:,i:12:end),1));
        spm_obs_std_jj_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_chl_ft(:,i:12:end),1));
        spm_obs_std_jj_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_chl_b_ft(:,i:12:end),1));
        spm_obs_std_jj_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_temp_ft(:,i:12:end),1));
        spm_obs_std_jj_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_temp_b_ft(:,i:12:end),1));
        spm_obs_std_jj_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_salt_ft(:,i:12:end),1));
        spm_obs_std_jj_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_salt_b_ft(:,i:12:end),1));
        spm_obs_std_jj_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_do_ft(:,i:12:end),1));
        spm_obs_std_jj_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_do_b_ft(:,i:12:end),1));

         %% spatial std
        tm_obs_std_gy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_no3_ft(:,i:12:end),2));
        tm_obs_std_gy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_no3_b_ft(:,i:12:end),2));
        tm_obs_std_gy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_nh4_ft(:,i:12:end),2));
        tm_obs_std_gy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_nh4_b_ft(:,i:12:end),2));
        tm_obs_std_gy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_chl_ft(:,i:12:end),2));
        tm_obs_std_gy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_chl_b_ft(:,i:12:end),2));
        tm_obs_std_gy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_temp_ft(:,i:12:end),2));
        tm_obs_std_gy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_temp_b_ft(:,i:12:end),2));
        tm_obs_std_gy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_salt_ft(:,i:12:end),2));
        tm_obs_std_gy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_salt_b_ft(:,i:12:end),2));
        tm_obs_std_gy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_do_ft(:,i:12:end),2));
        tm_obs_std_gy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_gy_do_b_ft(:,i:12:end),2));

        tm_obs_std_sgy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_no3_ft(:,i:12:end),2));
        tm_obs_std_sgy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_no3_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_nh4_ft(:,i:12:end),2));
        tm_obs_std_sgy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_nh4_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_chl_ft(:,i:12:end),2));
        tm_obs_std_sgy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_chl_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_temp_ft(:,i:12:end),2));
        tm_obs_std_sgy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_temp_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_salt_ft(:,i:12:end),2));
        tm_obs_std_sgy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_salt_b_ft(:,i:12:end),2));
        tm_obs_std_sgy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_do_ft(:,i:12:end),2));
        tm_obs_std_sgy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_sgy_do_b_ft(:,i:12:end),2));

        tm_obs_std_egy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_no3_ft(:,i:12:end),2));
        tm_obs_std_egy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_no3_b_ft(:,i:12:end),2));
        tm_obs_std_egy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_nh4_ft(:,i:12:end),2));
        tm_obs_std_egy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_nh4_b_ft(:,i:12:end),2));
        tm_obs_std_egy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_chl_ft(:,i:12:end),2));
        tm_obs_std_egy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_chl_b_ft(:,i:12:end),2));
        tm_obs_std_egy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_temp_ft(:,i:12:end),2));
        tm_obs_std_egy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_temp_b_ft(:,i:12:end),2));
        tm_obs_std_egy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_salt_ft(:,i:12:end),2));
        tm_obs_std_egy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_salt_b_ft(:,i:12:end),2));
        tm_obs_std_egy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_do_ft(:,i:12:end),2));
        tm_obs_std_egy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_egy_do_b_ft(:,i:12:end),2));

        tm_obs_std_jj_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_no3_ft(:,i:12:end),2));
        tm_obs_std_jj_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_no3_b_ft(:,i:12:end),2));
        tm_obs_std_jj_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_nh4_ft(:,i:12:end),2));
        tm_obs_std_jj_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_nh4_b_ft(:,i:12:end),2));
        tm_obs_std_jj_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_chl_ft(:,i:12:end),2));
        tm_obs_std_jj_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_chl_b_ft(:,i:12:end),2));
        tm_obs_std_jj_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_temp_ft(:,i:12:end),2));
        tm_obs_std_jj_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_temp_b_ft(:,i:12:end),2));
        tm_obs_std_jj_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_salt_ft(:,i:12:end),2));
        tm_obs_std_jj_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_salt_b_ft(:,i:12:end),2));
        tm_obs_std_jj_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_do_ft(:,i:12:end),2));
        tm_obs_std_jj_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanstd(nanmean(obs_jj_do_b_ft(:,i:12:end),2));

        
        %% mean
        obm_gy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_no3_ft(:,i:12:end)));
        obm_gy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_no3_b_ft(:,i:12:end)));
        obm_gy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_nh4_ft(:,i:12:end)));
        obm_gy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_nh4_b_ft(:,i:12:end)));
        obm_gy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_chl_ft(:,i:12:end)));
        obm_gy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_chl_b_ft(:,i:12:end)));
        obm_gy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_temp_ft(:,i:12:end)));
        obm_gy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_temp_b_ft(:,i:12:end)));
        obm_gy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_salt_ft(:,i:12:end)));
        obm_gy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_salt_b_ft(:,i:12:end)));
        obm_gy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_do_ft(:,i:12:end)));
        obm_gy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_gy_do_b_ft(:,i:12:end)));

        obm_sgy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_no3_ft(:,i:12:end)));
        obm_sgy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_no3_b_ft(:,i:12:end)));
        obm_sgy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_nh4_ft(:,i:12:end)));
        obm_sgy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_nh4_b_ft(:,i:12:end)));
        obm_sgy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_chl_ft(:,i:12:end)));
        obm_sgy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_chl_b_ft(:,i:12:end)));
        obm_sgy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_temp_ft(:,i:12:end)));
        obm_sgy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_temp_b_ft(:,i:12:end)));
        obm_sgy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_salt_ft(:,i:12:end)));
        obm_sgy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_salt_b_ft(:,i:12:end)));
        obm_sgy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_do_ft(:,i:12:end)));
        obm_sgy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_sgy_do_b_ft(:,i:12:end)));

        obm_egy_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_no3_ft(:,i:12:end)));
        obm_egy_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_no3_b_ft(:,i:12:end)));
        obm_egy_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_nh4_ft(:,i:12:end)));
        obm_egy_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_nh4_b_ft(:,i:12:end)));
        obm_egy_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_chl_ft(:,i:12:end)));
        obm_egy_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_chl_b_ft(:,i:12:end)));
        obm_egy_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_temp_ft(:,i:12:end)));
        obm_egy_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_temp_b_ft(:,i:12:end)));
        obm_egy_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_salt_ft(:,i:12:end)));
        obm_egy_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_salt_b_ft(:,i:12:end)));
        obm_egy_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_do_ft(:,i:12:end)));
        obm_egy_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_egy_do_b_ft(:,i:12:end)));

        obm_jj_no3_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_no3_ft(:,i:12:end)));
        obm_jj_no3_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_no3_b_ft(:,i:12:end)));
        obm_jj_nh4_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_nh4_ft(:,i:12:end)));
        obm_jj_nh4_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_nh4_b_ft(:,i:12:end)));
        obm_jj_chl_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_chl_ft(:,i:12:end)));
        obm_jj_chl_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_chl_b_ft(:,i:12:end)));
        obm_jj_temp_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_temp_ft(:,i:12:end)));
        obm_jj_temp_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_temp_b_ft(:,i:12:end)));
        obm_jj_salt_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_salt_ft(:,i:12:end)));
        obm_jj_salt_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_salt_b_ft(:,i:12:end)));
        obm_jj_do_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_do_ft(:,i:12:end)));
        obm_jj_do_b_ft(eom_d_each(1,i-1)+1:eom_d_each(1,i))=nanmean(nanmean(obs_jj_do_b_ft(:,i:12:end)));
    end
end

% save('koem_climate_3regime_v2_3sig_gy_to_2003_12.mat','obs*','obm*'); %default data: estimate 3sigma from spatially gy pt. (1st regime)
% save('koem_climate_3regime_v2_2sig_gy_to_2003_12.mat','obs*','obm*'); %default data: estimate 3sigma from spatially gy pt. (1st regime)
save(['koem_climate_3regime_v2_',num2str(sig),'sig_gy_to_2003_12.mat'],'obs*','obm*','tm_obs*','spm_obs*');
end

%% GY plot

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_no3_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_no3_ft./14 + obs_std_gy_no3_ft./14;
        lower_bound_plt = obm_gy_no3_ft./14 - obs_std_gy_no3_ft./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_no3_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_no3_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_no3_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 55])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4,'b','linew',2);

        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_nh4_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_nh4_ft./14 + obs_std_gy_nh4_ft./14;
        lower_bound_plt = obm_gy_nh4_ft./14 - obs_std_gy_nh4_ft./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_nh4_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_nh4_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_nh4_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 12])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);

        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_chl_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_chl_ft + obs_std_gy_chl_ft;
        lower_bound_plt = obm_gy_chl_ft - obs_std_gy_chl_ft;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_chl_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_chl_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);

        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_temp_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_temp_ft + obs_std_gy_temp_ft;
        lower_bound_plt = obm_gy_temp_ft - obs_std_gy_temp_ft;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_temp_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_temp_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_temp_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt,'b','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_salt_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_salt_ft + obs_std_gy_salt_ft;
        lower_bound_plt = obm_gy_salt_ft - obs_std_gy_salt_ft;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_salt_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_salt_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_salt_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency

        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3_b,'b','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_no3_b_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_no3_b_ft./14 + obs_std_gy_no3_b_ft./14;
        lower_bound_plt = obm_gy_no3_b_ft./14 - obs_std_gy_no3_b_ft./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_no3_b_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_no3_b_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_no3_b_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 55])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4_b,'b','linew',2);
       
                clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_nh4_b_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_nh4_b_ft./14 + obs_std_gy_nh4_b_ft./14;
        lower_bound_plt = obm_gy_nh4_b_ft./14 - obs_std_gy_nh4_b_ft./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_nh4_b_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_nh4_b_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_nh4_b_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 12])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl_b,'b','linew',2);
        
                clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_chl_b_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_chl_b_ft + obs_std_gy_chl_b_ft;
        lower_bound_plt = obm_gy_chl_b_ft - obs_std_gy_chl_b_ft;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_chl_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl_b_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_chl_b_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp_b,'b','linew',2);

                clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_temp_b_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_temp_b_ft + obs_std_gy_temp_b_ft;
        lower_bound_plt = obm_gy_temp_b_ft - obs_std_gy_temp_b_ft;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_temp_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_temp_b_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_temp_b_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt_b,'b','linew',2);
       
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_salt_b_ft)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_salt_b_ft + obs_std_gy_salt_b_ft;
        lower_bound_plt = obm_gy_salt_b_ft - obs_std_gy_salt_b_ft;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_salt_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_salt_b_ft),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_salt_b_ft),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot'),'-dpng')    
        
%%south gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        plot(obm_sgy_no3_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3_ft),obm_sgy_no3_ft./14 + obs_std_sgy_no3_ft./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3_ft),obm_sgy_no3_ft./14 - obs_std_sgy_no3_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_south_gy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4,'b','linew',2);
        plot(obm_sgy_nh4_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4_ft),obm_sgy_nh4_ft./14 + obs_std_sgy_nh4_ft./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4_ft),obm_sgy_nh4_ft./14 - obs_std_sgy_nh4_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_south_gy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);
        plot(obm_sgy_chl_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl_ft),obm_sgy_chl_ft + obs_std_sgy_chl_ft,'m-','linew',2);
        plot(1:length(obm_sgy_chl_ft),obm_sgy_chl_ft - obs_std_sgy_chl_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_south_gy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);
        plot(obm_sgy_temp_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp_ft),obm_sgy_temp_ft + obs_std_sgy_temp_ft,'m-','linew',2);
        plot(1:length(obm_sgy_temp_ft),obm_sgy_temp_ft - obs_std_sgy_temp_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_south_gy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt,'b','linew',2);
        plot(obm_sgy_salt_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt_ft),obm_sgy_salt_ft + obs_std_sgy_salt_ft,'m-','linew',2);
        plot(1:length(obm_sgy_salt_ft),obm_sgy_salt_ft - obs_std_sgy_salt_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_south_gy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3_b,'b','linew',2);
        plot(obm_sgy_no3_b_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3_b_ft),obm_sgy_no3_b_ft./14 + obs_std_sgy_no3_b_ft./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3_b_ft),obm_sgy_no3_b_ft./14 - obs_std_sgy_no3_b_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_south_gy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4_b,'b','linew',2);
        plot(obm_sgy_nh4_b_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4_b_ft),obm_sgy_nh4_b_ft./14 + obs_std_sgy_nh4_b_ft./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4_b_ft),obm_sgy_nh4_b_ft./14 - obs_std_sgy_nh4_b_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_south_gy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl_b,'b','linew',2);
        plot(obm_sgy_chl_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl_b_ft),obm_sgy_chl_b_ft + obs_std_sgy_chl_b_ft,'m-','linew',2);
        plot(1:length(obm_sgy_chl_b_ft),obm_sgy_chl_b_ft - obs_std_sgy_chl_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_south_gy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp_b,'b','linew',2);
        plot(obm_sgy_temp_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp_b_ft),obm_sgy_temp_b_ft + obs_std_sgy_temp_b_ft,'m-','linew',2);
        plot(1:length(obm_sgy_temp_b_ft),obm_sgy_temp_b_ft - obs_std_sgy_temp_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_south_gy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt_b,'b','linew',2);
        plot(obm_sgy_salt_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt_b_ft),obm_sgy_salt_b_ft + obs_std_sgy_salt_b_ft,'m-','linew',2);
        plot(1:length(obm_sgy_salt_b_ft),obm_sgy_salt_b_ft - obs_std_sgy_salt_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_south_gy'),'-dpng')    
        
%%east gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        plot(obm_egy_no3_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3_ft),obm_egy_no3_ft./14 + obs_std_egy_no3_ft./14,'m-','linew',2);
        plot(1:length(obm_egy_no3_ft),obm_egy_no3_ft./14 - obs_std_egy_no3_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4,'b','linew',2);
        plot(obm_egy_nh4_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4_ft),obm_egy_nh4_ft./14 + obs_std_egy_nh4_ft./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4_ft),obm_egy_nh4_ft./14 - obs_std_egy_nh4_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);
        plot(obm_egy_chl_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl_ft),obm_egy_chl_ft + obs_std_egy_chl_ft,'m-','linew',2);
        plot(1:length(obm_egy_chl_ft),obm_egy_chl_ft - obs_std_egy_chl_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);
        plot(obm_egy_temp_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp_ft),obm_egy_temp_ft + obs_std_egy_temp_ft,'m-','linew',2);
        plot(1:length(obm_egy_temp_ft),obm_egy_temp_ft - obs_std_egy_temp_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt,'b','linew',2);
        plot(obm_egy_salt_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt_ft),obm_egy_salt_ft + obs_std_egy_salt_ft,'m-','linew',2);
        plot(1:length(obm_egy_salt_ft),obm_egy_salt_ft - obs_std_egy_salt_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3_b,'b','linew',2);
        plot(obm_egy_no3_b_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3_b_ft),obm_egy_no3_b_ft./14 + obs_std_egy_no3_b_ft./14,'m-','linew',2);
        plot(1:length(obm_egy_no3_b_ft),obm_egy_no3_b_ft./14 - obs_std_egy_no3_b_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4_b,'b','linew',2);
        plot(obm_egy_nh4_b_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4_b_ft),obm_egy_nh4_b_ft./14 + obs_std_egy_nh4_b_ft./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4_b_ft),obm_egy_nh4_b_ft./14 - obs_std_egy_nh4_b_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl_b,'b','linew',2);
        plot(obm_egy_chl_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl_b_ft),obm_egy_chl_b_ft + obs_std_egy_chl_b_ft,'m-','linew',2);
        plot(1:length(obm_egy_chl_b_ft),obm_egy_chl_b_ft - obs_std_egy_chl_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp_b,'b','linew',2);
        plot(obm_egy_temp_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp_b_ft),obm_egy_temp_b_ft + obs_std_egy_temp_b_ft,'m-','linew',2);
        plot(1:length(obm_egy_temp_b_ft),obm_egy_temp_b_ft - obs_std_egy_temp_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt_b,'b','linew',2);
        plot(obm_egy_salt_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt_b_ft),obm_egy_salt_b_ft + obs_std_egy_salt_b_ft,'m-','linew',2);
        plot(1:length(obm_egy_salt_b_ft),obm_egy_salt_b_ft - obs_std_egy_salt_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot'),'-dpng')    
        
%% jinju

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        plot(obm_jj_no3_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3_ft),obm_jj_no3_ft./14 + obs_std_jj_no3_ft./14,'m-','linew',2);
        plot(1:length(obm_jj_no3_ft),obm_jj_no3_ft./14 - obs_std_jj_no3_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4,'b','linew',2);
        plot(obm_jj_nh4_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4_ft),obm_jj_nh4_ft./14 + obs_std_jj_nh4_ft./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4_ft),obm_jj_nh4_ft./14 - obs_std_jj_nh4_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);
        plot(obm_jj_chl_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl_ft),obm_jj_chl_ft + obs_std_jj_chl_ft,'m-','linew',2);
        plot(1:length(obm_jj_chl_ft),obm_jj_chl_ft - obs_std_jj_chl_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);
        plot(obm_jj_temp_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp_ft),obm_jj_temp_ft + obs_std_jj_temp_ft,'m-','linew',2);
        plot(1:length(obm_jj_temp_ft),obm_jj_temp_ft - obs_std_jj_temp_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt,'b','linew',2);
        plot(obm_jj_salt_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt_ft),obm_jj_salt_ft + obs_std_jj_salt_ft,'m-','linew',2);
        plot(1:length(obm_jj_salt_ft),obm_jj_salt_ft - obs_std_jj_salt_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3_b,'b','linew',2);
        plot(obm_jj_no3_b_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3_b_ft),obm_jj_no3_b_ft./14 + obs_std_jj_no3_b_ft./14,'m-','linew',2);
        plot(1:length(obm_jj_no3_b_ft),obm_jj_no3_b_ft./14 - obs_std_jj_no3_b_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4_b,'b','linew',2);
        plot(obm_jj_nh4_b_ft./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4_b_ft),obm_jj_nh4_b_ft./14 + obs_std_jj_nh4_b_ft./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4_b_ft),obm_jj_nh4_b_ft./14 - obs_std_jj_nh4_b_ft./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl_b,'b','linew',2);
        plot(obm_jj_chl_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl_b_ft),obm_jj_chl_b_ft + obs_std_jj_chl_b_ft,'m-','linew',2);
        plot(1:length(obm_jj_chl_b_ft),obm_jj_chl_b_ft - obs_std_jj_chl_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp_b,'b','linew',2);
        plot(obm_jj_temp_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp_b_ft),obm_jj_temp_b_ft + obs_std_jj_temp_b_ft,'m-','linew',2);
        plot(1:length(obm_jj_temp_b_ft),obm_jj_temp_b_ft - obs_std_jj_temp_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt_b,'b','linew',2);
        plot(obm_jj_salt_b_ft,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt_b_ft),obm_jj_salt_b_ft + obs_std_jj_salt_b_ft,'m-','linew',2);
        plot(1:length(obm_jj_salt_b_ft),obm_jj_salt_b_ft - obs_std_jj_salt_b_ft,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot'),'-dpng')    

        
        return
% %%
% %% GY plot
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3,'b','linew',2);
%         plot(obm_gy_no3./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_no3),obm_gy_no3./14 + obs_std_gy_no3./14,'m-','linew',2);
%         plot(1:length(obm_gy_no3),obm_gy_no3./14 - obs_std_gy_no3./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng')
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_nh4,'b','linew',2);
%         plot(obm_gy_nh4./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_nh4),obm_gy_nh4./14 + obs_std_gy_nh4./14,'m-','linew',2);
%         plot(1:length(obm_gy_nh4),obm_gy_nh4./14 - obs_std_gy_nh4./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4'),'-dpng')
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_chl,'b','linew',2);
%         plot(obm_gy_chl,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
%         plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
%         
%         
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_temp,'b','linew',2);
%         plot(obm_gy_temp,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_temp),obm_gy_temp + obs_std_gy_temp,'m-','linew',2);
%         plot(1:length(obm_gy_temp),obm_gy_temp - obs_std_gy_temp,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (o^C)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp'),'-dpng')    
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_salt,'b','linew',2);
%         plot(obm_gy_salt,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_salt),obm_gy_salt + obs_std_gy_salt,'m-','linew',2);
%         plot(1:length(obm_gy_salt),obm_gy_salt - obs_std_gy_salt,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt'),'-dpng')            
%  %% bot
%  
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3_b,'b','linew',2);
%         plot(obm_gy_no3_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_no3_b),obm_gy_no3_b./14 + obs_std_gy_no3_b./14,'m-','linew',2);
%         plot(1:length(obm_gy_no3_b),obm_gy_no3_b./14 - obs_std_gy_no3_b./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot'),'-dpng')
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_nh4_b,'b','linew',2);
%         plot(obm_gy_nh4_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_nh4_b),obm_gy_nh4_b./14 + obs_std_gy_nh4_b./14,'m-','linew',2);
%         plot(1:length(obm_gy_nh4_b),obm_gy_nh4_b./14 - obs_std_gy_nh4_b./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot'),'-dpng')
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_chl_b,'b','linew',2);
%         plot(obm_gy_chl_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_chl_b),obm_gy_chl_b + obs_std_gy_chl_b,'m-','linew',2);
%         plot(1:length(obm_gy_chl_b),obm_gy_chl_b - obs_std_gy_chl_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot'),'-dpng') 
%         
%         
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_temp_b,'b','linew',2);
%         plot(obm_gy_temp_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_temp_b),obm_gy_temp_b + obs_std_gy_temp_b,'m-','linew',2);
%         plot(1:length(obm_gy_temp_b),obm_gy_temp_b - obs_std_gy_temp_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (o^C)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot'),'-dpng')    
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_salt_b,'b','linew',2);
%         plot(obm_gy_salt_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_gy_salt_b),obm_gy_salt_b + obs_std_gy_salt_b,'m-','linew',2);
%         plot(1:length(obm_gy_salt_b),obm_gy_salt_b - obs_std_gy_salt_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot'),'-dpng')    
%         
% %%south gy plot
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3,'b','linew',2);
%         plot(obm_sgy_no3./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_no3),obm_sgy_no3./14 + obs_std_sgy_no3./14,'m-','linew',2);
%         plot(1:length(obm_sgy_no3),obm_sgy_no3./14 - obs_std_sgy_no3./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_south_gy'),'-dpng')
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_nh4,'b','linew',2);
%         plot(obm_sgy_nh4./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 + obs_std_sgy_nh4./14,'m-','linew',2);
%         plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 - obs_std_sgy_nh4./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_south_gy'),'-dpng')
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_chl,'b','linew',2);
%         plot(obm_sgy_chl,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_chl),obm_sgy_chl + obs_std_sgy_chl,'m-','linew',2);
%         plot(1:length(obm_sgy_chl),obm_sgy_chl - obs_std_sgy_chl,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_south_gy'),'-dpng') 
%         
%         
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_temp,'b','linew',2);
%         plot(obm_sgy_temp,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_temp),obm_sgy_temp + obs_std_sgy_temp,'m-','linew',2);
%         plot(1:length(obm_sgy_temp),obm_sgy_temp - obs_std_sgy_temp,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (o^C)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_south_gy'),'-dpng')    
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_salt,'b','linew',2);
%         plot(obm_sgy_salt,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_salt),obm_sgy_salt + obs_std_sgy_salt,'m-','linew',2);
%         plot(1:length(obm_sgy_salt),obm_sgy_salt - obs_std_sgy_salt,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_south_gy'),'-dpng')            
%  %% bot
%  
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3_b,'b','linew',2);
%         plot(obm_sgy_no3_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 + obs_std_sgy_no3_b./14,'m-','linew',2);
%         plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 - obs_std_sgy_no3_b./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_south_gy'),'-dpng')
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_nh4_b,'b','linew',2);
%         plot(obm_sgy_nh4_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 + obs_std_sgy_nh4_b./14,'m-','linew',2);
%         plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 - obs_std_sgy_nh4_b./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_south_gy'),'-dpng')
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_chl_b,'b','linew',2);
%         plot(obm_sgy_chl_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b + obs_std_sgy_chl_b,'m-','linew',2);
%         plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b - obs_std_sgy_chl_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_south_gy'),'-dpng') 
%         
%         
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_temp_b,'b','linew',2);
%         plot(obm_sgy_temp_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b + obs_std_sgy_temp_b,'m-','linew',2);
%         plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b - obs_std_sgy_temp_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (o^C)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_south_gy'),'-dpng')    
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_salt_b,'b','linew',2);
%         plot(obm_sgy_salt_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b + obs_std_sgy_salt_b,'m-','linew',2);
%         plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b - obs_std_sgy_salt_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_south_gy'),'-dpng')    
%         
% %%east gy plot
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3,'b','linew',2);
%         plot(obm_egy_no3./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_no3),obm_egy_no3./14 + obs_std_egy_no3./14,'m-','linew',2);
%         plot(1:length(obm_egy_no3),obm_egy_no3./14 - obs_std_egy_no3./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng')
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_nh4,'b','linew',2);
%         plot(obm_egy_nh4./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_nh4),obm_egy_nh4./14 + obs_std_egy_nh4./14,'m-','linew',2);
%         plot(1:length(obm_egy_nh4),obm_egy_nh4./14 - obs_std_egy_nh4./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4'),'-dpng')
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_chl,'b','linew',2);
%         plot(obm_egy_chl,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_chl),obm_egy_chl + obs_std_egy_chl,'m-','linew',2);
%         plot(1:length(obm_egy_chl),obm_egy_chl - obs_std_egy_chl,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
%         
%         
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_temp,'b','linew',2);
%         plot(obm_egy_temp,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_temp),obm_egy_temp + obs_std_egy_temp,'m-','linew',2);
%         plot(1:length(obm_egy_temp),obm_egy_temp - obs_std_egy_temp,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (o^C)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp'),'-dpng')    
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_salt,'b','linew',2);
%         plot(obm_egy_salt,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_salt),obm_egy_salt + obs_std_egy_salt,'m-','linew',2);
%         plot(1:length(obm_egy_salt),obm_egy_salt - obs_std_egy_salt,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt'),'-dpng')            
%  %% bot
%  
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3_b,'b','linew',2);
%         plot(obm_egy_no3_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 + obs_std_egy_no3_b./14,'m-','linew',2);
%         plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 - obs_std_egy_no3_b./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot'),'-dpng')
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_nh4_b,'b','linew',2);
%         plot(obm_egy_nh4_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 + obs_std_egy_nh4_b./14,'m-','linew',2);
%         plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 - obs_std_egy_nh4_b./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot'),'-dpng')
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_chl_b,'b','linew',2);
%         plot(obm_egy_chl_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_chl_b),obm_egy_chl_b + obs_std_egy_chl_b,'m-','linew',2);
%         plot(1:length(obm_egy_chl_b),obm_egy_chl_b - obs_std_egy_chl_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot'),'-dpng') 
%         
%         
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_temp_b,'b','linew',2);
%         plot(obm_egy_temp_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_temp_b),obm_egy_temp_b + obs_std_egy_temp_b,'m-','linew',2);
%         plot(1:length(obm_egy_temp_b),obm_egy_temp_b - obs_std_egy_temp_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (o^C)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot'),'-dpng')    
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_salt_b,'b','linew',2);
%         plot(obm_egy_salt_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_salt_b),obm_egy_salt_b + obs_std_egy_salt_b,'m-','linew',2);
%         plot(1:length(obm_egy_salt_b),obm_egy_salt_b - obs_std_egy_salt_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot'),'-dpng')    
%         
% %% jinju
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3,'b','linew',2);
%         plot(obm_jj_no3./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_no3),obm_jj_no3./14 + obs_std_jj_no3./14,'m-','linew',2);
%         plot(1:length(obm_jj_no3),obm_jj_no3./14 - obs_std_jj_no3./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng')
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_nh4,'b','linew',2);
%         plot(obm_jj_nh4./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_nh4),obm_jj_nh4./14 + obs_std_jj_nh4./14,'m-','linew',2);
%         plot(1:length(obm_jj_nh4),obm_jj_nh4./14 - obs_std_jj_nh4./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4'),'-dpng')
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_chl,'b','linew',2);
%         plot(obm_jj_chl,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_chl),obm_jj_chl + obs_std_jj_chl,'m-','linew',2);
%         plot(1:length(obm_jj_chl),obm_jj_chl - obs_std_jj_chl,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
%         
%         
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_temp,'b','linew',2);
%         plot(obm_jj_temp,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_temp),obm_jj_temp + obs_std_jj_temp,'m-','linew',2);
%         plot(1:length(obm_jj_temp),obm_jj_temp - obs_std_jj_temp,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (o^C)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp'),'-dpng')    
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_salt,'b','linew',2);
%         plot(obm_jj_salt,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_salt),obm_jj_salt + obs_std_jj_salt,'m-','linew',2);
%         plot(1:length(obm_jj_salt),obm_jj_salt - obs_std_jj_salt,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt'),'-dpng')            
%  %% bot
%  
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3_b,'b','linew',2);
%         plot(obm_jj_no3_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 + obs_std_jj_no3_b./14,'m-','linew',2);
%         plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 - obs_std_jj_no3_b./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot'),'-dpng')
% 
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_nh4_b,'b','linew',2);
%         plot(obm_jj_nh4_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 + obs_std_jj_nh4_b./14,'m-','linew',2);
%         plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 - obs_std_jj_nh4_b./14,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot'),'-dpng')
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_chl_b,'b','linew',2);
%         plot(obm_jj_chl_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_chl_b),obm_jj_chl_b + obs_std_jj_chl_b,'m-','linew',2);
%         plot(1:length(obm_jj_chl_b),obm_jj_chl_b - obs_std_jj_chl_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot'),'-dpng') 
%         
%         
%  fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_temp_b,'b','linew',2);
%         plot(obm_jj_temp_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_temp_b),obm_jj_temp_b + obs_std_jj_temp_b,'m-','linew',2);
%         plot(1:length(obm_jj_temp_b),obm_jj_temp_b - obs_std_jj_temp_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (o^C)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot'),'-dpng')    
%         
%         
% fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_salt_b,'b','linew',2);
%         plot(obm_jj_salt_b,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_salt_b),obm_jj_salt_b + obs_std_jj_salt_b,'m-','linew',2);
%         plot(1:length(obm_jj_salt_b),obm_jj_salt_b - obs_std_jj_salt_b,'m-','linew',2);
%         title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot'),'-dpng')    
%         
%                 
%         
%         
%         
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(gy_no3,'b','linew',2);
%         plot(obs_gy_no3(i,:)./14,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NO3']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%         close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_nh4(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_nh4(i,:)./14,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NH4']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_NH4_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off; close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_do(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_do(i,:).*0.7.*44.661,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL DO']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('DO (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_DO_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_temp(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_temp(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL temp']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('temp (^oC)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_temp_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_salt(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_salt(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL salt']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_salt_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_chl(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_chl(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL chla']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_chl_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%         close;
%    end
% end
