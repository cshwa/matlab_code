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

% return
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
find_t_ind2 = {'2010-01'};
find_date2=find(strcmp(ref_date,find_t_ind2) == 1);
yr10_filter(1:find_date) = NaN; % pass to 2003-12
yr10_filter(find_date+1:find_date2 -1) = 1; % remove 2004~2009
yr10_filter(find_date2:length(ref_date)) = NaN; % remove 2010~2018


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

% % save('koem_climate_3regime_v2_to_2003_12.mat','obs*','obm*'); % default data: estimate 3sigma from spatially all pt.
% % save('koem_climate_3regime_v2_3sig_gy_to_2003_12.mat','obs*','obm*'); %default data: estimate 3sigma from spatially gy pt. (1st regime)
% save('koem_climate_3regime_v2_3sig_gy_04to09.mat','obs*','obm*'); %default data: estimate 3sigma from spatially gy pt. (2nd regime)
% save('koem_climate_3regime_v2_2sig_gy_04to09.mat','obs*','obm*'); %default data: estimate 3sigma from spatially gy pt. (2nd regime)
save(['koem_climate_3regime_v2_',num2str(sig),'sig_gy_04to09.mat'],'obs*','obm*','spm_obs*','tm_obs*');
end
