close all; clear; clc;

sigsig=[1;];%2;3;
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

%% 5yr cutting
% ref_date{36} = 1997-01~1999-12'  1st cut
% ref_date{37}:ref_date{96} =2000-01~2004-12'  2nd cut
% ref_date{97}:ref_date{156} =2005-01~2009-12'  3rd cut

%% 10yr cutting
% ref_date{48} = 1997-01~2000-12'  1st cut
% ref_date{49}:ref_date{168} =2001-01~2010-12'  2nd cut
% ref_date{169}:length(ref_date) =2011-01~2018-12'  extra

find_t_ind = {'2000-12'};
find_date=find(strcmp(ref_date,find_t_ind) == 1);
find_t_ind2 = {'2010-12'};
find_date2=find(strcmp(ref_date,find_t_ind2) == 1);

yr10_filter_1(1:find_date) = 1; % pass to ~1999-12'  1st cut
yr10_filter_1(find_date+1:find_date2) = NaN; % remove 2000-01~2009-12'  2nd cut
yr10_filter_1(find_date2+1:length(ref_date)) = NaN; % remove 2010-01~2018-12'  extra

yr10_filter_2(1:find_date) = NaN; % pass to ~1999-12'  1st cut
yr10_filter_2(find_date+1:find_date2) = 1; % remove 2000-01~2009-12'  2nd cut
yr10_filter_2(find_date2+1:length(ref_date)) = NaN; % remove 2010-01~2018-12'  extra


for i = 1:65
regime_no3_1(i,:) = regime_no3(i,:) .* yr10_filter_1; 
regime_nh4_1(i,:) = regime_nh4(i,:) .* yr10_filter_1; 
regime_do_1(i,:) = regime_do(i,:) .* yr10_filter_1; 
regime_chl_1(i,:) = regime_chl(i,:) .* yr10_filter_1; 
regime_temp_1(i,:) = regime_temp(i,:) .* yr10_filter_1; 
regime_salt_1(i,:) = regime_salt(i,:) .* yr10_filter_1; 
regime_no3_b_1(i,:) = regime_no3_b(i,:) .* yr10_filter_1; 
regime_nh4_b_1(i,:) = regime_nh4_b(i,:) .* yr10_filter_1; 
regime_do_b_1(i,:) = regime_do_b(i,:) .* yr10_filter_1; 
regime_chl_b_1(i,:) = regime_chl_b(i,:) .* yr10_filter_1; 
regime_temp_b_1(i,:) = regime_temp_b(i,:) .* yr10_filter_1; 
regime_salt_b_1(i,:) = regime_salt_b(i,:) .* yr10_filter_1;
regime_po4_1(i,:) = regime_po4(i,:) .* yr10_filter_1; 
regime_po4_b_1(i,:) = regime_po4_b(i,:) .* yr10_filter_1;

regime_no3_2(i,:) = regime_no3(i,:) .* yr10_filter_2; 
regime_nh4_2(i,:) = regime_nh4(i,:) .* yr10_filter_2; 
regime_do_2(i,:) = regime_do(i,:) .* yr10_filter_2; 
regime_chl_2(i,:) = regime_chl(i,:) .* yr10_filter_2; 
regime_temp_2(i,:) = regime_temp(i,:) .* yr10_filter_2; 
regime_salt_2(i,:) = regime_salt(i,:) .* yr10_filter_2; 
regime_no3_b_2(i,:) = regime_no3_b(i,:) .* yr10_filter_2; 
regime_nh4_b_2(i,:) = regime_nh4_b(i,:) .* yr10_filter_2; 
regime_do_b_2(i,:) = regime_do_b(i,:) .* yr10_filter_2; 
regime_chl_b_2(i,:) = regime_chl_b(i,:) .* yr10_filter_2; 
regime_temp_b_2(i,:) = regime_temp_b(i,:) .* yr10_filter_2; 
regime_salt_b_2(i,:) = regime_salt_b(i,:) .* yr10_filter_2;
regime_po4_2(i,:) = regime_po4(i,:) .* yr10_filter_2;
regime_po4_b_2(i,:) = regime_po4_b(i,:) .* yr10_filter_2;

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

%sp_gy 1st
    obs_gy_no3_ft_1=squeeze(regime_no3_1(sp_gy,:));
    obs_gy_no3_b_ft_1=squeeze(regime_no3_b_1(sp_gy,:));
    obs_gy_nh4_ft_1=squeeze(regime_nh4_1(sp_gy,:));
    obs_gy_nh4_b_ft_1=squeeze(regime_nh4_b_1(sp_gy,:));
    obs_gy_chl_ft_1=squeeze(regime_chl_1(sp_gy,:));
    obs_gy_chl_b_ft_1=squeeze(regime_chl_b_1(sp_gy,:));
    obs_gy_temp_ft_1=squeeze(regime_temp_1(sp_gy,:));
    obs_gy_temp_b_ft_1=squeeze(regime_temp_b_1(sp_gy,:));
    obs_gy_salt_ft_1=squeeze(regime_salt_1(sp_gy,:));
    obs_gy_salt_b_ft_1=squeeze(regime_salt_b_1(sp_gy,:));
    obs_gy_do_ft_1=squeeze(regime_do_1(sp_gy,:));
    obs_gy_do_b_ft_1=squeeze(regime_do_b_1(sp_gy,:));
    obs_gy_po4_ft_1=squeeze(regime_po4_1(sp_gy,:));
    obs_gy_po4_b_ft_1=squeeze(regime_po4_b_1(sp_gy,:));

%sp_gy 2nd
    obs_gy_no3_ft_2=squeeze(regime_no3_2(sp_gy,:));
    obs_gy_no3_b_ft_2=squeeze(regime_no3_b_2(sp_gy,:));
    obs_gy_nh4_ft_2=squeeze(regime_nh4_2(sp_gy,:));
    obs_gy_nh4_b_ft_2=squeeze(regime_nh4_b_2(sp_gy,:));
    obs_gy_chl_ft_2=squeeze(regime_chl_2(sp_gy,:));
    obs_gy_chl_b_ft_2=squeeze(regime_chl_b_2(sp_gy,:));
    obs_gy_temp_ft_2=squeeze(regime_temp_2(sp_gy,:));
    obs_gy_temp_b_ft_2=squeeze(regime_temp_b_2(sp_gy,:));
    obs_gy_salt_ft_2=squeeze(regime_salt_2(sp_gy,:));
    obs_gy_salt_b_ft_2=squeeze(regime_salt_b_2(sp_gy,:));
    obs_gy_po4_ft_2=squeeze(regime_po4_2(sp_gy,:));
    obs_gy_po4_b_ft_2=squeeze(regime_po4_b_2(sp_gy,:));
    obs_gy_do_ft_2=squeeze(regime_do_2(sp_gy,:));
    obs_gy_do_b_ft_2=squeeze(regime_do_b_2(sp_gy,:));  

for i = 1:12
     if i==1
         clearvars t_eom_interval
         t_eom_interval = 1:eom_d_each(1,i);
     else    
         clearvars t_eom_interval
         t_eom_interval = eom_d_each(1,i-1)+1:eom_d_each(1,i);
     end
        %% time & spatial std 
        clearvars len_temp
        len_temp=size(obs_gy_no3_ft_1(:,i:12:end),1)*size(obs_gy_no3_ft_1(:,i:12:end),2);
        obs_std_gy_no3_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_no3_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_no3_b_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_no3_b_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_nh4_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_nh4_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_nh4_b_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_nh4_b_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_chl_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_chl_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_chl_b_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_chl_b_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_temp_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_temp_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_temp_b_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_temp_b_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_salt_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_salt_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_salt_b_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_salt_b_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_do_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_do_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_do_b_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_do_b_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_po4_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_po4_ft_1(:,i:12:end),1,len_temp));
        obs_std_gy_po4_b_ft_1(t_eom_interval)=nanstd(reshape(obs_gy_po4_b_ft_1(:,i:12:end),1,len_temp));

        
        clearvars len_temp
        len_temp=size(obs_gy_no3_ft_2(:,i:12:end),1)*size(obs_gy_no3_ft_2(:,i:12:end),2);
        obs_std_gy_no3_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_no3_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_no3_b_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_no3_b_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_nh4_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_nh4_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_nh4_b_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_nh4_b_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_chl_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_chl_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_chl_b_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_chl_b_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_temp_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_temp_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_temp_b_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_temp_b_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_salt_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_salt_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_salt_b_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_salt_b_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_do_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_do_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_do_b_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_do_b_ft_2(:,i:12:end),1,len_temp));        
        obs_std_gy_po4_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_po4_ft_2(:,i:12:end),1,len_temp));
        obs_std_gy_po4_b_ft_2(t_eom_interval)=nanstd(reshape(obs_gy_po4_b_ft_2(:,i:12:end),1,len_temp));
        
       %% time std
        spm_obs_std_gy_no3_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_no3_ft_1(:,i:12:end),1));
        spm_obs_std_gy_no3_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_no3_b_ft_1(:,i:12:end),1));
        spm_obs_std_gy_nh4_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_nh4_ft_1(:,i:12:end),1));
        spm_obs_std_gy_nh4_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_nh4_b_ft_1(:,i:12:end),1));
        spm_obs_std_gy_chl_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_chl_ft_1(:,i:12:end),1));
        spm_obs_std_gy_chl_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_chl_b_ft_1(:,i:12:end),1));
        spm_obs_std_gy_temp_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_temp_ft_1(:,i:12:end),1));
        spm_obs_std_gy_temp_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_temp_b_ft_1(:,i:12:end),1));
        spm_obs_std_gy_salt_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_salt_ft_1(:,i:12:end),1));
        spm_obs_std_gy_salt_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_salt_b_ft_1(:,i:12:end),1));
        spm_obs_std_gy_do_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_do_ft_1(:,i:12:end),1));
        spm_obs_std_gy_do_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_do_b_ft_1(:,i:12:end),1));
        spm_obs_std_gy_po4_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_po4_ft_1(:,i:12:end),1));
        spm_obs_std_gy_po4_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_po4_b_ft_1(:,i:12:end),1));

        spm_obs_std_gy_no3_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_no3_ft_2(:,i:12:end),1));
        spm_obs_std_gy_no3_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_no3_b_ft_2(:,i:12:end),1));
        spm_obs_std_gy_nh4_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_nh4_ft_2(:,i:12:end),1));
        spm_obs_std_gy_nh4_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_nh4_b_ft_2(:,i:12:end),1));
        spm_obs_std_gy_chl_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_chl_ft_2(:,i:12:end),1));
        spm_obs_std_gy_chl_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_chl_b_ft_2(:,i:12:end),1));
        spm_obs_std_gy_temp_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_temp_ft_2(:,i:12:end),1));
        spm_obs_std_gy_temp_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_temp_b_ft_2(:,i:12:end),1));
        spm_obs_std_gy_salt_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_salt_ft_2(:,i:12:end),1));
        spm_obs_std_gy_salt_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_salt_b_ft_2(:,i:12:end),1));
        spm_obs_std_gy_do_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_do_ft_2(:,i:12:end),1));
        spm_obs_std_gy_do_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_do_b_ft_2(:,i:12:end),1));
        spm_obs_std_gy_po4_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_po4_ft_2(:,i:12:end),1));
        spm_obs_std_gy_po4_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_po4_b_ft_2(:,i:12:end),1));
        
       %% spatial std
        tm_obs_std_gy_no3_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_no3_ft_1(:,i:12:end),2));
        tm_obs_std_gy_no3_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_no3_b_ft_1(:,i:12:end),2));
        tm_obs_std_gy_nh4_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_nh4_ft_1(:,i:12:end),2));
        tm_obs_std_gy_nh4_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_nh4_b_ft_1(:,i:12:end),2));
        tm_obs_std_gy_chl_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_chl_ft_1(:,i:12:end),2));
        tm_obs_std_gy_chl_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_chl_b_ft_1(:,i:12:end),2));
        tm_obs_std_gy_temp_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_temp_ft_1(:,i:12:end),2));
        tm_obs_std_gy_temp_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_temp_b_ft_1(:,i:12:end),2));
        tm_obs_std_gy_salt_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_salt_ft_1(:,i:12:end),2));
        tm_obs_std_gy_salt_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_salt_b_ft_1(:,i:12:end),2));
        tm_obs_std_gy_do_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_do_ft_1(:,i:12:end),2));
        tm_obs_std_gy_do_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_do_b_ft_1(:,i:12:end),2));
        tm_obs_std_gy_po4_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_po4_ft_1(:,i:12:end),2));
        tm_obs_std_gy_po4_b_ft_1(t_eom_interval)=nanstd(nanmean(obs_gy_po4_b_ft_1(:,i:12:end),2));
        
        tm_obs_std_gy_no3_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_no3_ft_2(:,i:12:end),2));
        tm_obs_std_gy_no3_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_no3_b_ft_2(:,i:12:end),2));
        tm_obs_std_gy_nh4_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_nh4_ft_2(:,i:12:end),2));
        tm_obs_std_gy_nh4_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_nh4_b_ft_2(:,i:12:end),2));
        tm_obs_std_gy_chl_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_chl_ft_2(:,i:12:end),2));
        tm_obs_std_gy_chl_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_chl_b_ft_2(:,i:12:end),2));
        tm_obs_std_gy_temp_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_temp_ft_2(:,i:12:end),2));
        tm_obs_std_gy_temp_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_temp_b_ft_2(:,i:12:end),2));
        tm_obs_std_gy_salt_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_salt_ft_2(:,i:12:end),2));
        tm_obs_std_gy_salt_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_salt_b_ft_2(:,i:12:end),2));
        tm_obs_std_gy_do_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_do_ft_2(:,i:12:end),2));
        tm_obs_std_gy_do_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_do_b_ft_2(:,i:12:end),2));
        tm_obs_std_gy_po4_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_po4_ft_2(:,i:12:end),2));
        tm_obs_std_gy_po4_b_ft_2(t_eom_interval)=nanstd(nanmean(obs_gy_po4_b_ft_2(:,i:12:end),2));
        
        %% mean
        obm_gy_no3_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_no3_ft_1(:,i:12:end)));
        obm_gy_no3_b_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_no3_b_ft_1(:,i:12:end)));
        obm_gy_nh4_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_nh4_ft_1(:,i:12:end)));
        obm_gy_nh4_b_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_nh4_b_ft_1(:,i:12:end)));
        obm_gy_chl_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_chl_ft_1(:,i:12:end)));
        obm_gy_chl_b_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_chl_b_ft_1(:,i:12:end)));
        obm_gy_temp_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_temp_ft_1(:,i:12:end)));
        obm_gy_temp_b_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_temp_b_ft_1(:,i:12:end)));
        obm_gy_salt_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_salt_ft_1(:,i:12:end)));
        obm_gy_salt_b_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_salt_b_ft_1(:,i:12:end)));
        obm_gy_do_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_do_ft_1(:,i:12:end)));
        obm_gy_do_b_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_do_b_ft_1(:,i:12:end)));
        obm_gy_po4_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_po4_ft_1(:,i:12:end)));
        obm_gy_po4_b_ft_1(t_eom_interval)=nanmean(nanmean(obs_gy_po4_b_ft_1(:,i:12:end)));

        obm_gy_no3_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_no3_ft_2(:,i:12:end)));
        obm_gy_no3_b_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_no3_b_ft_2(:,i:12:end)));
        obm_gy_nh4_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_nh4_ft_2(:,i:12:end)));
        obm_gy_nh4_b_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_nh4_b_ft_2(:,i:12:end)));
        obm_gy_chl_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_chl_ft_2(:,i:12:end)));
        obm_gy_chl_b_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_chl_b_ft_2(:,i:12:end)));
        obm_gy_temp_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_temp_ft_2(:,i:12:end)));
        obm_gy_temp_b_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_temp_b_ft_2(:,i:12:end)));
        obm_gy_salt_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_salt_ft_2(:,i:12:end)));
        obm_gy_salt_b_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_salt_b_ft_2(:,i:12:end)));
        obm_gy_do_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_do_ft_2(:,i:12:end)));
        obm_gy_do_b_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_do_b_ft_2(:,i:12:end)));
        obm_gy_po4_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_po4_ft_2(:,i:12:end)));
        obm_gy_po4_b_ft_2(t_eom_interval)=nanmean(nanmean(obs_gy_po4_b_ft_2(:,i:12:end)));
        
    
end
   save(['koem_climate_10yr_v2_',num2str(sig),'sig_gy.mat'],'obs*','obm*','spm_obs*','tm_obs*'); 
end

