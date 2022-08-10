close all; clear;clc

case_num = [1:3];
sigsig=[2;3;];
% 1: spatial & time std.
% 2: time std.
% 3: spatial std.
 dir_name ={'all_std','t_std','sp_std'};


cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
ltme =load('LTME_observation_data_v2.mat'); %LTME obs.

% cd D:\장기생태\Dynamic\KOEM\gy_2001\gy_2001_koem_daily\spmean
% load('koem_result_processed_std_spmean.mat'); % KOEM obs.
% % mpdata= load('mpdata_result.mat');

cd D:\장기생태\Dynamic\result\1997
c1997 = load('1997_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\1998
c1998 = load('1998_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\1999
c1999 = load('1999_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2000
c2000 = load('2000_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2001
c2001 = load('2001_mp_p_sewer_det_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2002
c2002 = load('2002_mp_p_sewer_det_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2003
c2003 = load('2003_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

c1997_eom = [1 (c1997.eom_d_each(1:end-1) + 1)];
c1998_eom = [(c1997.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(1:end-1) + 1)];
c1999_eom = [(c1997.eom_d_each(end)+c1998.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(1:end-1) + 1)];
c2000_eom = [(c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end) +c2000.eom_d_each(1:end-1)+ 1)];
c2001_eom = [(c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end)+c2000.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end) +c2000.eom_d_each(end)+c2001.eom_d_each(1:end-1)+ 1)];
c2002_eom = [(c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end)+c2000.eom_d_each(end)+c2001.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end) +c2000.eom_d_each(end)+c2001.eom_d_each(end)+c2002.eom_d_each(1:end-1) + 1)];
c2003_eom = [(c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end)+c2000.eom_d_each(end)+c2001.eom_d_each(end)+c2002.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end)+c2000.eom_d_each(end)+c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(1:end-1) + 1)];

com_eom = [c1997_eom c1998_eom c1999_eom c2000_eom c2001_eom c2002_eom c2003_eom];

c1997.sp_gy_no3(:,60)=[]; c1997.sp_gy_chl(:,60)=[]; c1997.sp_gy_nh4(:,60)=[]; c1997.sp_gy_temp(:,60)=[]; c1997.sp_gy_salt(:,60)=[];
c1998.sp_gy_no3(:,60)=[]; c1998.sp_gy_chl(:,60)=[]; c1998.sp_gy_nh4(:,60)=[]; c1998.sp_gy_temp(:,60)=[]; c1998.sp_gy_salt(:,60)=[];
c1999.sp_gy_no3(:,60)=[]; c1999.sp_gy_chl(:,60)=[]; c1999.sp_gy_nh4(:,60)=[]; c1999.sp_gy_temp(:,60)=[]; c1999.sp_gy_salt(:,60)=[];
c2000.sp_gy_no3(:,60)=[]; c2000.sp_gy_chl(:,60)=[]; c2000.sp_gy_nh4(:,60)=[]; c2000.sp_gy_temp(:,60)=[]; c2000.sp_gy_salt(:,60)=[];

c1997.sp_gy_no3_b(:,60)=[]; c1997.sp_gy_chl_b(:,60)=[]; c1997.sp_gy_nh4_b(:,60)=[]; c1997.sp_gy_temp_b(:,60)=[]; c1997.sp_gy_salt_b(:,60)=[];
c1998.sp_gy_no3_b(:,60)=[]; c1998.sp_gy_chl_b(:,60)=[]; c1998.sp_gy_nh4_b(:,60)=[]; c1998.sp_gy_temp_b(:,60)=[]; c1998.sp_gy_salt_b(:,60)=[];
c1999.sp_gy_no3_b(:,60)=[]; c1999.sp_gy_chl_b(:,60)=[]; c1999.sp_gy_nh4_b(:,60)=[]; c1999.sp_gy_temp_b(:,60)=[]; c1999.sp_gy_salt_b(:,60)=[];
c2000.sp_gy_no3_b(:,60)=[]; c2000.sp_gy_chl_b(:,60)=[]; c2000.sp_gy_nh4_b(:,60)=[]; c2000.sp_gy_temp_b(:,60)=[]; c2000.sp_gy_salt_b(:,60)=[];

c1997.gy_no3(60)=[]; c1997.gy_chl(60)=[]; c1997.gy_nh4(60)=[]; c1997.gy_temp(60)=[]; c1997.gy_salt(60)=[];
c1998.gy_no3(60)=[]; c1998.gy_chl(60)=[]; c1998.gy_nh4(60)=[]; c1998.gy_temp(60)=[]; c1998.gy_salt(60)=[];
c1999.gy_no3(60)=[]; c1999.gy_chl(60)=[]; c1999.gy_nh4(60)=[]; c1999.gy_temp(60)=[]; c1999.gy_salt(60)=[];
c2000.gy_no3(60)=[]; c2000.gy_chl(60)=[]; c2000.gy_nh4(60)=[]; c2000.gy_temp(60)=[]; c2000.gy_salt(60)=[];

c1997.gy_no3_b(60)=[]; c1997.gy_chl_b(60)=[]; c1997.gy_nh4_b(60)=[]; c1997.gy_temp_b(60)=[]; c1997.gy_salt_b(60)=[];
c1998.gy_no3_b(60)=[]; c1998.gy_chl_b(60)=[]; c1998.gy_nh4_b(60)=[]; c1998.gy_temp_b(60)=[]; c1998.gy_salt_b(60)=[];
c1999.gy_no3_b(60)=[]; c1999.gy_chl_b(60)=[]; c1999.gy_nh4_b(60)=[]; c1999.gy_temp_b(60)=[]; c1999.gy_salt_b(60)=[];
c2000.gy_no3_b(60)=[]; c2000.gy_chl_b(60)=[]; c2000.gy_nh4_b(60)=[]; c2000.gy_temp_b(60)=[]; c2000.gy_salt_b(60)=[];

        
        
mer_no3 = [c1997.gy_no3; c1998.gy_no3; c1999.gy_no3; c2000.gy_no3; c2001.gy_no3; ...
    c2002.gy_no3; c2003.gy_no3;];
mer_nh4 = [c1997.gy_nh4; c1998.gy_nh4; c1999.gy_nh4; c2000.gy_nh4; c2001.gy_nh4; ...
    c2002.gy_nh4; c2003.gy_nh4; ];
mer_chl = [c1997.gy_chl; c1998.gy_chl; c1999.gy_chl; c2000.gy_chl; c2001.gy_chl; ...
    c2002.gy_chl; c2003.gy_chl;];
mer_temp = [c1997.gy_temp; c1998.gy_temp; c1999.gy_temp; c2000.gy_temp; c2001.gy_temp; ...
    c2002.gy_temp; c2003.gy_temp;];
mer_salt = [c1997.gy_salt; c1998.gy_salt; c1999.gy_salt; c2000.gy_salt; c2001.gy_salt; ...
    c2002.gy_salt; c2003.gy_salt; ];

%bot merged
mer_no3_b = [c1997.gy_no3_b; c1998.gy_no3_b; c1999.gy_no3_b; c2000.gy_no3_b; c2001.gy_no3_b; ...
    c2002.gy_no3_b; c2003.gy_no3_b; ];
mer_nh4_b = [c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; ...
    c2002.gy_nh4_b; c2003.gy_nh4_b;];
mer_chl_b = [c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; ...
    c2002.gy_chl_b; c2003.gy_chl_b; ];
mer_temp_b = [c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; ...
    c2002.gy_temp_b; c2003.gy_temp_b; ];
mer_salt_b = [c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; ...
    c2002.gy_salt_b; c2003.gy_salt_b;];
        
% sptatial value
%surf merged
mer_sp_no3 = [c1997.sp_gy_no3 c1998.sp_gy_no3 c1999.sp_gy_no3 c2000.sp_gy_no3 c2001.sp_gy_no3 ...
    c2002.sp_gy_no3 c2003.sp_gy_no3 ];
mer_sp_nh4 = [c1997.sp_gy_nh4 c1998.sp_gy_nh4 c1999.sp_gy_nh4 c2000.sp_gy_nh4 c2001.sp_gy_nh4 ...
    c2002.sp_gy_nh4 c2003.sp_gy_nh4 ];
mer_sp_chl = [c1997.sp_gy_chl c1998.sp_gy_chl c1999.sp_gy_chl c2000.sp_gy_chl c2001.sp_gy_chl ...
    c2002.sp_gy_chl c2003.sp_gy_chl ];
mer_sp_temp = [c1997.sp_gy_temp c1998.sp_gy_temp c1999.sp_gy_temp c2000.sp_gy_temp c2001.sp_gy_temp ...
    c2002.sp_gy_temp c2003.sp_gy_temp ];
mer_sp_salt = [c1997.sp_gy_salt c1998.sp_gy_salt c1999.sp_gy_salt c2000.sp_gy_salt c2001.sp_gy_salt ...
    c2002.sp_gy_salt c2003.sp_gy_salt ];
%bot merged
mer_sp_no3_b = [c1997.sp_gy_no3_b c1998.sp_gy_no3_b c1999.sp_gy_no3_b c2000.sp_gy_no3_b c2001.sp_gy_no3_b ...
    c2002.sp_gy_no3_b c2003.sp_gy_no3_b ];
mer_sp_nh4_b = [c1997.sp_gy_nh4_b c1998.sp_gy_nh4_b c1999.sp_gy_nh4_b c2000.sp_gy_nh4_b c2001.sp_gy_nh4_b ...
    c2002.sp_gy_nh4_b c2003.sp_gy_nh4_b ];
mer_sp_chl_b = [c1997.sp_gy_chl_b c1998.sp_gy_chl_b c1999.sp_gy_chl_b c2000.sp_gy_chl_b c2001.sp_gy_chl_b ...
    c2002.sp_gy_chl_b c2003.sp_gy_chl_b ];
mer_sp_temp_b = [c1997.sp_gy_temp_b c1998.sp_gy_temp_b c1999.sp_gy_temp_b c2000.sp_gy_temp_b c2001.sp_gy_temp_b ...
    c2002.sp_gy_temp_b c2003.sp_gy_temp_b ];
mer_sp_salt_b = [c1997.sp_gy_salt_b c1998.sp_gy_salt_b c1999.sp_gy_salt_b c2000.sp_gy_salt_b c2001.sp_gy_salt_b ...
    c2002.sp_gy_salt_b c2003.sp_gy_salt_b ];

for case_num = 1:3
    for ixx=1:2
        clearvars -except ixx case_num sigsig c1997 c1998 c1999 c2000 c2001 c2002 c2003 ltme *_eom dir_name mer_*
        
        sig=sigsig(ixx);

cd D:\장기생태\Dynamic\koem
% clim1=load('koem_climate_3regime_v2_2sig_gy_to_2003_12.mat'); % from 'plot_KOEM_st_GY_model_compare_horizontal_plt_climate_3regime_v2_gy_2sig.m'
% clim1=load('koem_climate_3regime_v2_to_2003_12.mat');  % from 'plot_KOEM_st_GY_model_compare_horizontal_plt_seasonal_climate_3regime_v2.m'
% clim1=load('koem_climate_3regime_v2_3sig_gy_to_2003_12.mat');  % from 'plot_KOEM_st_seasonal_climate_3regime_v3_spatial_sigma_extract.m'
clim1=load(['koem_climate_3regime_v2_',num2str(sig),'sig_gy_to_2003_12.mat']);

% 1997~2018 : 22yr
% 1st regime : 1997 ~ 2003.12 (from koem nh4 shift)
% 2nd regime : 2004.01 ~ 2009.12 (from koem no3 shift)
% 3rd regime : 2010.01 ~ 2018.12 (from koem no3 shift)


%make time tick
for i = 1:12
t_tIc{i} = {num2str(i,'%02d')};
if i == 1
    t_tIc{i} = {['2001.' num2str(i,'%02d')]};
end
end

for i = 1:12
t_tIc2{i} = {num2str(i,'%02d')};
if i == 1
    t_tIc2{i} = {['2002.' num2str(i,'%02d')]};
end
end


clearvars com_t_tic
% com_t_tic = {'2001.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2002.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2003.01','02','03','04','05','06','07','08','09','10','11','12',...}
%             '2004.01','02','03','04','05','06','07','08','09','10','11','12'}
com_t_tic = {'1997','1998','1999','2000','2001', ...
            '2002', ...
            '2003'};
        


for i =  1:365
%% model 
%std
    % time & spatial std
    clearvars size2
    size2=size(mer_sp_no3(:,i:365:end),2);
    regm_std_no3(i)=std(reshape(mer_sp_no3(:,i:365:end),1,9*size2));
    regm_std_nh4(i)=std(reshape(mer_sp_nh4(:,i:365:end),1,9*size2));
    regm_std_chl(i)=std(reshape(mer_sp_chl(:,i:365:end),1,9*size2));
    regm_std_temp(i)=std(reshape(mer_sp_temp(:,i:365:end),1,9*size2));
    regm_std_salt(i)=std(reshape(mer_sp_salt(:,i:365:end),1,9*size2));
    
    % time std
    spm_regm_std_no3(i)=std(nanmean(mer_sp_no3(:,i:365:end),1));
    spm_regm_std_nh4(i)=std(nanmean(mer_sp_nh4(:,i:365:end),1));
    spm_regm_std_chl(i)=std(nanmean(mer_sp_chl(:,i:365:end),1));
    spm_regm_std_temp(i)=std(nanmean(mer_sp_temp(:,i:365:end),1));
    spm_regm_std_salt(i)=std(nanmean(mer_sp_salt(:,i:365:end),1));
    
    % spatial std
    tm_regm_std_no3(i)=std(nanmean(mer_sp_no3(:,i:365:end),2));
    tm_regm_std_nh4(i)=std(nanmean(mer_sp_nh4(:,i:365:end),2));
    tm_regm_std_chl(i)=std(nanmean(mer_sp_chl(:,i:365:end),2));
    tm_regm_std_temp(i)=std(nanmean(mer_sp_temp(:,i:365:end),2));
    tm_regm_std_salt(i)=std(nanmean(mer_sp_salt(:,i:365:end),2));
    
  
    % time & spatial std
    clearvars size2
    size2=size(mer_sp_no3(:,i:365:end),2);
    regm_std_no3_b(i)=std(reshape(mer_sp_no3_b(:,i:365:end),1,9*size2));
    regm_std_nh4_b(i)=std(reshape(mer_sp_nh4_b(:,i:365:end),1,9*size2));
    regm_std_chl_b(i)=std(reshape(mer_sp_chl_b(:,i:365:end),1,9*size2));
    regm_std_temp_b(i)=std(reshape(mer_sp_temp_b(:,i:365:end),1,9*size2));
    regm_std_salt_b(i)=std(reshape(mer_sp_salt_b(:,i:365:end),1,9*size2));
    
    % time std
    spm_regm_std_no3_b(i)=std(nanmean(mer_sp_no3_b(:,i:365:end),1));
    spm_regm_std_nh4_b(i)=std(nanmean(mer_sp_nh4_b(:,i:365:end),1));
    spm_regm_std_chl_b(i)=std(nanmean(mer_sp_chl_b(:,i:365:end),1));
    spm_regm_std_temp_b(i)=std(nanmean(mer_sp_temp_b(:,i:365:end),1));
    spm_regm_std_salt_b(i)=std(nanmean(mer_sp_salt_b(:,i:365:end),1));
    
    % spatial std
    tm_regm_std_no3_b(i)=std(nanmean(mer_sp_no3_b(:,i:365:end),2));
    tm_regm_std_nh4_b(i)=std(nanmean(mer_sp_nh4_b(:,i:365:end),2));
    tm_regm_std_chl_b(i)=std(nanmean(mer_sp_chl_b(:,i:365:end),2));
    tm_regm_std_temp_b(i)=std(nanmean(mer_sp_temp_b(:,i:365:end),2));
    tm_regm_std_salt_b(i)=std(nanmean(mer_sp_salt_b(:,i:365:end),2));

end  

% %          dir_name ={'all_std','t_std','sp_std'};
 cd(['D:\장기생태\Dynamic\KOEM\compare_sewer\1st_3clim\1st_3clim_',num2str(sig),'sig_',dir_name{case_num}]);      
   
if case_num == 1 
   
    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines
    %data=NaN(1,365);        
    %data(1:length(c2004.gy_no3))=c2004.gy_no3; 
    %         plot(nanmean([c1997.gy_no3; c1998.gy_no3; c1999.gy_no3; c2000.gy_no3; c2001.gy_no3; c2002.gy_no3; c2003.gy_no3;],1),'b','linew',2);
            plot(nanmean(mer_no3,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_no3,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_no3;
            lower_bound_plt = temp_d - regm_std_no3;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_ft./14 + clim1.obs_std_gy_no3_ft./14;
            lower_bound_plt = clim1.obm_gy_no3_ft./14 - clim1.obs_std_gy_no3_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL NO3']);
            xlabel('time','fontsize',13)
            ylabel('NO3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4))=c2004.gy_nh4; 
            plot(nanmean(mer_nh4,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_nh4,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_nh4;
            lower_bound_plt = temp_d - regm_std_nh4;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4; c1998.gy_nh4; c1999.gy_nh4; c2000.gy_nh4; c2001.gy_nh4; c2002.gy_nh4; c2003.gy_nh4;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_ft./14 + clim1.obs_std_gy_nh4_ft./14;
            lower_bound_plt = clim1.obm_gy_nh4_ft./14 - clim1.obs_std_gy_nh4_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    % data=NaN(1,366);        
    % data(1:length(c2004.gy_chl))=c2004.gy_chl; 
            plot(nanmean(mer_chl,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_chl,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_chl;
            lower_bound_plt = temp_d - regm_std_chl;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl; c1998.gy_chl; c1999.gy_chl; c2000.gy_chl; c2001.gy_chl; c2002.gy_chl; c2003.gy_chl;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_ft + clim1.obs_std_gy_chl_ft;
            lower_bound_plt = clim1.obm_gy_chl_ft - clim1.obs_std_gy_chl_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp))=c2004.gy_temp;
            plot(nanmean(mer_temp,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_temp,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_temp;
            lower_bound_plt = temp_d - regm_std_temp;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp; c1998.gy_temp; c1999.gy_temp; c2000.gy_temp; c2001.gy_temp; c2002.gy_temp; c2003.gy_temp;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_ft + clim1.obs_std_gy_temp_ft;
            lower_bound_plt = clim1.obm_gy_temp_ft - clim1.obs_std_gy_temp_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt))=c2004.gy_salt; 
            plot(nanmean(mer_salt,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_d
            salt_d=nanmean(mer_salt,1);
            nonan_data_plt = find(isnan(salt_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_d + regm_std_salt;
            lower_bound_plt = salt_d - regm_std_salt;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt; c1998.gy_salt; c1999.gy_salt; c2000.gy_salt; c2001.gy_salt; c2002.gy_salt; c2003.gy_salt;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_ft + clim1.obs_std_gy_salt_ft;
            lower_bound_plt = clim1.obm_gy_salt_ft - clim1.obs_std_gy_salt_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')            
     %% bot

     fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_no3_b))=c2004.gy_no3_b; 
            plot(nanmean(mer_no3_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p no3_b_d
            no3_b_d=nanmean(mer_no3_b,1);
            nonan_data_plt = find(isnan(no3_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = no3_b_d + regm_std_no3_b;
            lower_bound_plt = no3_b_d - regm_std_no3_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_no3_b; c1998.gy_no3_b; c1999.gy_no3_b; c2000.gy_no3_b; c2001.gy_no3_b; c2002.gy_no3_b; c2003.gy_no3_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_b_ft./14 + clim1.obs_std_gy_no3_b_ft./14;
            lower_bound_plt = clim1.obm_gy_no3_b_ft./14 - clim1.obs_std_gy_no3_b_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_b_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_b_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
            xlabel('time','fontsize',13)
            ylabel('no3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4_b))=c2004.gy_nh4_b; 
            plot(nanmean(mer_nh4_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p nh4_b_d
            nh4_b_d=nanmean(mer_nh4_b,1);
            nonan_data_plt = find(isnan(nh4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = nh4_b_d + regm_std_nh4_b;
            lower_bound_plt = nh4_b_d - regm_std_nh4_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_b_ft./14 + clim1.obs_std_gy_nh4_b_ft./14;
            lower_bound_plt = clim1.obm_gy_nh4_b_ft./14 - clim1.obs_std_gy_nh4_b_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_b_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_b_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_chl_b))=c2004.gy_chl_b; 
            plot(nanmean(mer_chl_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p chl_b_d
            chl_b_d=nanmean(mer_chl_b,1);
            nonan_data_plt = find(isnan(chl_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = chl_b_d + regm_std_chl_b;
            lower_bound_plt = chl_b_d - regm_std_chl_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_b_ft + clim1.obs_std_gy_chl_b_ft;
            lower_bound_plt = clim1.obm_gy_chl_b_ft - clim1.obs_std_gy_chl_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl bot.']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp_b))=c2004.gy_temp_b; 
            plot(nanmean(mer_temp_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
            temp_b_d=nanmean(mer_temp_b,1);
            nonan_data_plt = find(isnan(temp_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_b_d + regm_std_temp_b;
            lower_bound_plt = temp_b_d - regm_std_temp_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_b_ft + clim1.obs_std_gy_temp_b_ft;
            lower_bound_plt = clim1.obm_gy_temp_b_ft - clim1.obs_std_gy_temp_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt_b))=c2004.gy_salt_b; 
            plot(nanmean(mer_salt_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
            salt_b_d=nanmean(mer_salt_b,1);
            nonan_data_plt = find(isnan(salt_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_b_d + regm_std_salt_b;
            lower_bound_plt = salt_b_d - regm_std_salt_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_b_ft + clim1.obs_std_gy_salt_b_ft;
            lower_bound_plt = clim1.obm_gy_salt_b_ft - clim1.obs_std_gy_salt_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng') 
            
            close all;
elseif case_num==2

        fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    %data(1:length(c2004.gy_no3))=c2004.gy_no3; 
    %         plot(nanmean([c1997.gy_no3; c1998.gy_no3; c1999.gy_no3; c2000.gy_no3; c2001.gy_no3; c2002.gy_no3; c2003.gy_no3;],1),'b','linew',2);
            plot(nanmean(mer_no3,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_no3,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_no3;
            lower_bound_plt = temp_d - spm_regm_std_no3;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_ft./14 + clim1.spm_obs_std_gy_no3_ft./14;
            lower_bound_plt = clim1.obm_gy_no3_ft./14 - clim1.spm_obs_std_gy_no3_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL NO3']);
            xlabel('time','fontsize',13)
            ylabel('NO3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4))=c2004.gy_nh4; 
            plot(nanmean(mer_nh4,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_nh4,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_nh4;
            lower_bound_plt = temp_d - spm_regm_std_nh4;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4; c1998.gy_nh4; c1999.gy_nh4; c2000.gy_nh4; c2001.gy_nh4; c2002.gy_nh4; c2003.gy_nh4;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_ft./14 + clim1.spm_obs_std_gy_nh4_ft./14;
            lower_bound_plt = clim1.obm_gy_nh4_ft./14 - clim1.spm_obs_std_gy_nh4_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    % data=NaN(1,366);        
    % data(1:length(c2004.gy_chl))=c2004.gy_chl; 
            plot(nanmean(mer_chl,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_chl,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_chl;
            lower_bound_plt = temp_d - spm_regm_std_chl;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl; c1998.gy_chl; c1999.gy_chl; c2000.gy_chl; c2001.gy_chl; c2002.gy_chl; c2003.gy_chl;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_ft + clim1.spm_obs_std_gy_chl_ft;
            lower_bound_plt = clim1.obm_gy_chl_ft - clim1.spm_obs_std_gy_chl_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp))=c2004.gy_temp;
            plot(nanmean(mer_temp,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_temp,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_temp;
            lower_bound_plt = temp_d - spm_regm_std_temp;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp; c1998.gy_temp; c1999.gy_temp; c2000.gy_temp; c2001.gy_temp; c2002.gy_temp; c2003.gy_temp;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_ft + clim1.spm_obs_std_gy_temp_ft;
            lower_bound_plt = clim1.obm_gy_temp_ft - clim1.spm_obs_std_gy_temp_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt))=c2004.gy_salt; 
            plot(nanmean(mer_salt,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_d
            salt_d=nanmean(mer_salt,1);
            nonan_data_plt = find(isnan(salt_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_d + spm_regm_std_salt;
            lower_bound_plt = salt_d - spm_regm_std_salt;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt; c1998.gy_salt; c1999.gy_salt; c2000.gy_salt; c2001.gy_salt; c2002.gy_salt; c2003.gy_salt;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_ft + clim1.spm_obs_std_gy_salt_ft;
            lower_bound_plt = clim1.obm_gy_salt_ft - clim1.spm_obs_std_gy_salt_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')            
     %% bot

     fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_no3_b))=c2004.gy_no3_b; 
            plot(nanmean(mer_no3_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p no3_b_d
            no3_b_d=nanmean(mer_no3_b,1);
            nonan_data_plt = find(isnan(no3_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = no3_b_d + spm_regm_std_no3_b;
            lower_bound_plt = no3_b_d - spm_regm_std_no3_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_no3_b; c1998.gy_no3_b; c1999.gy_no3_b; c2000.gy_no3_b; c2001.gy_no3_b; c2002.gy_no3_b; c2003.gy_no3_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_b_ft./14 + clim1.spm_obs_std_gy_no3_b_ft./14;
            lower_bound_plt = clim1.obm_gy_no3_b_ft./14 - clim1.spm_obs_std_gy_no3_b_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_b_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_b_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
            xlabel('time','fontsize',13)
            ylabel('no3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4_b))=c2004.gy_nh4_b; 
            plot(nanmean(mer_nh4_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p nh4_b_d
            nh4_b_d=nanmean(mer_nh4_b,1);
            nonan_data_plt = find(isnan(nh4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = nh4_b_d + spm_regm_std_nh4_b;
            lower_bound_plt = nh4_b_d - spm_regm_std_nh4_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_b_ft./14 + clim1.spm_obs_std_gy_nh4_b_ft./14;
            lower_bound_plt = clim1.obm_gy_nh4_b_ft./14 - clim1.spm_obs_std_gy_nh4_b_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_b_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_b_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_chl_b))=c2004.gy_chl_b; 
            plot(nanmean(mer_chl_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p chl_b_d
            chl_b_d=nanmean(mer_chl_b,1);
            nonan_data_plt = find(isnan(chl_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = chl_b_d + spm_regm_std_chl_b;
            lower_bound_plt = chl_b_d - spm_regm_std_chl_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_b_ft + clim1.spm_obs_std_gy_chl_b_ft;
            lower_bound_plt = clim1.obm_gy_chl_b_ft - clim1.spm_obs_std_gy_chl_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl bot.']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp_b))=c2004.gy_temp_b; 
            plot(nanmean(mer_temp_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
            temp_b_d=nanmean(mer_temp_b,1);
            nonan_data_plt = find(isnan(temp_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_b_d + spm_regm_std_temp_b;
            lower_bound_plt = temp_b_d - spm_regm_std_temp_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_b_ft + clim1.spm_obs_std_gy_temp_b_ft;
            lower_bound_plt = clim1.obm_gy_temp_b_ft - clim1.spm_obs_std_gy_temp_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt_b))=c2004.gy_salt_b; 
            plot(nanmean(mer_salt_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
            salt_b_d=nanmean(mer_salt_b,1);
            nonan_data_plt = find(isnan(salt_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_b_d + spm_regm_std_salt_b;
            lower_bound_plt = salt_b_d - spm_regm_std_salt_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_b_ft + clim1.spm_obs_std_gy_salt_b_ft;
            lower_bound_plt = clim1.obm_gy_salt_b_ft - clim1.spm_obs_std_gy_salt_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng')
            
            close all;
elseif case_num == 3
        fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    %data(1:length(c2004.gy_no3))=c2004.gy_no3; 
    %         plot(nanmean([c1997.gy_no3; c1998.gy_no3; c1999.gy_no3; c2000.gy_no3; c2001.gy_no3; c2002.gy_no3; c2003.gy_no3;],1),'b','linew',2);
            plot(nanmean(mer_no3,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_no3,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_no3;
            lower_bound_plt = temp_d - tm_regm_std_no3;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_ft./14 + clim1.tm_obs_std_gy_no3_ft./14;
            lower_bound_plt = clim1.obm_gy_no3_ft./14 - clim1.tm_obs_std_gy_no3_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL NO3']);
            xlabel('time','fontsize',13)
            ylabel('NO3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4))=c2004.gy_nh4; 
            plot(nanmean(mer_nh4,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_nh4,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_nh4;
            lower_bound_plt = temp_d - tm_regm_std_nh4;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4; c1998.gy_nh4; c1999.gy_nh4; c2000.gy_nh4; c2001.gy_nh4; c2002.gy_nh4; c2003.gy_nh4;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_ft./14 + clim1.tm_obs_std_gy_nh4_ft./14;
            lower_bound_plt = clim1.obm_gy_nh4_ft./14 - clim1.tm_obs_std_gy_nh4_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    % data=NaN(1,366);        
    % data(1:length(c2004.gy_chl))=c2004.gy_chl; 
            plot(nanmean(mer_chl,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_chl,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_chl;
            lower_bound_plt = temp_d - tm_regm_std_chl;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl; c1998.gy_chl; c1999.gy_chl; c2000.gy_chl; c2001.gy_chl; c2002.gy_chl; c2003.gy_chl;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_ft + clim1.tm_obs_std_gy_chl_ft;
            lower_bound_plt = clim1.obm_gy_chl_ft - clim1.tm_obs_std_gy_chl_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp))=c2004.gy_temp;
            plot(nanmean(mer_temp,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_temp,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_temp;
            lower_bound_plt = temp_d - tm_regm_std_temp;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp; c1998.gy_temp; c1999.gy_temp; c2000.gy_temp; c2001.gy_temp; c2002.gy_temp; c2003.gy_temp;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_ft + clim1.tm_obs_std_gy_temp_ft;
            lower_bound_plt = clim1.obm_gy_temp_ft - clim1.tm_obs_std_gy_temp_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt))=c2004.gy_salt; 
            plot(nanmean(mer_salt,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_d
            salt_d=nanmean(mer_salt,1);
            nonan_data_plt = find(isnan(salt_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_d + tm_regm_std_salt;
            lower_bound_plt = salt_d - tm_regm_std_salt;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt; c1998.gy_salt; c1999.gy_salt; c2000.gy_salt; c2001.gy_salt; c2002.gy_salt; c2003.gy_salt;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_ft + clim1.tm_obs_std_gy_salt_ft;
            lower_bound_plt = clim1.obm_gy_salt_ft - clim1.tm_obs_std_gy_salt_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')            
     %% bot

     fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_no3_b))=c2004.gy_no3_b; 
            plot(nanmean(mer_no3_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p no3_b_d
            no3_b_d=nanmean(mer_no3_b,1);
            nonan_data_plt = find(isnan(no3_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = no3_b_d + tm_regm_std_no3_b;
            lower_bound_plt = no3_b_d - tm_regm_std_no3_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_no3_b; c1998.gy_no3_b; c1999.gy_no3_b; c2000.gy_no3_b; c2001.gy_no3_b; c2002.gy_no3_b; c2003.gy_no3_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_b_ft./14 + clim1.tm_obs_std_gy_no3_b_ft./14;
            lower_bound_plt = clim1.obm_gy_no3_b_ft./14 - clim1.tm_obs_std_gy_no3_b_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_b_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_b_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
            xlabel('time','fontsize',13)
            ylabel('no3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4_b))=c2004.gy_nh4_b; 
            plot(nanmean(mer_nh4_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p nh4_b_d
            nh4_b_d=nanmean(mer_nh4_b,1);
            nonan_data_plt = find(isnan(nh4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = nh4_b_d + tm_regm_std_nh4_b;
            lower_bound_plt = nh4_b_d - tm_regm_std_nh4_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_b_ft./14 + clim1.tm_obs_std_gy_nh4_b_ft./14;
            lower_bound_plt = clim1.obm_gy_nh4_b_ft./14 - clim1.tm_obs_std_gy_nh4_b_ft./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_b_ft./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_b_ft),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_chl_b))=c2004.gy_chl_b; 
            plot(nanmean(mer_chl_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p chl_b_d
            chl_b_d=nanmean(mer_chl_b,1);
            nonan_data_plt = find(isnan(chl_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = chl_b_d + tm_regm_std_chl_b;
            lower_bound_plt = chl_b_d - tm_regm_std_chl_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_b_ft + clim1.tm_obs_std_gy_chl_b_ft;
            lower_bound_plt = clim1.obm_gy_chl_b_ft - clim1.tm_obs_std_gy_chl_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl bot.']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp_b))=c2004.gy_temp_b; 
            plot(nanmean(mer_temp_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
            temp_b_d=nanmean(mer_temp_b,1);
            nonan_data_plt = find(isnan(temp_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_b_d + tm_regm_std_temp_b;
            lower_bound_plt = temp_b_d - tm_regm_std_temp_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_b_ft + clim1.tm_obs_std_gy_temp_b_ft;
            lower_bound_plt = clim1.obm_gy_temp_b_ft - clim1.tm_obs_std_gy_temp_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt_b))=c2004.gy_salt_b; 
            plot(nanmean(mer_salt_b,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
            salt_b_d=nanmean(mer_salt_b,1);
            nonan_data_plt = find(isnan(salt_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_b_d + tm_regm_std_salt_b;
            lower_bound_plt = salt_b_d - tm_regm_std_salt_b;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_b_ft)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_b_ft + clim1.tm_obs_std_gy_salt_b_ft;
            lower_bound_plt = clim1.obm_gy_salt_b_ft - clim1.tm_obs_std_gy_salt_b_ft;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_b_ft,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_b_ft),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_b_ft),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng') 
            
            close all;
end
    end
end