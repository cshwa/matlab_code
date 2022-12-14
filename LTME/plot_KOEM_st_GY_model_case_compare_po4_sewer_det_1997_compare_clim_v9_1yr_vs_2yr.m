close all; clear;clc

sig = 3;

cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
ltme =load('LTME_observation_data_v2.mat'); %LTME obs.

cd D:\장기생태\Dynamic\koem
% clim1=load(['koem_climate_10yr_v2_3sig_gy.mat']);
clim1=load(['koem_climate_10yr_v5_',num2str(sig),'sig_gy_each_st.mat']); % full (9point) data (extract_sigma_for_each_st)
% clim1=load(['koem_climate_10yr_v5_2sig_gy_8point.mat'],'obs*','obm*','spm_obs*','tm_obs*'); %extract gwangyang-port point

% cd D:\장기생태\Dynamic\KOEM\gy_2001\gy_2001_koem_daily\spmean
% load('koem_result_processed_std_spmean.mat'); % KOEM obs.
% % mpdata= load('mpdata_result.mat');

% cd D:\장기생태\Dynamic\result\1997
% defau = load('1997_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
% tunn1=load('1997_tunn1.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
% tunn2 = load('1997_tunn2.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
% tunn3 = load('1997_tunn3.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
% tunn4=load('1997_tunn4.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
% tunn5 = load('1997_tunn5.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
% tunn6 = load('1997_tunn6.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
% tunn7=load('1997_tunn7.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
cd F:\장기생태_output_임시
% tunn1=load('1997_sewer_re_1yr.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
% tunn4=load('1997_sewer_re_2yr.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');
defau=load('1997_sewer_re_v9grid_t1.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*'); %1yr
% tunn5=load('1997_sewer_re_v9_1yr.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*'); % origin
tunn4=load('1997_sewer_re_v9grid_t1_2yr.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*'); %2yr
% tunn6=load('1997_sewer_re_v9_2_1yr.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*');

% cd D:\장기생태\Dynamic\KOEM\gy_2001\phymin_half_saturation
% phymin= load('phymin_result.mat');
% 
% cd D:\장기생태\Dynamic\KOEM\gy_2001\mp_p_ws_half_saturation
% mp_p_ws= load('mp_p_ws_result.mat');

% % mpdata.gy_temp
% cd D:\장기생태\Dynamic\KOEM\compare_sewer\1997
if sig ==2
% cd(['F:\장기생태_output_임시\1997\',num2str(sig),'sig'])  %origin
else
%     cd(['F:\장기생태_output_임시\1997\']) %origin
cd F:\장기생태_output_임시\1997\1yr_vs_2yr
end

%% GY plot zoo
% fig = figure; hold on;
%         plot(defau.gy_zoo,'b','linew',2);
%         %plot(tunn6.gy_zoo(2:end),'r','linew',2);
%          plot(defau.gy_phy,'m','linew',2);
%         %plot(tunn6.gy_phy(2:end),'g','linew',2);
%         legend('zoo','zoo_w_t_8','phy','phy_w_t_8')
%         xlim([1 365]);
%         title(['GY 1997 daily KOEM OBS vs. MODEL Zoo & Phy surface']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('concentration (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
% 
%         
% fig = figure; hold on;      
%         plot(defau.gy_zoo_b,'b','linew',2);
%         %plot(tunn6.gy_zoo_b(2:end),'r','linew',2);
%         plot(defau.gy_phy_b,'m','linew',2);
%         %plot(tunn6.gy_phy_b(2:end),'g','linew',2);
%         legend('zoo','zoo_w_t_8','phy','phy_w_t_8')
%         xlim([1 365]);
%         title(['GY 1997 daily KOEM OBS vs. MODEL Zoo & Phy bottom']);
%         xlabel('time(days on 2001)','fontsize',13)
%         ylabel('concentration (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)              
        
%% GY plot
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_po4,'b','linew',2);
        %plot(tunn1.gy_po4,'k','linew',2);
%         plot(tunn2.gy_po4,'g','linew',2);
        plot(tunn4.gy_po4,'m','linew',2);
        %plot(tunn5.gy_po4,'c','linew',2);
        %plot(tunn6.gy_po4,'color',[.5 .5 .5],'linew',2);
%         %plot(tunn7.gy_po4,'g','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_po4_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_po4_ft_1./30.973762 + clim1.obs_std_gy_po4_ft_1./30.973762;
        lower_bound_plt = clim1.obm_gy_po4_ft_1./30.973762 - clim1.obs_std_gy_po4_ft_1./30.973762;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_po4_ft_1./30.973762,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_po4_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_po4_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency   
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        
        title(['GY KOEM OBS vs. daily MODEL po4']);
        xlabel('time','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 3.5])
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_po4'),'-dpng')




fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_no3,'b','linew',2);
        %plot(tunn1.gy_no3,'k','linew',2);
%         plot(tunn2.gy_no3,'g','linew',2);
        plot(tunn4.gy_no3,'m','linew',2);
        %plot(tunn5.gy_no3,'c','linew',2);
        %plot(tunn6.gy_no3,'color',[.5 .5 .5],'linew',2);
%         %plot(tunn7.gy_no3,'g','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_no3_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_no3_ft_1./14 + clim1.obs_std_gy_no3_ft_1./14;
        lower_bound_plt = clim1.obm_gy_no3_ft_1./14 - clim1.obs_std_gy_no3_ft_1./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_no3_ft_1./14,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_no3_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_no3_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency   
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        
        title(['GY KOEM OBS vs. daily MODEL NO3']);
        xlabel('time','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 70])
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_nh4,'b','linew',2);
        %plot(tunn1.gy_nh4,'k','linew',2);
%         plot(tunn2.gy_nh4,'g','linew',2);
        plot(tunn4.gy_nh4,'m','linew',2);
        %plot(tunn5.gy_nh4,'c','linew',2);
        %plot(tunn6.gy_nh4,'color',[.5 .5 .5],'linew',2);
%         %plot(tunn7.gy_nh4,'g','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_nh4_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_nh4_ft_1./14 + clim1.obs_std_gy_nh4_ft_1./14;
        lower_bound_plt = clim1.obm_gy_nh4_ft_1./14 - clim1.obs_std_gy_nh4_ft_1./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_nh4_ft_1./14,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_nh4_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_nh4_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency    
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        
        title(['GY KOEM OBS vs. daily MODEL nh4']);
        xlabel('time','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 14])
        yticks(1:1:14)
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_chl,'b','linew',2);
        %plot(tunn1.gy_chl,'k','linew',2);
%         plot(tunn2.gy_chl,'g','linew',2);
        plot(tunn4.gy_chl,'m','linew',2);
        %plot(tunn5.gy_chl,'c','linew',2);
        %plot(tunn6.gy_chl,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_chl,'g','linew',2);
             
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_chl_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_chl_ft_1 + clim1.obs_std_gy_chl_ft_1;
        lower_bound_plt = clim1.obm_gy_chl_ft_1 - clim1.obs_std_gy_chl_ft_1;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_chl_ft_1,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_chl_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_chl_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency 
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        
        title(['GY KOEM OBS vs. daily MODEL chl']);
        xlabel('time','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        yticks(1:1:15)
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_temp,'b','linew',2);
        %plot(tunn1.gy_temp,'k','linew',2);
%         plot(tunn2.gy_temp,'g','linew',2);
        plot(tunn4.gy_temp,'m','linew',2);
        %plot(tunn5.gy_temp,'c','linew',2);
        %plot(tunn6.gy_temp,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_temp,'g','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_temp_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_temp_ft_1 + clim1.obs_std_gy_temp_ft_1;
        lower_bound_plt = clim1.obm_gy_temp_ft_1 - clim1.obs_std_gy_temp_ft_1;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_temp_ft_1,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_temp_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_temp_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency    
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        
        title(['GY KOEM OBS vs. daily MODEL temp']);
        xlabel('time','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        ylim([0 35])
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_salt,'b','linew',2);
        %plot(tunn1.gy_salt,'k','linew',2);
%         plot(tunn2.gy_salt,'g','linew',2);
        plot(tunn4.gy_salt,'m','linew',2);
        %plot(tunn5.gy_salt,'c','linew',2);
        %plot(tunn6.gy_salt,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_salt,'g','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_salt_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_salt_ft_1 + clim1.obs_std_gy_salt_ft_1;
        lower_bound_plt = clim1.obm_gy_salt_ft_1 - clim1.obs_std_gy_salt_ft_1;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_salt_ft_1,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_salt_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_salt_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency    
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        title(['GY KOEM OBS vs. daily MODEL salt']);
        xlabel('time','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        ylim([0 35])
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')            
 %% bot
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_po4_b,'b','linew',2);
        %plot(tunn1.gy_po4_b,'k','linew',2);
%         plot(tunn2.gy_po4,'g','linew',2);
        plot(tunn4.gy_po4_b,'m','linew',2);
        %plot(tunn5.gy_po4_b,'c','linew',2);
        %plot(tunn6.gy_po4_b,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_po4_b,'g','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_po4_b_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 + clim1.obs_std_gy_po4_b_ft_1./30.973762;
        lower_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 - clim1.obs_std_gy_po4_b_ft_1./30.973762;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_po4_b_ft_1./30.973762,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_po4_b_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_po4_b_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        title(['GY KOEM OBS vs. daily MODEL po4 bot.']);
        xlabel('time','fontsize',13)
        ylabel('po4 bot (umol/m^3)','fontsize',13)
        grid on
        ylim([0 3.5])
        set(gca,'fontsize',13)
        xlim([1 365])
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_po4_bot'),'-dpng')
 
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_no3_b,'b','linew',2);
        %plot(tunn1.gy_no3_b,'k','linew',2);
%         plot(tunn2.gy_no3,'g','linew',2);
        plot(tunn4.gy_no3_b,'m','linew',2);
        %plot(tunn5.gy_no3_b,'c','linew',2);
        %plot(tunn6.gy_no3_b,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_no3_b,'g','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_no3_b_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_no3_b_ft_1./14 + clim1.obs_std_gy_no3_b_ft_1./14;
        lower_bound_plt = clim1.obm_gy_no3_b_ft_1./14 - clim1.obs_std_gy_no3_b_ft_1./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_no3_b_ft_1./14,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_no3_b_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_no3_b_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
        xlabel('time','fontsize',13)
        ylabel('no3 bot (umol/m^3)','fontsize',13)
        grid on
        ylim([0 70])
        set(gca,'fontsize',13)
        xlim([1 365])
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_nh4_b,'b','linew',2);
        %plot(tunn1.gy_nh4_b,'k','linew',2);
%         plot(tunn2.gy_nh4,'g','linew',2);
        plot(tunn4.gy_nh4_b,'m','linew',2);
        %plot(tunn5.gy_nh4_b,'c','linew',2);
        %plot(tunn6.gy_nh4_b,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_nh4_b,'g','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_nh4_b_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_nh4_b_ft_1./14 + clim1.obs_std_gy_nh4_b_ft_1./14;
        lower_bound_plt = clim1.obm_gy_nh4_b_ft_1./14 - clim1.obs_std_gy_nh4_b_ft_1./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_nh4_b_ft_1./14,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_nh4_b_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_nh4_b_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
        xlabel('time','fontsize',13)
        ylabel('nh4 bot (umol/m^3)','fontsize',13)
        grid on
        ylim([0 14])
        yticks(1:1:14)
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_chl_b,'b','linew',2);
        %plot(tunn1.gy_chl_b,'k','linew',2);
%         plot(tunn2.gy_chl,'g','linew',2);
        plot(tunn4.gy_chl_b,'m','linew',2);
        %plot(tunn5.gy_chl_b,'c','linew',2);
        %plot(tunn6.gy_chl_b,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_chl_b,'g','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_chl_b_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_chl_b_ft_1 + clim1.obs_std_gy_chl_b_ft_1;
        lower_bound_plt = clim1.obm_gy_chl_b_ft_1 - clim1.obs_std_gy_chl_b_ft_1;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_chl_b_ft_1,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_chl_b_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_chl_b_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency   
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        title(['GY KOEM OBS vs. daily MODEL chl bot.']);
        xlabel('time','fontsize',13)
        ylabel('chl bot (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        xlim([1 365])
        yticks(1:1:15)
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_temp_b,'b','linew',2);
        %plot(tunn1.gy_temp_b,'k','linew',2);
%         plot(tunn2.gy_temp,'g','linew',2);
        plot(tunn4.gy_temp_b,'m','linew',2);
        %plot(tunn5.gy_temp_b,'c','linew',2);
        %plot(tunn6.gy_temp_b,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_temp_b,'g','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_temp_b_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_temp_b_ft_1 + clim1.obs_std_gy_temp_b_ft_1;
        lower_bound_plt = clim1.obm_gy_temp_b_ft_1 - clim1.obs_std_gy_temp_b_ft_1;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_temp_b_ft_1,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_temp_b_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_temp_b_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency        
        legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        title(['GY KOEM OBS vs. daily MODEL temp bot.']);
        xlabel('time','fontsize',13)
        ylabel('temp bot (o^C)','fontsize',13)
        grid on
        ylim([0 35])
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%% 2001
        plot(defau.gy_salt_b,'b','linew',2);
        %plot(tunn1.gy_salt_b,'k','linew',2);
%         plot(tunn2.gy_salt,'g','linew',2);
        plot(tunn4.gy_salt_b,'m','linew',2);
        %plot(tunn5.gy_salt_b,'c','linew',2);
        %plot(tunn6.gy_salt_b,'color',[.5 .5 .5],'linew',2);
        %plot(tunn7.gy_salt_b,'g','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(clim1.obm_gy_salt_b_ft_1)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = clim1.obm_gy_salt_b_ft_1 + clim1.obs_std_gy_salt_b_ft_1;
        lower_bound_plt = clim1.obm_gy_salt_b_ft_1 - clim1.obs_std_gy_salt_b_ft_1;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(clim1.obm_gy_salt_b_ft_1,'r-','linew',2); 
        plot(1:length(clim1.obm_gy_salt_b_ft_1),upper_bound_plt,'r-','linew',2);
        plot(1:length(clim1.obm_gy_salt_b_ft_1),lower_bound_plt,'r-','linew',2);
        
        alpha(0.3) %transparency
                legend('contr.','t4'); %legend('contr.','t4','t5'); %legend('contr.','t1','t4','t5','t6'); %legend('contr.','t1','t4','t5','t6','t7'); 
        title(['GY KOEM OBS vs. daily MODEL salt bot.']);
        xlabel('time','fontsize',13)
        ylabel('salt bot (psu)','fontsize',13)
        grid on
        ylim([0 35])
        xlim([1 365])
        set(gca,'fontsize',13)
        print(fig,strcat('1997_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng')            
return        
        
%%south gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_no3,'b','linew',2);
        plot(pcase.sgy_no3,'k','linew',2);
        plot(mp_p_ws.sgy_no3,'g','linew',2); 
        plot(obm_sgy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 + obs_std_sgy_no3./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 - obs_std_sgy_no3./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_nh4,'b','linew',2);
        plot(pcase.sgy_nh4,'k','linew',2);
        plot(mp_p_ws.sgy_nh4,'g','linew',2);
        plot(obm_sgy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 + obs_std_sgy_nh4./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 - obs_std_sgy_nh4./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_south_sgy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_chl,'b','linew',2);
        plot(pcase.sgy_chl,'k','linew',2);
        plot(mp_p_ws.sgy_chl,'g','linew',2);
        plot(obm_sgy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl),obm_sgy_chl + obs_std_sgy_chl,'m-','linew',2);
        plot(1:length(obm_sgy_chl),obm_sgy_chl - obs_std_sgy_chl,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_temp,'b','linew',2);
        plot(pcase.sgy_temp,'k','linew',2);
        plot(mp_p_ws.sgy_temp,'g','linew',2);
        plot(obm_sgy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp),obm_sgy_temp + obs_std_sgy_temp,'m-','linew',2);
        plot(1:length(obm_sgy_temp),obm_sgy_temp - obs_std_sgy_temp,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_salt,'b','linew',2);
        plot(pcase.sgy_salt,'k','linew',2);
        plot(mp_p_ws.sgy_salt,'g','linew',2);   
        plot(obm_sgy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt),obm_sgy_salt + obs_std_sgy_salt,'m-','linew',2);
        plot(1:length(obm_sgy_salt),obm_sgy_salt - obs_std_sgy_salt,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_south_sgy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_no3_b,'b','linew',2);
        plot(pcase.sgy_no3_b,'k','linew',2);
        plot(mp_p_ws.sgy_no3_b,'g','linew',2);
        plot(obm_sgy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 + obs_std_sgy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 - obs_std_sgy_no3_b./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_nh4_b,'b','linew',2);
        plot(pcase.sgy_nh4_b,'k','linew',2);
        plot(mp_p_ws.sgy_nh4_b,'g','linew',2); 
        plot(obm_sgy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 + obs_std_sgy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 - obs_std_sgy_nh4_b./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_south_sgy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_chl_b,'b','linew',2);
        plot(pcase.sgy_chl_b,'k','linew',2);
        plot(mp_p_ws.sgy_chl_b,'g','linew',2);  
        plot(obm_sgy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b + obs_std_sgy_chl_b,'m-','linew',2);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b - obs_std_sgy_chl_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_temp_b,'b','linew',2);
        plot(pcase.sgy_temp_b,'k','linew',2);
        plot(mp_p_ws.sgy_temp_b,'g','linew',2);        
        plot(obm_sgy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b + obs_std_sgy_temp_b,'m-','linew',2);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b - obs_std_sgy_temp_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_salt_b,'b','linew',2);
        plot(pcase.sgy_salt_b,'k','linew',2);
        plot(mp_p_ws.sgy_salt_b,'g','linew',2);    
        plot(obm_sgy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b + obs_std_sgy_salt_b,'m-','linew',2);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b - obs_std_sgy_salt_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_south_sgy'),'-dpng')    
        
%%east gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_no3,'b','linew',2);
        plot(pcase.egy_no3,'k','linew',2);
        plot(mp_p_ws.egy_no3,'g','linew',2);  
        plot(obm_egy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 + obs_std_egy_no3./14,'m-','linew',2);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 - obs_std_egy_no3./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines       
        plot(mpdata.egy_nh4,'b','linew',2);
        plot(pcase.egy_nh4,'k','linew',2);
        plot(mp_p_ws.egy_nh4,'g','linew',2); 
        plot(obm_egy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 + obs_std_egy_nh4./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 - obs_std_egy_nh4./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_egy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_chl,'b','linew',2);
        plot(pcase.egy_chl,'k','linew',2);
        plot(mp_p_ws.egy_chl,'g','linew',2); 
        plot(obm_egy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl),obm_egy_chl + obs_std_egy_chl,'m-','linew',2);
        plot(1:length(obm_egy_chl),obm_egy_chl - obs_std_egy_chl,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_temp,'b','linew',2);
        plot(pcase.egy_temp,'k','linew',2);
        plot(mp_p_ws.egy_temp,'g','linew',2); 
        plot(obm_egy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp),obm_egy_temp + obs_std_egy_temp,'m-','linew',2);
        plot(1:length(obm_egy_temp),obm_egy_temp - obs_std_egy_temp,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_salt,'b','linew',2);
        plot(pcase.egy_salt,'k','linew',2);
        plot(mp_p_ws.egy_salt,'g','linew',2); 
        plot(obm_egy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt),obm_egy_salt + obs_std_egy_salt,'m-','linew',2);
        plot(1:length(obm_egy_salt),obm_egy_salt - obs_std_egy_salt,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_egy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_no3_b,'b','linew',2);
        plot(pcase.egy_no3_b,'k','linew',2);
        plot(mp_p_ws.egy_no3_b,'g','linew',2); 
        plot(obm_egy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 + obs_std_egy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 - obs_std_egy_no3_b./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_nh4_b,'b','linew',2);
        plot(pcase.egy_nh4_b,'k','linew',2);
        plot(mp_p_ws.egy_nh4_b,'g','linew',2); 
        plot(obm_egy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 + obs_std_egy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 - obs_std_egy_nh4_b./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_egy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_chl_b,'b','linew',2);
        plot(pcase.egy_chl_b,'k','linew',2);
        plot(mp_p_ws.egy_chl_b,'g','linew',2); 
        plot(obm_egy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b + obs_std_egy_chl_b,'m-','linew',2);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b - obs_std_egy_chl_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_temp_b,'b','linew',2);
        plot(pcase.egy_temp_b,'k','linew',2);
        plot(mp_p_ws.egy_temp_b,'g','linew',2); 
        plot(obm_egy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b + obs_std_egy_temp_b,'m-','linew',2);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b - obs_std_egy_temp_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_salt_b,'b','linew',2);
        plot(pcase.egy_salt_b,'k','linew',2);
        plot(mp_p_ws.egy_salt_b,'g','linew',2); 
        plot(obm_egy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b + obs_std_egy_salt_b,'m-','linew',2);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b - obs_std_egy_salt_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_egy'),'-dpng')    
        
%% jinju

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_no3,'b','linew',2);
        plot(pcase.jj_no3,'k','linew',2);
        plot(mp_p_ws.jj_no3,'g','linew',2); 
        plot(obm_jj_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 + obs_std_jj_no3./14,'m-','linew',2);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 - obs_std_jj_no3./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_nh4,'b','linew',2);
        plot(pcase.jj_nh4,'k','linew',2);
        plot(mp_p_ws.jj_nh4,'g','linew',2); 
        plot(obm_jj_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 + obs_std_jj_nh4./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 - obs_std_jj_nh4./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_chl,'b','linew',2);
        plot(pcase.jj_chl,'k','linew',2);
        plot(mp_p_ws.jj_chl,'g','linew',2); 
        plot(obm_jj_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl),obm_jj_chl + obs_std_jj_chl,'m-','linew',2);
        plot(1:length(obm_jj_chl),obm_jj_chl - obs_std_jj_chl,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_temp,'b','linew',2);
        plot(pcase.jj_temp,'k','linew',2);
        plot(mp_p_ws.jj_temp,'g','linew',2); 
        plot(obm_jj_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp),obm_jj_temp + obs_std_jj_temp,'m-','linew',2);
        plot(1:length(obm_jj_temp),obm_jj_temp - obs_std_jj_temp,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_salt,'b','linew',2);
        plot(pcase.jj_salt,'k','linew',2);
        plot(mp_p_ws.jj_salt,'g','linew',2); 
        plot(obm_jj_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt),obm_jj_salt + obs_std_jj_salt,'m-','linew',2);
        plot(1:length(obm_jj_salt),obm_jj_salt - obs_std_jj_salt,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_jj'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_no3_b,'b','linew',2);
        plot(pcase.jj_no3_b,'k','linew',2);
        plot(mp_p_ws.jj_no3_b,'g','linew',2); 
        plot(obm_jj_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 + obs_std_jj_no3_b./14,'m-','linew',2);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 - obs_std_jj_no3_b./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_nh4_b,'b','linew',2);
        plot(pcase.jj_nh4_b,'k','linew',2);
        plot(mp_p_ws.jj_nh4_b,'g','linew',2); 
        plot(obm_jj_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 + obs_std_jj_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 - obs_std_jj_nh4_b./14,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_chl_b,'b','linew',2);
        plot(pcase.jj_chl_b,'k','linew',2);
        plot(mp_p_ws.jj_chl_b,'g','linew',2); 
        plot(obm_jj_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b + obs_std_jj_chl_b,'m-','linew',2);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b - obs_std_jj_chl_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_temp_b,'b','linew',2);
        plot(pcase.jj_temp_b,'k','linew',2);
        plot(mp_p_ws.jj_temp_b,'g','linew',2); 
        plot(obm_jj_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b + obs_std_jj_temp_b,'m-','linew',2);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b - obs_std_jj_temp_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_salt_b,'b','linew',2);
        plot(pcase.jj_salt_b,'k','linew',2);
        plot(mp_p_ws.jj_salt_b,'g','linew',2); 
        plot(obm_jj_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b + obs_std_jj_salt_b,'m-','linew',2);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b - obs_std_jj_salt_b,'m-','linew',2);
        title(['GY 1997 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_jj'),'-dpng')    


