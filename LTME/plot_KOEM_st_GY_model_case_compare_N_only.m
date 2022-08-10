close all; clear;clc

% cd D:\장기생태\Dynamic\KOEM\gy_2001\detclo_p
% pcase=load('detclo_p_result.mat');

cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
load LTME_observation_data_v2.mat %LTME obs.

cd D:\장기생태\Dynamic\KOEM\gy_2001\gy_2001_koem_daily\spmean
koem=load('koem_result_processed_std_spmean.mat'); % KOEM obs.
% mpdata= load('mpdata_result.mat');

% cd D:\장기생태\Dynamic\KOEM\compare_sewer\det
% mp_p_sewer_det = load('mp_p_sewer_det_result.mat');

cd D:\장기생태\Dynamic\KOEM\gy_2001\gy_2001_koem_daily\spmean
mpdata= load('mpdata_result.mat');

% cd D:\장기생태\Dynamic\KOEM\compare_sewer\
% mp_p_sewer= load('mp_p_sewer_result.mat');

% cd D:\장기생태\Dynamic\KOEM\gy_2001\phymin_half_saturation
% phymin= load('phymin_result.mat');
% 
% cd D:\장기생태\Dynamic\KOEM\gy_2001\mp_p_ws_half_saturation
% mp_p_ws= load('mp_p_ws_result.mat');

% % mpdata.gy_temp
cd D:\장기생태\Dynamic\KOEM\compare_sewer\N_only


fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_no3,'b','linew',2);
%         plot(pcase.gy_no3,'k','linew',2);
%         plot(mp_p_sewer.gy_no3,'g','linew',2);
%         plot(mp_p_sewer_det.gy_no3,'m','linew',2);
        xlim([1 365]);
      clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_no3./14 + koem.obs_std_gy_no3./14;
        lower_bound_plt = koem.obm_gy_no3./14 - koem.obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_no3./14,'r-','linew',2); 
        plot(1:length(koem.obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_no3),lower_bound_plt,'r-','linew',2);        
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 75])
        alpha(0.3)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_gy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_nh4,'b','linew',2);
%         plot(pcase.gy_nh4,'k','linew',2);
%         plot(mp_p_sewer.gy_nh4,'g','linew',2);
%         plot(mp_p_sewer_det.gy_nh4,'m','linew',2);
         xlim([1 365]);
       clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_nh4./14 + koem.obs_std_gy_nh4./14;
        lower_bound_plt = koem.obm_gy_nh4./14 - koem.obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_nh4./14,'r-','linew',2); 
        plot(1:length(koem.obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_nh4),lower_bound_plt,'r-','linew',2); 
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 14])
        alpha(0.3) %transparency
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_gy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_chl,'b','linew',2);
%         plot(pcase.gy_chl,'k','linew',2);
%         plot(mp_p_sewer.gy_chl,'g','linew',2);
%         plot(mp_p_sewer_det.gy_chl,'m','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_chl + koem.obs_std_gy_chl;
        lower_bound_plt = koem.obm_gy_chl - koem.obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_chl,'r-','linew',2); 
        plot(1:length(koem.obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_chl),lower_bound_plt,'r-','linew',2);  
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        xlim([1 365])
        ylim([0 30])
        yticks(1:2:30)
        alpha(0.3) %transparency  
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_gy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_temp,'b','linew',2);
%         plot(pcase.gy_temp,'k','linew',2);
%         plot(mp_p_sewer.gy_temp,'g','linew',2);  
%         plot(mp_p_sewer_det.gy_temp,'m','linew',2);
         xlim([1 365]);
 clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_temp + koem.obs_std_gy_temp;
        lower_bound_plt = koem.obm_gy_temp - koem.obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_temp,'r-','linew',2); 
        plot(1:length(koem.obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_temp),lower_bound_plt,'r-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        alpha(0.3) %transparency
        ylim([0 35])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_gy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_salt,'b','linew',2);
%         plot(pcase.gy_salt,'k','linew',2);
%         plot(mp_p_sewer.gy_salt,'g','linew',2);   
%         plot(mp_p_sewer_det.gy_salt,'m','linew',2);
         xlim([1 365]);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_salt + koem.obs_std_gy_salt;
        lower_bound_plt = koem.obm_gy_salt - koem.obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_salt,'r-','linew',2); 
        plot(1:length(koem.obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_salt),lower_bound_plt,'r-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
         alpha(0.3) %transparency 
        ylim([0 35])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_gy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_no3_b,'b','linew',2);
%         plot(pcase.gy_no3_b,'k','linew',2);
%         plot(mp_p_sewer.gy_no3_b,'g','linew',2);
%         plot(mp_p_sewer_det.gy_no3_b,'m','linew',2);
        xlim([1 365]);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_no3_b./14 + koem.obs_std_gy_no3_b./14;
        lower_bound_plt = koem.obm_gy_no3_b./14 - koem.obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_no3_b./14,'r-','linew',2); 
        plot(1:length(koem.obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_no3_b),lower_bound_plt,'r-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        alpha(0.3) %transparency    
        ylim([0 75])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_gy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_nh4_b,'b','linew',2);
%         plot(pcase.gy_nh4_b,'k','linew',2);
%         plot(mp_p_sewer.gy_nh4_b,'g','linew',2);
%         plot(mp_p_sewer_det.gy_nh4_b,'m','linew',2);
        xlim([1 365]);
         clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_nh4_b./14 + koem.obs_std_gy_nh4_b./14;
        lower_bound_plt = koem.obm_gy_nh4_b./14 - koem.obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_nh4_b./14,'r-','linew',2); 
        plot(1:length(koem.obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_nh4_b),lower_bound_plt,'r-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 14])
        alpha(0.3) %transparency  
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_gy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_chl_b,'b','linew',2);
%         plot(pcase.gy_chl_b,'k','linew',2);
%         plot(mp_p_sewer.gy_chl_b,'g','linew',2);
%         plot(mp_p_sewer_det.gy_chl_b,'m','linew',2);
         xlim([1 365]);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_chl_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_chl_b + koem.obs_std_gy_chl_b;
        lower_bound_plt = koem.obm_gy_chl_b - koem.obs_std_gy_chl_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_chl_b,'r-','linew',2); 
        plot(1:length(koem.obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_chl_b),lower_bound_plt,'r-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 30])
        yticks(1:2:30)
        alpha(0.3) %transparency  
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_gy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_temp_b,'b','linew',2);
%         plot(pcase.gy_temp_b,'k','linew',2);
%         plot(mp_p_sewer.gy_temp_b,'g','linew',2);    
%         plot(mp_p_sewer_det.gy_temp_b,'m','linew',2);
         xlim([1 365]);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_temp_b + koem.obs_std_gy_temp_b;
        lower_bound_plt = koem.obm_gy_temp_b - koem.obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_temp_b,'r-','linew',2); 
        plot(1:length(koem.obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_temp_b),lower_bound_plt,'r-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
         ylim([0 35])
         alpha(0.3) %transparency   
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_gy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.gy_salt_b,'b','linew',2);
%         plot(pcase.gy_salt_b,'k','linew',2);
%         plot(mp_p_sewer.gy_salt_b,'g','linew',2);
%         plot(mp_p_sewer_det.gy_salt_b,'m','linew',2);
         xlim([1 365]);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(koem.obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = koem.obm_gy_salt_b + koem.obs_std_gy_salt_b;
        lower_bound_plt = koem.obm_gy_salt_b - koem.obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(koem.obm_gy_salt_b,'r-','linew',2); 
        plot(1:length(koem.obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(koem.obm_gy_salt_b),lower_bound_plt,'r-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        alpha(0.3) %transparency
        ylim([0 35])
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_gy'),'-dpng')    
    
        return
%%south gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_no3,'b','linew',2);
        plot(pcase.sgy_no3,'k','linew',2);
        plot(mp_p_sewer.sgy_no3,'g','linew',2); 
        plot(obm_sgy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 + obs_std_sgy_no3./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 - obs_std_sgy_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
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
        plot(mp_p_sewer.sgy_nh4,'g','linew',2);
        plot(obm_sgy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 + obs_std_sgy_nh4./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 - obs_std_sgy_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
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
        plot(mp_p_sewer.sgy_chl,'g','linew',2);
        plot(obm_sgy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl),obm_sgy_chl + obs_std_sgy_chl,'m-','linew',2);
        plot(1:length(obm_sgy_chl),obm_sgy_chl - obs_std_sgy_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
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
        plot(mp_p_sewer.sgy_temp,'g','linew',2);
        plot(obm_sgy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp),obm_sgy_temp + obs_std_sgy_temp,'m-','linew',2);
        plot(1:length(obm_sgy_temp),obm_sgy_temp - obs_std_sgy_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_salt,'b','linew',2);
        plot(pcase.sgy_salt,'k','linew',2);
        plot(mp_p_sewer.sgy_salt,'g','linew',2);   
        plot(obm_sgy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt),obm_sgy_salt + obs_std_sgy_salt,'m-','linew',2);
        plot(1:length(obm_sgy_salt),obm_sgy_salt - obs_std_sgy_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
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
        plot(mp_p_sewer.sgy_no3_b,'g','linew',2);
        plot(obm_sgy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 + obs_std_sgy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 - obs_std_sgy_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_nh4_b,'b','linew',2);
        plot(pcase.sgy_nh4_b,'k','linew',2);
        plot(mp_p_sewer.sgy_nh4_b,'g','linew',2); 
        plot(obm_sgy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 + obs_std_sgy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 - obs_std_sgy_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
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
        plot(mp_p_sewer.sgy_chl_b,'g','linew',2);  
        plot(obm_sgy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b + obs_std_sgy_chl_b,'m-','linew',2);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b - obs_std_sgy_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_temp_b,'b','linew',2);
        plot(pcase.sgy_temp_b,'k','linew',2);
        plot(mp_p_sewer.sgy_temp_b,'g','linew',2);        
        plot(obm_sgy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b + obs_std_sgy_temp_b,'m-','linew',2);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b - obs_std_sgy_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_salt_b,'b','linew',2);
        plot(pcase.sgy_salt_b,'k','linew',2);
        plot(mp_p_sewer.sgy_salt_b,'g','linew',2);    
        plot(obm_sgy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b + obs_std_sgy_salt_b,'m-','linew',2);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b - obs_std_sgy_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
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
        plot(mp_p_sewer.egy_no3,'g','linew',2);  
        plot(obm_egy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 + obs_std_egy_no3./14,'m-','linew',2);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 - obs_std_egy_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines       
        plot(mpdata.egy_nh4,'b','linew',2);
        plot(pcase.egy_nh4,'k','linew',2);
        plot(mp_p_sewer.egy_nh4,'g','linew',2); 
        plot(obm_egy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 + obs_std_egy_nh4./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 - obs_std_egy_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
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
        plot(mp_p_sewer.egy_chl,'g','linew',2); 
        plot(obm_egy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl),obm_egy_chl + obs_std_egy_chl,'m-','linew',2);
        plot(1:length(obm_egy_chl),obm_egy_chl - obs_std_egy_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_temp,'b','linew',2);
        plot(pcase.egy_temp,'k','linew',2);
        plot(mp_p_sewer.egy_temp,'g','linew',2); 
        plot(obm_egy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp),obm_egy_temp + obs_std_egy_temp,'m-','linew',2);
        plot(1:length(obm_egy_temp),obm_egy_temp - obs_std_egy_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_salt,'b','linew',2);
        plot(pcase.egy_salt,'k','linew',2);
        plot(mp_p_sewer.egy_salt,'g','linew',2); 
        plot(obm_egy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt),obm_egy_salt + obs_std_egy_salt,'m-','linew',2);
        plot(1:length(obm_egy_salt),obm_egy_salt - obs_std_egy_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
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
        plot(mp_p_sewer.egy_no3_b,'g','linew',2); 
        plot(obm_egy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 + obs_std_egy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 - obs_std_egy_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_nh4_b,'b','linew',2);
        plot(pcase.egy_nh4_b,'k','linew',2);
        plot(mp_p_sewer.egy_nh4_b,'g','linew',2); 
        plot(obm_egy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 + obs_std_egy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 - obs_std_egy_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
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
        plot(mp_p_sewer.egy_chl_b,'g','linew',2); 
        plot(obm_egy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b + obs_std_egy_chl_b,'m-','linew',2);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b - obs_std_egy_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_temp_b,'b','linew',2);
        plot(pcase.egy_temp_b,'k','linew',2);
        plot(mp_p_sewer.egy_temp_b,'g','linew',2); 
        plot(obm_egy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b + obs_std_egy_temp_b,'m-','linew',2);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b - obs_std_egy_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_salt_b,'b','linew',2);
        plot(pcase.egy_salt_b,'k','linew',2);
        plot(mp_p_sewer.egy_salt_b,'g','linew',2); 
        plot(obm_egy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b + obs_std_egy_salt_b,'m-','linew',2);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b - obs_std_egy_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
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
        plot(mp_p_sewer.jj_no3,'g','linew',2); 
        plot(obm_jj_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 + obs_std_jj_no3./14,'m-','linew',2);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 - obs_std_jj_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_nh4,'b','linew',2);
        plot(pcase.jj_nh4,'k','linew',2);
        plot(mp_p_sewer.jj_nh4,'g','linew',2); 
        plot(obm_jj_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 + obs_std_jj_nh4./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 - obs_std_jj_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
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
        plot(mp_p_sewer.jj_chl,'g','linew',2); 
        plot(obm_jj_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl),obm_jj_chl + obs_std_jj_chl,'m-','linew',2);
        plot(1:length(obm_jj_chl),obm_jj_chl - obs_std_jj_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_temp,'b','linew',2);
        plot(pcase.jj_temp,'k','linew',2);
        plot(mp_p_sewer.jj_temp,'g','linew',2); 
        plot(obm_jj_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp),obm_jj_temp + obs_std_jj_temp,'m-','linew',2);
        plot(1:length(obm_jj_temp),obm_jj_temp - obs_std_jj_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_salt,'b','linew',2);
        plot(pcase.jj_salt,'k','linew',2);
        plot(mp_p_sewer.jj_salt,'g','linew',2); 
        plot(obm_jj_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt),obm_jj_salt + obs_std_jj_salt,'m-','linew',2);
        plot(1:length(obm_jj_salt),obm_jj_salt - obs_std_jj_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
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
        plot(mp_p_sewer.jj_no3_b,'g','linew',2); 
        plot(obm_jj_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 + obs_std_jj_no3_b./14,'m-','linew',2);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 - obs_std_jj_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_nh4_b,'b','linew',2);
        plot(pcase.jj_nh4_b,'k','linew',2);
        plot(mp_p_sewer.jj_nh4_b,'g','linew',2); 
        plot(obm_jj_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 + obs_std_jj_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 - obs_std_jj_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_chl_b,'b','linew',2);
        plot(pcase.jj_chl_b,'k','linew',2);
        plot(mp_p_sewer.jj_chl_b,'g','linew',2); 
        plot(obm_jj_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b + obs_std_jj_chl_b,'m-','linew',2);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b - obs_std_jj_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_temp_b,'b','linew',2);
        plot(pcase.jj_temp_b,'k','linew',2);
        plot(mp_p_sewer.jj_temp_b,'g','linew',2); 
        plot(obm_jj_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b + obs_std_jj_temp_b,'m-','linew',2);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b - obs_std_jj_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_salt_b,'b','linew',2);
        plot(pcase.jj_salt_b,'k','linew',2);
        plot(mp_p_sewer.jj_salt_b,'g','linew',2); 
        plot(obm_jj_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b + obs_std_jj_salt_b,'m-','linew',2);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b - obs_std_jj_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_jj'),'-dpng')    


