close all; clc; clear;
clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_9points.mat');
regime1_kdc_zerosw=load('1regime_part_res_kdc_zerosw.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime1_noriv_chl=load('1regime_part_res_noriv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime1_riv_chl=load('1regime_part_res_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_kdc_zerosw=load('2regime_part_res_2regime_kdc_zerosw.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_noriv_chl=load('2regime_part_res_2regime_noriv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_riv_chl=load('2regime_part_res_2regime_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');

P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;

xlim_set = [1 365];

%% plot obs 1st and 2nd at the same time which is time axis
plt_filter_1{1} = 32:45;
plt_filter_2{1} = 46:59;
plt_filter_1{2} = 121:135.5;
plt_filter_2{2} = 136:151;
plt_filter_1{3} = 213:227.5;
plt_filter_2{3} = 228:243;
plt_filter_1{4} = 305:319;
plt_filter_2{4} = 320:334;

clearvars clim2*
clim2=clim1.obm_gy_po4_ft_1(:,3);
clim2_std = clim1.obs_std_gy_po4_ft_1(:,3);

clim2_2=clim1.obm_gy_po4_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_po4_ft_2(:,3);

% fig = figure; hold on;   
% %             plot(nanmean(mer_po4_1,1),'r','linew',2);
% %             clearvars *_bound_plt nonan_data_plt discon_p temp_d
% %             temp_d=nanmean(mer_po4_1,1);
% %             nonan_data_plt = find(isnan(temp_d)==0);
% %             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
% %             upper_bound_plt = temp_d + regm_std_po4_1;
% %             lower_bound_plt = temp_d - regm_std_po4_1;
% %             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
% 
%             clearvars *_bound_plt nonan_data_plt discon_p
%             nonan_data_plt = find(isnan(clim2)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             discon_p = find(nonan_diff ~= 1); % find discontinuous point
%             discon_p(4) = length(nonan_data_plt); %end point
%             upper_bound_plt = clim2./30.973762 + clim2_std./30.973762;
%             lower_bound_plt = clim2./30.973762 - clim2_std./30.973762;
%             for i = 1:4
%             if i == 1
%             patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
%             else
%             patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
%             end
%             end 
%             plot(clim2./30.973762,'r-','linew',2); 
%             plot(1:length(clim2),upper_bound_plt,'r-','linew',2);
%             plot(1:length(clim2),lower_bound_plt,'r-','linew',2);
% 
%             alpha(0.3) %transparency      
% 
%             title(['GY KOEM OBS vs. daily MODEL po4']);
%             xlabel('time','fontsize',13)
%             ylabel('po4 (umol/m^3)','fontsize',13)
%             grid on
%             ylim([0 5])
%             xlim(xlim_set)
%             set(gca,'fontsize',13)
%     %         %xticks(com_eom(1:12:end))
%     %         %xticklabels(com_t_tic)
%     %         %xtickangle(90)
%     print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4'),'-dpng')
   
    clearvars clim2*
clim2=clim1.obm_gy_no3_ft_1(:,3);
clim2_std = clim1.obs_std_gy_no3_ft_1(:,3);
clim2_2=clim1.obm_gy_no3_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_no3_ft_2(:,3);

    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_no3,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_no3,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_no3_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_no3_1;
%             lower_bound_plt = temp_d - regm_std_no3_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
% 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL NO3']);
            xlabel('time','fontsize',13)
            ylabel('NO3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')

    
        clearvars clim2*
clim2=clim1.obm_gy_nh4_ft_1(:,3);
clim2_std = clim1.obs_std_gy_nh4_ft_1(:,3);
clim2_2=clim1.obm_gy_nh4_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_nh4_ft_2(:,3);

    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_nh4,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_nh4,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_nh4_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_nh4_1;
%             lower_bound_plt = temp_d - regm_std_nh4_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
             % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')


    clearvars clim2*
clim2=clim1.obm_gy_chl_ft_1(:,3);
clim2_std = clim1.obs_std_gy_chl_ft_1(:,3);
clim2_2=clim1.obm_gy_chl_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_chl_ft_2(:,3);

    fig = figure; hold on;
            plot(regime1_kdc_zerosw.gy_chl_mid(20,:),'r','linew',2);
%             plot(nanmean(regime1_kdc_zerosw.gy_chl_mid,1),'r--','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_chl,1),'g','linew',2);
%             plot(nanmean(regime2_kdc_zerosw.gy_chl_mid,1),'g--','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_chl_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_chl_1;
%             lower_bound_plt = temp_d - regm_std_chl_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
             % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 30])
            xlim(xlim_set)
            set(gca,'fontsize',13)

    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


     clearvars clim2*
clim2=clim1.obm_gy_temp_ft_1(:,3);
clim2_std = clim1.obs_std_gy_temp_ft_1(:,3);
clim2_2=clim1.obm_gy_temp_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_temp_ft_2(:,3);

    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_temp,1),'r','linew',2);
             plot(nanmean(regime2_kdc_zerosw.gy_temp,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_temp_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_temp_1;
%             lower_bound_plt = temp_d - regm_std_temp_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp_1; c1998.gy_temp_1; c1999.gy_temp_1; c2000.gy_temp_1; c2001.gy_temp_1; c2002.gy_temp_1; c2003.gy_temp_1;],1),'r','linew',2);
                         % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')

         clearvars clim2*
clim2=clim1.obm_gy_salt_ft_1(:,3);
clim2_std = clim1.obs_std_gy_salt_ft_1(:,3);
clim2_2=clim1.obm_gy_salt_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_salt_ft_2(:,3);   

    fig = figure; hold on; 
            plot(nanmean(regime1_kdc_zerosw.gy_salt,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_salt,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p salt_d
%             salt_d=nanmean(mer_salt_1,1);
%             nonan_data_plt = find(isnan(salt_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = salt_d + regm_std_salt_1;
%             lower_bound_plt = salt_d - regm_std_salt_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_1; c1998.gy_salt_1; c1999.gy_salt_1; c2000.gy_salt_1; c2001.gy_salt_1; c2002.gy_salt_1; c2003.gy_salt_1;],1),'r','linew',2);
                 % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')
    
    
%     clearvars clim2*
% clim2=clim1.obm_gy_salt_ft_1(:,3);
% clim2_std = clim1.obs_std_gy_salt_ft_1(:,3);
%     
%     fig = figure; hold on; 
%             plot(nanmean(regime1_kdc_zerosw.gy_salt,1),'r','linew',2);
%     
%      %% bot
% fig = figure; hold on;
%     %         plot(zeros(365,1),'g--','linew',2); % zero lines
% 
%     %data=NaN(1,365);        
%     % data(1:length(c2004.gy_po4_b))=c2004.gy_po4_b_1; 
%             plot(nanmean(mer_po4_b_1,1),'r','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p po4_b_d
%             po4_b_d=nanmean(mer_po4_b_1,1);
%             nonan_data_plt = find(isnan(po4_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = po4_b_d + regm_std_po4_b_1;
%             lower_bound_plt = po4_b_d - regm_std_po4_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
%     %         plot(nanmean([c1997.gy_po4_b_1; c1998.gy_po4_b_1; c1999.gy_po4_b_1; c2000.gy_po4_b_1; c2001.gy_po4_b_1; c2002.gy_po4_b_1; c2003.gy_po4_b_1;],1),'r','linew',2);
%                    clearvars *_bound_plt nonan_data_plt discon_p
%             nonan_data_plt = find(isnan(clim1.obm_gy_po4_b_ft_1)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             discon_p = find(nonan_diff ~= 1); % find discontinuous point
%             discon_p(4) = length(nonan_data_plt); %end point
%             upper_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 + clim1.obs_std_gy_po4_b_ft_1./30.973762;
%             lower_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 - clim1.obs_std_gy_po4_b_ft_1./30.973762;
%             for i = 1:4
%             if i == 1
%             patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
%             else
%             patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
%             end
%             end 
%             plot(clim1.obm_gy_po4_b_ft_1./30.973762,'r-','linew',2); 
%             plot(1:length(clim1.obm_gy_po4_b_ft_1),upper_bound_plt,'r-','linew',2);
%             plot(1:length(clim1.obm_gy_po4_b_ft_1),lower_bound_plt,'r-','linew',2);
% 
% 
%             alpha(0.3) %transparency      
% 
%             title(['GY KOEM OBS vs. daily MODEL po4 bot.']);
%             xlabel('time','fontsize',13)
%             ylabel('po4 (umol/m^3)','fontsize',13)
%             grid on
%  ylim([0 5])
%             xlim(xlim_set)
%             set(gca,'fontsize',13)
%     %         %xticks(com_eom(1:12:end))
%     %         %xticklabels(com_t_tic)
%     %         %xtickangle(90)
%             print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4_bot'),'-dpng')

clearvars clim2*
clim2=clim1.obm_gy_no3_b_ft_1(:,3);
clim2_std = clim1.obs_std_gy_no3_b_ft_1(:,3);
clim2_2=clim1.obm_gy_no3_b_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_no3_b_ft_2(:,3);

     fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_no3_b,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_no3_b,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p no3_b_d
%             no3_b_d=nanmean(mer_no3_b_1,1);
%             nonan_data_plt = find(isnan(no3_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = no3_b_d + regm_std_no3_b_1;
%             lower_bound_plt = no3_b_d - regm_std_no3_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_no3_b_1; c1998.gy_no3_b_1; c1999.gy_no3_b_1; c2000.gy_no3_b_1; c2001.gy_no3_b_1; c2002.gy_no3_b_1; c2003.gy_no3_b_1;],1),'r','linew',2);
% 1regime
             clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
            xlabel('time','fontsize',13)
            ylabel('no3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

      clearvars clim2*
clim2=clim1.obm_gy_nh4_b_ft_1(:,3);
clim2_std = clim1.obs_std_gy_nh4_b_ft_1(:,3);
clim2_2=clim1.obm_gy_nh4_b_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_nh4_b_ft_2(:,3);

                
    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_nh4_b,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_nh4_b,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p nh4_b_d
%             nh4_b_d=nanmean(mer_nh4_b_1,1);
%             nonan_data_plt = find(isnan(nh4_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = nh4_b_d + regm_std_nh4_b_1;
%             lower_bound_plt = nh4_b_d - regm_std_nh4_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
%     %         plot(nanmean([c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b;],1),'r','linew',2);
            % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')

          clearvars clim2*
clim2=clim1.obm_gy_chl_b_ft_1(:,3);
clim2_std = clim1.obs_std_gy_chl_b_ft_1(:,3);
clim2_2=clim1.obm_gy_chl_b_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_chl_b_ft_2(:,3);
                
    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_chl_b,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_chl_b,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p chl_b_d
%             chl_b_d=nanmean(mer_chl_b_1,1);
%             nonan_data_plt = find(isnan(chl_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = chl_b_d + regm_std_chl_b_1;
%             lower_bound_plt = chl_b_d - regm_std_chl_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b;],1),'r','linew',2);
                        % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL chl bot.']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')


            clearvars clim2*
clim2=clim1.obm_gy_temp_b_ft_1(:,3);
clim2_std = clim1.obs_std_gy_temp_b_ft_1(:,3);
clim2_2=clim1.obm_gy_temp_b_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_temp_b_ft_2(:,3);  
                
    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_temp_b,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_temp_b,1),'g','linew',2);  
%             clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
%             temp_b_d=nanmean(mer_temp_b_1,1);
%             nonan_data_plt = find(isnan(temp_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_b_d + regm_std_temp_b_1;
%             lower_bound_plt = temp_b_d - regm_std_temp_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'r','linew',2);
                     % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')

            clearvars clim2*
clim2=clim1.obm_gy_salt_b_ft_1(:,3);
clim2_std = clim1.obs_std_gy_salt_b_ft_1(:,3);
clim2_2=clim1.obm_gy_salt_b_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_salt_b_ft_2(:,3);   
                
    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_salt_b,1),'r','linew',2);  
            plot(nanmean(regime2_kdc_zerosw.gy_salt_b,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
%             salt_b_d=nanmean(mer_salt_b_1,1);
%             nonan_data_plt = find(isnan(salt_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = salt_b_d + regm_std_salt_b_1;
%             lower_bound_plt = salt_b_d - regm_std_salt_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'r','linew',2);
              % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'r');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'r-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'r-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'g');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'g-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'g-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng') 
            
            close all;
