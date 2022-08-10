close all; clc; clear;
clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_9points.mat');
% regime1_kdc
% regime1_kdc_zerosw=load('regime1_kdc_zerosw.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime1_kdc_zerosw=load('part_kdc_nosew_1reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_noriv_chl=load('1regime_part_res_noriv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_riv_chl=load('1regime_part_res_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc=
regime2_kdc_zerosw=load('part_kdc_nosew_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_noriv_chl=load('2regime_part_res_2regime_noriv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_riv_chl=load('2regime_part_res_2regime_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');

kd_1reg=ncread('C:\Users\user\Desktop\GY_메타\day_clim_ERA5_1regime_AttFac.nc','AttFac');
kd_2reg=ncread('C:\Users\user\Desktop\GY_메타\day_clim_ERA5_2regime_AttFac.nc','AttFac');
% path_nc='C:\Users\user\Desktop\GY_메타';
% swrad=ncread([path_nc,'\day_clim_ERA5_1regime_swrad.nc'],'swrad');
% swrad_2reg=ncread([path_nc,'\day_clim_ERA5_2regime_swrad.nc'],'swrad');

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

clearvars temp_d sp_gy_*_mm
for i = 1:9 %st
    for j = 1: size(regime1_kdc_zerosw.sp_gy_chl_mid,3) %time    
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_chl_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_chl_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_no3_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_no3_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_nh4_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_nh4_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_mm(i,:,j)=temp_d;
        % 2 regime
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_chl_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_chl_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_no3_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_no3_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_nh4_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_nh4_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_2reg_mm(i,:,j)=temp_d;
        % 3 regime
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_chl_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_chl_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_no3_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_no3_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_nh4_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_nh4_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_2reg_mm(i,:,j)=temp_d;
        
    end
end
% plot(squeeze(mean(sp_gy_chl_mm(:,3,:),1)))
% squeeze(regime1_kdc_zerosw.sp_gy_z(:,1,1)) % 9st's depth

% clearvars clim2*
% clim2=clim1.obm_gy_po4_ft_1(:,3);
% clim2_std = clim1.obs_std_gy_po4_ft_1(:,3);
% 
% clim2_2=clim1.obm_gy_po4_ft_2(:,3);
% clim2_2_std = clim1.obs_std_gy_po4_ft_2(:,3);

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

%% Eppley maximum growth rate for temp.
mu_0 = 1.0; %Eppley temperature-limited growth param.
mu_max = mu_0 .* 0.59 .* 1.066.^(regime1_kdc_zerosw.sp_gy_temp_mid); %maximum growth rate, zcor= sp_gy_temp_mm
mu_max_2reg = mu_0 .* 0.59 .* 1.066.^(regime2_kdc_zerosw.sp_gy_temp_mid); %maximum growth rate  zcor = sp_gy_temp_2reg_mm

%% growth rate with radiation

z_ref = regime1_kdc_zerosw.sp_gy_zw;
z_ref_2reg = regime2_kdc_zerosw.sp_gy_zw;
I_0 = regime1_kdc_zerosw.sp_gy_swrad; %radiation [watt/m^2]
I_0_2reg = regime2_kdc_zerosw.sp_gy_swrad; %radiation
par = 0.43; % Fraction of shortwave radiation that is photosynthetically active 
Kchl=0.025; 
rho0 = 1025.00; % Mean density [kg/m^3]
Cp=3985.00; % spcific heat for seawater [Joules/kg/degc]
alpha = 0.025; % P-I initial slop
% PAR_sur = par .* I_0 .* rho0 .* Cp;
PAR_sur = par .* I_0 ;
I_up = PAR_sur(:,1:size(z_ref,3));
I=NaN(size(z_ref,1),size(z_ref,3),size(z_ref,2)-1);
for i=20:-1:2
    clearvars att
att = squeeze((z_ref(:,i,:)-z_ref(:,i-1,:)).*(kd_1reg(1:size(z_ref,1),1,1:size(z_ref,3)) + (Kchl .* regime1_kdc_zerosw.sp_gy_chl_mid(:,i,:))));
if i==20
    I(:,:,i) = I_up .* (1 - exp(-1.* att))./att;
    I_up = I_up  .* exp(-1.* att);
else 
    I(:,:,i) = I_up .* (1 - exp(-1.* att))./att;
    I_up = I_up  .* exp(-1.* att);
end
end
I=permute(I,[1 3 2]);
f_I = (alpha.*I) ./ sqrt(mu_max.^2 + alpha^2 .* I.^2);

% PAR_sur_2reg = par .* I_0_2reg .* rho0 .* Cp;
PAR_sur_2reg = par .* I_0_2reg ;
I_up_2reg = PAR_sur_2reg(:,1:size(z_ref_2reg,3));
I_2reg=NaN(size(z_ref_2reg,1),size(z_ref_2reg,3),size(z_ref_2reg,2)-1);
for i=20:-1:2
    clearvars att
att = squeeze((z_ref_2reg(:,i,:)-z_ref_2reg(:,i-1,:)).*(kd_2reg(1:size(z_ref_2reg,1),1,1:size(z_ref_2reg,3)) + (Kchl .* regime2_kdc_zerosw.sp_gy_chl_mid(:,i,:))));
if i==20
    I_2reg(:,:,i) = I_up_2reg .* (1 - exp(-1.* att))./att;
    I_up_2reg = I_up_2reg  .* exp(-1.* att);
else 
    I_2reg(:,:,i) = I_up_2reg .* (1 - exp(-1.* att))./att;
    I_up_2reg = I_up_2reg  .* exp(-1.* att);
end
end
I_2reg=permute(I_2reg,[1 3 2]);
f_I_2reg = (alpha.*I_2reg) ./ sqrt(mu_max_2reg.^2 + alpha^2 .* I_2reg.^2);

%% growth rate with nutrients
% inverse half saturation
% Kno3=(10)^-1;
% Knh4=(10)^-1;
Kno3=(0.1)^-1; % inverse half saturation
Knh4=(0.1)^-1; % inverse half saturation
Lno3 = (regime1_kdc_zerosw.sp_gy_no3_mid./((Kno3^-1) + regime1_kdc_zerosw.sp_gy_no3_mid)) .* (1./(1+regime1_kdc_zerosw.sp_gy_nh4_mid./(Knh4^-1)));
LNH4 = (regime1_kdc_zerosw.sp_gy_nh4_mid./((Knh4^-1) + regime1_kdc_zerosw.sp_gy_nh4_mid));
nut_term = (Lno3 + LNH4);

Lno3_2reg = (regime2_kdc_zerosw.sp_gy_no3_mid./((Kno3^-1) + regime2_kdc_zerosw.sp_gy_no3_mid)) .* (1./(1+regime2_kdc_zerosw.sp_gy_nh4_mid./(Knh4^-1)));
LNH4_2reg = (regime2_kdc_zerosw.sp_gy_nh4_mid./((Knh4^-1) + regime2_kdc_zerosw.sp_gy_nh4_mid));
nut_term_2reg = (Lno3_2reg + LNH4_2reg);

%% growth rate
mu = mu_max .* f_I .* nut_term;
mu_2reg = mu_max_2reg .* f_I_2reg .* nut_term_2reg;

fig = figure; hold on;
plot(squeeze(nanmean(mu(:,18,:),1)),'r','linew',2);
plot(squeeze(nanmean(mu_2reg(:,18,:),1)),'g','linew',2);
% alpha(0.3) %transparency      
title(['total growth rate with with u_m_a_x * f(I) * (Lno3 + Lnh4)']);
xlabel('time','fontsize',13)
ylabel('total growth rate','fontsize',13)
grid on
ylim([-inf inf])
xlim(xlim_set)
set(gca,'fontsize',13)
print(fig,strcat('total_growth_rate'),'-dpng')


fig = figure; hold on;
plot(squeeze(nanmean(nut_term(:,18,:),1)),'r','linew',2);
plot(squeeze(nanmean(nut_term_2reg(:,18,:),1)),'g','linew',2);
% alpha(0.3) %transparency      
% title(['growth rate with nutrients']);
xlabel('time','fontsize',13)
ylabel('growth rate with nutrients','fontsize',13)
grid on
ylim([-inf 1.1])
xlim(xlim_set)
set(gca,'fontsize',13)
print(fig,strcat('growth_rate_with_nutrients'),'-dpng')


fig = figure; hold on;
plot(squeeze(nanmean(mu_max(:,18,:),1)),'r','linew',2);
plot(squeeze(nanmean(mu_max_2reg(:,18,:),1)),'g','linew',2);
% alpha(0.3) %transparency      
% title(['growth rate with nutrients']);
xlabel('time','fontsize',13)
ylabel('Eppley maximum growth rate for temp.','fontsize',13)
grid on
% ylim([0 70])
xlim(xlim_set)
set(gca,'fontsize',13)
print(fig,strcat('Eppley_maximum_growth_rate_for_temp'),'-dpng')


fig = figure; hold on;
plot(squeeze(nanmean(f_I(:,18,:),1)),'r','linew',2);
plot(squeeze(nanmean(f_I_2reg(:,18,:),1)),'g','linew',2);
% alpha(0.3) %transparency      
% title(['growth rate with nutrients']);
xlabel('time','fontsize',13)
ylabel('growth rate with irr..','fontsize',13)
grid on
% ylim([0 70])
xlim(xlim_set)
set(gca,'fontsize',13)
print(fig,strcat('growth_rate_with_irr'),'-dpng')


clearvars clim2*
clim2=clim1.obm_gy_din_ft_1(:,3);
clim2_std = clim1.obs_std_gy_din_ft_1(:,3);
clim2_2=clim1.obm_gy_din_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_din_ft_2(:,3);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)) + squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)),'r','linew',2);
%             plot((squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)).*N_MW + squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)).*N_MW)./N_MW,'k--','linew',2);
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm(:,3,:),1)) + squeeze(nanmean(sp_gy_nh4_2reg_mm(:,3,:),1)),'g','linew',2);
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

            title(['GY KOEM OBS vs. daily MODEL DIN']);
            xlabel('time','fontsize',13)
            ylabel('DIN (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_DIN'),'-dpng')

    clearvars clim2*
clim2=clim1.obm_gy_no3_ft_1(:,3);
clim2_std = clim1.obs_std_gy_no3_ft_1(:,3);
clim2_2=clim1.obm_gy_no3_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_no3_ft_2(:,3);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm(:,3,:),1)),'g','linew',2);
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
            plot(squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_nh4_2reg_mm(:,3,:),1)),'g','linew',2);
            
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
            plot(squeeze(nanmean(sp_gy_chl_mm(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_chl_2reg_mm(:,3,:),1)),'g','linew',2);
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
            ylim([0 21])
            xlim(xlim_set)
            set(gca,'fontsize',13)

    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


     clearvars clim2*
clim2=clim1.obm_gy_temp_ft_1(:,3);
clim2_std = clim1.obs_std_gy_temp_ft_1(:,3);
clim2_2=clim1.obm_gy_temp_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_temp_ft_2(:,3);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_temp_mm(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_temp_2reg_mm(:,3,:),1)),'g','linew',2);
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
            plot(squeeze(nanmean(sp_gy_salt_mm(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_salt_2reg_mm(:,3,:),1)),'g','linew',2);
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
%             plot(nanmean(regime1_noriv_chl.gy_salt,1),'r','linew',2);
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
clim2=clim1.obm_gy_din_b_ft_1(:,3);
clim2_std = clim1.obs_std_gy_din_b_ft_1(:,3);
clim2_2=clim1.obm_gy_din_b_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_din_b_ft_2(:,3);

    fig = figure; hold on;
            plot(squeeze(nanmean(regime1_kdc_zerosw.sp_gy_no3_b,1)) + squeeze(nanmean(regime1_kdc_zerosw.sp_gy_nh4_b,1)),'r','linew',2);
%             plot((squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)).*N_MW + squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)).*N_MW)./N_MW,'k--','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw.sp_gy_no3_b,1)) + squeeze(nanmean(regime2_kdc_zerosw.sp_gy_nh4_b,1)),'g','linew',2);
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

            title(['GY KOEM OBS vs. daily MODEL DIN']);
            xlabel('time','fontsize',13)
            ylabel('DIN (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_DIN_bot'),'-dpng')


clearvars clim2*
clim2=clim1.obm_gy_no3_b_ft_1(:,3);
clim2_std = clim1.obs_std_gy_no3_b_ft_1(:,3);
clim2_2=clim1.obm_gy_no3_b_ft_2(:,3);
clim2_2_std = clim1.obs_std_gy_no3_b_ft_2(:,3);

     fig = figure; hold on;
            plot(squeeze(nanmean(regime1_kdc_zerosw.sp_gy_no3_b,1)),'r','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw.sp_gy_no3_b,1)),'g','linew',2);
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
            plot(squeeze(nanmean(regime1_kdc_zerosw.sp_gy_nh4_b,1)),'r','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw.sp_gy_nh4_b,1)),'g','linew',2);
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
            plot(squeeze(nanmean(regime1_kdc_zerosw.sp_gy_chl_b,1)),'r','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw.sp_gy_chl_b,1)),'g','linew',2);
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
            ylim([0 20])
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
