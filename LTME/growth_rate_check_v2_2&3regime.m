close all; clc; clear;
clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_9points.mat');
% regime1_kdc
% regime1_kdc_zerosw=load('kdc_nosew_1reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('part_nosew_t3_1reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_noriv_chl=load('1regime_part_res_noriv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_riv_chl=load('1regime_part_res_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc=
% regime2_kdc_zerosw=load('kdc_nosew_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_kdc_zerosw=load('part_nosew_t3_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_noriv_chl=load('2regime_part_res_2regime_noriv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_riv_chl=load('2regime_part_res_2regime_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime3_kdc_zerosw=load('kdc_nosew_3reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime3_kdc_zerosw=load('part_nosew_t3_3reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');

% kd_1reg=ncread('J:\장기생태_2021\Dynamic\result\3regime\day_clim_ERA5_1regime_AttFac.nc','AttFac');
kd_2reg=ncread('J:\장기생태_2021\Dynamic\result\3regime\day_clim_ERA5_2regime_AttFac.nc','AttFac');
kd_3reg=ncread('J:\장기생태_2021\Dynamic\result\3regime\day_clim_ERA5_3regime_AttFac.nc','AttFac');
% path_nc='J:\장기생태_2021\Dynamic\result\3regime';
% swrad=ncread([path_nc,'\day_clim_ERA5_1regime_swrad.nc'],'swrad');
% swrad_2reg=ncread([path_nc,'\day_clim_ERA5_2regime_swrad.nc'],'swrad');

case_name = 't4_pre'

%% Eppley maximum growth rate for temp.
mu_0 = 1.0; %Eppley temperature-limited growth param.
coeffi = 0.59 ; %default
% coeffi = 0.69
% mu_max = mu_0 .* coeffi .* 1.066.^(regime1_kdc_zerosw.sp_gy_temp_mid); %maximum growth rate, zcor= sp_gy_temp_mm
mu_max_2reg = mu_0 .* coeffi .* 1.066.^(regime2_kdc_zerosw.sp_gy_temp_mid); %maximum growth rate  zcor = sp_gy_temp_2reg_mm
mu_max_3reg = mu_0 .* coeffi .* 1.066.^(regime3_kdc_zerosw.sp_gy_temp_mid); %maximum growth rate  zcor = sp_gy_temp_2reg_mm

%% growth rate with radiation
% z_ref = regime1_kdc_zerosw.sp_gy_zw;
z_ref_2reg = regime2_kdc_zerosw.sp_gy_zw;
z_ref_3reg = regime3_kdc_zerosw.sp_gy_zw;
% I_0 = regime1_kdc_zerosw.sp_gy_swrad; %radiation [watt/m^2]
I_0_2reg = regime2_kdc_zerosw.sp_gy_swrad; %radiation
I_0_3reg = regime3_kdc_zerosw.sp_gy_swrad; %radiation
par = 0.43; % Fraction of shortwave radiation that is photosynthetically active 
Kchl=0.025; 
rho0 = 1025.00; % Mean density [kg/m^3]
Cp=3985.00; % spcific heat for seawater [Joules/kg/degc]
alpha = 0.025; % P-I initial slop default
% alpha = 0.09; % P-I initial slop
% PAR_sur = par .* I_0 .* rho0 .* Cp;
% PAR_sur = par .* I_0 ;
% I_up = PAR_sur(:,1:size(z_ref,3));
% I=NaN(size(z_ref,1),size(z_ref,3),size(z_ref,2)-1);
% for i=20:-1:2
%     clearvars att
% att = squeeze((z_ref(:,i,:)-z_ref(:,i-1,:)).*(kd_1reg(1:size(z_ref,1),1,1:size(z_ref,3)) + (Kchl .* regime1_kdc_zerosw.sp_gy_chl_mid(:,i,:))));
% if i==20
%     I(:,:,i) = I_up .* (1 - exp(-1.* att))./att;
%     I_up = I_up  .* exp(-1.* att);
% else 
%     I(:,:,i) = I_up .* (1 - exp(-1.* att))./att;
%     I_up = I_up  .* exp(-1.* att);
% end
% end
% I=permute(I,[1 3 2]);
% f_I = (alpha.*I) ./ sqrt(mu_max.^2 + alpha^2 .* I.^2);

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

% PAR_sur_3reg = par .* I_0_3reg .* rho0 .* Cp;
PAR_sur_3reg = par .* I_0_3reg ;
I_up_3reg = PAR_sur_3reg(:,1:size(z_ref_3reg,3));
I_3reg=NaN(size(z_ref_3reg,1),size(z_ref_3reg,3),size(z_ref_3reg,2)-1);
for i=20:-1:2
    clearvars att
att = squeeze((z_ref_3reg(:,i,:)-z_ref_3reg(:,i-1,:)).*(kd_3reg(1:size(z_ref_3reg,1),1,1:size(z_ref_3reg,3)) + (Kchl .* regime3_kdc_zerosw.sp_gy_chl_mid(:,i,:))));
if i==20
    I_3reg(:,:,i) = I_up_3reg .* (1 - exp(-1.* att))./att;
    I_up_3reg = I_up_3reg  .* exp(-1.* att);
else 
    I_3reg(:,:,i) = I_up_3reg .* (1 - exp(-1.* att))./att;
    I_up_3reg = I_up_3reg  .* exp(-1.* att);
end
end
I_3reg=permute(I_3reg,[1 3 2]);
f_I_3reg = (alpha.*I_3reg) ./ sqrt(mu_max_3reg.^2 + alpha^2 .* I_3reg.^2);

%% growth rate with nutrients
% inverse half saturation
% Kno3=(0.007);
% Knh4=(0.007);
% Kno3=(1.5);
% Knh4=(1.5);
% Kno3=(10);
% Knh4=(10);
% Kno3=(0.1);
% Knh4=(0.1);
%t =4
Kno3=(.6667);
Knh4=(.6667);
%t=5
% Kno3=(142);
% Knh4=(142);
% Kno3=(10)^-1;
% Knh4=(10)^-1;
% Kno3=(0.1)^-1; % inverse half saturation
% Knh4=(0.1)^-1; % inverse half saturation
% Lno3 = (regime1_kdc_zerosw.sp_gy_no3_mid .* Kno3 ) .* (1.0 ./ (1.0 + regime1_kdc_zerosw.sp_gy_nh4_mid .* Knh4)) ./ (1.0 + regime1_kdc_zerosw.sp_gy_no3_mid .* Kno3);
% LNH4 = (regime1_kdc_zerosw.sp_gy_nh4_mid.* Knh4 )./(1 + regime1_kdc_zerosw.sp_gy_nh4_mid .* Knh4);
% nut_term = (Lno3 + LNH4);

Lno3_2reg = (regime2_kdc_zerosw.sp_gy_no3_mid .* Kno3 ) .* (1.0 ./ (1.0 + regime2_kdc_zerosw.sp_gy_nh4_mid .* Knh4)) ./ (1.0 + regime2_kdc_zerosw.sp_gy_no3_mid .* Kno3);
LNH4_2reg = (regime2_kdc_zerosw.sp_gy_nh4_mid.* Knh4 )./(1 + regime2_kdc_zerosw.sp_gy_nh4_mid .* Knh4);
nut_term_2reg = (Lno3_2reg + LNH4_2reg);

Lno3_3reg = (regime3_kdc_zerosw.sp_gy_no3_mid .* Kno3 ) .* (1.0 ./ (1.0 + regime3_kdc_zerosw.sp_gy_nh4_mid .* Knh4)) ./ (1.0 + regime3_kdc_zerosw.sp_gy_no3_mid .* Kno3);
LNH4_3reg = (regime3_kdc_zerosw.sp_gy_nh4_mid.* Knh4 )./(1 + regime3_kdc_zerosw.sp_gy_nh4_mid .* Knh4);
nut_term_3reg = (Lno3_3reg + LNH4_3reg);

%% growth rate
% mu = mu_max .* f_I .* nut_term;
mu_2reg = mu_max_2reg .* f_I_2reg .* nut_term_2reg;
mu_3reg = mu_max_3reg .* f_I_3reg .* nut_term_3reg;

fig = figure; hold on;
% plot(squeeze(nanmean(mu(:,18,:),1)),'r','linew',2);
plot(squeeze(nanmean(mu_2reg(:,18,:),1)),'g','linew',2);
plot(squeeze(nanmean(mu_3reg(:,18,:),1)),'b','linew',2);
% alpha(0.3) %transparency      
title(['total growth rate with with u_m_a_x * f(I) * (Lno3 + Lnh4)']);
xlabel('time','fontsize',13)
ylabel('total growth rate','fontsize',13)
grid on
ylim([-inf inf])
% xlim(xlim_set)
set(gca,'fontsize',13)
print(fig,strcat(['total_growth_rate_',case_name]),'-dpng')


fig = figure; hold on;
% plot(squeeze(nanmean(nut_term(:,18,:),1)),'r--','linew',2);
plot(squeeze(nanmean(nut_term_2reg(:,18,:),1)),'g','linew',2);
plot(squeeze(nanmean(nut_term_3reg(:,18,:),1)),'b','linew',2);
% alpha(0.3) %transparency      
% title(['growth rate with nutrients']);
xlabel('time','fontsize',13)
ylabel('growth rate with nutrients','fontsize',13)
grid on
ylim([-inf 1.1])
% xlim(xlim_set)
set(gca,'fontsize',13)
print(fig,strcat(['growth_rate_with_nutrients',case_name]),'-dpng')


fig = figure; hold on;
% plot(squeeze(nanmean(mu_max(:,18,:),1)),'r','linew',2);
plot(squeeze(nanmean(mu_max_2reg(:,18,:),1)),'g','linew',2);
plot(squeeze(nanmean(mu_max_3reg(:,18,:),1)),'b','linew',2);
% alpha(0.3) %transparency      
% title(['growth rate with nutrients']);
xlabel('time','fontsize',13)
ylabel('Eppley maximum growth rate for temp.','fontsize',13)
grid on
% ylim([0 70])
% xlim(xlim_set)
set(gca,'fontsize',13)
print(fig,strcat(['Eppley_maximum_growth_rate_for_temp',case_name]),'-dpng')


fig = figure; hold on;
% plot(squeeze(nanmean(f_I(:,18,:),1)),'r','linew',2);
plot(squeeze(nanmean(f_I_2reg(:,18,:),1)),'g','linew',2);
plot(squeeze(nanmean(f_I_3reg(:,18,:),1)),'b','linew',2);
% alpha(0.3) %transparency      
% title(['growth rate with nutrients']);
xlabel('time','fontsize',13)
ylabel('growth rate with irr..','fontsize',13)
grid on
% ylim([0 70])
% xlim(xlim_set)
set(gca,'fontsize',13)
print(fig,strcat(['growth_rate_with_irr',case_name]),'-dpng')
