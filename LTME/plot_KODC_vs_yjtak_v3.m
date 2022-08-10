close all; clear; clc;
cd /home/cshwa/GY_compare
% load kodc_40016_2001_2010_just_mean_v2.mat
load kodc_20501_2001_2010_just_mean_v2.mat
% 

cd /data1/yjtak/auto_fennel/Output/ERA5_rain/control

lon=ncread('./2010/stdep1_avg_mon_2010_01.nc','lon');
lat=ncread('./2010/stdep1_avg_mon_2010_01.nc','lat');
dep_mod=ncread('./2010/stdep1_avg_mon_2010_01.nc','depth');
mask_pre=ncread('./2010/stdep1_avg_mon_2010_01.nc','temp');
[lat_m lon_m]=meshgrid(lat,lon);

clearvars lon lat

lon = lon_m; lat = lat_m;
lon_obs = lon_20501; lat_obs = lat_20501;
% lon_obs = lon_40016; lat_obs = lat_40016;
mask=squeeze(mask_pre(:,:,1));
mask(~isnan(mask)) = 1;

for i = 1:length(lon_obs)
clearvars temp_2d
    temp_2d = sqrt((lon - lon_obs(i)).^2 + (lat - lat_obs(i)).^2);
    temp_2d = temp_2d .* mask;
    near_point(i) =find(nanmin(nanmin(temp_2d))==temp_2d);
    dist_2d(:,:,i) = temp_2d;
end

for i = 1:length(near_point)
    clearvars row col
    [row,col]=ind2sub(size(lon),near_point(i)); %% 1d to 2d index
    near_point_2d(i,:) = [row, col]; 
end

% min(min(dist_2d))
% find(min(min(dist_2d))==dist_2d)
k = 0;
for i = 1:10 %year
     for j = 1:12 %month
        k= k+1;
        t_model(:,:,:,i,j) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'temp');
        s_model(:,:,:,i,j) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'salt');
%         do_model(:,:,:,k) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'oxygen');
        no3_model(:,:,:,i,j) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'NO3');
        nh4_model(:,:,:,i,j) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'NH4');
%         zeta(:,:,k) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'zeta');
     end
end

% t_mod_in
figure
plot(squeeze(t_model(near_point_2d(1,1),near_point_2d(1,2),end,:)));

clearvars z_mod
for i =1:size(t_model,5)
% monthly climate
        t_mod_in(:,i)=squeeze(mean(t_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
        s_mod_in(:,i)=squeeze(mean(s_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
%         do_mod_in(:,i)=squeeze(mth_do(near_point_2d(1,1),near_point_2d(1,2),:,i));
        nh4_mod_in(:,i)=squeeze(mean(nh4_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
        no3_mod_in(:,i)=squeeze(mean(no3_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
 
%monthly std        
        std_t_mod(:,i)=squeeze(nanstd(t_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),0,4));
        std_s_mod(:,i)=squeeze(nanstd(s_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),0,4)); % 0 =weight. 4 = dim.
        std_nh4_mod(:,i)=squeeze(nanstd(nh4_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),0,4));
        std_no3_mod(:,i)=squeeze(nanstd(no3_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),0,4));
end

% monthly mean for all years        
        t_mod_pre(:,:,i)=squeeze(mean(t_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
        s_mod_pre(:,:,i)=squeeze(mean(s_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
        nh4_mod_pre(:,:,i)=squeeze(mean(nh4_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
        no3_mod_pre(:,:,i)=squeeze(mean(no3_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
clearvars mod_clim_*_sur mod_clim_*_bot

mod_clim_t_sur = t_mod_in(1,:);
mod_clim_s_sur = s_mod_in(1,:);
mod_clim_nh4_sur = nh4_mod_in(1,:);
mod_clim_no3_sur = no3_mod_in(1,:);

mod_clim_t_bot = t_mod_in(3,:);
mod_clim_s_bot = s_mod_in(3,:);
mod_clim_nh4_bot = nh4_mod_in(3,:);
mod_clim_no3_bot = no3_mod_in(3,:);

% make monthly on the obs
% eom_366_end =eom_d_each(1,:) %366
% eom_366_start(1,1) =1 %366
% eom_366_start(1,1:12) = eom_366_end(1,1:11)+1

% for i = 1:12
%     obs_temp_sur(i) = nanmean(yp_w_temp_04(eom_366_start(i):eom_366_end(i)));
%     obs_temp_bot(i) = nanmean(yp_w_temp_04_b(eom_366_start(i):eom_366_end(i)));
%     obs_salt_sur(i) = nanmean(yp_w_salt_04(eom_366_start(i):eom_366_end(i)));
%     obs_salt_bot(i) = nanmean(yp_w_salt_04_b(eom_366_start(i):eom_366_end(i)));
%     obs_no3_sur(i) = nanmean(yp_w_no3_04(eom_366_start(i):eom_366_end(i)));
%     obs_no3_bot(i) = nanmean(yp_w_no3_04_b(eom_366_start(i):eom_366_end(i)));
% %     obs_nh4_sur(i) = nanmean(yp_w_nh4_04(eom_366_start(i):eom_366_end(i)));
% %     obs_nh4_bot(i) = nanmean(yp_w_nh4_04_b(eom_366_start(i):eom_366_end(i)));
% end

%surface
return
%% temp
cd /home/cshwa/GY_compare

figure;
scatter(1:12,reg_clim_temp);
hold on
plot(1:12, mod_clim_t_sur,'--','color','r','linew',2)
errorbar(1:12, mod_clim_t_sur,std_t_mod(1,:),'r', 'LineStyle', 'none')
errorbar(1:12, reg_clim_temp,std_temp,'b', 'LineStyle', 'none');
xlabel('time(month)','fontsize',13)
ylabel('temp (deg ^oC)','fontsize',13)
title('temp -polynomial(surf)','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 12])


figure;
scatter(1:12,reg_clim_temp_b);
hold on
plot(1:12, mod_clim_t_bot,'--','color','r','linew',2)
errorbar(1:12, mod_clim_t_bot,std_t_mod(3,:),'r', 'LineStyle', 'none')
errorbar(1:12, reg_clim_temp_b,std_temp_b,'b', 'LineStyle', 'none');
xlabel('time(month)','fontsize',13)
ylabel('temp (deg ^oC)','fontsize',13)
title('temp -polynomial(bot)','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 12])

%% salt
figure;
scatter(1:12,reg_clim_salt);
hold on
plot(1:12, mod_clim_s_sur,'--','color','r','linew',2)
errorbar(1:12, mod_clim_s_sur,std_s_mod(1,:),'r', 'LineStyle', 'none')
errorbar(1:12, reg_clim_salt(1:12),std_salt(1:12),'b', 'LineStyle', 'none');
xlabel('time(month)','fontsize',13)
ylabel('salt (psu)','fontsize',13)
title('salt -polynomial(surf)','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([30 inf])
xlim([1 12])

figure;
scatter(1:12,reg_clim_salt_b);
hold on
plot(1:12, mod_clim_s_bot,'--','color','r','linew',2)
errorbar(1:12, mod_clim_s_bot,std_s_mod(3,:),'r', 'LineStyle', 'none')
errorbar(1:12, reg_clim_salt_b(1:12),std_salt_b(1:12),'b', 'LineStyle', 'none');
xlabel('time(month)','fontsize',13)
ylabel('salt (psu)','fontsize',13)
title('salt -polynomial(surf)','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([30 inf])
xlim([1 12])


%% no3
figure;
scatter(1:12,reg_clim_no3);
hold on
plot(1:12, mod_clim_no3_sur.*14,'--','color','r','linew',2)
errorbar(1:12, mod_clim_no3_sur.*14,std_no3_mod(1,:).*14,'r', 'LineStyle', 'none')
errorbar(1:12, reg_clim_no3(1:12),std_no3(1:12),'b', 'LineStyle', 'none');
xlabel('time(month)','fontsize',13)
ylabel('no3 (ug/L)','fontsize',13)
title('no3 -polynomial(surf)','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 250])
xlim([1 12])

figure;
scatter(1:12,reg_clim_no3);
hold on
plot(1:12, mod_clim_no3_bot.*14,'--','color','r','linew',2)
errorbar(1:12, mod_clim_no3_bot.*14,std_no3_mod(3,:).*14,'r', 'LineStyle', 'none')
errorbar(1:12, reg_clim_no3(1:12),std_no3_b(1:12),'b', 'LineStyle', 'none');
xlabel('time(month)','fontsize',13)
ylabel('no3 (ug/L)','fontsize',13)
title('no3 -polynomial(bot)','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 250])
xlim([1 12])




ylim([2 18]) 
xlim([1 366])



%% no3
figure;
scatter(1:366,reg_clim_no3 + regm_no3);
hold on
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('��������(ǥ��)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% temp
figure;
scatter(1:366,reg_clim_temp + regm_temp);
hold on
plot(1:366, yp_w_temp_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('temp (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial temp.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])



%bottom
%% DO
figure;
scatter(1:366,reg_clim_do_b + regm_do_b);
hold on
plot(1:366, yp_w_do_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('��������(����)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([2 18]) 
xlim([1 366])

%% salt
figure;
scatter(1:366,reg_clim_salt_b + regm_salt_b);
hold on
plot(1:366, yp_w_salt_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('salta (mg/m^3)','fontsize',13)
title('sumjin(songjung)-polynomial salt.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3_b + regm_no3_b);
hold on
plot(1:366, yp_w_no3_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('��������(����)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% temp
figure;
scatter(1:366,reg_clim_temp_b + regm_temp_b);
hold on
plot(1:366, yp_w_temp_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('temp (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial temp.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% ALL
figure;
hold on;
% plot(1:366, yp_w_do_04_b,'--','color','b','linew',2)
% plot(1:366, yp_w_salt_04_b./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','b','linew',2)
plot(1:366, yp_w_no3_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (ug/L)','fontsize',13)
title('��������(ǥ��&����)-polynomial.','fontsize',13)
grid on;
legend('no3-sur','no3-bot');
set(gca,'fontsize',13)
xlim([1 366])


