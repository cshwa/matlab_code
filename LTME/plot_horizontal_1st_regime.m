clc;close all;clear all;

%% -- read nc_file
cd /data1/cshwa/ext_hdd/gwangyang/model/result/2001/input
% /data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/grid_sumjin_v1970_fix_3m_v3.nc
file_path = 'grid_gy_v11_s.nc';
% x = ncread('ocean_avg_0001.nc','xi_rho');
% y = ncread('ocean_avg_0001.nc','eta_rho');
lon = ncread(file_path,'lon_rho');
lat = ncread(file_path,'lat_rho');
lon_u = ncread(file_path,'lon_u');
lat_u = ncread(file_path,'lat_u');
lon_v = ncread(file_path,'lon_v');
lat_v = ncread(file_path,'lat_v');
mask_r = ncread(file_path,'mask_rho');
mask_u = ncread(file_path,'mask_u');
mask_v = ncread(file_path,'mask_v');

cd /data1/cshwa/ext_hdd/gwangyang/model/result

% load 1980_1986_mean_mth_GY.mat;
k=0
for refyear=1991:2000
for i = 1:12
   mod_file=['/data1/cshwa/ext_hdd/gwangyang/model/result/',num2str(refyear),'/','mp_p_sewer_det_f_monthly_',num2str(refyear),'_',num2str(i,'%02d'),'.nc'];
if exist(mod_file,'file') ~= 2
    mod_file=[]; 
    mod_file=['/data1/cshwa/ext_hdd/gwangyang/model/result/',num2str(refyear),'/','mp_p_sewer_det_monthly_',num2str(refyear),'_',num2str(i,'%02d'),'.nc'];
if exist(mod_file,'file') ~= 2
    mod_file=[]; 
    mod_file=['/data1/cshwa/ext_hdd/gwangyang/model/result/',num2str(refyear),'/','mp_p_sewer_det_f_2nd_monthly_',num2str(refyear),'_',num2str(i,'%02d'),'.nc'];
end
end
k=k+1;
mod_u(:,:,:,k)=ncread(mod_file, 'u');
mod_v(:,:,:,k)=ncread(mod_file, 'v');
mod_no3(:,:,:,k)=ncread(mod_file, 'NO3');
mod_nh4(:,:,:,k)=ncread(mod_file, 'NH4');
mod_chl(:,:,:,k)=ncread(mod_file, 'chlorophyll');
mod_do(:,:,:,k)=ncread(mod_file, 'oxygen');
mod_temp(:,:,:,k)=ncread(mod_file, 'temp');
mod_salt(:,:,:,k)=ncread(mod_file, 'salt');
mod_po4(:,:,:,k)=ncread(mod_file, 'tPO4');
mod_zoo(:,:,:,k)=ncread(mod_file, 'zooplankton');
mod_phy(:,:,:,k)=ncread(mod_file, 'phytoplankton');     
end
end

% missing value will be NaN;
mod_no3(mod_no3 > 10^25) =NaN; 
mod_nh4(mod_nh4 > 10^25) =NaN; 
mod_do(mod_do > 10^25) =NaN; 
mod_chl(mod_chl > 10^25) =NaN; 
mod_temp(mod_temp > 10^25) =NaN; 
mod_salt(mod_salt > 10^25) =NaN; 
mod_po4(mod_po4 > 10^25) =NaN; 
mod_zoo(mod_zoo > 10^25) =NaN; 
mod_phy(mod_phy > 10^25) =NaN; 

for i = 1:12
clim_mod_po4(:,:,:,i)=nanmean(mod_po4(:,:,:,i:12:end),4);
clim_mod_no3(:,:,:,i)=nanmean(mod_no3(:,:,:,i:12:end),4);
clim_mod_nh4(:,:,:,i)=nanmean(mod_nh4(:,:,:,i:12:end),4);
clim_mod_chl(:,:,:,i)=nanmean(mod_chl(:,:,:,i:12:end),4);
clim_mod_temp(:,:,:,i)=nanmean(mod_temp(:,:,:,i:12:end),4);
clim_mod_salt(:,:,:,i)=nanmean(mod_salt(:,:,:,i:12:end),4);
clim_mod_u(:,:,:,i)=nanmean(mod_u(:,:,:,i:12:end),4);
clim_mod_v(:,:,:,i)=nanmean(mod_v(:,:,:,i:12:end),4);
end

clearvars *_365

min_lon = min(min(lon));
max_lon =  max(max(lon));
min_lat = min(min(lat));
max_lat = max(max(lat));
% min_lonc = min(lonc);
% max_lonc = max(lonc);
% min_latc = min(latc);
% max_latc = max(latc);
min_temp = min(clim_mod_temp);
max_temp = max(clim_mod_salt);
min_min_temp = min(clim_mod_salt);
max_max_temp = max(clim_mod_salt);

% latlim = [min_lat, max_lat];
% lonlim = [min_lon, max_lon];

for i = 1:20
    for j = 1:12
    clearvars U_1 V_1
        U_1 = griddata(lon_u,lat_u, squeeze(clim_mod_u(:,:,i,j)).*(mask_u./mask_u), lon, lat);
        V_1 = griddata(lon_v,lat_v, squeeze(clim_mod_v(:,:,i,j)).*(mask_v./mask_v), lon, lat);
        u_r(:,:,i,j) = U_1;
        v_r(:,:,i,j) = V_1;
    end
end

save('gy_clim_mod_1st_regime.mat','clim_mod_*','*_r');

return

close all; clear;clc;
load gy_clim_mod_1st_regime.mat


file_path = 'D:\장기생태\Dynamic\input\cshwa_2000\grid_gy_v11_s.nc';
% x = ncread('ocean_avg_0001.nc','xi_rho');
% y = ncread('ocean_avg_0001.nc','eta_rho');
lon = ncread(file_path,'lon_rho');
lat = ncread(file_path,'lat_rho');
lon_u = ncread(file_path,'lon_u');
lat_u = ncread(file_path,'lat_u');
lon_v = ncread(file_path,'lon_v');
lat_v = ncread(file_path,'lat_v');
mask_r = ncread(file_path,'mask_rho');
mask_u = ncread(file_path,'mask_u');
mask_v = ncread(file_path,'mask_v');

min_lon = min(min(lon));
max_lon =  max(max(lon));
min_lat = min(min(lat));
max_lat = max(max(lat));
% min_lonc = min(lonc);
% max_lonc = max(lonc);
% min_latc = min(latc);
% max_latc = max(latc);
min_temp = min(clim_mod_temp);
max_temp = max(clim_mod_salt);
min_min_temp = min(clim_mod_salt);
max_max_temp = max(clim_mod_salt);


% cblim = [6 31];
% xtick=linspace( floor(lonlim(1)),ceil(lonlim(2)),11 );
% ytick=linspace( floor(latlim(1)),ceil(latlim(2)),11 );

temp_range=[0:1:31];
salt_range=[15:1:35]; salt_range_c=[0:.5:35];
chl_range=[0:.5:20]; chl_range_c=[0:.5:5];
no3_range=[0:5:150]; no3_range_c=[0:5:60];
nh4_range=[0:1:29]; nh4_range_c=[0:1:12]; 
po4_range=[0:.1:17]; po4_range_c=[0:.2:3];

% referenc vector
scale_vec = 0.2;


clearvars *_ref
u_ref = zeros(size(lon,1),size(lon,2));
v_ref = zeros(size(lon,1),size(lon,2));
u_ref(11,141) = 0.1 * scale_vec;
% u_ref(161,21) = 0.1;

% save('plot_mth_1980_1986.mat');
% close all; clear all; clc;
% load('plot_mth_1980_1986.mat');

cd D:\장기생태\Dynamic\result\clim_pic_gy\1st_1991-2000

interval = 5;
for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. surf.');
% pcolor(lon,lat,squeeze(clim_mod_temp(:,:,20,i)).*(mask_r./mask_r)); shading flat;
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
    squeeze(v_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
% m_grid('linestyle','none','box','fancy','tickdir','in');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','r');
text(lon(11,131),lat(11,131),'0.5m/s','fontsize',15,'Color','r');

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(v).png'),'png');
close; 
end


for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. temp.(^oC) surf.');
pcolor(lon,lat,squeeze(clim_mod_temp(:,:,20,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_temp(:,:,20,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',13,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(temp_range) max(temp_range)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(t).png'),'png');
close; 
end

for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. salt.(psu) surf.');
pcolor(lon,lat,squeeze(clim_mod_salt(:,:,20,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_salt(:,:,20,i)),salt_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(salt_range_c) max(salt_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(s).png'),'png');
close; 
end


for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. chl.(ug/L) surf.');
pcolor(lon,lat,squeeze(clim_mod_chl(:,:,20,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_chl(:,:,20,i)),chl_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(chl_range_c) max(chl_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(chl).png'),'png');
close; 
end

for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. no3.(mmol N m^-^3) surf.');
pcolor(lon,lat,squeeze(clim_mod_no3(:,:,20,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_no3(:,:,20,i)),no3_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(no3_range_c) max(no3_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(no3).png'),'png');
close; 
end


for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. nh4.(mmol N m^-^3) surf.');
pcolor(lon,lat,squeeze(clim_mod_nh4(:,:,20,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_nh4(:,:,20,i)),nh4_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(nh4_range_c) max(nh4_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(nh4).png'),'png');
close; 
end

for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. po4.(mmol P m^-^3) surf.');
pcolor(lon,lat,squeeze(clim_mod_po4(:,:,20,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_po4(:,:,20,i)),po4_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(po4_range_c) max(po4_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(po4).png'),'png');
close; 
end



%% bottom
interval = 5;
for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. bot.');
% pcolor(lon,lat,squeeze(clim_mod_temp(:,:,1,i)).*(mask_r./mask_r)); shading flat;
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,1,i) .* scale_vec), ...
    squeeze(v_r(1:interval:end,1:interval:end,1,i) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
% m_grid('linestyle','none','box','fancy','tickdir','in');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','r');
text(lon(11,131),lat(11,131),'0.5m/s','fontsize',15,'Color','r');

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(v).png'),'png');
close; 
end


for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. temp.(^oC) bot.');
pcolor(lon,lat,squeeze(clim_mod_temp(:,:,1,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_temp(:,:,1,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',13,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(temp_range) max(temp_range)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(t).png'),'png');
close; 
end

for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. salt.(psu) bot.');
pcolor(lon,lat,squeeze(clim_mod_salt(:,:,1,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_salt(:,:,1,i)),salt_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(salt_range_c) max(salt_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(s).png'),'png');
close; 
end


for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. chl.(ug/L) bot.');
pcolor(lon,lat,squeeze(clim_mod_chl(:,:,1,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_chl(:,:,1,i)),chl_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(chl_range_c) max(chl_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(chl).png'),'png');
close; 
end

for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. no3.(mmol N m^-^3) bot.');
pcolor(lon,lat,squeeze(clim_mod_no3(:,:,1,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_no3(:,:,1,i)),no3_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(no3_range_c) max(no3_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(no3).png'),'png');
close; 
end


for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. nh4.(mmol N m^-^3) bot.');
pcolor(lon,lat,squeeze(clim_mod_nh4(:,:,1,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_nh4(:,:,1,i)),nh4_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(nh4_range_c) max(nh4_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(nh4).png'),'png');
close; 
end

for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1991~2000.',num2str(i,'%02d'), ' mon. po4.(mmol P m^-^3) bot.');
pcolor(lon,lat,squeeze(clim_mod_po4(:,:,1,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=contour(lon,lat,squeeze(clim_mod_po4(:,:,1,i)),po4_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',11,'labelspacing',700,'fontweight','bold');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(po4_range_c) max(po4_range_c)]);

ylim([34.723 35])
xlim([min(min(lon)) 127.9])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])

saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(po4).png'),'png');
close; 
end