clc;close all;clear all;

%% -- read nc_file
cd /data1/cshwa/ext_hdd/yellow/result
% /data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/grid_sumjin_v1970_fix_3m_v3.nc
file_path = 'roms_grd_auto_rdrg2_new8_smooth.nc';
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


% load 1980_1986_mean_mth_GY.mat;
k=0
for refyear=1991:2000
for i = 1:12
   mod_file=['/data1/cshwa/ext_hdd/yellow/result/',num2str(refyear),'/','yellow_monthly_',num2str(refyear),'_',num2str(i,'%02d'),'.nc'];
if exist(mod_file,'file') ~= 2
    mod_file=[]; 
    mod_file=['/data1/cshwa/ext_hdd/yellow/result/',num2str(refyear),'/','spinup1_monthly',num2str(i,'%02d'),'.nc'];
end
k=k+1;
mod_u(:,:,:,k)=ncread(mod_file, 'u');
mod_v(:,:,:,k)=ncread(mod_file, 'v');
% mod_no3(:,:,:,i)=ncread(mod_file, 'NO3');
% mod_nh4(:,:,:,i)=ncread(mod_file, 'NH4');
% mod_chl(:,:,:,i)=ncread(mod_file, 'chlorophyll');
% mod_do(:,:,:,i)=ncread(mod_file, 'oxygen');
mod_temp(:,:,:,k)=ncread(mod_file, 'temp');
mod_salt(:,:,:,k)=ncread(mod_file, 'salt');
% mod_po4(:,:,:,i)=ncread(mod_file, 'tPO4');
% mod_zoo(:,:,:,i)=ncread(mod_file, 'zooplankton');
% mod_phy(:,:,:,i)=ncread(mod_file, 'phytoplankton');     
end
end

% missing value will be NaN;
% mod_no3(mod_no3 > 10^25) =NaN; 
% mod_nh4(mod_nh4 > 10^25) =NaN; 
% mod_do(mod_do > 10^25) =NaN; 
% mod_chl(mod_chl > 10^25) =NaN; 
mod_temp(mod_temp > 10^25) =NaN; 
mod_salt(mod_salt > 10^25) =NaN; 
% mod_po4(mod_po4 > 10^25) =NaN; 
% mod_zoo(mod_zoo > 10^25) =NaN; 
% mod_phy(mod_phy > 10^25) =NaN; 

for i = 1:12
% clim_mod_no3(:,:,:,i)=nanmean(mod_no3(:,:,:,i:12:end),4);
% clim_mod_nh4(:,:,:,i)=nanmean(mod_nh4(:,:,:,i:12:end),4);
% clim_mod_chl(:,:,:,i)=nanmean(mod_chl(:,:,:,i:12:end),4);
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

for i = 1:40
    for j = 1:12
    clearvars U_1 V_1
        U_1 = griddata(lon_u,lat_u, squeeze(clim_mod_u(:,:,i,j)).*(mask_u./mask_u), lon, lat);
        V_1 = griddata(lon_v,lat_v, squeeze(clim_mod_v(:,:,i,j)).*(mask_v./mask_v), lon, lat);
        u_r(:,:,i,j) = U_1;
        v_r(:,:,i,j) = V_1;
    end
end

save('yellow_clim_mod_1st_regime.mat','clim_mod_*','*_r');

return

close all;clear all; clc;
load yellow_clim_mod_1st_regime.mat


% /data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/grid_sumjin_v1970_fix_3m_v3.nc
file_path = '/data1/cshwa/ext_hdd/yellow/result/roms_grd_auto_rdrg2_new8_smooth.nc';

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

% 4seasons data
seasons4{1} = [1,2,12];
seasons4{2} = [3,4,5];
seasons4{3} = [6,7,8];
seasons4{4} = [9,10,11];

% temp_po4 = clim_mod_po4;
% temp_no3=clim_mod_no3;
% temp_nh4=clim_mod_nh4;
% temp_chl=clim_mod_chl;
temp_temp=clim_mod_temp;
temp_salt=clim_mod_salt;
temp_u=clim_mod_u;
temp_v=clim_mod_v;

clearvars clim_mod_*
for i = 1:4
% clim_mod_po4(:,:,:,i)=nanmean(temp_po4(:,:,:,seasons4{i}),4);
% clim_mod_no3(:,:,:,i)=nanmean(temp_no3(:,:,:,seasons4{i}),4);
% clim_mod_nh4(:,:,:,i)=nanmean(temp_nh4(:,:,:,seasons4{i}),4);
% clim_mod_chl(:,:,:,i)=nanmean(temp_chl(:,:,:,seasons4{i}),4);
clim_mod_temp(:,:,:,i)=nanmean(temp_temp(:,:,:,seasons4{i}),4);
clim_mod_salt(:,:,:,i)=nanmean(temp_salt(:,:,:,seasons4{i}),4);
clim_mod_u(:,:,:,i)=nanmean(temp_u(:,:,:,seasons4{i}),4);
clim_mod_v(:,:,:,i)=nanmean(temp_v(:,:,:,seasons4{i}),4);
end


temp_range=[0:1:31];
salt_range=[30:.5:35]; salt_range_c=[30:.5:35];
% nh4_range=[0:0.5:20];
% no3_range=[0:1:10];

% referenc vector
scale_vec = 4;

clearvars *_ref
u_ref = zeros(size(lon,1),size(lon,2));
v_ref = zeros(size(lon,1),size(lon,2));
u_ref(20,120) = 0.5 * scale_vec;
% u_ref(161,21) = 0.1;


% save('plot_mth_1980_1986.mat');
% close all; clear all; clc;
% load('plot_mth_1980_1986.mat');
cd /data1/cshwa/ext_hdd/yellow/result/clim_pic/1st_1991-2000_4season_fulls

i_season{1}=['DJF'];
i_season{2}=['MAM'];
i_season{3}=['JJA'];
i_season{4}=['SON'];


interval = 10;
for i = 1:4
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('yellow-',' 1991~2000.',i_season{i}, '  surf.');
% pcolor(lon,lat,squeeze(clim_mod_temp(:,:,40,i)).*(mask_r./mask_r)); shading flat;
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,40,i) .* scale_vec), ...
    squeeze(v_r(1:interval:end,1:interval:end,40,i) .* scale_vec), ...
    'AutoScale','off','LineWidth',1,'Color','k');

m_grid('linestyle','none','box','fancy','tickdir','in');
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);

m_quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','w');
m_text(lon(20,100),lat(20,100),'0.5m/s','fontsize',15,'Color','w');
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
% xlim([min_lon max_lon]);
% ylim([min_lat max_lat]);
saveas(gcf,strcat('monthly_mean_',i_season{i},'m_(v).png'),'png');
close; 
end


for i = 1:4
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('yellow-',' 1991~2000.',i_season{i}, '  surf.');
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_pcolor(lon,lat,squeeze(clim_mod_temp(:,:,40,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=m_contour(lon,lat,squeeze(clim_mod_temp(:,:,40,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',13,'labelspacing',700,'fontweight','bold');
m_grid('linestyle','none','box','fancy','tickdir','in');
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
%     squeeze(u_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
%     squeeze(v_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
%     'AutoScale','off','LineWidth',1,'Color','k');
% quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(temp_range) max(temp_range)]);
% text(lon(20,100),lat(20,100),'0.5m/s','fontsize',16);
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
saveas(gcf,strcat('monthly_mean_',i_season{i},'m_(t).png'),'png');
close; 
end


for i = 1:4
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('yellow-',' 1991~2000.',i_season{i}, '  surf.');
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_pcolor(lon,lat,squeeze(clim_mod_salt(:,:,40,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=m_contour(lon,lat,squeeze(clim_mod_salt(:,:,40,i)),salt_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',13,'labelspacing',700,'fontweight','bold');
m_grid('linestyle','none','box','fancy','tickdir','in');
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
%     squeeze(u_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
%     squeeze(v_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
%     'AutoScale','off','LineWidth',1,'Color','k');
% quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(salt_range_c) max(salt_range_c)]);
% text(lon(20,100),lat(20,100),'0.5m/s','fontsize',16);
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
saveas(gcf,strcat('monthly_mean_',i_season{i},'m_(s).png'),'png');
% close; 
end




% bottom

interval = 3;
for i = 1:4
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('yellow-',' 1991~2000.',i_season{i}, '  bot.');
% pcolor(lon,lat,squeeze(clim_mod_temp(:,:,1,i)).*(mask_r./mask_r)); shading flat;
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,1,i) .* scale_vec), ...
    squeeze(v_r(1:interval:end,1:interval:end,1,i) .* scale_vec), ...
    'AutoScale','off','LineWidth',1,'Color','k');
m_grid('linestyle','none','box','fancy','tickdir','in');
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);

m_quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','w');
m_text(lon(20,100),lat(20,100),'0.5m/s','fontsize',15,'Color','w');
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
% xlim([min_lon max_lon]);
% ylim([min_lat max_lat]);
saveas(gcf,strcat('monthly_mean_',i_season{i},'m_bot(v).png'),'png');
close; 
end


for i = 1:4
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('yellow-',' 1991~2000.',i_season{i}, '  bot.');
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_pcolor(lon,lat,squeeze(clim_mod_temp(:,:,1,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=m_contour(lon,lat,squeeze(clim_mod_temp(:,:,1,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',13,'labelspacing',700,'fontweight','bold');
m_grid('linestyle','none','box','fancy','tickdir','in');
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
%     squeeze(u_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
%     squeeze(v_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
%     'AutoScale','off','LineWidth',1,'Color','k');
% quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(temp_range) max(temp_range)]);
% text(lon(20,100),lat(20,100),'0.5m/s','fontsize',16);
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
saveas(gcf,strcat('monthly_mean_',i_season{i},'m_bot(t).png'),'png');
close; 
end


for i = 1:4
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('yellow-',' 1991~2000.',i_season{i}, '  bot.');
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_pcolor(lon,lat,squeeze(clim_mod_salt(:,:,1,i)).*(mask_r./mask_r)); shading flat;
[cs,h]=m_contour(lon,lat,squeeze(clim_mod_salt(:,:,1,i)),salt_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',13,'labelspacing',700,'fontweight','bold');
m_grid('linestyle','none','box','fancy','tickdir','in');
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
%     squeeze(u_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
%     squeeze(v_r(1:interval:end,1:interval:end,20,i) .* scale_vec), ...
%     'AutoScale','off','LineWidth',1,'Color','k');
% quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');

title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([min(salt_range_c) max(salt_range_c)]);
% text(lon(20,100),lat(20,100),'0.5m/s','fontsize',16);
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
saveas(gcf,strcat('monthly_mean_',i_season{i},'m_bot(s).png'),'png');
close; 
end