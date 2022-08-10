
close all; clear; clc;

cd D:\������\yellow_sea
% /data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/grid_sumjin_v1970_fix_3m_v3.nc
file_path = 'roms_grd_auto_rdrg2_new8_smooth.nc';

lon = ncread(file_path,'lon_rho');
lat = ncread(file_path,'lat_rho');
lon_u = ncread(file_path,'lon_u');
lat_u = ncread(file_path,'lat_u');
lon_v = ncread(file_path,'lon_v');
lat_v = ncread(file_path,'lat_v');
mask_r = ncread(file_path,'mask_rho');
mask_u = ncread(file_path,'mask_u');
mask_v = ncread(file_path,'mask_v');

load yellow_clim_mod_2nd_regime.mat


cd D:\������\yellow_sea\input\ERA5
u=ncread('auto_ERA5_clim01-10_Uwind_monthly.nc','Uwind');
v=ncread('auto_ERA5_clim01-10_Vwind_monthly.nc','Vwind');

min_lon = min(min(lon));
max_lon =  max(max(lon));
min_lat = min(min(lat));
max_lat = max(max(lat));
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
scale_vec = 0.2;

clearvars *_ref
u_ref = zeros(size(lon,1),size(lon,2));
v_ref = zeros(size(lon,1),size(lon,2));
u_ref(20,120) = 5 * scale_vec;
% u_ref(161,21) = 0.1;
% u_ref(161,21) = 0.1;

% save('plot_mth_1980_1986.mat');
% close all; clear all; clc;
% load('plot_mth_1980_1986.mat');
cd /data1/cshwa/ext_hdd/yellow/result/clim_pic/2nd_2001-2010_4season_fulls

i_season{1}=['DJF'];
i_season{2}=['MAM'];
i_season{3}=['JJA'];
i_season{4}=['SON'];

interval = 10;
for i = 1:4
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('yellow-',' 2001~2010.',i_season{i}, '  surf.');
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

interval = 15;
for i = 1:4
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('yellow-',' 2001~2010.',i_season{i}, '  surf.');
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_pcolor(lon,lat,squeeze(clim_mod_salt(:,:,40,i)).*(mask_r./mask_r)); shading flat;
m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u(1:interval:end,1:interval:end,i) .* scale_vec), ...
    squeeze(v(1:interval:end,1:interval:end,i) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');


[cs,h]=m_contour(lon,lat,squeeze(clim_mod_salt(:,:,40,i)),salt_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',13,'labelspacing',700,'fontweight','bold');
m_grid('linestyle','none','box','fancy','tickdir','in');
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);

m_quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','w');
m_text(lon(20,100),lat(20,100),'5 m/s','fontsize',15,'Color','w');
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
saveas(gcf,strcat('monthly_mean_',i_season{i},'m_(sv).png'),'png');
% close; 
end

