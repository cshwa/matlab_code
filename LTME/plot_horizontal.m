clc;close all;clear all;
%% -- read nc_file
% /data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/grid_sumjin_v1970_fix_3m_v3.nc
file_path = '/data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/grid_sumjin_v1970_fix_3m_v3.nc';
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
load 1980_1986_mean_mth_all_GY.mat

clearvars *_365

min_lon = min(min(lon));
max_lon =  max(max(lon));
min_lat = min(min(lat));
max_lat = max(max(lat));
% min_lonc = min(lonc);
% max_lonc = max(lonc);
% min_latc = min(latc);
% max_latc = max(latc);
min_temp = min(temp_mth);
max_temp = max(temp_mth);
min_min_temp = min(temp_mth);
max_max_temp = max(temp_mth);
% latlim = [min_lat, max_lat];
% lonlim = [min_lon, max_lon];

for i = 1:20
    for j = 1:12
    clearvars U_1 V_1
        U_1 = griddata(lon_u,lat_u, squeeze(u_mth(:,:,i,j)).*(mask_u./mask_u), lon, lat);
        V_1 = griddata(lon_v,lat_v, squeeze(v_mth(:,:,i,j)).*(mask_v./mask_v), lon, lat);
        u_r(:,:,i,j) = U_1;
        v_r(:,:,i,j) = V_1;
    end
end

% cblim = [6 31];
% xtick=linspace( floor(lonlim(1)),ceil(lonlim(2)),11 );
% ytick=linspace( floor(latlim(1)),ceil(latlim(2)),11 );

% referenc vector
clearvars *_ref
u_ref = zeros(size(lon,1),size(lon,2));
v_ref = zeros(size(lon,1),size(lon,2));
u_ref(11,161) = 0.1;
% u_ref(161,21) = 0.1;

% save('plot_mth_1980_1986.mat');
close all; clear all; clc;
load('plot_mth_1980_1986.mat');

interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 1980~1986.',num2str(i,'%02d'), ' mon. surf.');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in');
pcolor(lon,lat,squeeze(temp_mth(:,:,20,i)).*(mask_r./mask_r)); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,20,i)), ...
    squeeze(v_r(1:interval:end,1:interval:end,20,i)), ...
    'AutoScale','off','LineWidth',1,'Color','k');
quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 30]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(t,v).png'),'png');
% close; 
end

interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 1980~1986.',num2str(i,'%02d'), ' mon. surf. salt');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in');
pcolor(lon,lat,squeeze(salt_mth(:,:,20,i)).*(mask_r./mask_r)); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,20,i)), ...
    squeeze(v_r(1:interval:end,1:interval:end,20,i)), ...
    'AutoScale','off','LineWidth',1,'Color','k');
quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 35]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(s,v).png'),'png');
% close; 
end

interval = 10;
for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1980~1986.',num2str(i,'%02d'), ' mon. zeta(m) surf.');
pcolor(lon,lat,squeeze(zeta_mth(:,:,i)).*(mask_r./mask_r)); shading flat;
% quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
%     squeeze(u_r(1:interval:end,1:interval:end,20,i)), ...
%     squeeze(v_r(1:interval:end,1:interval:end,20,i)), ...
%     'AutoScale','off','LineWidth',1,'Color','k');
% quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-0.05 0.05]);colormap(jet);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
% text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
% saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(z,v).png'),'png');
saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(zeta).png'),'png');
% close; 
end

% bottom
interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 1980~1986.',num2str(i,'%02d'), ' mon. bot.');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in');
pcolor(lon,lat,squeeze(temp_mth(:,:,1,i)).*(mask_r./mask_r)); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,1,i)), ...
    squeeze(v_r(1:interval:end,1:interval:end,1,i)), ...
    'AutoScale','off','LineWidth',1,'Color','k');
quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 30]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(t,v).png'),'png');
% close; 
end

interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 1980~1986.',num2str(i,'%02d'), ' mon. bot. salt');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in');
pcolor(lon,lat,squeeze(salt_mth(:,:,1,i)).*(mask_r./mask_r)); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,1,i)), ...
    squeeze(v_r(1:interval:end,1:interval:end,1,i)), ...
    'AutoScale','off','LineWidth',1,'Color','k');
quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 35]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_bot(s,v).png'),'png');
% close; 
end





for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 1980~1986.',num2str(i,'%02d'), ' monthly bottom mean ROMS');
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_grid('box','fancy','tickdir','in'); 
m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), squeeze(u_mth(1:interval:end,1:interval:end,1,i)), squeeze(v_mth(1:interval:end,1:interval:end,1,i)), 4,'LineWidth',1,'Color','k');
m_pcolor(lon,lat,squeeze(temp_mth(:,:,1,i))); 
shading interp;
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
% caxis([0 36]);
saveas(gcf,strcat('monthly_mean_bot',num2str(i,'%02d'),'m_(t,v).png'),'png');
close; 
end


%%%%%%%%   sst only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;close all;clear all;
%% -- read nc_file

temp_file = strcat('monthly_spinup_*.nc');
d_filename = dir(temp_file); % make a list of nc in the current folder.
for i = 1:12
file_num = num2str(i);
temp = ncread(d_filename(i).name,'temp');
lon = ncread(d_filename(i).name,'lon_rho');
lat = ncread(d_filename(i).name,'lat_rho');
lon_u = ncread(d_filename(i).name,'lon_u');
lat_u = ncread(d_filename(i).name,'lat_u');
lon_v = ncread(d_filename(i).name,'lon_v');
lat_v = ncread(d_filename(i).name,'lat_v');
u = ncread(d_filename(i).name,'u');
v = ncread(d_filename(i).name,'v');
z = ncread(d_filename(i).name,'s_rho'); % depth

u = squeeze(u(:,:,30));
v = squeeze(v(:,:,30));
sst = temp(:,:,30); % diminishing the dimension
squeeze(sst);

index = find(sst >= 1000);
index_cold = find(sst <= -1);
reshape(sst,1988616,1);
sst(index) = NaN;
reshape(sst,1862,1068);

index_U = find(u >= 1000);
reshape(u,1987548,1);
u(index_U) = NaN;
reshape(u,1861,1068);

index_V = find(v >= 1000);
reshape(v,1986754,1);
v(index_V) = NaN;
reshape(v,1862,1067);

min_lon = min(lon);
max_lon = max(lon);
min_lat = min(lat);
max_lat = max(lat);
% min_lonc = min(lonc);
% max_lonc = max(lonc);
% min_latc = min(latc);
% max_latc = max(latc);
min_temp = min(sst);
max_temp = max(sst);
min_min_temp = min(sst);
max_max_temp = max(sst);
latlim = [min_lat, max_lat];
lonlim = [min_lon, max_lon];



% cblim = [6 31];
xtick=linspace( floor(lonlim(1)),ceil(lonlim(2)),11 );
ytick=linspace( floor(latlim(1)),ceil(latlim(2)),11 );

[lon1,lat1,v_1] = griddata(lon_v, lat_v, v, lon_u,lat_u,);

figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('NP-',' 2007.',num2str(file_num), 'm - 4yr monthly mean ROMS');
% m_proj('mercator','lon',[98 284], 'lat',[-20 65]);
m_proj('mercator','lon',[115 164], 'lat',[15 52]);
m_grid('box','fancy','tickdir','in');
m_quiver(lon_u(1:10:end,1:10:end),lat_u(1:10:end,1:10:end), u(1:10:end,1:10:end), v_1(1:10:end,1:10:end), 4,'LineWidth',1,'Color','k');
m_quiver(lon_u(1:3:end,1:3:end),lat_u(1:3:end,1:3:end), u(1:3:end,1:3:end), v_1(1:3:end,1:3:end), 4,'LineWidth',1,'Color','k');
m_pcolor(lon,lat,sst); 
shading interp;
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 36]);
saveas(gcf,strcat('NWP_monthly_mean_4yr_',num2str(file_num),'m_temp.png'),'png');
close; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot NP% %%%%%%%%
clc;close all;clear all;
%% -- read nc_file

temp_file = strcat('monthly_spinup_*.nc');
xscale = 75; %% size of image (inches)
yscale = 42.5; %% size of image (inches)
d_filename = dir(temp_file); % make a list of nc in the current folder.
for i = 1:12
file_num = num2str(i);
temp = ncread(d_filename(i).name,'temp');
lon = ncread(d_filename(i).name,'lon_rho');
lat = ncread(d_filename(i).name,'lat_rho');
lon_u = ncread(d_filename(i).name,'lon_u');
lat_u = ncread(d_filename(i).name,'lat_u');
lon_v = ncread(d_filename(i).name,'lon_v');
lat_v = ncread(d_filename(i).name,'lat_v');
u = ncread(d_filename(i).name,'u');
v = ncread(d_filename(i).name,'v');
z = ncread(d_filename(i).name,'s_rho'); % depth

u = squeeze(u(:,:,30));
v = squeeze(v(:,:,30));
sst = temp(:,:,30); % diminishing the dimension
squeeze(sst);

index = find(sst >= 1000);
index_cold = find(sst <= -1);
reshape(sst,1988616,1);
sst(index) = NaN;
reshape(sst,1862,1068);

index_U = find(u >= 1000);
reshape(u,1987548,1);
u(index_U) = NaN;
reshape(u,1861,1068);

index_V = find(v >= 1000);
reshape(v,1986754,1);
v(index_V) = NaN;
reshape(v,1862,1067);

min_lon = min(lon);
max_lon = max(lon);
min_lat = min(lat);
max_lat = max(lat);
% min_lonc = min(lonc);
% max_lonc = max(lonc);
% min_latc = min(latc);
% max_latc = max(latc);
min_temp = min(sst);
max_temp = max(sst);
min_min_temp = min(sst);
max_max_temp = max(sst);
latlim = [min_lat, max_lat];
lonlim = [min_lon, max_lon];



% cblim = [6 31];
xtick=linspace( floor(lonlim(1)),ceil(lonlim(2)),11 );
ytick=linspace( floor(latlim(1)),ceil(latlim(2)),11 );

[lon1,lat1,v_1] = griddata(lon_v, lat_v, v, lon_u,lat_u);

figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('NP-',' 2007.',num2str(file_num), 'm - 4yr monthly mean ROMS');
m_proj('mercator','lon',[115 164], 'lat',[15 52]);
m_grid('box','fancy','tickdir','in');
% m_quiver(lon_u(1:10:end,1:10:end),lat_u(1:10:end,1:10:end),
% u(1:10:end,1:10:end), v_1(1:10:end,1:10:end), 2,'LineWidth',0.5,'Color','k');
m_quiver(lon_u(1:3:end,1:3:end),lat_u(1:3:end,1:3:end), u(1:3:end,1:3:end), v_1(1:3:end,1:3:end), 4,'LineWidth',1,'Color','k');
m_pcolor(lon,lat,sst); 
shading interp;
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 36]);
set(gcf,'PaperPosition',[0 0 xscale yscale]);
saveas(gcf,strcat('NWP_monthly_mean_4yr_',num2str(file_num),'m_temp.png'),'png');
close; 
end

%%%%%%%% check temp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;close all;clear all;
%% -- read nc_file
for filenum = 1:12
    
temp_file = strcat(num2str(filenum),'month.mat');

load(temp_file);

index = find(monthly_mean >= 1000);
index_cold = find(monthly_mean <= -1);
reshape(monthly_mean,1988616,1);
monthly_mean(index) = NaN;
reshape(monthly_mean,1862,1068);

min_temp = min(monthly_mean);
max_temp = max(monthly_mean);
filenum
min_min_temp = min(min_temp)
max_max_temp = max(max_temp)
end