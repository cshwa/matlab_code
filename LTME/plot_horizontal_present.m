clc;close all;clear all;
%% -- read nc_file
% /data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/grid_sumjin_v1970_fix_3m_v3.nc
file_path = '/home/cshwa/GY_compare/grid_gy_v11_s.nc';
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

% lv = ncread(file_path,'s_rho'); % depth
% temp_file = strcat(num2str(filenum),'month.mat');
% u_file = strcat('U_',num2str(filenum),'month.mat');
% v_file = strcat('V_',num2str(filenum),'month.mat');
load 2011_2016_uv_mean_mth_GY.mat
load 2011_2016_TS_mean_mth_GY.mat
% load(u_file);
% load(v_file);
% U = U_monthly_mean;
% V = V_monthly_mean;
% temp = ncread('ocean_avg_1476.nc','temp'); 
% t = ncread('ocean_avg_1476.nc','ocean_time');
% temp = temp(:,:,30); % diminishing the dimension
% index = find(monthly_mean >= 1000);
% index_cold = find(monthly_mean <= -1);
% reshape(monthly_mean,1988616,1);
% monthly_mean(index) = NaN;
% reshape(monthly_mean,1862,1068);
% 
% index_U = find(U >= 1000);
% reshape(U,1987548,1);
% U(index_U) = NaN;
% reshape(U,1861,1068);
% 
% index_V = find(V >= 1000);
% reshape(V,1986754,1);
% V(index_V) = NaN;
% reshape(V,1862,1067);

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

save('plot_mth_2011_2016.mat');

close all; clear all; clc;
load plot_mth_2011_2016.mat


interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 2011~2016.',num2str(i,'%02d'), ' mon. surf.');
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
saveas(gcf,strcat('monthly_mean_2011_',num2str(i,'%02d'),'m_(t,v).png'),'png');
close; 
end

interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 2011~2016.',num2str(i,'%02d'), ' mon. surf. salt');
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
saveas(gcf,strcat('monthly_mean_2011_',num2str(i,'%02d'),'m_(s,v).png'),'png');
close; 
end

interval = 10;
for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 2011~2016.',num2str(i,'%02d'), ' mon. zeta(m) surf.');
pcolor(lon,lat,squeeze(zeta_mth(:,:,i)).*(mask_r./mask_r)); shading flat;
% quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
%     squeeze(u_r(1:interval:end,1:interval:end,20,i)), ...
%     squeeze(v_r(1:interval:end,1:interval:end,20,i)), ...
%     'AutoScale','off','LineWidth',1,'Color','k');
% quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-0.05 0.05]); colormap(jet);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
% text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
% saveas(gcf,strcat('monthly_mean_2011_',num2str(i,'%02d'),'m_(z,v).png'),'png');
saveas(gcf,strcat('monthly_mean_2011_',num2str(i,'%02d'),'m_(zeta).png'),'png');
% close; 
end

% bottom
interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 2011~2016.',num2str(i,'%02d'), ' mon. bot.');
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
saveas(gcf,strcat('monthly_mean_2011_',num2str(i,'%02d'),'m_bot(t,v).png'),'png');
close; 
end

interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' 2011~2016.',num2str(i,'%02d'), ' mon. bot. salt');
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
saveas(gcf,strcat('monthly_mean_2011_',num2str(i,'%02d'),'m_bot(s,v).png'),'png');
close; 
end


% 
% for i = 1:12
% figure; hold on;
% %     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
% title_name = strcat('GY-',' diff.',num2str(i,'%02d'), ' monthly bottom mean ROMS');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in'); 
% m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), squeeze(u_mth(1:interval:end,1:interval:end,1,i)), squeeze(v_mth(1:interval:end,1:interval:end,1,i)), 4,'LineWidth',1,'Color','k');
% m_pcolor(lon,lat,squeeze(temp_mth(:,:,1,i))); 
% shading interp;
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% title(title_name,'fontsize',30);set(gca,'FontSize',18);
% xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
% colorbar;
% % caxis([0 36]);
% saveas(gcf,strcat('monthly_mean_bot_2011_',num2str(i,'%02d'),'m_(t,v).png'),'png');
% close; 
% end
% 

%%%%%%%%   diff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;close all;clear all;
%% -- read file

load('plot_mth_1980_1986.mat');
temp_mth_80=temp_mth;
salt_mth_80=salt_mth;
u_mth_80=u_mth;
v_mth_80=v_mth;
u_r_80 = u_r;
v_r_80 = v_r;
zeta_80 = zeta_mth;
mask_r_80 = mask_r;
mask_u_80 = mask_u;
mask_v_80 = mask_v;

clearvars *_mth u_r v_r
load('plot_mth_2011_2016.mat');

temp_temp=temp_mth;
salt_temp=salt_mth;
u_temp=u_r;
v_temp=v_r;
zeta_temp = zeta_mth;

for i  = 1:20
   for j = 1:12
    temp_mth(:,:,i,j) = squeeze(temp_temp(:,:,i,j)) .*(mask_r./mask_r);
    salt_mth(:,:,i,j) = squeeze(salt_temp(:,:,i,j)) .*(mask_r./mask_r);
    u_r(:,:,i,j) = squeeze(u_temp(:,:,i,j)) .*(mask_r./mask_r);
    v_r(:,:,i,j) = squeeze(v_temp(:,:,i,j)) .*(mask_r./mask_r);
   end
end

clearvars zeta_mth
for i = 1:12
    zeta_mth(:,:,i) = squeeze(zeta_temp(:,:,i)) .* (mask_r./mask_r);
end
temp_mth(find(isnan(temp_mth) == 1)) = 0;
salt_mth(find(isnan(salt_mth) == 1)) = 0;
u_r(find(isnan(u_r) == 1)) = 0;
v_r(find(isnan(v_r) == 1)) = 0;

temp_diff = temp_mth - temp_mth_80;
salt_diff = salt_mth - salt_mth_80;
u_diff = u_r - u_r_80;
v_diff = v_r - v_r_80;
zeta_diff = zeta_mth - zeta_80;


interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' diff.',num2str(i,'%02d'), ' mon. surf.');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in');
pcolor(lon,lat,squeeze(temp_diff(:,:,20,i)).*(mask_r_80./mask_r_80)); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,20,i)), ...
    squeeze(v_r(1:interval:end,1:interval:end,20,i)), ...
    'AutoScale','off','LineWidth',1,'Color','k');
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r_80(1:interval:end,1:interval:end,20,i)), ...
    squeeze(v_r_80(1:interval:end,1:interval:end,20,i)), ...
    'AutoScale','off','LineWidth',1,'Color','r');
quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-2 2]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
saveas(gcf,strcat('monthly_mean_diff2_',num2str(i,'%02d'),'m_(t,v).png'),'png');
% close; 
end

interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' diff.',num2str(i,'%02d'), ' mon. surf. salt');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in');
pcolor(lon,lat,squeeze(salt_diff(:,:,20,i)).*(mask_r_80./mask_r_80)); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,20,i)), ...
    squeeze(v_r(1:interval:end,1:interval:end,20,i)), ...
    'AutoScale','off','LineWidth',1,'Color','k');
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r_80(1:interval:end,1:interval:end,20,i)), ...
    squeeze(v_r_80(1:interval:end,1:interval:end,20,i)), ...
    'AutoScale','off','LineWidth',1,'Color','r');
quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-2 2]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
saveas(gcf,strcat('monthly_mean_diff2_',num2str(i,'%02d'),'m_(s,v).png'),'png');
% close; 
end

interval = 10;
for i = 1:12
figure; hold on;
title_name = strcat('GY-',' diff.',num2str(i,'%02d'), ' mon. zeta(m) surf.');
pcolor(lon,lat,squeeze(zeta_diff(:,:,i)).*(mask_r_80./mask_r_80)); shading flat;
% quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
%     squeeze(u_r(1:interval:end,1:interval:end,20,i)), ...
%     squeeze(v_r(1:interval:end,1:interval:end,20,i)), ...
%     'AutoScale','off','LineWidth',1,'Color','k');
% quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-0.05 0.05]); colormap(jet);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
% text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
% saveas(gcf,strcat('monthly_mean_diff2_',num2str(i,'%02d'),'m_(z,v).png'),'png');
saveas(gcf,strcat('monthly_mean_diff2_',num2str(i,'%02d'),'m_(zeta).png'),'png');
% close; 
end


% bottom
interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' diff.',num2str(i,'%02d'), ' mon. bot.');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in');
pcolor(lon,lat,squeeze(temp_diff(:,:,1,i)).*(mask_r_80./mask_r_80)); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,1,i)), ...
    squeeze(v_r(1:interval:end,1:interval:end,1,i)), ...
    'AutoScale','off','LineWidth',1,'Color','k');
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r_80(1:interval:end,1:interval:end,1,i)), ...
    squeeze(v_r_80(1:interval:end,1:interval:end,1,i)), ...
    'AutoScale','off','LineWidth',1,'Color','r');
quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-2 2]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
saveas(gcf,strcat('monthly_mean_diff2_',num2str(i,'%02d'),'m_bot(t,v).png'),'png');
% close; 
end

interval = 10;
for i = 1:12
figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
title_name = strcat('GY-',' diff.',num2str(i,'%02d'), ' mon. bot. salt');
% m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
% m_grid('box','fancy','tickdir','in');
pcolor(lon,lat,squeeze(salt_diff(:,:,1,i)).*(mask_r_80./mask_r_80)); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r(1:interval:end,1:interval:end,1,i)), ...
    squeeze(v_r(1:interval:end,1:interval:end,1,i)), ...
    'AutoScale','off','LineWidth',1,'Color','k');
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(u_r_80(1:interval:end,1:interval:end,1,i)), ...
    squeeze(v_r_80(1:interval:end,1:interval:end,1,i)), ...
    'AutoScale','off','LineWidth',1,'Color','r');
quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','k');
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-2 2]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
text(lon(11,151),lat(11,151),'0.1m/s','fontsize',16);
saveas(gcf,strcat('monthly_mean_diff2_',num2str(i,'%02d'),'m_bot(s,v).png'),'png');
% close; 
end


% speed plot
speed=sqrt((u_r.*u_r)+(v_r.*v_r));
speed_80=sqrt((u_r_80.*u_r_80)+(v_r_80.*v_r_80));
speed_diff = speed - speed_80;



interval = 10;
for i = 1:12
figure; hold on;
title_name = strcat('GY-','speed diff.',num2str(i,'%02d'), ' mon. surf.');
pcolor(lon,lat,squeeze(speed_diff(:,:,20,i)).*(mask_r_80./mask_r_80)); shading flat;
% contour(lon,lat,squeeze(speed_diff(:,:,20,i)),'k','showtext','on');
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-.2 .2]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
colormap(jet)
saveas(gcf,strcat('monthly_mean_diff2_',num2str(i,'%02d'),'m_(speed).png'),'png');
% close; 
end

% bottom
interval = 10;
for i = 1:12
figure; hold on;
title_name = strcat('GY-','speed diff.',num2str(i,'%02d'), ' mon. bot.');
pcolor(lon,lat,squeeze(speed_diff(:,:,1,i)).*(mask_r_80./mask_r_80)); shading flat;
% contour(lon,lat,squeeze(speed_diff(:,:,20,i)),'k','showtext','on');
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-.2 .2]);
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
colormap(jet)
saveas(gcf,strcat('monthly_mean_diff2_',num2str(i,'%02d'),'m_bot(speed).png'),'png');
% close; 
end

