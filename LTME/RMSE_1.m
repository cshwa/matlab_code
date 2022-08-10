clear all; clc; close all; 

temp_file = strcat('monthly_spinup_*.nc');
d_mat = dir('2007*')
d_filename = dir(temp_file); % make a list of nc in the current folder.
load location.mat

lon_ob = ncread('avhrr-only-v2.20100101.nc','lon');
lat_ob = ncread('avhrr-only-v2.20100101.nc','lat');

for i = 1:12
file_num = num2str(i);
temp = ncread(d_filename(i).name,'temp');
lon = ncread(d_filename(i).name,'lon_rho');
lat = ncread(d_filename(i).name,'lat_rho');
z = ncread(d_filename(i).name,'s_rho'); % depth
load(d_mat(i).name);

sst = temp(:,:,30); % diminishing the dimension
squeeze(sst);

index = find(sst >= 1000);
index_cold = find(sst <= -1);
reshape(sst,1988616,1);
sst(index) = NaN;
reshape(sst,1862,1068);

index1 = find(monthly_mean <= -9);
monthly_mean(index1) = NaN;

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

[lat1,lon1] = meshgrid(lat_ob,lon_ob);

[lon2,lat2,mth_mean] = griddata(lon1,lat1, monthly_mean, lon, lat);
% [lon2,lat2,V_1] = griddata(lon_v,lat_v, V, lon, lat);

index_nan = find(isnan(mth_mean) == 1);

% mth_mean(index_nan) = -999.9;

RMSE = sqrt((mth_mean - sst).^2);

% index_error = find(RMSE <= 1.5);

figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% 지저분함
title_name = strcat(' 2007.',num2str(file_num), 'm - 4yr mth mean ROMS RMSE');
m_proj('mercator','lon',[98 284], 'lat',[-20 65]);
m_grid('box','fancy','tickdir','in'); 
m_pcolor(lon2,lat2,RMSE); 
shading interp;
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 1.5]);
saveas(gcf,strcat('monthly_mean_4yr_',num2str(file_num),'m_RMSE.png'),'png');
close; 


% index_error = find(RMSE <= 1.5);
% RMSE(index_error) = NaN;
% figure; hold on;
% %     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% 지저분함
% title_name = strcat(' 2007.',num2str(file_num), 'm - 1yr mth mean ROMS RMSE');
% m_proj('mercator','lon',[98 284], 'lat',[-20 65]);
% m_grid('box','fancy','tickdir','in'); 
% m_pcolor(lon2,lat2,RMSE); 
% shading interp;
% m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
% title(title_name,'fontsize',30);set(gca,'FontSize',18);
% xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
% colorbar;
% caxis([0 25]);
% saveas(gcf,strcat('monthly_mean_1yr_',num2str(file_num),'m_RMSE_part.png'),'png');
% close; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPACE MEAN RMSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all; 

temp_file = strcat('monthly_spinup_*.nc');
d_mat = dir('2007*')
d_filename = dir(temp_file); % make a list of nc in the current folder.
load location.mat

lon_ob = ncread('avhrr-only-v2.20070101.nc','lon');
lat_ob = ncread('avhrr-only-v2.20070101.nc','lat');

for i = 1:12
file_num = num2str(i);
temp = ncread(d_filename(i).name,'temp');
lon = ncread(d_filename(i).name,'lon_rho');
lat = ncread(d_filename(i).name,'lat_rho');
z = ncread(d_filename(i).name,'s_rho'); % depth
load(d_mat(i).name);

sst = temp(:,:,30); % diminishing the dimension
squeeze(sst);

index = find(sst >= 1000);
index_cold = find(sst <= -1);
reshape(sst,1988616,1);
sst(index) = NaN;
reshape(sst,1862,1068);

index1 = find(monthly_mean <= -9);
monthly_mean(index1) = NaN;

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

[lat1,lon1] = meshgrid(lat_ob,lon_ob);

[lon2,lat2,mth_mean] = griddata(lon1,lat1, monthly_mean, lon, lat);
% [lon2,lat2,V_1] = griddata(lon_v,lat_v, V, lon, lat);

index_nan = find(isnan(mth_mean) == 1);

% mth_mean(index_nan) = -999.9;

RMSE = sqrt((mth_mean - sst).^2);

RMSE_mean = reshape(RMSE,1988616,1);
index_nan1 = find(isnan(RMSE_mean) == 1);
RMSE_mean(index_nan1) = 0;
result(i) = mean(RMSE_mean)
end
save('result.mat','result');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RMSE            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all; 

temp_file = strcat('monthly_spinup_*.nc');
d_mat = dir('2007*')
d_filename = dir(temp_file); % make a list of nc in the current folder.
load location.mat

lon_ob = ncread('avhrr-only-v2.20070101.nc','lon');
lat_ob = ncread('avhrr-only-v2.20070101.nc','lat');

for i = 1:12
file_num = num2str(i);
clearvars monthly_mean sst climate_sst
load(['2001-2010_avhrr_sst_climate_',num2str(i),'mon.mat']);

sst = squeeze(t_mod_in(:,:,i)); % diminishing the dimension
sst(mask2)=NaN;

monthly_mean=climate_sst;
% index = find(sst >= 1000);
% index_cold = find(sst <= -1);
% reshape(sst,1988616,1);
% sst(index) = NaN;
% reshape(sst,1862,1068);

index1 = find(monthly_mean <= -9);
monthly_mean(index1) = NaN;

min_lon = min(lon_m);
max_lon = max(lon_m);
min_lat = min(lat_m);
max_lat = max(lat_m);
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

[lat1,lon1] = meshgrid(double(lat_ob),double(lon_ob));

[lon2,lat2,mth_mean] = griddata(lon1,lat1, monthly_mean, lon_m, lat_m);
% [lon2,lat2,V_1] = griddata(lon_v,lat_v, V, lon, lat);

index_nan = find(isnan(mth_mean) == 1);

% mth_mean(index_nan) = -999.9;
bias = mth_mean - sst;
minus_1 = (mth_mean - sst).^2;
minus_11 = reshape(minus_1,size(minus_1,1)*size(minus_1,2),1);
index_nan1 = find(isnan(minus_11) == 1);
minus_11(index_nan1) = 0;
RM = mean(minus_11);
RMSE = sqrt(RM);

result(i) = mean(RMSE)
end
save('result_fix.mat','result');

plot(result) % rmse = 0.8318, bias = 
