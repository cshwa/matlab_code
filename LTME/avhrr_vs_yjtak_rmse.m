clear all; clc; close all; 
cd D:\장기생태\Dynamic
lon_ob = ncread('amsr-avhrr-v2.20100101.nc','lon');
lat_ob = ncread('amsr-avhrr-v2.20100101.nc','lat');

cd D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary

lon=ncread('./2010/stdep1_avg_mon_2010_01.nc','lon');
lat=ncread('./2010/stdep1_avg_mon_2010_01.nc','lat');
dep_mod=ncread('./2010/stdep1_avg_mon_2010_01.nc','depth');
mask_pre=ncread('./2010/stdep1_avg_mon_2010_01.nc','temp');
[lat_m lon_m]=meshgrid(lat,lon);

clearvars lon lat


mask=squeeze(mask_pre(:,:,1));
mask(~isnan(mask)) = 1;

k = 0;
for i = 1:10 %year
     for j = 1:12 %month
        k= k+1;
        t_model(:,:,:,i,j) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'temp');
%         s_model(:,:,:,i,j) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'salt');
% %         do_model(:,:,:,k) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'oxygen');
%         no3_model(:,:,:,i,j) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'NO3');
%         nh4_model(:,:,:,i,j) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'NH4');
%         zeta(:,:,k) = ncread(['./',num2str(2000+i),'/stdep1_avg_mon_',num2str(2000+i),'_',num2str(j,'%02d'),'.nc'],'zeta');
     end
end

test=squeeze(t_model(:,:,1,1,1));
mask2=find(test > 1000);
test(mask2)=NaN;

clearvars z_mod t_mod_in
for i =1:size(t_model,5)
% monthly climate
        t_mod_in(:,:,i)=double(squeeze(mean(t_model(:,:,1,:,i),4)));
%         s_mod_in(:,i)=squeeze(mean(s_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
% %         do_mod_in(:,i)=squeeze(mth_do(near_point_2d(1,1),near_point_2d(1,2),:,i));
%         nh4_mod_in(:,i)=squeeze(mean(nh4_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
%         no3_mod_in(:,i)=squeeze(mean(no3_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),4));
 
%monthly std        
        std_t_mod(:,:,i)=double(squeeze(nanstd(t_model(:,:,1,:,i),0,4)));
%         std_s_mod(:,i)=squeeze(nanstd(s_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),0,4)); % 0 =weight. 4 = dim.
%         std_nh4_mod(:,i)=squeeze(nanstd(nh4_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),0,4));
%         std_no3_mod(:,i)=squeeze(nanstd(no3_model(near_point_2d(1,1),near_point_2d(1,2),:,:,i),0,4));
end

for i = 1:12
file_num = num2str(i);
clearvars monthly_mean sst climate_sst
load(['2001-2010_avhrr_sst_climate_',num2str(i),'mon.mat']);

sst = squeeze(t_mod_in(:,:,i)); % diminishing the dimension
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

RMSE = sqrt((mth_mean - sst).^2);

% index_error = find(RMSE <= 1.5);

figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% 지저분함
title_name = strcat(' 2001-2010.',num2str(file_num), 'm - mth mean ROMS RMSE');
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_grid('box','fancy','tickdir','in'); 
m_pcolor(lon2,lat2,RMSE); 
shading interp;
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 1.5]);
saveas(gcf,strcat('monthly_mean_yjtak_',num2str(file_num),'m_RMSE.png'),'png');
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
bias_f(i) = nanmean(nanmean(bias));
minus_1 = (mth_mean - sst).^2;
minus_11 = reshape(minus_1,size(minus_1,1)*size(minus_1,2),1);
index_nan1 = find(isnan(minus_11) == 1);
minus_11(index_nan1) = NaN;
RM = nanmean(minus_11);
RMSE = sqrt(RM);

result(i) = mean(RMSE)
end
save('result_fix.mat','result');

figure;
plot(result) 
title_name = strcat(' 2001-2010 - mth mean ROMS RMSE');
title(title_name,'fontsize',15);
xlabel('time (month)','fontsize', 15);ylabel('RMSE','fontsize',15);
set(gca,'fontsize',15);
xlim([1 12]); grid on;

figure;
plot(bias_f) 
title_name = strcat(' 2001-2010 - mth mean ROMS Bias');
title(title_name,'fontsize',15);
xlabel('time (month)','fontsize', 15);ylabel('Bias','fontsize',15);
set(gca,'fontsize',15);
xlim([1 12]); grid on;


% rmse = 1.0971, bias = -0.5665

figure; hold on;
%     m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20); %% 지저분함
title_name = strcat(' 2001-2010.',num2str(file_num), 'm - mth mean ROMS RMSE');
m_proj('mercator','lon',[min_lon max_lon], 'lat',[min_lat max_lat]);
m_grid('box','fancy','tickdir','in'); 
m_pcolor(lon2,lat2,RMSE); 
shading interp;
m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([0 1.5]);
saveas(gcf,strcat('monthly_mean_yjtak_',num2str(file_num),'m_RMSE.png'),'png');
close; 
