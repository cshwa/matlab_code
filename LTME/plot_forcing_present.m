%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;

file_path = 'D:\장기생태\Dynamic\Input\grid_gy_v11_s.nc';
lon = ncread(file_path,'lon_rho');
lat = ncread(file_path,'lat_rho');
lon_u = ncread(file_path,'lon_u');
lat_u = ncread(file_path,'lat_u');
lon_v = ncread(file_path,'lon_v');
lat_v = ncread(file_path,'lat_v');
mask_r = ncread(file_path,'mask_rho');
mask_u = ncread(file_path,'mask_u');
mask_v = ncread(file_path,'mask_v');


k = [2011:2016];
for i = 1:6
%      for j =  1:365
        h=k(i);
        if leapyear(h) == 1
            clearvars temp_*
        temp_temp = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_',num2str(h,'%04d'),'Tair_366.nc'],'Tair'));
        temp_temp(:,:,end) = []; a_temp(:,:,i,:) = temp_temp;
        temp_p = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_',num2str(h,'%04d'),'Pair_366.nc'],'Pair'));
        temp_p(:,:,end) = []; a_p(:,:,i,:) = temp_p;
        temp_u = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_',num2str(h,'%04d'),'Uwind_366.nc'],'Uwind'));
        temp_u(:,:,end) = []; a_u(:,:,i,:) = temp_u;
        temp_v = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_',num2str(h,'%04d'),'Vwind_366.nc'],'Vwind'));
        temp_v(:,:,end) = []; a_v(:,:,i,:) = temp_v;
%         temp_trans = double(ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/river_',num2str(h,'%04d'),'_realts_biofennel_Gwangyang.nc'],'river_transport'));      
%         temp_trans(:,end) = []; r_trans(:,:,i) = temp_trans;
        elseif leapyear(h) == 0 && h~=2013
            a_temp(:,:,i,:) = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_',num2str(h,'%04d'),'Tair_365.nc'],'Tair'));
            a_p(:,:,i,:) = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_',num2str(h,'%04d'),'Pair_365.nc'],'Pair'));
            a_u(:,:,i,:) = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_',num2str(h,'%04d'),'Uwind_365.nc'],'Uwind'));
            a_v(:,:,i,:) = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_',num2str(h,'%04d'),'Vwind_365.nc'],'Vwind'));
    %         r_trans(:,:,i) = double(ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/river_',num2str(h,'%04d'),'_realts_biofennel_Gwangyang.nc'],'river_transport'));
        elseif h == 2013
            a_temp(:,:,i,:) = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_sumjin_',num2str(h,'%04d'),'Tair_365.nc'],'Tair'));
            a_p(:,:,i,:) = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_sumjin_',num2str(h,'%04d'),'Pair_365.nc'],'Pair'));
            a_u(:,:,i,:) = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_sumjin_',num2str(h,'%04d'),'Uwind_365.nc'],'Uwind'));
            a_v(:,:,i,:) = double(ncread(['D:\장기생태\Dynamic\05_airforcing\result\',num2str(h),'\frc_ecmwf_sumjin_',num2str(h,'%04d'),'Vwind_365.nc'],'Vwind'));
        end
end

cd G:\장기생태\Dynamic\06_river
load('climate_nam_2011to2016.mat'); %'nam_transp_365','nam_clim_trans_mth');
cd G:\장기생태\Dynamic\06_river\data\sj_00_16_year
load('climate_songjung_2011to2016.mat'); %'clim_song_11_16_d','clim_trans_mth');

r_trans_mth = [clim_trans_mth; nam_clim_trans_mth;];

clearvars *_365
for i = 1:365
    a_temp_365(:,:,i)=nanmean(squeeze(a_temp(:,:,:,i)),3);
    a_p_365(:,:,i)=nanmean(squeeze(a_p(:,:,:,i)),3);
    a_u_365(:,:,i)=nanmean(squeeze(a_u(:,:,:,i)),3);
    a_v_365(:,:,i)=nanmean(squeeze(a_v(:,:,:,i)),3);
%     r_trans_365(:,i)=nanmean(squeeze(r_trans(:,i,:)),2);
end

save('2011_2016_mean_forcing_GY.mat','-v7.3');

save('2011_2016_forcing_GY.mat','*_365','-v7.3');


% make 1980~present
k=0
for i = 1981:1981
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end


clearvars em em_seanson k_f kk
k_f=[1:12];
for i = 1:length(k_f)
    em{i} = sum(eom_d(1:k_f(i)-1))+1:sum(eom_d(1:k_f(i)));
end

clearvars *_mth
for i = 1:12
    a_temp_mth(:,:,i)=nanmean(a_temp_365(:,:,em{i}),3);
    a_p_mth(:,:,i)=nanmean(a_p_365(:,:,em{i}),3);
    a_u_mth(:,:,i)=nanmean(a_u_365(:,:,em{i}),3);
    a_v_mth(:,:,i)=nanmean(a_v_365(:,:,em{i}),3);
%     r_trans_mth(:,i)=nanmean(squeeze(r_trans_365(:,em{i})),2);
end

min_lon = min(min(lon));
max_lon =  max(max(lon));
min_lat = min(min(lat));
max_lat = max(max(lat));


clearvars *_ref
u_ref = zeros(size(lon,1),size(lon,2));
v_ref = zeros(size(lon,1),size(lon,2));
u_ref(11,161) = 10*0.01;


interval = 10;
for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1980~1986.',num2str(i,'%02d'), ' mon. air temp.');
pcolor(lon,lat,squeeze(a_temp_mth(:,:,i))); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(a_u_mth(1:interval:end,1:interval:end,i)).*0.01, ...
    squeeze(a_v_mth(1:interval:end,1:interval:end,i)).*0.01, ...
    'AutoScale','off','LineWidth',1,'Color','k');
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([-1 26]);
colormap(jet)
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
if i == 5 || i == 10
    quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','r');
    text(lon(11,151),lat(11,151),'10m/s','fontsize',16,'color','r');
elseif i == 1 || i == 2 || i == 3 || i == 4 || i == 6 || i == 7 || i == 8 || i == 9 || i == 11 || i == 12
    quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color',[1 1 1]);
    text(lon(11,151),lat(11,151),'10m/s','fontsize',16,'color',[1 1 1]);
end
saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(air_temp).png'),'png');
end


interval = 10;
for i = 1:12
figure; hold on;
title_name = strcat('GY-',' 1980~1986.',num2str(i,'%02d'), ' mon. air pressure.');
pcolor(lon,lat,squeeze(a_p_mth(:,:,i))); shading flat;
quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end), ...
    squeeze(a_u_mth(1:interval:end,1:interval:end,i)).*0.01, ...
    squeeze(a_v_mth(1:interval:end,1:interval:end,i)).*0.01, ...
    'AutoScale','off','LineWidth',1,'Color','k');
title(title_name,'fontsize',30);set(gca,'FontSize',18);
xlabel('Lon','fontsize', 20);ylabel('Lat','fontsize',20);
colorbar;
caxis([1000 1030]);
colormap(jet)
xlim([min_lon max_lon]);
ylim([min_lat max_lat]);
if i == 5 || i == 10
    quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','r');
    text(lon(11,151),lat(11,151),'10m/s','fontsize',16,'color','r');
elseif i == 1 || i == 2 || i == 3 || i == 4 || i == 6 || i == 7 || i == 8 || i == 9 || i == 11 || i == 12
    quiver(lon,lat,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color',[1 1 1]);
    text(lon(11,151),lat(11,151),'10m/s','fontsize',16,'color',[1 1 1]);
end
saveas(gcf,strcat('monthly_mean_',num2str(i,'%02d'),'m_(air_p).png'),'png');
end

% figure; hold on;
% for i = 1:7
%     plot(squeeze(r_trans(1,:,i)))
% end
% xlim([1 365])

figure; hold on;
for i = 1:2
    plot(squeeze(r_trans_mth(i,:)))
end
xlim([1 12])




save('2011_2016_mean_mth_all_GY.mat','-v7.3');
save('2011_2016_mean_mth_GY.mat','temp_mth','salt_mth','u_mth','v_mth','-v7.3');