close all; clear all; clc;

tair_list=dir(['auto_ERA5_*_Tair.nc']);

od_i = 0
for i = 1982:2020
    od_i = od_i+1;
    clearvars temp_t
    temp_t=mean(ncread(['auto_ERA5_',num2str(i),'_Tair.nc'],'Tair'),'all');
    temp_air(od_i) =  temp_t;
end

lon=ncread('/data2/cshwa/make_frc_ecmwf_m/SFC_ERA5/grid_gy_v11_s.nc','lon_rho');
lat=ncread('/data2/cshwa/make_frc_ecmwf_m/SFC_ERA5/grid_gy_v11_s.nc','lat_rho');
area_mask=NaN(size(lon,1),size(lon,2));
area_mask(find(lon >= 127.593 & lon <= 127.82 & 34.67<= lat & 34.945 >= lat))=1;

od_i = 0
for i = 1982:2020
    od_i = od_i+1;
    clearvars temp_t
    temp_t=ncread(['auto_ERA5_',num2str(i),'_Tair.nc'],'Tair');
    for j = 1:size(temp_t,3)
    temp_t2(:,:,j) = temp_t(:,:,j) .* area_mask;
    end
    temp_t3 =  nanmean(temp_t2, 'all');
    temp_air(od_i) =  temp_t3;
end

plot(temp_air)
save('Tair_ERA5_yearly_gy.mat')



close all; clear; clc;

load Tair_ERA5_aug_gy.mat

figure;
plot(temp_aug,'linew',2)
xlim([1 10]); xticks(1:10);
xticklabels({'1982-1986','1987-1992','1993-2004','2007-2015','2016-2020'});
xtickangle(45); grid on;
ylabel('기온(^oC)')
title('8월 평균 기온');
xlim([1 5])
set(gca,'fontsize',15,'fontweight','bold');

save('Tair_ERA5_yearly_gy.mat')