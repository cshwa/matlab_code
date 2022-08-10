close all; clear all; clc;

cd /home/cshwa/cmip5/CMIP5_bio
% GY -RCM v3 load data, wind climate fixed(deg, speed), river temp fixed
% but boundary delta is not added (made mistake) 
list_gy_v3_pre1=dir('ncra_*to*_2yr.nc');
list_gy_v3_pre2=dir('ncra_*_cir_2yr.nc');


list_gy_v3 = [list_gy_v3_pre1; list_gy_v3_pre2;];

clearvars *_gy_all_v3

t_gy_all_v3 = NaN(20,10);
s_gy_all_v3 = NaN(20,10);
chl_gy_all_v3 = NaN(20,10);
no3_gy_all_v3 = NaN(20,10);
nh4_gy_all_v3 = NaN(20,10);
o2_gy_all_v3 = NaN(20,10);

lon=ncread('grid_gy_v11_s.nc','lon_rho');
lat=ncread('grid_gy_v11_s.nc','lat_rho');
% latLim = [34.67 34.945];
% lonLim = [127.593 127.82];
% area_mask=NaN(size(lon,1),size(lon,2));
% area_mask(find(lon >= 127.593 & lon <= 127.82 & 34.67<= lat & 34.945 >= lat))=1;
area_mask=ones(size(lon,1),size(lon,2));

% v3_i=[1:2,4:10];
v3_i = [1:10];
for i = 1:length(list_gy_v3)
    for j=1:20 %s_rho
    clearvars t_pre s_pre chl_pre nh4_pre no3_pre mask
t_pre=ncread(list_gy_v3(i).name,'temp');
s_pre=ncread(list_gy_v3(i).name,'salt');
chl_pre=ncread(list_gy_v3(i).name,'chlorophyll');
nh4_pre=ncread(list_gy_v3(i).name,'NH4');
no3_pre=ncread(list_gy_v3(i).name,'NO3');
o2_pre=ncread(list_gy_v3(i).name,'oxygen');
mask=ncread(list_gy_v3(i).name,'mask_rho');

% t_gy_all_v3(j,i) = squeeze(nanmean(nanmean(squeeze(t_pre(:,:,j)).*(mask./mask),1),2));
% s_gy_all_v3(j,i) = squeeze(nanmean(nanmean(squeeze(s_pre(:,:,j)).*(mask./mask),1),2));
% chl_gy_all_v3(j,i) = squeeze(nanmean(nanmean(squeeze(chl_pre(:,:,j)).*(mask./mask),1),2));
% no3_gy_all_v3(j,i) = squeeze(nanmean(nanmean(squeeze(no3_pre(:,:,j)).*(mask./mask),1),2));
% nh4_gy_all_v3(j,i) = squeeze(nanmean(nanmean(squeeze(nh4_pre(:,:,j)).*(mask./mask),1),2));
% o2_gy_all_v3(j,i) = squeeze(nanmean(nanmean(squeeze(o2_pre(:,:,j)).*(mask./mask),1),2));
t_gy_all_v3(j,v3_i(i)) = squeeze(nanmean(nanmean(squeeze(t_pre(:,:,j)).*(mask./mask) .* area_mask,1),2));
s_gy_all_v3(j,v3_i(i)) = squeeze(nanmean(nanmean(squeeze(s_pre(:,:,j)).*(mask./mask).* area_mask,1),2));
chl_gy_all_v3(j,v3_i(i)) = squeeze(nanmean(nanmean(squeeze(chl_pre(:,:,j)).*(mask./mask).* area_mask,1),2));
no3_gy_all_v3(j,v3_i(i)) = squeeze(nanmean(nanmean(squeeze(no3_pre(:,:,j)).*(mask./mask).* area_mask,1),2));
nh4_gy_all_v3(j,v3_i(i)) = squeeze(nanmean(nanmean(squeeze(nh4_pre(:,:,j)).*(mask./mask).* area_mask,1),2));
o2_gy_all_v3(j,v3_i(i)) = squeeze(nanmean(nanmean(squeeze(o2_pre(:,:,j)).*(mask./mask).* area_mask,1),2));
    end
end


fig_loc = ['/data1/cshwa/cmip5/'];


figure; hold on;
plot(t_gy_all_v3(end,:)','r','linew',2);
plot(t_gy_all_v3(1,:)','r--','linew',2);
xlim([1 10]); xticks(1:10);
xticklabels({'1982-1986','1987-1992','1993-2004','2007-2015','2016-2020'});
xtickangle(45); grid on;
ylabel('temp(^oC)')
title('GY model  plot temp');
xlim([1 5])
set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(s_gy_all_v3(end,:)','r','linew',2);
plot(s_gy_all_v3(1,:)','r--','linew',2);
xlim([1 10]); xticks(1:10);
xticklabels({'1982-1986','1987-1992','1993-2004','2007-2015','2016-2020'});

xtickangle(45); grid on;
ylabel('g/kg')
title('GY model  plot salt');
xlim([1 5])
set(gca,'fontsize',15,'fontweight','bold');


figure; hold on;
plot(no3_gy_all_v3(end,:)','r');
plot(no3_gy_all_v3(1,:)','r--');

xlim([1 94]); xticks(1:5:94);
xticklabels(2007:5:2100)
xtickangle(45); grid on;

xlim([1 10]); xticks(1:10);
xticklabels({'1982-1986','1987-1992','1993-2004','2007-2015','2016-2020'});

xtickangle(45); grid on;
ylabel('NO3 (mmol / m^3)')
title('CMIP5 IPSL-CM5A-MR raw. RCP plot no3 0m & 35m');
% saveas(gca,[fig_loc,'IPSL-CM5A-MR_interp_RCP_yrplot_NO3_0m_&_30m.png'],'png')
xlim([1 5])


figure; hold on;
plot(nh4_gy_all_v3(end,:)','r');
plot(nh4_gy_all_v3(1,:)','r--');
xlim([1 94]); xticks(1:5:94);
xticklabels(2007:5:2100)
xtickangle(45); grid on;

xlim([1 10]); xticks(1:10);
xticklabels({'1982-1986','1987-1992','1993-2004','2007-2015','2016-2020'});

xtickangle(45); grid on;
ylabel('NH4 (mmol / m^3)')
legend('GCM','GY','GY-v3');
title('CMIP5 IPSL-CM5A-MR raw. RCP plot NH4 0m & 35m');
% saveas(gca,[fig_loc,'IPSL-CM5A-MR_interp_RCP_yrplot_NO3_0m_&_30m.png'],'png')
xlim([1 5])


figure; hold on;
plot((no3_gy_all_v3(end,:)' + nh4_gy_all_v3(end,:)'),'r','linew',2);
plot((no3_gy_all_v3(1,:)' + nh4_gy_all_v3(1,:)'),'r--','linew',2);

xlim([1 94]); xticks(1:5:94);
xticklabels(2007:5:2100)
xtickangle(45); grid on;
xlim([1 10]); xticks(1:10);
xticklabels({'1982-1986','1987-1992','1993-2004','2007-2015','2016-2020'});
xtickangle(45); grid on;
ylabel('DIN (mmol / m^3)')
title('GY model  plot DIN');
xlim([1 5])
set(gca,'fontsize',15,'fontweight','bold');

figure; hold on;
plot(o2_gy_all_v3(end,:)','r','linew',2);
plot(o2_gy_all_v3(1,:)','r--','linew',2);
xlim([1 94]); xticks(1:5:94);
xticklabels(2007:5:2100)
xtickangle(45); grid on;
 yticks(250:5:270); ylim([250 270])
xlim([1 10]); xticks(1:10);
xticklabels({'1982-1986','1987-1992','1993-2004','2007-2015','2016-2020'});
xtickangle(45); grid on;
ylabel('o2 (mmol / m^3)')
title('GY model  plot oxygen');
xlim([1 5])
set(gca,'fontsize',15,'fontweight','bold');

% figure; hold on;
% plot(mod3_chl_rcp85_all_sp(1,:)'.*10^6,'r');
% 
% plot(mod3_chl_rcp85_all_sp(4,:)'.*10^6,'r--');
% 
% xlim([1 94]); xticks(1:5:94);
% xticklabels(2007:5:2100)
% xtickangle(45); grid on;
% ylabel('chl (ug/L)')
% legend('RCP26','RCP45','RCP85');
% title('CMIP5 IPSL-CM5A-MR interp. RCP plot chl 0m & 30m');
% saveas(gca,[fig_loc,'IPSL-CM5A-MR_interp_RCP_yrplot_chl_0m_&_30m.png'],'png')

figure; hold on;
plot(chl_gy_all_v3(end,:)','r','linew',2);
plot(chl_gy_all_v3(1,:)','r--','linew',2);

xlim([1 94]); xticks(1:5:94);
xticklabels(2007:5:2100)
xtickangle(45); grid on;

xlim([1 10]); xticks(1:10);
xticklabels({'1982-1986','1987-1992','1993-2004','2007-2015','2016-2020'});

xtickangle(45); grid on;
ylabel('chl (ug/L)')

title('GY model  plot Chl.a');
xlim([1 5])
set(gca,'fontsize',15,'fontweight','bold');


