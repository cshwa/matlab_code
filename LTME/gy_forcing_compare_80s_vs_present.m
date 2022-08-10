%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% air temp 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% present
close all; clear all; clc;
cd D:\장기생태\Dynamic\03_initial_TS
lon=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m_v3.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m_v3.nc','h');

% figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
% grid on; caxis([0 40]);

yy = 2012:2017;
isleap = leapyear(yy);
tair = NaN(size(lon,1),size(lon,2),366,6);
for i = 1:6
% 2012, 2016 = 366
cd(['G:\장기생태\NPZD\' num2str(yy(i)) '\forcing'])
if isleap(i) == 1
    tair(:,:,:,i) = ncread(['frc_ecmwf_' num2str(yy(i)) 'Tair_366.nc'],'Tair');
else
    tair(:,:,1:365,i) = ncread(['frc_ecmwf_' num2str(yy(i)) 'Tair_365.nc'],'Tair');
end
end

% pcolor(lon,lat, squeeze(tair(:,:,1,1)).*(mask./mask)); hh=colorbar; shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'temp (^oC)');
% grid on;

% make calandr
k=0
for i = 2012:2012
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

kk=[1:12];
for i = 1:length(kk)
    if i == 1
        em{i} = 1:sum(eom_d(1:kk(i)));
    else
        em{i} = sum(eom_d(1:kk(i)-1))+1:sum(eom_d(1:kk(i)));
    end
end

%make monthly mean
for i = 1:12 %mth
        t_a_clim(:,:,i)= squeeze(nanmean(nanmean(tair(:,:,em{i},:),4),3));
    for j = 1:6 %year
        t_a_mth(:,:,i,j)= nanmean(tair(:,:,em{i},j),3);
    end
end

for i = 1:12
    figure;
pcolor(lon,lat, squeeze(t_a_clim(:,:,i)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'temp (^oC)');
grid on;
caxis([-1 25]);
end

% detect time axes error on frc_ecmwf_2017Tair_365_ori.nc 
plot(squeeze(nanmean(nanmean(tair(:,:,:,5),2),1)));hold on;
plot(squeeze(nanmean(nanmean(tair(:,:,:,4),2),1)),'r')
plot(squeeze(nanmean(nanmean(tair(:,:,:,6),2),1)),'k')
cd D:\장기생태\Dynamic
save('present_atemp.mat','-v7.3');

%% past 1980s
close all; clear all; clc;
cd D:\장기생태\Dynamic\03_initial_TS
lon=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m_v3.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m_v3.nc','h');

% figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
% grid on; caxis([0 40]);

cd D:\장기생태\Dynamic\Input\cshwa_1980s_new
yy = 1980:1985;
isleap = leapyear(yy);
tair = NaN(size(lon,1),size(lon,2),366,6);
for i = 1:6
% 1980, 1984 = 366
if isleap(i) == 1
    tair(:,:,:,i) = ncread(['Gwangyang_' num2str(yy(i)) '_Tair.nc'],'Tair');
else
    tair(:,:,1:365,i) = ncread(['Gwangyang_' num2str(yy(i)) '_Tair.nc'],'Tair');
end
end

% make calandr
k=0
for i = 1980:1980
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

kk=[1:12];
for i = 1:length(kk)
    if i == 1
        em{i} = 1:sum(eom_d(1:kk(i)));
    else
        em{i} = sum(eom_d(1:kk(i)-1))+1:sum(eom_d(1:kk(i)));
    end
end

%make monthly mean
for i = 1:12 %mth
        t_a_clim(:,:,i)= squeeze(nanmean(nanmean(tair(:,:,em{i},:),4),3));
    for j = 1:6 %year
        t_a_mth(:,:,i,j)= nanmean(tair(:,:,em{i},j),3);
    end
end

for i = 1:12
    figure;
pcolor(lon,lat, squeeze(t_a_clim(:,:,i)).*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'temp (^oC)');
grid on;
caxis([-1 25]);
end

cd D:\장기생태\Dynamic
save('past_atemp.mat','-v7.3');

%% compare
close all; clear all; clc;
cd D:\장기생태\Dynamic\03_initial_TS
lon=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m_v3.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m_v3.nc','h');

cd D:\장기생태\Dynamic

load past_atemp.mat
t_a_clim_80 = t_a_clim;
t_a_mth_80 = t_a_mth;

clearvars -except t_a_clim_80 t_a_mth_80 lon lat mask
load present_atemp.mat

% make spatial mean
t_a_sm_80 = squeeze(nanmean(nanmean(t_a_clim_80,2),1));
t_a_sm = squeeze(nanmean(nanmean(t_a_clim,2),1));

t_a_sm_mth_80 = squeeze(nanmean(nanmean(t_a_mth_80,2),1));
t_a_sm_mth = squeeze(nanmean(nanmean(t_a_mth,2),1));

figure; hold on;
plot(t_a_sm_80,'b','linew',2); plot(t_a_sm,'r','linew',2); 
plot(1:12, repmat(mean(t_a_sm_80),12,1),'--','color','b','linew',2);
plot(1:12, repmat(mean(t_a_sm),12,1),'--','color','r','linew',2);
legend('temp_1_9_8_0_~_1_9_8_5', 'temp_2_0_1_2_~_2_0_1_7');
ylabel(gca,'Air temperature (^oC)','fontsize',15,'fontweight','bold'); xlabel(gca,'month','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 12]);

t_a_sm_mth_80_f = [t_a_sm_mth_80(:,1); t_a_sm_mth_80(:,2); t_a_sm_mth_80(:,3); t_a_sm_mth_80(:,4); t_a_sm_mth_80(:,5); t_a_sm_mth_80(:,6)];
t_a_sm_mth_f = [t_a_sm_mth(:,1); t_a_sm_mth(:,2); t_a_sm_mth(:,3); t_a_sm_mth(:,4); t_a_sm_mth(:,5); t_a_sm_mth(:,6)];

figure; hold on;
plot(t_a_sm_mth_80_f,'b','linew',2); plot(t_a_sm_mth_f,'r','linew',2); 
plot(1:72, repmat(mean(t_a_sm_mth_80_f),72,1),'--','color','b','linew',2);
plot(1:72, repmat(mean(t_a_sm_mth_f),72,1),'--','color','r','linew',2);
le=legend('temp_1_9_8_0_~_1_9_8_5', 'temp_2_0_1_2_~_2_0_1_7');
set(le,'fontsize',15,'fontweight','bold')
ylabel(gca,'Air temperature (^oC)','fontsize',15,'fontweight','bold'); xlabel(gca,'time(year)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 12*6]); ylim([-5 35]);
set(gca,'xtick',[1:12:72]);
set(gca,'xticklabel',1:1:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% River
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% present
close all; clear all; clc;
cd D:\장기생태\Dynamic\03_initial_TS
lon=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m_v3.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m_v3.nc','h');

% figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
% grid on; caxis([0 40]);

yy = 2012:2017;
isleap = leapyear(yy);
riv_trans = NaN(2,366,6);
riv_temp = NaN(2,20,366,6);
% riv_salt ; present = 1
isleap(5)=0
for i = 1:6
% 2012, 2016 = 366  (it's 2016 365 time on river file
cd(['G:\장기생태\NPZD\' num2str(yy(i)) '\forcing'])
if isleap(i) == 1
    riv_trans(:,:,i) = ncread(['ocean_rivers' num2str(yy(i)) '_ts.nc'],'river_transport');
    riv_temp(:,:,:,i) = ncread(['ocean_rivers' num2str(yy(i)) '_ts.nc'],'river_temp');
else
    riv_trans(:,1:365,i) = ncread(['ocean_rivers' num2str(yy(i)) '_ts.nc'],'river_transport');
    riv_temp(:,:,1:365,i) = ncread(['ocean_rivers' num2str(yy(i)) '_ts.nc'],'river_temp');
end
end

riv_nh4 = NaN(2,20,366,3);
riv_no3 = NaN(2,20,366,3);
for i = 4:6
% 2012, 2016 = 366  (it's 2016 365 time on river file
cd(['G:\장기생태\NPZD\' num2str(yy(i)) '\forcing'])
if isleap(i) == 1
    riv_nh4(:,:,:,i-3) = ncread(['ocean_rivers' num2str(yy(i)) '_ts.nc'],'river_NH4');
    riv_no3(:,:,:,i-3) = ncread(['ocean_rivers' num2str(yy(i)) '_ts.nc'],'river_NO3');
else
    riv_nh4(:,:,1:365,i-3) = ncread(['ocean_rivers' num2str(yy(i)) '_ts.nc'],'river_NH4');
    riv_no3(:,:,1:365,i-3) = ncread(['ocean_rivers' num2str(yy(i)) '_ts.nc'],'river_NO3');
end
end

% make calandr
k=0
for i = 2012:2012
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

kk=[1:12];
for i = 1:length(kk)
    if i == 1
        em{i} = 1:sum(eom_d(1:kk(i)));
    else
        em{i} = sum(eom_d(1:kk(i)-1))+1:sum(eom_d(1:kk(i)));
    end
end

%make daily mean
 r_trans_clim(:,:)= squeeze(nanmean(riv_trans(:,:,:),3));
 r_temp_clim(:,:) = squeeze(nanmean(riv_temp(:,20,:,:),4));
 r_no3_clim(:,:) = squeeze(nanmean(riv_no3(:,20,:,:),4));
 r_nh4_clim(:,:) = squeeze(nanmean(riv_nh4(:,20,:,:),4));  
  
% %  check that it's really same with variation of the depth and time
% for i = 1:3
% plot(squeeze(riv_NO3(2,1,:,i))); hold on;
% end 
% for i = 1:3
% plot(squeeze(riv_NH4(2,1,:,i))); hold on;
% end
  
% for i = 1:12
figure; hold on;
plot(r_trans_clim(1,:),'b','linew',2); plot(r_trans_clim(2,:),'r','linew',2); 
legend('Sumjin_2_0_1_2_~_2_0_1_7', 'Nam_2_0_1_2_~_2_0_1_7');
ylabel(gca,'discharge (m^3/s)','fontsize',15,'fontweight','bold'); xlabel(gca,'days','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 366]); 

cd D:\장기생태\Dynamic
save('present_river.mat','-v7.3');

%% past 1980s
close all; clear all; clc;
cd D:\장기생태\Dynamic\03_initial_TS
lon=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m_v3.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m_v3.nc','h');

cd D:\장기생태\Dynamic\Input\cshwa_1980s_new
yy = 1980:1985;
isleap = leapyear(yy);
riv_trans = NaN(2,366,6);
riv_temp = NaN(2,20,366,6);
% riv_salt ; present = 1
for i = 1:6
%  it's all 365 time on river file
if isleap(i) == 1
    riv_trans(:,:,i) = ncread(['river_' num2str(yy(i)) '_realts_biofennel_Gwangyang.nc'],'river_transport');
    riv_temp(:,:,:,i) = ncread(['river_' num2str(yy(i)) '_realts_biofennel_Gwangyang.nc'],'river_temp');
else
    riv_trans(:,1:365,i) = ncread(['river_' num2str(yy(i)) '_realts_biofennel_Gwangyang.nc'],'river_transport');
    riv_temp(:,:,1:365,i) = ncread(['river_' num2str(yy(i)) '_realts_biofennel_Gwangyang.nc'],'river_temp');
end
end

riv_nh4 = NaN(2,20,366,6);
riv_no3 = NaN(2,20,366,6);
for i = 1:6
%  it's all 365 time on river file
if isleap(i) == 1
    riv_nh4(:,:,:,i) = ncread(['river_' num2str(yy(i)) '_realts_biofennel_Gwangyang.nc'],'river_NH4');
    riv_no3(:,:,:,i) = ncread(['river_' num2str(yy(i)) '_realts_biofennel_Gwangyang.nc'],'river_NO3');
else
    riv_nh4(:,:,1:365,i) = ncread(['river_' num2str(yy(i)) '_realts_biofennel_Gwangyang.nc'],'river_NH4');
    riv_no3(:,:,1:365,i) = ncread(['river_' num2str(yy(i)) '_realts_biofennel_Gwangyang.nc'],'river_NO3');
end
end

% make calandr
k=0
for i = 2012:2012
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

kk=[1:12];
for i = 1:length(kk)
    if i == 1
        em{i} = 1:sum(eom_d(1:kk(i)));
    else
        em{i} = sum(eom_d(1:kk(i)-1))+1:sum(eom_d(1:kk(i)));
    end
end

%make daily mean
 r_trans_clim(:,:)= squeeze(nanmean(riv_trans(:,:,:),3));
 r_temp_clim(:,:) = squeeze(nanmean(riv_temp(:,20,:,:),4)); % same with variation of the depth
 r_nh4_clim(:,:) = squeeze(riv_nh4(:,20,:,1)); % same with variation of the depth and time
 r_no3_clim(:,:) = squeeze(riv_no3(:,20,:,1)); % same with variation of the depth and time
 
% %  check that it's really same with variation of the depth and time
%  for i = 1:6
%      plot(squeeze(riv_nh4(1,1,:,i))); hold on;
%  end
%  
%  for i = 1:6
%      plot(squeeze(riv_no3(1,1,:,i))./14); hold on;
%  end
 
              
% for i = 1:12
figure; hold on;
plot(r_trans_clim(1,:),'b','linew',2); plot(r_trans_clim(2,:),'r','linew',2); 
legend('Sumjin_1_9_8_0_~_1_9_8_5', 'Nam_1_9_8_0_~_1_9_8_5');
ylabel(gca,'discharge (m^3/s)','fontsize',15,'fontweight','bold'); xlabel(gca,'days','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 366]); 

cd D:\장기생태\Dynamic
save('past_river.mat','-v7.3');


%% compare
close all; clear all; clc;

cd D:\장기생태\Dynamic
load past_river.mat

r_nh4_clim_80 = r_nh4_clim;
r_no3_clim_80 = r_no3_clim;
r_temp_clim_80 = r_temp_clim;
r_trans_clim_80 = r_trans_clim;

riv_nh4_80 = riv_nh4;
riv_no3_80 = riv_no3;
riv_temp_80 = riv_temp;
riv_trans_80 = riv_trans;

clearvars -except *_80
load present_river.mat

% transp
figure; hold on;
plot(r_trans_clim_80(1,:),'b','linew',2); plot(r_trans_clim(1,:),'r','linew',2); 
plot(1:366, repmat(mean(r_trans_clim_80(1,:)),366,1),'--','color','b','linew',2);
plot(1:366, repmat(mean(r_trans_clim(1,:)),366,1),'--','color','r','linew',2);
legend('sumjin_1_9_8_0_~_1_9_8_5', 'sumjin_2_0_1_2_~_2_0_1_7');
ylabel(gca,'River discharge (m^3/s)','fontsize',15,'fontweight','bold'); xlabel(gca,'days','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 366]);ylim([0 1000])

figure; hold on;
plot(r_trans_clim_80(2,:),'b','linew',2); plot(r_trans_clim(2,:),'r','linew',2); 
plot(1:366, repmat(mean(r_trans_clim_80(2,:)),366,1),'--','color','b','linew',2);
plot(1:366, repmat(mean(r_trans_clim(2,:)),366,1),'--','color','r','linew',2);
legend('Nam_1_9_8_0_~_1_9_8_5', 'Nam_2_0_1_2_~_2_0_1_7');
ylabel(gca,'River discharge (m^3/s)','fontsize',15,'fontweight','bold'); xlabel(gca,'days','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 366]); ylim([0 1000]);

clearvars *_f
riv_trans_f(:,1) = [squeeze(riv_trans(1,:,1))'; squeeze(riv_trans(1,:,2))'; squeeze(riv_trans(1,:,3))'; squeeze(riv_trans(1,:,4))'; squeeze(riv_trans(1,:,5))'; squeeze(riv_trans(1,:,6))'];
riv_trans_f(:,2) = [squeeze(riv_trans(2,:,1))'; squeeze(riv_trans(2,:,2))'; squeeze(riv_trans(2,:,3))'; squeeze(riv_trans(2,:,4))'; squeeze(riv_trans(2,:,5))'; squeeze(riv_trans(2,:,6))'];
riv_trans_80_f(:,1) = [squeeze(riv_trans_80(1,:,1))'; squeeze(riv_trans_80(1,:,2))'; squeeze(riv_trans_80(1,:,3))'; squeeze(riv_trans_80(1,:,4))'; squeeze(riv_trans_80(1,:,5))'; squeeze(riv_trans_80(1,:,6))'];
riv_trans_80_f(:,2) = [squeeze(riv_trans_80(2,:,1))'; squeeze(riv_trans_80(2,:,2))'; squeeze(riv_trans_80(2,:,3))'; squeeze(riv_trans_80(2,:,4))'; squeeze(riv_trans_80(2,:,5))'; squeeze(riv_trans_80(2,:,6))'];

figure; hold on;
plot(riv_trans_80_f(:,1),'b','linew',2); plot(riv_trans_f(:,1),'r','linew',2); 
% plot(1:length(riv_trans_80_f(:,1)), repmat(nanmean(riv_trans_80_f(:,1)),length(riv_trans_f(:,1)),1),'--','color','b','linew',2);
% plot(1:length(riv_trans_f(:,1)), repmat(nanmean(riv_trans_f(:,1)),length(riv_trans_f(:,1)),1),'--','color','r','linew',2);
legend('sumjin_1_9_8_0_~_1_9_8_5', 'sumjin_2_0_1_2_~_2_0_1_7');
ylabel(gca,'River discharge (m^3/s)','fontsize',15,'fontweight','bold'); xlabel(gca,'time(year)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 length(riv_trans_f(:,1))]); ylim([0 3650]);
set(gca,'xtick',[1:366:length(riv_trans_f(:,1))]);
set(gca,'xticklabel',1:1:6);

figure; hold on;
plot(riv_trans_80_f(:,2),'b','linew',2); plot(riv_trans_f(:,2),'r','linew',2); 
% plot(1:length(riv_trans_80_f(:,1)), repmat(nanmean(riv_trans_80_f(:,1)),length(riv_trans_f(:,1)),1),'--','color','b','linew',2);
% plot(1:length(riv_trans_f(:,1)), repmat(nanmean(riv_trans_f(:,1)),length(riv_trans_f(:,1)),1),'--','color','r','linew',2);
legend('Nam_1_9_8_0_~_1_9_8_5', 'Nam_2_0_1_2_~_2_0_1_7');
ylabel(gca,'River discharge (m^3/s)','fontsize',15,'fontweight','bold'); xlabel(gca,'time(year)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 length(riv_trans_f(:,2))]); ylim([0 3650]);
set(gca,'xtick',[1:366:length(riv_trans_f(:,2))]);
set(gca,'xticklabel',1:1:6);

% no3
figure; hold on;
plot(r_no3_clim_80(1,:),'b','linew',2); plot(r_no3_clim(1,:),'r','linew',2); 
plot(1:366, repmat(mean(r_no3_clim_80(1,:)),366,1),'--','color','b','linew',2);
plot(1:366, repmat(nanmean(r_no3_clim(1,:)),366,1),'--','color','r','linew',2);
legend('sumjin_1_9_8_0_s', 'sumjin_2_0_1_5_~_2_0_1_7');
ylabel(gca,'NO3 (uM/L)','fontsize',15,'fontweight','bold'); xlabel(gca,'days','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 366]); ylim([0 140])

figure; hold on;
plot(r_no3_clim_80(2,:),'b','linew',2); plot(r_no3_clim(2,:),'r','linew',2); 
plot(1:366, repmat(mean(r_no3_clim_80(2,:)),366,1),'--','color','b','linew',2);
plot(1:366, repmat(nanmean(r_no3_clim(2,:)),366,1),'--','color','r','linew',2);
legend('Nam_1_9_8_0_s', 'Nam_2_0_1_5_~_2_0_1_7');
ylabel(gca,'NO3 (uM/L)','fontsize',15,'fontweight','bold'); xlabel(gca,'days','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 366]); ylim([0 140])

riv_no3_f(:,1) = [squeeze(riv_no3(1,20,:,1)); squeeze(riv_no3(1,20,:,2)); squeeze(riv_no3(1,20,:,3));];
riv_no3_f(:,2) = [squeeze(riv_no3(2,20,:,1)); squeeze(riv_no3(2,20,:,2)); squeeze(riv_no3(2,20,:,3));];
riv_no3_80_f(:,1) = [squeeze(riv_no3_80(1,20,:,1)); squeeze(riv_no3_80(1,20,:,2)); squeeze(riv_no3_80(1,20,:,3));];
riv_no3_80_f(:,2) = [squeeze(riv_no3_80(2,20,:,1)); squeeze(riv_no3_80(2,20,:,2)); squeeze(riv_no3_80(2,20,:,3));];

figure; hold on;
plot(riv_no3_80_f(:,1),'b','linew',2); plot(riv_no3_f(:,1),'r','linew',2); 
% plot(1:length(riv_no3_80_f(:,1)), repmat(nanmean(riv_no3_80_f(:,1)),length(riv_no3_f(:,1)),1),'--','color','b','linew',2);
% plot(1:length(riv_no3_f(:,1)), repmat(nanmean(riv_no3_f(:,1)),length(riv_no3_f(:,1)),1),'--','color','r','linew',2);
legend('sumjin_1_9_8_2_~_1_9_8_4', 'sumjin_2_0_1_5_~_2_0_1_7');
ylabel(gca,'NO3 (uM/L)','fontsize',15,'fontweight','bold'); xlabel(gca,'time(year)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 length(riv_no3_f(:,1))]); ylim([0 181]);
set(gca,'xtick',[1:366:length(riv_no3_f(:,1))]);
set(gca,'xticklabel',1:1:6);

figure; hold on;
plot(riv_no3_80_f(:,2),'b','linew',2); plot(riv_no3_f(:,2),'r','linew',2); 
% plot(1:length(riv_no3_80_f(:,1)), repmat(nanmean(riv_no3_80_f(:,1)),length(riv_no3_f(:,1)),1),'--','color','b','linew',2);
% plot(1:length(riv_no3_f(:,1)), repmat(nanmean(riv_no3_f(:,1)),length(riv_no3_f(:,1)),1),'--','color','r','linew',2);
legend('Nam_1_9_8_2_~_1_9_8_4', 'Nam_2_0_1_5_~_2_0_1_7');
ylabel(gca,'NO3 (uM/L)','fontsize',15,'fontweight','bold'); xlabel(gca,'time(year)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 length(riv_no3_f(:,2))]); ylim([0 181]);
set(gca,'xtick',[1:366:length(riv_no3_f(:,2))]);
set(gca,'xticklabel',1:1:6);

% nh4
figure; hold on;
plot(r_nh4_clim_80(1,:),'b','linew',2); plot(r_nh4_clim(1,:),'r','linew',2); 
plot(1:366, repmat(mean(r_nh4_clim_80(1,:)),366,1),'--','color','b','linew',2);
plot(1:366, repmat(nanmean(r_nh4_clim(1,:)),366,1),'--','color','r','linew',2);
legend('sumjin_1_9_8_0_s', 'sumjin_2_0_1_5_~_2_0_1_7');
ylabel(gca,'NH4 (uM/L)','fontsize',15,'fontweight','bold'); xlabel(gca,'days','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 366]); ylim([0 14.5]); % ylim([-inf 14.5])

figure; hold on;
plot(r_nh4_clim_80(2,:),'b','linew',2); plot(r_nh4_clim(2,:),'r','linew',2); 
plot(1:366, repmat(mean(r_nh4_clim_80(2,:)),366,1),'--','color','b','linew',2);
plot(1:366, repmat(nanmean(r_nh4_clim(2,:)),366,1),'--','color','r','linew',2);
legend('Nam_1_9_8_0_s', 'Nam_2_0_1_5_~_2_0_1_7');
ylabel(gca,'NH4 (uM/L)','fontsize',15,'fontweight','bold'); xlabel(gca,'days','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 366]); ylim([0 14.5])

riv_nh4_f(:,1) = [squeeze(riv_nh4(1,20,:,1)); squeeze(riv_nh4(1,20,:,2)); squeeze(riv_nh4(1,20,:,3));];
riv_nh4_f(:,2) = [squeeze(riv_nh4(2,20,:,1)); squeeze(riv_nh4(2,20,:,2)); squeeze(riv_nh4(2,20,:,3));];
riv_nh4_80_f(:,1) = [squeeze(riv_nh4_80(1,20,:,1)); squeeze(riv_nh4_80(1,20,:,2)); squeeze(riv_nh4_80(1,20,:,3));];
riv_nh4_80_f(:,2) = [squeeze(riv_nh4_80(2,20,:,1)); squeeze(riv_nh4_80(2,20,:,2)); squeeze(riv_nh4_80(2,20,:,3));];

figure; hold on;
plot(riv_nh4_80_f(:,1),'b','linew',2); plot(riv_nh4_f(:,1),'r','linew',2); 
% plot(1:length(riv_nh4_80_f(:,1)), repmat(nanmean(riv_nh4_80_f(:,1)),length(riv_nh4_f(:,1)),1),'--','color','b','linew',2);
% plot(1:length(riv_nh4_f(:,1)), repmat(nanmean(riv_nh4_f(:,1)),length(riv_nh4_f(:,1)),1),'--','color','r','linew',2);
legend('sumjin_1_9_8_2_~_1_9_8_4', 'sumjin_2_0_1_5_~_2_0_1_7');
ylabel(gca,'NH4 (uM/L)','fontsize',15,'fontweight','bold'); xlabel(gca,'time(year)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 length(riv_nh4_f(:,1))]); ylim([0 28]);
set(gca,'xtick',[1:366:length(riv_nh4_f(:,1))]);
set(gca,'xticklabel',1:1:6);

figure; hold on;
plot(riv_nh4_80_f(:,2),'b','linew',2); plot(riv_nh4_f(:,2),'r','linew',2); 
% plot(1:length(riv_nh4_80_f(:,1)), repmat(nanmean(riv_nh4_80_f(:,1)),length(riv_nh4_f(:,1)),1),'--','color','b','linew',2);
% plot(1:length(riv_nh4_f(:,1)), repmat(nanmean(riv_nh4_f(:,1)),length(riv_nh4_f(:,1)),1),'--','color','r','linew',2);
legend('Nam_1_9_8_2_~_1_9_8_4', 'Nam_2_0_1_5_~_2_0_1_7');
ylabel(gca,'NH4 (uM/L)','fontsize',15,'fontweight','bold'); xlabel(gca,'time(year)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
grid on; xlim([1 length(riv_nh4_f(:,2))]); ylim([0 28]);
set(gca,'xtick',[1:366:length(riv_nh4_f(:,2))]);
set(gca,'xticklabel',1:1:6);





