cd F:\ROMS\roms_tools\Run
start
cd D:\장기생태\Dynamic\06_river
close all; clear; clc; 
trans=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_transport');
no3=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_NO3');
nh4=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_NH4');
chl=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_chlorophyll');

figure; 
plot(trans(1,:),'r'); hold on; plot(trans(2,:),'b'); legend('songjung','Nam');

figure; 
plot(squeeze(no3(1,20,:)),'r'); hold on; plot(squeeze(no3(2,20,:)),'b')
xlabel('time (days)','fontsize',13); xlim([1 366]);
ylabel('NO3 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily MODEL river NO3']);


figure; 
plot(squeeze(nh4(1,20,:)),'r'); hold on; plot(squeeze(nh4(2,20,:)),'b')
xlabel('time (days)','fontsize',13); xlim([1 366]);
ylabel('NH4 (mmol N / m^3)','fontsize',13)
grid on;
title(['GY 2001 daily MODEL river NH4']);

figure; 
plot(squeeze(chl(1,20,:)),'r'); hold on; plot(squeeze(chl(2,20,:)),'b')
xlabel('time (days)','fontsize',13); xlim([1 366]);
ylabel('Chla (ug/L)','fontsize',13)
grid on;
title(['GY 2001 daily MODEL river Chla']);



close all; clear; clc; 
trans=ncread('river_2001_realts_biofennel_GY_jinwall.nc','river_transport');
no3=ncread('river_2001_realts_biofennel_GY_jinwall.nc','river_NO3');
nh4=ncread('river_2001_realts_biofennel_GY_jinwall.nc','river_NH4');
chl=ncread('river_2001_realts_biofennel_GY_jinwall.nc','river_chlorophyll');

figure; 
plot(trans(1,:),'r'); hold on; plot(trans(2,:),'b'); legend('songjung','Nam');
xlim([1 365]);
ylim([0 5000])

figure; 
plot(squeeze(no3(1,20,:)),'r'); hold on; plot(squeeze(no3(2,20,:)),'b')

figure; 
plot(squeeze(nh4(1,20,:)),'r'); hold on; plot(squeeze(nh4(2,20,:)),'b')

figure; 
plot(squeeze(chl(1,20,:)),'r'); hold on; plot(squeeze(chl(2,20,:)),'b')

load('jinwall_reconst_climate_3sig.mat');

figure; hold on;
plot(squeeze(nh4(1,20,:)),'ro','linew',2); 
plot(yp_est_nh4,'linew',2); % 04~ NO3
plot(yp_est_nh4_adv,'g','linew',2); % ~04 NO3
xlim([1 366]);
xlabel('time(days)','fontsize',13)
ylabel('NH4-N (mmol N /m^3)','fontsize',13)
grid on;
set(gca,'fontsize',13)
legend('input','jinwall-re','jinwall-re(past)');


figure; hold on;
plot(squeeze(no3(1,20,:)),'ro','linew',2); 
plot(yp_est_no3,'linew',2); % 04~ NO3
plot(yp_est_no3_adv,'g','linew',2); % ~04 NO3
xlim([1 366]);
xlabel('time(days)','fontsize',13)
ylabel('NO3-N (mmol N /m^3)','fontsize',13)
grid on;
set(gca,'fontsize',13)
legend('input','jinwall-re','jinwall-re(past)');

