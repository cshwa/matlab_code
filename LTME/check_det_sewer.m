close all; clear; clc;
t_year = 2001
Fname =['river_' num2str(t_year) '_realts_biofennel_GY_sewer_det.nc'];

dis=ncread(Fname,'river_transport');
Lden=ncread(Fname,'river_LDeN');
sden=ncread(Fname,'river_SDeN');

figure;plot(Lden(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');


figure;plot(squeeze(Lden(1:2,20,:))'./1000*14); ylim([0 30])
figure;plot(squeeze(sden(:,20,:))'); ylim([0 30])

figure;plot((dis.*squeeze(Lden(:,20,:)))');