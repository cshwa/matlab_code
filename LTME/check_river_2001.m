temp=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_temp');
dis=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');





temp=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_Oxyg');
dis=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');


temp=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_chlorophyll');
dis=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');


temp=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_NH4');
dis=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)).*14./1000);hold on; plot(squeeze(temp(2,20,:)).*14./1000,'r');



temp=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_NO3');
dis=ncread('river_2001_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)).*14./1000);hold on; plot(squeeze(temp(2,20,:)).*14./1000,'r');

