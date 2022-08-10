close all; clear; clc;
riv=load('datong_1961_2017.dat');
dis=riv(:,2);

dis_1=dis(362:481,1); %1991~2000
dis_2=dis(482:601,1); %2001~2010

yt=[1:10]';

yt_s  = num2str(yt);

figure;
plot(dis(362:481,1),'b'); hold on;  plot(dis(482:601,1),'r');
xlim([1 120]); alpha(0.3); grid on; ylabel('datong discharge (m^3/s)'); xlabel('year');
set(gca,'fontsize',13); xticks(1:12:120); xticklabels(gca,yt_s)

for i = 1:12
dis_1_clim(i)= mean(dis_1(i:12:end));
dis_2_clim(i)= mean(dis_2(i:12:end));
end

figure;
plot(dis_1_clim,'b'); hold on; plot(dis_2_clim,'r'); 
text(2,50000,['1991~2000 : ',num2str(mean(dis_1_clim),'%0.2f'),' m^3/s'],'color','b','fontweight','bold');
text(2,45000,['2001~2010 : ',num2str(mean(dis_2_clim),'%0.2f'),' m^3/s'],'color','r','fontweight','bold');
xlim([1 12]); alpha(0.3); grid on; ylabel('datong discharge (m^3/s)'); xlabel('month');
set(gca,'fontsize',13); xticks(1:12)