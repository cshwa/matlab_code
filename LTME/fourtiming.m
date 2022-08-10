% high 시기 low 시기 ctd -  idronaut 
clc; clear all; close all;

% load CTD high data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\1high';
filename = 'h_den.xlsx';
filename = fullfile(directory,filename);
h_den =  xlsread(filename);

filename = 'h_sal.xlsx';
filename = fullfile(directory,filename);
h_sal =  xlsread(filename);

filename = 'h_temp.xlsx';
filename = fullfile(directory,filename);
h_temp =  xlsread(filename);

% load CTD ebb data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\ctd\ebb';
filename = 'i_den.xlsx';
filename = fullfile(directory,filename);
e_den =  xlsread(filename);

filename = 'i_sal.xlsx';
filename = fullfile(directory,filename);
e_sal =  xlsread(filename);

filename = 'i_temp.xlsx';
filename = fullfile(directory,filename);
e_temp =  xlsread(filename);

% load CTD low data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\CTD_idronaut\3low';
filename = 'l_den.xlsx';
filename = fullfile(directory,filename);
l_den =  xlsread(filename);

filename = 'l_sal.xlsx';
filename = fullfile(directory,filename);
l_sal =  xlsread(filename);

filename = 'l_temp.xlsx';
filename = fullfile(directory,filename);
l_temp =  xlsread(filename);

% load CTD flood data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\ctd\flood';
filename = 'i_den.xlsx';
filename = fullfile(directory,filename);
f_den =  xlsread(filename);

filename = 'i_sal.xlsx';
filename = fullfile(directory,filename);
f_sal =  xlsread(filename);

filename = 'i_temp.xlsx';
filename = fullfile(directory,filename);
f_temp =  xlsread(filename);

%%
hfig = figure()
set(hfig, 'Position', [100 100 1400 600])
% high
subplot(341)
contourf(flipud(h_temp));
shading flat;
colorbar;
title('High','fontsize',14);
ylabel('depth (m)');
caxis([5 10]);
text(2,3,'temperature');
set(gca, 'ylim',[0 12],'ytick',[0:2:12],'yticklabel',[12:-2:0]');
set(gca, 'xtick',[1:3:14],'xticklabel',[1:3:14]);
subplot(345)
contourf(flipud(h_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,3,'salinity');
subplot(349)
contourf(flipud(h_den));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,3,'density');

% ebb
subplot(3,4,2)
contourf(flipud(e_temp));
shading flat;
title('Ebb','fontsize',14);
colorbar;caxis([5 10]);
text(2,3,'temperature');
subplot(3,4,6)
contourf(flipud(e_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,3,'salinity');
subplot(3,4,10)
contourf(flipud(e_den));
shading flat;
colorbar;colorbar;caxis([1000 1030]);
text(2,3,'density');
xlabel('station');

% low
subplot(343)
contourf(flipud(l_temp));
shading flat;
title('Low','fontsize',14);
colorbar;caxis([5 10]);
text(2,3,'temperature');
subplot(347)
contourf(flipud(l_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,3,'salinity');
subplot(3,4,11)
contourf(flipud(l_den));
shading flat;
colorbar;colorbar;caxis([1000 1030]);
text(2,3,'density');
xlabel('station');

% flood
subplot(344)
contourf(flipud(f_temp));
shading flat;
title('Flood','fontsize',14);
colorbar;caxis([5 10]);
text(2,3,'temperature');
subplot(348)
contourf(flipud(f_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,3,'salinity');
subplot(3,4,12)
contourf(flipud(f_den));
shading flat;
colorbar;colorbar;caxis([1000 1030]);
text(2,3,'density');
xlabel('station');

%% depth mean density
figure()
en = 10; % 수심 cutting
subplot(121)
[a b] = size(f_den(1:en,4:end));
plot(1:b, nanmean(f_den(1:en,4:end)),'r.','markersize',20);
hold on;
x = [1:b];
fd = nanmean(f_den(1:en,4:end));
p = polyfit(x,fd,1);
y = p(1)*x+p(2);
plot(x,y,'r-');
[a b] = size(e_den(1:en,4:end));
plot(1:b, nanmean(e_den(:,4:end,1:en)),'b.','markersize',20);
x = [1:b];
ed = nanmean(e_den(1:en,4:end));
p = polyfit(x,ed,1);
y = p(1)*x+p(2);
plot(x,y,'b-');

legend('flood','ebb');

subplot(122)
[a b] = size(h_den(1:en,4:end));
plot(1:b, nanmean(h_den(1:en,4:end)),'r.','markersize',20);
hold on;
x = [1:b];
hd = nanmean(h_den(1:en,4:end));
p = polyfit(x,hd,1);
y = p(1)*x+p(2);
plot(x,y,'r-');
[a b] = size(l_den(1:en,4:end));
plot(1:b, nanmean(l_den(1:en,4:end)),'b.','markersize',20);
x = [1:b];
ld = nanmean(l_den(1:en,4:end));
p = polyfit(x,ld,1);
y = p(1)*x+p(2);
plot(x,y,'b-');
legend('high','low');
plot(6.5,1010:0.1:1024,'k');
% ylim([1010 1024]);
xlim([0 13]);




%% strarification paramete
h = (max(h_sal)-min(h_sal))./nanmean(h_sal);
e = (max(e_sal)-min(e_sal))./nanmean(e_sal);
l = (max(l_sal)-min(l_sal))./nanmean(l_sal);
f = (max(f_sal)-min(f_sal))./nanmean(f_sal);

figure()
plot(h,'m.');
hold on;
plot(e,'b.');
plot(l,'g.');
plot(f,'r.');
legend('high','ebb','low','flood');
plot(1:0.01:15,0.32,'k');
plot(1:0.01:15,0.15,'k');

sp(3,:) = f;
sp(1,:) = e;
sp(2,:) = l;
save('sp20150226.dat','sp','-ascii');

%% horizontal density gradient
% diff_den(1) = nanmean(h_den(:,1)) - nanmean(h_den(:,end));
% diff_den(2) = nanmean(e_den(:,1)) - nanmean(e_den(:,end));
% diff_den(3) = nanmean(l_den(:,1)) - nanmean(l_den(:,end));
% 
% distance = 11.90*1000;
% NdiffDen = [7.96404415584425,7.81162261904751,5.86997285714301];
% NdiffDen = NdiffDen/distance;
% SdiffDen = diff_den/distance;
% 
% figure()
% plot(NdiffDen,'b.','markersize',20);
% hold on;
% plot(SdiffDen,'r.','markersize',20);
% legend('Neap','Spring');
% ylabel('\delta\rho/\deltax','fontsize',14);
% set(gca,'xtick',[1:1:3],'xlim',[0 4],'xticklabel',{'high','ebb','low'},'fontsize',14);

% save density gradient for comparing to other obseved data
hdg(4,:) = nanmean(h_den);
hdg(1,:) = nanmean(e_den);
hdg(2,:) = nanmean(l_den);
hdg(3,:) = nanmean(f_den);
save('hdg20150226.dat','hdg','-ascii');

distance = [0 1.158887931 2.308527646 3.547673354 4.606953311 ...
    5.609994106 6.536190309 7.509171283 8.732269443 9.889092071 ...
    10.78271688 11.60272967 12.55500728 13.37585421];
save('distance20131105.dat','distance','-ascii');


