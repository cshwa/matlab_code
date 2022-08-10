% flood 시기ctd 비교 - castaway, idronaut 
clc; clear all; close all;

% load Castaway data
directory = 'D:\SilverStar\data\01_SumJin\2015_02_26\ctd\flood';
filename = 'ca_den.xlsx';
filename = fullfile(directory,filename);
ca_den =  xlsread(filename);

filename = 'ca_sal.xlsx';
filename = fullfile(directory,filename);
ca_sal =  xlsread(filename);

filename = 'ca_temp.xlsx';
filename = fullfile(directory,filename);
ca_temp =  xlsread(filename);

filename = 'i_den.xlsx';
filename = fullfile(directory,filename);
i_den =  xlsread(filename);

filename = 'i_sal.xlsx';
filename = fullfile(directory,filename);
i_sal =  xlsread(filename);

filename = 'i_temp.xlsx';
filename = fullfile(directory,filename);
i_temp =  xlsread(filename);

for i = 1: floor(38/3)
    ca_salm(i,:) = nanmean(ca_sal((i-1)*3+1:i*3,:));
    ca_denm(i,:) = nanmean(ca_den((i-1)*3+1:i*3,:));
    ca_tempm(i,:) = nanmean(ca_temp((i-1)*3+1:i*3,:));
end
% ca_salm(i+1,:) = nanmean(ca_sal(37:38,:));
% ca_tenm(i+1,:) = nanmean(ca_den(37:38,:));
% ca_tempm(i+1,:) = nanmean(ca_temp(37:38,:));

hfig = figure()
set(hfig, 'Position', [100 100 1200 200])
for i = 1:14
subplot(1,14,i)
hold on;
plot(ca_salm(:,i+1),'r','linewidth',2);
plot(i_sal(:,i),'b','linewidth',2);
% legend('ca','i');
az = 90;
el = 90;
view(az, el);
ylim([10 35]);
xlim([0 15]);
if i == 1
    xlabel('depth');
    ylabel('salinity');
end
end

hfig = figure()
set(hfig, 'Position', [100 100 1200 200])
for i = 1:14
subplot(1,14,i)
hold on;
plot(ca_denm(:,i+1),'r','linewidth',2);
plot(i_den(:,i),'b','linewidth',2);
% legend('ca','i');
az = 90;
el = 90;
view(az, el);
ylim([1010 1030]);
xlim([0 15]);
if i == 1
    xlabel('depth');
    ylabel('density');
end
end

hfig = figure()
set(hfig, 'Position', [100 100 1200 200])
for i = 1:14
subplot(1,14,i)
hold on;
plot(ca_tempm(:,i+1),'r','linewidth',2);
plot(i_temp(:,i),'b','linewidth',2);
% legend('ca','i');
az = 90;
el = 90;
view(az, el);
ylim([6 10]);
xlim([0 15]);
if i == 1
    xlabel('depth');
    ylabel('temperature');
end
end

diff_temp = i_temp - ca_tempm(:,2:end);
max(max(diff_temp))
min(min(diff_temp))


diff_sal = i_sal - ca_salm(:,2:end);
max(max(diff_sal))
min(min(diff_sal))


diff_den = i_den - ca_denm(:,2:end);
max(max(diff_den))
min(min(diff_den))

figure()
plot(i_temp(1,:),'b.')
hold on;
plot(ca_tempm(1,2:end),'r.')
legend('Idronaut','Castaway');
title('surface temperature','fontsize',14)
ylabel('temperature','fontsize',14);
xlabel('seaward                   station                    landward','fontsize',14);

figure()
plot(i_sal(1,:),'b.')
hold on;
plot(ca_salm(1,2:end),'r.')
legend('Idronaut','Castaway');
title('surface salinity','fontsize',14)
ylabel('salinity','fontsize',14);
xlabel('seaward                   station                    landward','fontsize',14);

hfig = figure()
set(hfig, 'Position', [100 100 1000 600])
subplot(321)
contourf(flipud(ca_tempm(:,2:end)));
shading flat;
colorbar;
title('Castaway','fontsize',14);
ylabel('depth (m)');
caxis([5 10]);
text(2,3,'temperature');
subplot(323)
contourf(flipud(ca_salm(:,2:end)));
shading flat;
colorbar;caxis([10 35]);
text(2,3,'salinity');
subplot(325)
contourf(flipud(ca_denm(:,2:end)));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,3,'density');
subplot(322)
contourf(flipud(i_temp));
shading flat;
title('Idronaut','fontsize',14);
colorbar;caxis([5 10]);
text(2,3,'temperature');
subplot(324)
contourf(flipud(i_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,3,'salinity');
subplot(326)
contourf(flipud(i_den));
shading flat;
colorbar;colorbar;caxis([1000 1030]);
text(2,3,'density');
xlabel('station');

%% test

hfig = figure()
set(hfig, 'Position', [100 100 1000 600])
subplot(321)
contourf(flipud(ca_temp(:,2:end)));
shading flat;
colorbar;
title('Castaway','fontsize',14);
ylabel('depth (m)');
caxis([5 10]);
text(2,3,'temperature');
subplot(323)
contourf(flipud(ca_sal(:,2:end)));
shading flat;
colorbar;caxis([10 35]);
text(2,3,'salinity');
subplot(325)
contourf(flipud(ca_den(:,2:end)));
shading flat;
xlabel('station');
colorbar;caxis([1000 1030]);
text(2,3,'density');
subplot(322)
contourf(flipud(i_temp));
shading flat;
title('Idronaut','fontsize',14);
colorbar;caxis([5 10]);
text(2,3,'temperature');
subplot(324)
contourf(flipud(i_sal));
shading flat;
colorbar;caxis([10 35]);
text(2,3,'salinity');
subplot(326)
contourf(flipud(i_den));
shading flat;
colorbar;colorbar;caxis([1000 1030]);
text(2,3,'density');
xlabel('station');

