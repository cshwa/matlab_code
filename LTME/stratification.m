% flood ½Ã stratification
clc; clear all; close all;


% wn = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\am_flood\neap_flood_sal.dat');
% ws =load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\am_flood\spring_flood_sal.dat');
% ebb 
wn = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\pm_ebb\spring_ebb_sal.dat');
ws = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\pm_ebb\spring_ebb_sal.dat');

sn = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\pm_flood\neap_flood_sal.dat');
ss = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\am_flood\spring_flood_sal.dat');

sun = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150827_neap\pm_flood\neap_flood_sal.dat');
sus = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_summer_ys\CTD\150816_spring\am_flood\neap_flood_sal.dat');
% wn = fliplr(wn);
% ws = fliplr(ws);
% sn = fliplr(sn);
% ss = fliplr(ss);
%%
FigHandle = figure;
set(FigHandle, 'Position', [100, 30, 850, 600])
subplot(311)
[a b]= size(wn);
x = 6:30;
for i = 1:b
y(i) = abs(nanmax(wn(3:end,i))-wn(2,i))/nanmean(wn(2:end,i));
end
mean(y(7:8))

plot(x,y,'b.','markersize',20);
hold on;
[a b]= size(ws);
x = 6:30;
for i = 1:b
y(i) = abs(nanmax(ws(3:end,i))-ws(2,i))/nanmean(ws(2:end,i));
end
mean(y(7:8))

plot(x,y,'r.','markersize',20);
plot(0:0.1:30,0.15,'k');
plot(0:0.1:30, 0.32,'k');
ylabel('\deltaS/<S>','fontsize',14);
ylim([0 2]);

subplot(312)
hold on;
[a b]= size(sn);
x = 1:25;
for i = 1:b
y(i) = abs(nanmax(sn(3:end,i))-sn(2,i))/nanmean(sn(2:end,i));
end
mean(y(12:13))
plot(x,y,'b.','Markersize',20);

hold on;
[a b]= size(ss);
x = 1:30;
for i = 1:b
y(i) = abs(nanmax(ss(3:end,i))-ss(2,i))/nanmean(ss(2:end,i));
end
mean(y(12:13))
plot(x,y,'r.','Markersize',20');
plot(0:0.1:30,0.15,'k');
plot(0:0.1:30, 0.32,'k');
ylabel('\deltaS/<S>','fontsize',14);
ylim([0 2]);

subplot(313)
[a b]= size(sun);
x = 1:30;
for i = 1:b
y(i) = abs(nanmax(sun(3:end,i))-sun(2,i))/nanmean(sun(2:end,i));
end
mean(y(12:13))
plot(x,y,'b.','markersize',20);
% plot(wn(1,:),y,'b.');
ylabel('\deltaS/<S>','fontsize',14);

hold on;
[a b]= size(sus);
x = 1:30;
for i = 1:b
y(i) = abs(nanmax(sus(3:end,i))-sus(2,i))/nanmean(sus(2:end,i));
end
mean(y(12:13))
plot(x,y,'r.','markersize',20);
plot(0:0.1:30,0.15,'k');
plot(0:0.1:30, 0.32,'k');
ylim([0 2]);
legend('Neap','Spring');
ylabel('\deltaS/<S>','fontsize',14);
xlabel('Station','fontsize',20);
text(7,-0.2,'Landward','fontsize',14);
text(22,-0.2,'Seaward','fontsize',14);
text(20,1,'> 0.32 Stratified','fontsize',14);
text(20,0.6,'< 0.15 Well-mixed','fontsize',14);

%% depth mean salinity
%{
figure()

subplot(311)
temp = wn;
x = temp(1,:);
mean_sal = nanmean(temp(2:end,:));
plot(6:30,mean_sal,'b.','markersize',15);
hold on;
temp = ws;
x = temp(1,:);
mean_sal = nanmean(temp(2:end,:));
plot(6:30,mean_sal,'r.','markersize',15);
xlim([0 30]);
ylabel('Salinity','fontsize',14);
text(1,32,'Winter','fontsize',14);
legend('Neap','Spring');
title('Depth mean salinity','fontsize',14);


subplot(312)
temp = sn;
x = temp(1,:);
mean_sal = nanmean(temp(2:end,:));
plot(1:24,mean_sal,'b.','markersize',15);
hold on;
temp = ss;
x = temp(1,:);
mean_sal = nanmean(temp(2:end,:));
plot(1:30,mean_sal,'r.','markersize',15);
ylabel('Salinity','fontsize',14);
text(1,32,'Spring','fontsize',14);


subplot(313)
temp = sun;
x = temp(1,:);
mean_sal = nanmean(temp(2:end,:));
plot(1:30,mean_sal,'b.','markersize',15);
hold on;
temp = sus;
x = temp(1,:);
mean_sal = nanmean(temp(2:end,:));
plot(1:30,mean_sal,'r.','markersize',15);
xlabel('Station','fontsize',14);
ylabel('Salinity','fontsize',14);
text(1,32,'Summer','fontsize',14);

%}



