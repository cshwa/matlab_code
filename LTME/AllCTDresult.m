% 2015 spring Neap / Spring CTD 결과 비교
clc; clear all; close all;

load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\am_ebb\neap_ebb_temp.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\am_ebb\neap_ebb_den.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\am_ebb\neap_ebb_sal.dat');

load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\pm_flood\neap_flood_temp.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\pm_flood\neap_flood_den.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150513_CTD_neap\pm_flood\neap_flood_sal.dat');

load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\pm_ebb\spring_ebb_temp.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\pm_ebb\spring_ebb_sal.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\pm_ebb\spring_ebb_den.dat');

load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\am_flood\spring_flood_temp.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\am_flood\spring_flood_sal.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_spring_ys\CTD\150519_CTD_spring\am_flood\spring_flood_den.dat');

%% 수심 정리 32m 까지만 사용
dep_en = 33;
neap_flood_den = neap_flood_den(1:dep_en,:);
neap_flood_sal = neap_flood_sal(1:dep_en,:);
neap_flood_temp = neap_flood_temp(1:dep_en,:);

neap_ebb_den = neap_ebb_den(1:dep_en,:);
neap_ebb_sal = neap_ebb_sal(1:dep_en,:);
neap_ebb_temp = neap_ebb_temp(1:dep_en,:);

spring_ebb_den = spring_ebb_den(1:dep_en,:);
spring_ebb_sal = spring_ebb_sal(1:dep_en,:);
spring_ebb_temp = spring_ebb_temp(1:dep_en,:);

spring_flood_den = spring_flood_den(1:dep_en,:);
spring_flood_sal = spring_flood_sal(1:dep_en,:);
spring_flood_temp = spring_flood_temp(1:dep_en,:);

%% caxis 경계값 설정
sal_st = 12;
sal_en = 20;
%% 전체 비교 vs 상류비교
figure()
% neap ebb
subplot(341)
contourf(flipud(neap_ebb_temp(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([sal_st sal_en]);
title('Neap ebb','fontsize',14);
ylabel('Depth (m)','fontsize',14);
text(2, 6, 'Temperature','fontsize',14);

subplot(345)
contourf(flipud(neap_ebb_sal(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([0 30]);
text(2, 6, 'Salinity','fontsize',14);

subplot(349)
contourf(flipud(neap_ebb_den(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([1000 1027]);
xlabel('Station','fontsize',14);
text(2, 6, 'Density','fontsize',14);


% neap flood
subplot(342)
contourf(flipud(neap_flood_temp(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([sal_st sal_en]);
title('Neap flood','fontsize',14);

subplot(346)
contourf(flipud(neap_flood_sal(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([0 30]);

subplot(3,4,10)
contourf(flipud(neap_flood_den(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([1000 1027]);

% spring ebb
subplot(343)
contourf(flipud(spring_ebb_temp(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([sal_st sal_en]);
title('Spring ebb','fontsize',14);

subplot(347)
contourf(flipud(spring_ebb_sal(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([0 30]);

subplot(3,4,11)
contourf(flipud(spring_ebb_den(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([1000 1027]);

% spring flood
subplot(344)
contourf(flipud(spring_flood_temp(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([sal_st sal_en]);
title('Neap flood','fontsize',14);

subplot(348)
contourf(flipud(spring_flood_sal(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([0 30]);

subplot(3,4,12)
contourf(flipud(spring_flood_den(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[2 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 30]);
caxis([1000 1027]);


%% 입구에서 상류 비교
y_en = 22;  % 전체 y_en = 2;
x_en = 15;  % 전체 x_en = 30;

figure()
% neap ebb
subplot(341)
contourf(flipud(neap_ebb_temp(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([sal_st sal_en]);
title('Neap ebb','fontsize',14);
ylabel('Depth (m)','fontsize',14);
text(2, y_en+2, 'Temperature','fontsize',14);

subplot(345)
contourf(flipud(neap_ebb_sal(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([0 30]);
text(2, y_en+2, 'Salinity','fontsize',14);

subplot(349)
contourf(flipud(neap_ebb_den(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([1000 1027]);
xlabel('Station','fontsize',14);
text(2, y_en+2, 'Density','fontsize',14);


% neap flood
subplot(342)
contourf(flipud(neap_flood_temp(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([sal_st sal_en]);
title('Neap flood','fontsize',14);

subplot(346)
contourf(flipud(neap_flood_sal(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([0 30]);

subplot(3,4,10)
contourf(flipud(neap_flood_den(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([1000 1027]);

% spring ebb
subplot(343)
contourf(flipud(spring_ebb_temp(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([sal_st sal_en]);
title('Spring ebb','fontsize',14);

subplot(347)
contourf(flipud(spring_ebb_sal(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([0 30]);

subplot(3,4,11)
contourf(flipud(spring_ebb_den(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([1000 1027]);

% spring flood
subplot(344)
contourf(flipud(spring_flood_temp(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([sal_st sal_en]);
title('Neap flood','fontsize',14);

subplot(348)
contourf(flipud(spring_flood_sal(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([0 30]);

subplot(3,4,12)
contourf(flipud(spring_flood_den(2:end,:))); shading flat;
set(gca,'ytick',[2:5:32],'yticklabel',[30:-5:0],'ylim',[y_en 32]);
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[1 x_en]);
caxis([1000 1027]);

%% density 비교 with errorbar
en = 30;
st_en = 15;
den_st = 990;
den_en = 1030;

figure()
subplot(221)
data = neap_flood_den(2:en,1:st_en);
y =  nanmean(data); 
x = 1:size(data,2); 
e = nanstd(data); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[0 st_en]);
set(gca,'ylim',[den_st den_en]);
title('Neap flood','fontsize',14);
ylabel('<\rho>','fontsize',14);

subplot(222)
data = neap_ebb_den(2:en,1:st_en);
y =  nanmean(data); 
x = 1:size(data,2); 
e = nanstd(data); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[0 st_en]);
set(gca,'ylim',[den_st den_en]);
title('Neap ebb','fontsize',14);


subplot(223)
data = spring_flood_den(2:en,1:st_en);
y =  nanmean(data); 
x = 1:size(data,2); 
e = nanstd(data); 
errorbar(x, y,e,'r.','markersize',20');
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[0 st_en]);
set(gca,'ylim',[den_st den_en]);
title('Spring flood','fontsize',14);
xlabel('Station','fontsize',14);

subplot(224)
data = spring_ebb_den(2:en,1:st_en);
y =  nanmean(data); 
x = 1:size(data,2); 
e = nanstd(data); 
errorbar(x, y,e,'r.','markersize',20');
set(gca,'xtick',[1:3:30],'xticklabel',[1:3:30],'xlim',[0 st_en]);
set(gca,'ylim',[den_st den_en]);
title('Spring ebb','fontsize',14);



%% density 비교
st_en = 15;
figure()
plot(nanmean(neap_ebb_den(2:end,1:st_en)),'b.-','markersize',20);
hold on;
plot(nanmean(neap_flood_den(2:end,1:st_en)),'b.-','markersize',20);
plot(nanmean(spring_ebb_den(2:end,1:st_en)),'r.-','markersize',20);
plot(nanmean(spring_flood_den(2:end,1:st_en)),'r.-','markersize',20);
xlim([0 st_en]);
ylabel('Density','fontsize',14);
xlabel('Stiation','fontsize',14);

%% density 비교 - 깊은골 무시
depth_en = 7;
st_st = 1;
st_en = 24;
figure()
plot(nanmean(neap_ebb_den(2:depth_en,st_st:st_en)),'b.-','markersize',20);
hold on;
plot(nanmean(neap_flood_den(2:depth_en,st_st:st_en)),'b.-','markersize',20);
plot(nanmean(spring_ebb_den(2:depth_en,st_st:st_en)),'r.-','markersize',20);
plot(nanmean(spring_flood_den(2:depth_en,st_st:st_en)),'r.-','markersize',20);
xlim([0 st_en]);
ylabel('Density','fontsize',14);
xlabel('Stiation','fontsize',14);