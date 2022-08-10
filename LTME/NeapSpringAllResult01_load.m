% depth mean density - spring/neap 비교
clc; clear all; close all;

load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\am_flood\neap_flood_den.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\am_flood\neap_flood_temp.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\am_flood\neap_flood_sal.dat');

load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\pm_ebb\neap_ebb_den.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\pm_ebb\neap_ebb_temp.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\pm_ebb\neap_ebb_sal.dat');

load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\am_flood\spring_flood_den.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\am_flood\spring_flood_temp.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\am_flood\spring_flood_sal.dat');

load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\pm_ebb\spring_ebb_den.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\pm_ebb\spring_ebb_temp.dat');
load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\pm_ebb\spring_ebb_sal.dat');

% 수심 정리
% neap_flood_den = neap_flood_den(1:33,:);
% neap_flood_sal = neap_flood_sal(1:33,:);
% neap_flood_temp = neap_flood_temp(1:33,:);
% 
% spring_ebb_den = spring_ebb_den(1:33,:);
% spring_ebb_sal = spring_ebb_sal(1:33,:);
% spring_ebb_temp = spring_ebb_temp(1:33,:);
% 
% spring_flood_den = spring_flood_den(1:33,:);
% spring_flood_sal = spring_flood_sal(1:33,:);
% spring_flood_temp = spring_flood_temp(1:33,:);

%% 전체 비교 vs 상류비교
% ctd_en = 36;
% st = 06; en = 30;
% figure()
% % neap ebb
% subplot(341)
% contourf(flipud(neap_ebb_temp(2:ctd_en,:))); shading flat;
% % colorbar; 
% caxis([7 10]);
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% title('Neap ebb','fontsize',14);
% ylabel('Depth (m)','fontsize',14);
% 
% subplot(345)
% contourf(flipud(neap_ebb_sal(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([0 30]);
% 
% subplot(349)
% contourf(flipud(neap_ebb_den(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([1000 1027]);
% xlabel('Station','fontsize',14);
% 
% % neap flood
% subplot(342)
% contourf(flipud(neap_flood_temp(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([6.5 9.5]);
% title('Neap flood','fontsize',14);
% 
% subplot(346)
% contourf(flipud(neap_flood_sal(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([0 30]);
% 
% % spring ebb
% subplot(3,4,10)
% contourf(flipud(neap_flood_den(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([1000 1027]);
% 
% subplot(343)
% contourf(flipud(spring_ebb_temp(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([6.5 9.5]);
% title('Spring ebb','fontsize',14);
% 
% subplot(347)
% contourf(flipud(spring_ebb_sal(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([0 30]);
% 
% subplot(3,4,11)
% contourf(flipud(spring_ebb_den(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([1000 1027]);
% 
% subplot(344)
% contourf(flipud(spring_flood_temp(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([6.5 9.5]);
% title('Neap flood','fontsize',14);
% 
% subplot(348)
% contourf(flipud(spring_flood_sal(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([0 30]);
% 
% subplot(3,4,12)
% contourf(flipud(spring_flood_den(2:ctd_en,:))); shading flat;
% set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% caxis([1000 1027]);


%% 입구에서 상류 비교
% en = 15;
% figure()
% subplot(341)
% contourf(flipud(neap_ebb_temp(2:en,1:10))); shading flat;
% caxis([6.5 9.5]);
% title('Neap ebb','fontsize',14);
% subplot(345)
% contourf(flipud(neap_ebb_sal(2:en,1:10))); shading flat;
% caxis([0 30]);
% subplot(349)
% contourf(flipud(neap_ebb_den(2:en,1:10))); shading flat;
% caxis([1000 1027]);
% 
% subplot(342)
% contourf(flipud(neap_flood_temp(2:en,1:10))); shading flat;
% caxis([6.5 9.5]);
% title('Neap flood','fontsize',14);
% subplot(346)
% contourf(flipud(neap_flood_sal(2:en,1:10))); shading flat;
% caxis([0 30]);
% subplot(3,4,10)
% contourf(flipud(neap_flood_den(2:en,1:10))); shading flat;
% caxis([1000 1027]);
% 
% 
% subplot(343)
% contourf(flipud(spring_ebb_temp(2:en,1:10))); shading flat;
% caxis([6.5 9.5]);
% title('Spring ebb','fontsize',14);
% subplot(347)
% contourf(flipud(spring_ebb_sal(2:en,1:10))); shading flat;
% caxis([0 30]);
% subplot(3,4,11)
% contourf(flipud(spring_ebb_den(2:en,1:10))); shading flat;
% caxis([1000 1027]);
% 
% subplot(344)
% contourf(flipud(spring_flood_temp(2:en,1:10))); shading flat;
% caxis([6.5 9.5]); 
% title('Neap flood','fontsize',14);
% subplot(348)
% contourf(flipud(spring_flood_sal(2:en,1:10))); shading flat;
% caxis([0 30]);
% subplot(3,4,12)
% contourf(flipud(spring_flood_den(2:en,1:10))); shading flat;
% caxis([1000 1027]);

%% density 비교 with errorbar
% figure()
% subplot(221)
% data = neap_flood_den(2:en,1:10);
% y =  nanmean(data); 
% x = 1:size(data,2); 
% e = nanstd(data); 
% errorbar(x, y,e,'b.','markersize',20');
% set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 11]);
% set(gca,'ylim',[1000 1030]);
% title('Density (Neap flood)','fontsize',14);
% ylabel('Mean T','fontsize',14);
% 
% subplot(222)
% data = neap_ebb_den(2:en,1:10);
% y =  nanmean(data); 
% x = 1:size(data,2); 
% e = nanstd(data); 
% errorbar(x, y,e,'b.','markersize',20');
% set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 11]);
% set(gca,'ylim',[1000 1030]);
% title('Density (Neap ebb)','fontsize',14);
% ylabel('Mean T','fontsize',14);
% 
% subplot(223)
% data = spring_flood_den(2:en,1:10);
% y =  nanmean(data); 
% x = 1:size(data,2); 
% e = nanstd(data); 
% errorbar(x, y,e,'r.','markersize',20');
% set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 11]);
% set(gca,'ylim',[1000 1030]);
% title('Density (Spring flood)','fontsize',14);
% ylabel('Mean T','fontsize',14);
% 
% subplot(224)
% data = spring_ebb_den(2:en,1:10);
% y =  nanmean(data); 
% x = 1:size(data,2); 
% e = nanstd(data); 
% errorbar(x, y,e,'r.','markersize',20');
% set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 11]);
% set(gca,'ylim',[1000 1030]);
% title('Density (Spring ebb)','fontsize',14);
% ylabel('Mean T','fontsize',14);



%% density 비교
% st = 6; en = 30;
% st_en = 22;
% figure()
% plot(nanmean(neap_ebb_den(2:end,1:st_en)),'b.-','markersize',20);
% hold on;
% plot(nanmean(neap_flood_den(2:end,1:st_en)),'b.-','markersize',20);
% plot(nanmean(spring_ebb_den(2:end,1:st_en)),'r.-','markersize',20);
% plot(nanmean(spring_flood_den(2:end,1:st_en)),'r.-','markersize',20);
% set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
% xlim([0 st_en]);
% ylabel('Density','fontsize',14);
% xlabel('Stiation','fontsize',14);

%% density 비교 - 깊은골 무시
% depth_en = 7;
% st_st = 10;
% st_en = 15;
% figure()
% plot(nanmean(neap_ebb_den(2:depth_en,st_st:st_en)),'b.-','markersize',20);
% hold on;
% plot(nanmean(neap_flood_den(2:depth_en,st_st:st_en)),'b.-','markersize',20);
% plot(nanmean(spring_ebb_den(2:depth_en,st_st:st_en)),'r.-','markersize',20);
% plot(nanmean(spring_flood_den(2:depth_en,st_st:st_en)),'r.-','markersize',20);
% xlim([0 st_en]);
% ylabel('Density','fontsize',14);
% xlabel('Stiation','fontsize',14);