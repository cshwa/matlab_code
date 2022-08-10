
ctd_en = 38;

sal = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\pm_ebb\neap_ebb_sal.dat');
station = 22;
y =  nanmean(sal(2:ctd_en,:)); 
x = 1:station; 
e = nanstd(sal(2:ctd_en,:)); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
set(gca,'ylim',[0 40]);
title('Depth mean salinity ','fontsize',14);
ylabel('¢¶','fontsize',14);

sal = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150226_ys_neap tide\am_flood\neap_flood_sal.dat');
station = 25;
hold on;
y =  nanmean(sal(2:ctd_en,:)); 
x = 1:station; 
e = nanstd(sal(2:ctd_en,:)); 
errorbar(x, y,e,'c.','markersize',20');
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
set(gca,'ylim',[-2 40]);
title('Depth mean salinity ','fontsize',14);
ylabel('¢¶','fontsize',14);

sal = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\pm_ebb\spring_ebb_sal.dat');
station = 25;
y =  nanmean(sal(2:ctd_en,:)); 
x = 1:station; 
e = nanstd(sal(2:ctd_en,:)); 
errorbar(x, y,e,'r.','markersize',20');
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
set(gca,'ylim',[-2 40]);
title('Depth mean salinity ','fontsize',14);
ylabel('¢¶','fontsize',14);

sal = load('D:\SilverStar\data\01_Sumjin_yeosu_camp\2015_winter_ys\20150305 ys_spring tide\am_flood\spring_flood_sal.dat');
y =  nanmean(sal(2:ctd_en,:)); 
x = 1:station; 
e = nanstd(sal(2:ctd_en,:)); 
errorbar(x, y,e,'m.','markersize',20');
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
set(gca,'ylim',[-2 40]);
title('Depth mean salinity ','fontsize',14);
ylabel('¢¶','fontsize',14);

