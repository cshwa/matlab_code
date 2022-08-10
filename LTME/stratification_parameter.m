%% salinity stratifiacation parameter
% clc; clear all; close all;

load neap_ebb_sal.dat
load neap_flood_sal.dat
load spring_ebb_sal.dat
load spring_flood_sal.dat

ctd_en = 36;
st = 01;
en = 30;

%%
figure()
subplot(221)
sal = neap_flood_sal;
contourf(flipud(sal(3:ctd_en,:))); shading flat;
colormap jet
colorbar;  caxis([0 30]);
text(3,5,'Salinity','fontsize',12);
ylabel('Depth (m)','fontsize',12);
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[20 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);

subplot(222)
sal = spring_flood_sal;
contourf(flipud(sal(3:ctd_en,:))); shading flat;
colorbar; colorbar; caxis([0 30]);
text(st+1,5,'Salinity','fontsize',12);
ylabel('Depth (m)','fontsize',12);
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[20 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);

subplot(223)
data = neap_flood_sal(3:end,:);
sp = abs(nanmax(data)-nanmin(data))./nanmean(data);
plot(sp,'b.');
ylim([0 3]);
% plot(0:0.5:30, 0.15,'k.','markersize',2);
% plot(0:0.5:30, 0.32,'k.','markersize',2);
% title('Stratification parameter (Spring flood)','fontsize',14);
xlabel('Station number','fontsize',14);
ylabel('\deltaS/<S>','fontsize',14);
subplot(224)
data = spring_flood_sal(3:end,:);
sp = abs(nanmax(data)-nanmin(data))./nanmean(data);
plot(sp,'r.');
ylim([0 3]);
% plot(0:0.5:30, 0.15,'k.','markersize',2);
% plot(0:0.5:30, 0.32,'k.','markersize',2);
% title('Stratification parameter (Spring flood)','fontsize',14);
xlabel('Station number','fontsize',14);

%%
figure()
subplot(221)
sal = neap_ebb_sal;
contourf(flipud(sal(3:ctd_en,:))); shading flat;
colorbar; colorbar; caxis([0 30]);
t = text(st+1,5,'Neap','fontsize',12);
t.FontSize = 14;
ylabel('Depth (m)','fontsize',12);
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);

subplot(222)
sal = spring_ebb_sal;
contourf(flipud(sal(3:ctd_en,:))); shading flat;
colorbar; colorbar; caxis([0 30]);
t = text(st+1,5,'Spring','fontsize',12);
t.FontSize = 14;
ylabel('Depth (m)','fontsize',12);
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);

subplot(223)
data = neap_ebb_sal(3:end,:);
sp = abs(nanmax(data)-nanmin(data))./nanmean(data);
plot(sp,'b.');
ylim([0 3]);
% plot(0:0.5:30, 0.15,'k.','markersize',2);
% plot(0:0.5:30, 0.32,'k.','markersize',2);
% title('Stratification parameter (Spring flood)','fontsize',14);
% xlabel('Station number','fontsize',14);

subplot(224)
data = spring_ebb_sal(3:end,:);
sp = abs(nanmax(data)-nanmin(data))./nanmean(data);
plot(sp,'r.');
ylim([0 3]);
% plot(0:0.5:30, 0.15,'k.','markersize',2);
% plot(0:0.5:30, 0.32,'k.','markersize',2);
% title('Stratification parameter (Spring flood)','fontsize',14);
% xlabel('Station number','fontsize',14);
