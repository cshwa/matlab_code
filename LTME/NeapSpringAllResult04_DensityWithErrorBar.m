%% density 비교 with errorbar
figure()
subplot(221)
data = neap_flood_den(2:en,1:10);
y =  nanmean(data); 
x = 1:size(data,2); 
e = nanstd(data); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 11]);
set(gca,'ylim',[1000 1030]);
title('Density (Neap flood)','fontsize',14);
ylabel('Mean T','fontsize',14);

subplot(222)
data = neap_ebb_den(2:en,1:10);
y =  nanmean(data); 
x = 1:size(data,2); 
e = nanstd(data); 
errorbar(x, y,e,'b.','markersize',20');
set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 11]);
set(gca,'ylim',[1000 1030]);
title('Density (Neap ebb)','fontsize',14);
ylabel('Mean T','fontsize',14);

subplot(223)
data = spring_flood_den(2:en,1:10);
y =  nanmean(data); 
x = 1:size(data,2); 
e = nanstd(data); 
errorbar(x, y,e,'r.','markersize',20');
set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 11]);
set(gca,'ylim',[1000 1030]);
title('Density (Spring flood)','fontsize',14);
ylabel('Mean T','fontsize',14);

subplot(224)
data = spring_ebb_den(2:en,1:10);
y =  nanmean(data); 
x = 1:size(data,2); 
e = nanstd(data); 
errorbar(x, y,e,'r.','markersize',20');
set(gca,'xtick',[1:3:25],'xticklabel',[6:3:35],'xlim',[0 11]);
set(gca,'ylim',[1000 1030]);
title('Density (Spring ebb)','fontsize',14);
ylabel('Mean T','fontsize',14);

%% density 비교
st = 6; en = 30;
st_en = 22;
figure()
plot(nanmean(neap_ebb_den(2:end,1:st_en)),'b.-','markersize',20);
hold on;
plot(nanmean(neap_flood_den(2:end,1:st_en)),'b.-','markersize',20);
plot(nanmean(spring_ebb_den(2:end,1:st_en)),'r.-','markersize',20);
plot(nanmean(spring_flood_den(2:end,1:st_en)),'r.-','markersize',20);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
xlim([0 st_en]);
ylabel('Density','fontsize',14);
xlabel('Stiation','fontsize',14);

