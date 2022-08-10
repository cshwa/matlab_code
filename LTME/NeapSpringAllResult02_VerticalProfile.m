%% 전체 비교 vs 상류비교
ctd_en = 36;
st = 06; en = 30;
figure()
% neap ebb
subplot(341)
contourf(flipud(neap_ebb_temp(2:ctd_en,:))); shading flat;
% colorbar; 
caxis([7 10]);
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
title('Neap ebb','fontsize',14);
ylabel('Depth (m)','fontsize',14);

subplot(345)
contourf(flipud(neap_ebb_sal(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([0 30]);

subplot(349)
contourf(flipud(neap_ebb_den(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([1000 1027]);
xlabel('Station','fontsize',14);

% neap flood
subplot(342)
contourf(flipud(neap_flood_temp(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([6.5 9.5]);
title('Neap flood','fontsize',14);

subplot(346)
contourf(flipud(neap_flood_sal(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([0 30]);

% spring ebb
subplot(3,4,10)
contourf(flipud(neap_flood_den(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([1000 1027]);

subplot(343)
contourf(flipud(spring_ebb_temp(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([6.5 9.5]);
title('Spring ebb','fontsize',14);

subplot(347)
contourf(flipud(spring_ebb_sal(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([0 30]);

subplot(3,4,11)
contourf(flipud(spring_ebb_den(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([1000 1027]);

subplot(344)
contourf(flipud(spring_flood_temp(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([6.5 9.5]);
title('Neap flood','fontsize',14);

subplot(348)
contourf(flipud(spring_flood_sal(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([0 30]);

subplot(3,4,12)
contourf(flipud(spring_flood_den(2:ctd_en,:))); shading flat;
set(gca,'ytick',[0:5:35],'yticklabel',[35:-5:0],'ylim',[0 35]);
set(gca,'xtick',[1-st:5:en],'xticklabel',[0:5:en],'xlim',[1-st en-st+1]);
caxis([1000 1027]);