pcolor(lon5,lat5,(vort5mean(:,:))');
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis([-1.1e-5,1.1e-5]);
c=colorbar;
c.Label.String= 'Vorticity';
c.Label.FontSize= 12;
xlabel('Longitude (^oE)','fontsize',27);
ylabel('Latitude (^oN)','fontsize',27);
set(gcf,'PaperPosition',[0 0 20 15]);
axis equal;
set(gca,'fontsize',15);
colormap(bwrmap);
qv=quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
squeeze(u2mean(1:3:120,1:3:110,1,1))'*2,squeeze(v2mean(1:3:120,1:3:110,1,1))'*2,'k','AutoScale','off');
set(qv,'Color',[0 0 0]);  
hold off;
tex=text(135.9,35.2,'0.2m/s');
set(tex,'fontsize',10);
fig=gcf;
fig.InvertHardcopy='off';
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
saveas(gcf,filename,'tiff');