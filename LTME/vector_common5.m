set(gca,...
    'XTickLabel',{'128.0^oE','128.5^oE','129.0^oE','129.5^oE'},...
    'XTick',[128.0 128.5 129.0 129.5]);
set(gca,...
    'YTickLabel',{'-160m','-120m','-80m','-40m'},...
        'YTick',[-160 -120 -80 -40]);
set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
set(gca,'fontsize',12);
set(gca, 'color', [0.8,0.8,0.8]);
fig=gcf;
fig.InvertHardcopy='off';
saveas(gcf,filename,'tiff');