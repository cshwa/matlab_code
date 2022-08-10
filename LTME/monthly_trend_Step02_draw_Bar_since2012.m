
% with Step04_pop_climitology.m
map = [0.2 0.0 0.6
    0.2549    0.4118    0.8824
    0.1333    0.5451    0.1333
    1.0000    0.8431         0
    1.0000    0.3804    0.0118 
    0.2510    0.8784    0.8157];


figure('Position', [100, 100, 900, 700])
subplot(311)
data = GD(1:9,:);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';     
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:9],'xticklabel',fname(1:9),'ylim',[-1 5]);  %since 2012
title('Since 2012','fontsize',14);

subplot(312)
data = GD(10:18,:);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';   
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:9],'xticklabel',fname(10:18),'ylim',[-1 1]); %since 2012

subplot(313)
data = GD(18:26,:);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';   
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:9],'xticklabel',fname(18:26),'ylim',[-1 1]); %since 2012

% subplot(414)
% data = GD(18:end,:);
% h = bar(data);
% colormap(map)
% grid on
% l = cell(1,5);
% % l = fname;
% l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';   
% legend(h,l,'location','eastoutside');
% set(gca,'xtick',[1:10],'xticklabel',fname(18:27),'ylim',[-1 1], ...
%     'xlim',[0 11]);  %since 2012

%--- save figure --------------------
outpath = '.\figure\';
% out_vari = fname(i);
out=[outpath,'bar_since2012_',num2str(p_dep)];
set(gcf,'renderer','painter');
set(gcf, 'PaperUnits', 'inches');
x_width=9;
y_width=7;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,out,'tif');
%------------------------------------