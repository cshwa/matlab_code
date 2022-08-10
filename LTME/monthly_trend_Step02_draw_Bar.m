
% with Step04_pop_climitology.m
map = [0.2 0.0 0.6
    0.2549    0.4118    0.8824
    0.1333    0.5451    0.1333
    1.0000    0.8431         0
    1.0000    0.3804    0.0118 
    0.2510    0.8784    0.8157];


figure('Position', [100, 100, 900, 700])
subplot(411)
data = GD(1:7,:);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';     
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:7],'xticklabel',fname(1:7),'ylim',[-0.2 0.2]);

subplot(412)
data = GD(8:13,:);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';   
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:6],'xticklabel',fname(8:13),'ylim',[-0.2 0.2]);

subplot(413)
data = GD(14:17,:);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';   
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:4],'xticklabel',fname(14:17),'ylim',[-0.2 0.2]);

subplot(414)
data = GD(18:end,:);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';   
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:10],'xticklabel',fname(18:27),'ylim',[-0.2 0.2], ...
    'xlim',[0 11]);

%--- save figure --------------------
% outpath = '.\figure\';
% % out_vari = fname(i);
% out=[outpath,'bar_since1984_',num2str(p_dep)];
% set(gcf,'renderer','painter');
% set(gcf, 'PaperUnits', 'inches');
% x_width=9;
% y_width=7;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% saveas(gcf,out,'tif');
%------------------------------------