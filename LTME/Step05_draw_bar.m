
% load GD

% with Step04_pop_climitology.m
map = [0.2 0.0 0.6
    0.2549    0.4118    0.8824
    1.0000    0.8431         0
    1.0000    0.3804    0.0118 ];
%     0.2510    0.8784    0.8157];


figure('Position', [100, 100, 1100, 700])
subplot(311)
data = GD(1:9,[1 2 4 5 ]);
h = bar(data);
colormap(map)
grid on
set(gca,'xtick',[1:9],'xticklabel',fname(1:9),'ylim',[-0.05 0.15]);
ylabel('C^{\circ}/year','fontsize',14);
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='8'; l{4}='10'; l{5}='12'; l{6}='12';    
% l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';    
legend(h,l,'location','eastoutside');
% legend(h,l);


subplot(312)
data = GD(10:18,[1 2 4 5 ]);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='8'; l{4}='10'; l{5}='12'; l{6}='12'; 
% l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';    
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:9],'xticklabel',fname(10:18),'ylim',[-0.05 0.15]);
ylabel('C^{\circ}/year','fontsize',14);

subplot(313)
data = GD(19:27,[1 2 4 5 ]);
h = bar(data);
colormap(map)
grid on
l = cell(1,5);
% l = fname;
l{1}='2'; l{2}='4'; l{3}='8'; l{4}='10'; l{5}='12'; l{6}='12';   
% l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';    
legend(h,l,'location','eastoutside');
set(gca,'xtick',[1:9],'xticklabel',fname(19:27),'ylim',[-0.05 0.15]);
ylabel('C^{\circ}/year','fontsize',14);

% subplot(414)
% data = GD(18:end,[1 2 4 5]);
% h = bar(data);
% colormap(map)
% grid on
% l = cell(1,5);
% % l = fname;
% % l{1}='2'; l{2}='4'; l{3}='6'; l{4}='8'; l{5}='10'; l{6}='12';   
% % legend(h,l,'location','eastoutside');
% set(gca,'xtick',[1:10],'xticklabel',fname(18:27),'ylim',[-0.2 0.2], ...
%     'xlim',[0 11]);

%--- save figure --------------------
outpath = '.\figure_season\';
% out_vari = fname(i);
out=[outpath,'bar_since1984_',num2str(p_dep)];
set(gcf,'renderer','painter');
set(gcf, 'PaperUnits', 'inches');
x_width=11;
y_width=5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,out,'tif');
%------------------------------------