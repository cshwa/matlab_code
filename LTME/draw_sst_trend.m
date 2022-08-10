clc;clear all;
close all;

% ��Ⱓ�� yearly-mean SST trend �ڷᰡ �ִ� ��� �׸� �׸���

%--- color map ------------------------------------------------------------
cmap = [255 85 0; 255 170 0;255 255 0;170 255 85;85 255 170; 0 255 255; ...
        0 170 255; 0 85 255;0 0 255];
cmap = flipud(cmap)/255;

%--- ����� �ڷ� load -----------------------------------------------------
load E:\11.���\������_3�ܰ�\draw_figure\data_from_yukyung\all_trend.txt
data = all_trend(:,3);
%--- �Ʒ� �ڷ�� climate�� ���ܵ� 30�Ⱓ�� daily mean ���� �ִ� ���̹Ƿ� trend
% ���ϴ� ó���� �ʿ���
% load E:\11.���\������_3�ܰ�\draw_figure\data_from_yukyung\daily_trend_tongyeong.mat
% temp = daily_trend_tongyeong;
% y = nanmean(temp);
% x = 1:length(y);
% scatter(x,y);
% h = lsline;
% p1 = polyfit(get(h,'xdata'),get(h,'ydata'),1)
%%
%--- m_map ���� coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 20; % marker size
c_lim = 0.04;
for i = 1:length(cmap)-1    
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > 0.04);
    m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
    m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
%--- ���� �����ϱ� -------------------------------------
title('NFRDI SST trend (1984-2014)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output �̸� �����ϱ� -------------------------
    out=[outpath,'NFRDI_yearly_trend'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%     saveas(gcf,out,'tif');
%------------------------------------------------------


