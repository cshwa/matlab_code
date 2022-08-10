close all;clear all;clc


% 장기간의 yearly-mean SST trend 자료가 있는 경우 그림 그리기

%--- color map ------------------------------------------------------------
cmap = [255 85 0; 255 170 0;255 255 0;170 255 85;85 255 170; 0 255 255; ...
        0 170 255; 0 85 255;0 0 255];
cmap = flipud(cmap)/255;

%--- 사용할 자료 load -----------------------------------------------------
% load E:\11.사업\장기생태_3단계\draw_figure\data_from_yukyung\all_trend.txt
% data = all_trend(:,3);
%--- 아래 자료는 climate가 제외된 30년간의 daily mean 값이 있는 것이므로 trend
% 구하는 처리가 필요함
% load E:\11.사업\장기생태_3단계\draw_figure\data_from_yukyung\daily_trend_tongyeong.mat
% temp = daily_trend_tongyeong;
% y = nanmean(temp);
% x = 1:length(y);
% scatter(x,y);
% h = lsline;
% p1 = polyfit(get(h,'xdata'),get(h,'ydata'),1)
clearvars lon lat
load('kodc_203_2001_2010_just_mean_yearly', 'coeff_*','lon_*','lat_*');
lon = lon_203; lat=lat_203; 
coef_s = coeff_salt; coef_t = coeff_temp; coef_n = coeff_no3; 

clearvars coeff_salt coeff_temp coeff_no3 lon_* lat_*
load('kodc_204_2001_2010_just_mean_yearly', 'coeff_*','lon_*','lat_*');
lon = cat(2,lon,lon_204); lat=cat(2,lat,lat_204);
coef_s = cat(1,coef_s,coeff_salt); coef_t =cat(1,coef_t, coeff_temp); coef_n = cat(1,coef_n,coeff_no3); 

clearvars coeff_salt coeff_temp coeff_no3 lon_* lat_*
load('kodc_205_2001_2010_just_mean_yearly', 'coeff_*','lon_*','lat_*');
lon = cat(2,lon,lon_205); lat=cat(2,lat,lat_205);
coef_s = cat(1,coef_s,coeff_salt); coef_t =cat(1,coef_t, coeff_temp); coef_n = cat(1,coef_n,coeff_no3);

clearvars coeff_salt coeff_temp coeff_no3 lon_* lat_*
load('kodc_206_2001_2010_just_mean_yearly', 'coeff_*','lon_*','lat_*');
lon = cat(2,lon,lon_206); lat=cat(2,lat,lat_206);
coef_s = cat(1,coef_s,coeff_salt); coef_t =cat(1,coef_t, coeff_temp); coef_n = cat(1,coef_n,coeff_no3);

clearvars coeff_salt coeff_temp coeff_no3 lon_* lat_*
load('kodc_400_2001_2010_just_mean_yearly', 'coeff_*','lon_*','lat_*');
lon = cat(2,lon,lon_400); lat=cat(2,lat,lat_400);
coef_s = cat(1,coef_s,coeff_salt); coef_t =cat(1,coef_t, coeff_temp); coef_n = cat(1,coef_n,coeff_no3);

return

%temp : - 0.04 ~ 0.04 
max(coef_t(:,1))
min(coef_t(:,1))

%salt : -0.015 ~ 0.015
max(coef_s(:,1))
min(coef_s(:,1))

%no3 : 0~1
max(coef_n(:,1))
min(coef_n(:,1))



%%temp
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size
c_lim = 0.04;
clearvars data
data = coef_t(:,1);
for i = 1:length(cmap)-1    
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > 0.04);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','^oC/year')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SST trend (1980-2019)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_trend_temp'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------

%% salt
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size
c_lim = 0.015;
clearvars data
data = coef_s(:,1);
for i = 1:length(cmap)-1    
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > 0.015);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','psu/year')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SSS trend (1980-2019)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_trend_salt'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------


%% no3
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size
c_lim = 1;
clearvars data
data = coef_n(:,1);
for i = 1:length(cmap)-1
    c_int = [0:c_lim/8:c_lim];
    [a b] = find(data > 1);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','ug/L /year')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SSNO3 trend (1988-2019)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_trend_no3'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------

