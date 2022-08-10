clc;clear all;
close all;

%--- 정선 관측 자료와 해당 위치 위성자료의 위치정보과 트렌드 값 --------------
% load E:\11.사업\장기생태_3단계\Data\위성자료와수과원자료비교\NOAA_m_NFRDI\SST_trend.mat
% NIFS Serial Oceanographic observation

load E:\11.사업\장기생태_3단계\1차년도\draw_figure\nfrdi_trends.txt
nfrdi_trend=nfrdi_trends;
data = nfrdi_trend(:,4);
lon = nfrdi_trend(:,3);
lat = nfrdi_trend(:,2);
%--- jet color ---------------------------
% cmap = [255 85 0; 255 170 0;255 255 0;170 255 85;85 255 170; 0 255 255; ...
%         0 170 255; 0 85 255;0 0 255];
cmap = [0 0 255; 0 170 255; 255 255 0; 255 170 0; 255 85 0];
cmap = (cmap)/255;
%--- geoshow 용 coastal line data load ------------------------------------
%{
dum2=load('fine_coast2_closed_SS.dat');
coa_lon=dum2(:,1); coa_lat=dum2(:,2);
lat_lim=[33.125 35.125]; lon_lim=[126.125 129.125];

figure('Position', [100, 100, 900, 300]); 
subaxis(1,2,1,'SpacingVert',0,'MR',0); 
line_gap=0.5;
axesm('MapProjection','mercator','MapParallels',[],...
 'MapLatLimit',lat_lim,'MapLonLimit',lon_lim,...
 'MLabelLocation',line_gap,'MLineLocation',line_gap,'FontSize', 15, ...
 'PLabelLocation',line_gap,'PLineLocation',line_gap,'FontSize', 15, ...
 'GColor',.5*[1 1 1],'GLinestyle','-',...
 'MLabelParallel','south'); 
framem;
gridm('off');
mlabel;
plabel; 

geoshow(coa_lat, coa_lon,'Color','black')
% ,'DisplayType','polygon','facecolor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
axis tight; 
hold on; 
ms = 20; % marker size
c_lim = 0.04;
data = NFRDI_t;
for i = 1:length(cmap)-1    
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > 0.04);
    geoshow(lat(a),lon(a),'DisplayType','point','MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker', 'o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
    geoshow(lat(a),lon(a),'DisplayType','point','MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'MarkerSize',ms,'Marker', 'o');
hold on; 
end
colormap(cmap)
cbh = colorbar;
% set(cbh,'ytick',[0:1/9:1],'YTicklabel',c_int);

subaxis(1,2,2,'SpacingVert',0,'MR',0); 
line_gap=0.5;
axesm('MapProjection','mercator','MapParallels',[],...
 'MapLatLimit',lat_lim,'MapLonLimit',lon_lim,...
 'MLabelLocation',line_gap,'MLineLocation',line_gap,'FontSize', 15, ...
 'PLabelLocation',line_gap,'PLineLocation',line_gap,'FontSize', 15, ...
 'GColor',.5*[1 1 1],'GLinestyle','-',...
 'MLabelParallel','south'); 
framem;
gridm('off');
mlabel;
plabel; 

geoshow(coa_lat, coa_lon,'Color','black')
% ,'DisplayType','polygon','facecolor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
axis tight; 
hold on; 

data = AVHRR_t;
for i = 1:length(cmap)-1
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
    geoshow(lat(a),lon(a),'DisplayType','point','MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'MarkerSize',ms,'Marker', 'o');
hold on; 
end
colormap(cmap)
cbh = colorbar;
set(cbh,'ytick',[0:1/9:1],'YTicklabel',c_int);

%}
%%
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
% m_plot(lon, lat,'k.','markersize',15);
ms = 18; % marker size
c_lim = 0.03;
c_min = -0.02;
for i = 1:length(cmap)-1    
    c_int = [c_min:(c_lim-(c_min))/5:c_lim];
    [a b] = find(data > 0.03);
    m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(5,:), ...
            'MarkerFaceColor',cmap(5,:),'MarkerSize',ms,'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
    m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
hold on; 
end

%--- 연안정지 관측 자료도 같이 나타내기 -----------------------------------
load E:\11.사업\장기생태_3단계\1차년도\draw_figure\data_from_yukyung\all_trend.txt
data = all_trend(:,3);
lon = all_trend(:,2);
lat = all_trend(:,1);

for i = 1:length(cmap)-1    
    c_int = [c_min:(c_lim-(c_min))/5:c_lim];
    [a b] = find(data > 0.03);
    m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(5,:), ...
            'MarkerFaceColor',cmap(5,:),'MarkerSize',ms,'Marker','s');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
    m_plot(lon(a),lat(a),'Marker', 's','MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
hold on; 
end





colormap(cmap);
cbh = colorbar;
set(cbh,'ytick',[0:1/5:1],'YTicklabel',c_int);
title('SST trend (1984-2014)','fontsize',14);

%--- save figure --------------------------------------
%     outpath = '.\figure\';
%     out=[outpath,'NFRDI_KOHA_yearly_trend'];
%     set(gcf,'renderer','painter');
%     set(gcf, 'PaperUnits', 'inches');
%     x_width=5;
%     y_width=3;
%     set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%     saveas(gcf,out,'tif');
%------------------------------------------------------


