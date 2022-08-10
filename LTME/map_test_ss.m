clc;clear all;close all;

dum2=load('fine_coast2_closed_SS.dat');
coa_lon=dum2(:,1); coa_lat=dum2(:,2);

lat_lim=[33.125 35.125]; lon_lim=[126.125 129.125];
% lat_lim=[31 39]; lon_lim=[124 133];
figure; 
line_gap=0.5;
axesm('MapProjection','mercator','MapParallels',[],...
 'MapLatLimit',lat_lim,'MapLonLimit',lon_lim,...
 'MLabelLocation',line_gap,'MLineLocation',line_gap,'FontSize', 15, ...
 'PLabelLocation',line_gap,'PLineLocation',line_gap,'FontSize', 15, ...
 'GColor',.5*[1 1 1],'GLinestyle','-',...
 'MLabelParallel','south'); framem;gridm;mlabel;plabel; 

geoshow(coa_lat, coa_lon)
% ,'DisplayType','polygon','facecolor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
axis tight; 
hold on; 
plotm(35,132,'k.','MarkerSize',8);

%--- 정선별 정점 정보 가져오기 ---------------------------------------------
n_info=load('E:\11.사업\장기생태_3단계\Data\국립수산과학원(KODC)_정선관측\position_silver.txt');
n_lat = n_info(:,2);
n_lon = n_info(:,3);
plotm(n_lat,n_lon,'r.','MarkerSize',8)

%%
%--- m_map 으로 coastal line load -----------------------------------------

figure1 = figure('Position', [100, 100, 900, 300]); 
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
m_plot(n_lat, n_lon,'r.','markersize',15);
