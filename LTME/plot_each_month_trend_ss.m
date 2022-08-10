clc;clear all;close all;

%load file name
% name = dir('avhrr*');
ss = 1;
for mon = 1:1:12
mon=[num2char(mon,2)];

tn = 2016-1982+1;
% tn = 3;
for yr = 1:tn
        filepath1 = 'monthly\';   %monthly folder
        name11 = ['avhrr_monthly' num2str(yr+1981) '_' mon '.nc'];
        file1 = strcat(filepath1, name11);
        nc = netcdf(file1); 
        temp2(yr,:,:) = nc{'temp'}(:);
      
        lon = nc{'long'}(:);
        lat = nc{'lat'}(:);
end
lim = [126 129 33.4 35];
[a b]=find(abs((lon(1,:)-lim(1)))==0.1250);
min_lon = min(b);
[a b]=find(abs((lon(1,:)-lim(2)))==0.1250);
max_lon = max(b);
[min_lat b]=find(abs((lat(:,1)-lim(3)))<=0.0250);
[a b]=find(abs((lat(:,1)-lim(4)))<=0.1250);
max_lat = max(a);

temp2 = temp2(:,min_lat:max_lat,min_lon:max_lon);

%trend data;
x = [1:tn]';
% trend_data = zeros(720,1440)+NaN;
% trend_idx = find(isnan(squeeze(temp2(1,:)))==0);
%--- edited by silver 원하는 영역만 계산하기 ---------
nu_lon = max_lon-min_lon+1;
nu_lat = max_lat-min_lat+1;
trend_data = zeros(nu_lat,nu_lon)+NaN;
trend_idx = find(isnan(squeeze(temp2(1,:)))==0);
%---------------------------------------------------
for i=trend_idx
    data = temp2(:,i);
    [a,b] = polyfit(x,data,1);
    trend_data(i) = a(1);
end

T_trend_data(ss,1:8,1:14) = trend_data;

%%
figure1 = figure;
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126 129],'lat',[33.4 35]);
m_gshhs_l('color','k');
%m_gshhs_l('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
% m_pcolor(lon,lat,trend_data);
m_pcolor(lon(min_lat:max_lat,min_lon:max_lon),lat(min_lat:max_lat,min_lon:max_lon),trend_data);

shading flat;
set(gca,'FontSize',14);
tt = strcat(mon, ' SST trend(1982-2016)');
title(tt,'fontsize',20,'fontweight','bold');  
colorbar;colormap jet;
caxis([-0,0.05]);
% saveas(gcf,strcat(tt,'.png'),'png');

set(gcf,'renderer','painter');
set(gcf, 'PaperUnits', 'inches');
x_width=5;
y_width=3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,strcat(tt),'tif')

ss = ss+1;
clear temp2
close all;
end
