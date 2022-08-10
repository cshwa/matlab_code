clc;clear all;close all;

%load file name
name = dir('avhrr*');

%load nc file
name2 = name(1).name;
nc = netcdf(name2);
lon = nc{'long'}(:);
lat = nc{'lat'}(:);

temp2 = zeros(35,720,1440);
for i = 1:35
    name2 = name(i).name;
    nc = netcdf(name2);
    temp2(i,:,:) = nc{'temp'}(:);
end
%trend data;
x = [1:35]';
trend_data = zeros(720,1440)+NaN;
trend_idx = find(isnan(squeeze(temp2(1,:)))==0);

for i=trend_idx
    data = temp2(:,i);
    [a,b] = polyfit(x,data,1);
    trend_data(i) = a(1);
end


figure1 = figure;
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 2 1]*6);

m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_gshhs_l('color','k');
%m_gshhs_l('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
m_pcolor(lon,lat,trend_data);shading flat;
set(gca,'FontSize',14);
tt = strcat('AVHRR(1982-2016)');
title(tt,'fontsize',20,'fontweight','bold');  
colorbar;colormap jet;
caxis([-0.05,0.05]);
saveas(gcf,strcat(tt,'.png'),'png');
close all;
