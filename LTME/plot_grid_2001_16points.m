close all; clear; clc;
cd D:\장기생태\Dynamic\KOEM
lat_st = [34.8625,34.85139,34.88389,34.90222,34.92083,34.83194,34.86944,...
    34.90278,34.89167,34.9025,34.86361,34.91,34.94,34.73556,34.76278,34.76444,...
    34.69056, 34.62];
lon_st = [127.7403,127.6797,127.6506,127.6822,127.8233,127.8011,127.7889,...
    127.7989,127.7597,127.7239,127.7103,127.7,127.8264,127.7661,127.7614,...
    127.8053,127.845,127.8069];

grd_file='grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');
h=ncread(grd_file, 'h');

clearvars name_tag
name_tag{1}={'광양항01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['광양만',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['여수연안',num2str(k,'%02d')]};
end

%plot only in the domain
figure; hold on; pcolor(x,y,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
plot(x(87,:),y(87,:),'r--','linew',2); %vertical diagonal
% plot(x(:,76),y(:,76),'--','color',[.5 .5 .5],'linew',2); %vertical meridional
plot(x(:,74),y(:,74),'--','color','b','linew',2); %vertical meridional
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 15
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','r','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','r'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','r','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','r'); %plot only 97
    end
end
% ylim([34.723 35])
% xlim([min(min(x)) 127.9])
% xlim([min(min(x)) max(max(x))])
% ylim([min(min(y)) max(max(y))])
ylim([min(min(y)) 35])
xlim([127.5753 127.85])
plot_google_map('MapType','terrain','Scale',2,'Resize',3)   % overlay google map
% ylim([34.723 35])
% xlim([min(min(x)) 127.9])
% xlim([min(min(x)) max(max(x))])
% ylim([min(min(y)) max(max(y))])
ylim([min(min(y)) 35])
xlim([127.5753 127.85])



