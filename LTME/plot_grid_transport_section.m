close all; clear; clc;
cd D:\장기생태\Dynamic\KOEM

grd_file='grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');
h=ncread(grd_file, 'h');

cd D:\장기생태\Dynamic
load picked_transport_cal_points

%plot only in the domain
figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;
out_name = ['all_section_transport_RCP85'];

figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos) 
hold on; pcolor(x,y,(mask./mask));  hold on; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
plot([x(87,117) x(83,98)], [y(87,117) y(83,98)],'g','linew',2); %vertical diagonal
plot([x(72,113) x(83,121)], [y(72,113) y(83,121)],'b','linew',2); %vertical diagonal
plot([x(31,102) x(42,104)], [y(31,102) y(42,104)],'k','linew',2); %vertical diagonal
plot([x(78,90) x(76,66)], [y(78,90) y(76,66)],'c','linew',2); %vertical diagonal
plot([x(77,66) x(103,68)], [y(77,66) y(103,68)],'r','linew',2); %vertical diagonal
plot([x(131,122) x(132,105)], [y(131,122) y(132,105)],'m','linew',2); %vertical diagonal
ylim([34.835 35])
xlim([min(min(x)) 127.89])
print(gcf,[out_name,'_grid','.png'],'-dpng','-r200');

%%
plot(x(87,117), y(87,117),'g.','linew',4); %vertical diagonal
plot(x(83,98), y(83,98),'g.','linew',4); %vertical diagonal

plot(x(72,113), y(72,113),'b.','linew',4); %vertical diagonal
plot(x(83,121), y(83,121),'b.','linew',4); %vertical diagonal

plot(x(31,102), y(31,102),'k.','linew',4); %vertical diagonal
plot(x(42,102), y(42,102),'k.','linew',4); %vertical diagonal

plot(x(78,90), y(78,90),'c.','linew',4); %vertical diagonal
plot(x(76,66), y(76,66),'c.','linew',4); %vertical diagonal

plot(x(77,66), y(77,66),'r.','linew',4); %vertical diagonal
plot(x(103,66), y(103,66),'r.','linew',4); %vertical diagonal

plot(x(131,122), y(131,122),'m.','linew',4); %vertical diagonal
plot(x(132,105), y(132,105),'m.','linew',4); %vertical diagonal



% plot(x(73,111), y(73,111),'m.','linew',2); %vertical diagonal
% plot(x(73,105), y(73,105),'m.','linew',2); %vertical diagonal

% plot(x(:,76),y(:,76),'--','color',[.5 .5 .5],'linew',2); %vertical meridional
plot(x(:,74),y(:,74),'--','color',[.5 .5 .5],'linew',2); %vertical meridional
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



