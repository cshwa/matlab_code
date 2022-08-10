close all; clear; clc;
cd D:\장기생태\Dynamic\KOEM
load KOEM_name_tag.mat
load KOEM_st_info_(decimal_deg).mat
lat_koem = [lat_1; lat_2; lat_3;];
lon_koem = [lon_1; lon_2; lon_3;];
grd_file='grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');
h=ncread(grd_file, 'h');
for i  = 1:size(name_tag,1); name_tag1{i,1}=name_tag{i,1}; end;
for i = 1:size(name_tag,1); tag_mask{i,1}=name_tag{i,2}; end
clearvars name_tag;
name_tag = name_tag1;

%plot only in the domain
figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if tag_mask{i} == 1
        text(lon_koem(i),lat_koem(i), char(name_tag{i}),'color','r'); %plot only 97
        plot(lon_koem(i),lat_koem(i),'.','color','r'); %plot only 97
    end
end
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

%sewer site
sewer_ind(1,:) =  [112,68];
sewer_ind(2,:) =  [99,29];
sewer_ind(3,:) =  [106,3];
sewer_ind(4,:) =  [121,49];
sewer_ind(5,:) =  [63,56];
sewer_ind(6,:) =  [84,5];
sewer_ind(7,:) =  [56,24];
sewer_ind(8,:) =  [23,28];
sewer_ind(9,:) =  [126,61];

figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:size(sewer_ind,1)
plot(x(sewer_ind(i,2),sewer_ind(i,1)),y(sewer_ind(i,2),sewer_ind(i,1)),'.','color','c'); 
text(x(sewer_ind(i,2),sewer_ind(i,1)),y(sewer_ind(i,2),sewer_ind(i,1))-0.003,[num2str(i),'st'],'color','c'); 
end
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

%% TOPO
figure; h_plt=pcolor(x,y,h.*(mask./mask)); hh=colorbar; shading flat; hold on
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
for i=1:size(sewer_ind,1)
plot(x(sewer_ind(i,2),sewer_ind(i,1)),y(sewer_ind(i,2),sewer_ind(i,1)),'.','color','c'); 
text(x(sewer_ind(i,2),sewer_ind(i,1)),y(sewer_ind(i,2),sewer_ind(i,1))-0.003,[num2str(i),'st'],'color','c'); 
end

grd_file='D:\장기생태\Dynamic\06_river\하수종말처리장\grid_gy_v11_s_9_pollution.nc';
x_tak=ncread(grd_file, 'lon_rho');
y_tak=ncread(grd_file, 'lat_rho');
mask_tak=ncread(grd_file, 'mask_rho');
h_tak=ncread(grd_file, 'h');

figure; h_plt=pcolor(x_tak,y_tak,h_tak.*(mask_tak./mask_tak)); hh=colorbar; shading flat; hold on
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
for i=1:size(sewer_ind,1)
plot(x(sewer_ind(i,2),sewer_ind(i,1)),y(sewer_ind(i,2),sewer_ind(i,1)),'.','color','c'); 
text(x(sewer_ind(i,2),sewer_ind(i,1)),y(sewer_ind(i,2),sewer_ind(i,1))-0.003,[num2str(i),'st'],'color','c'); 
end

% diff
mask_re =mask_tak - mask;
figure; h_plt=pcolor(x,y,h.*(mask_re./mask_re)); hh=colorbar; shading flat; hold on
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;
% plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
for i=1:size(sewer_ind,1)
plot(x(sewer_ind(i,2),sewer_ind(i,1)),y(sewer_ind(i,2),sewer_ind(i,1)),'.','color','c'); 
text(x(sewer_ind(i,2),sewer_ind(i,1)),y(sewer_ind(i,2),sewer_ind(i,1))-0.003,[num2str(i),'st'],'color','c'); 
end

% optional
sp_gy = [4,22,28,29,30,32,33,34,35];
for i=1:length(sp_gy)
    i_plt=sp_gy(i);
%         text(lon_koem(i_plt),lat_koem(i_plt), char(name_tag{i_plt}),'color','r'); %plot only 97
        st_plt=plot(lon_koem(i_plt),lat_koem(i_plt),'+','color','r'); %plot only 97
end
set(get(get(h_plt,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(st_plt(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('시계열 자료 생성 위치')


figure; pcolor(x,y,h.*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;
plot_google_map('MapType','satellite')  % overlay google map
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])