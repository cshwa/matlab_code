close all; clear; clc; 
lon=ncread('grid_sumjin_v1970_fix_3m.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m.nc','h');

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;

lon2=ncread('grid_sumjin_v1970_fix_3m_v2.nc','lon_rho');
lat2=ncread('grid_sumjin_v1970_fix_3m_v2.nc','lat_rho');
mask2=ncread('grid_sumjin_v1970_fix_3m_v2.nc','mask_rho');
h2=ncread('grid_sumjin_v1970_fix_3m_v2.nc','h');

figure; pcolor(lon2,lat2,h2.*(mask2./mask2)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;

lons=ncread('grid_gy_v11_s_0613_depth.nc','lon_rho');
lats=ncread('grid_gy_v11_s_0613_depth.nc','lat_rho');
masks=ncread('grid_gy_v11_s_0613_depth.nc','mask_rho');
hs=ncread('grid_gy_v11_s_0613_depth.nc','h');

figure; pcolor(lons,lats,hs.*(masks./masks)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;


lon3=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lon_rho');
lat3=ncread('grid_sumjin_v1970_fix_3m_v3.nc','lat_rho');
mask3=ncread('grid_sumjin_v1970_fix_3m_v3.nc','mask_rho');
h3=ncread('grid_sumjin_v1970_fix_3m_v3.nc','h');

figure; pcolor(lon3,lat3,h3.*(mask3./mask3)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;




%Nolyang strait
h1 = h; 
h1(find(lon  >= 127.85 & lon  <= 127.95 & lat >= 34.94 & lat <= 34.96))=NaN
figure; pcolor(lon,lat,h1); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;

lon3=ncread('grid_sumjin_v1970_fix_3m.nc','lon_rho');
lat3=ncread('grid_sumjin_v1970_fix_3m.nc','lat_rho');
mask3=ncread('grid_sumjin_v1970_fix_3m.nc','mask_rho');
h3=ncread('grid_sumjin_v1970_fix_3m.nc','h');

figure; pcolor(lon3,lat3,h3.*(mask3./mask3)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
caxis([0 50])
grid on;


clear; clc; 
load merge_1970_depth_gjb_fix.mat
figure; pcolor(X,Y,h_merge_1); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading interp;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;


close all; clear; clc;   % -v2
lon=ncread('grid_sumjin_v1970_fix.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix.nc','h');

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading flat;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;

figure; pcolor(lon,lat,h.*(mask./mask)); hh=colorbar; shading interp;
set(gca,'fontsize',15,'fontweight','bold'); ylabel(hh,'depth (m)');
grid on;


close all; clear; clc;   % -v3
lon=ncread('grid_sumjin_v1970_fix_3m.nc','lon_rho');
lat=ncread('grid_sumjin_v1970_fix_3m.nc','lat_rho');
mask=ncread('grid_sumjin_v1970_fix_3m.nc','mask_rho');
h=ncread('grid_sumjin_v1970_fix_3m.nc','h');
load KOEM_st_info_(decimal_deg).mat
lat_koem = [lat_1; lat_2; lat_3;];
lon_koem = [lon_1; lon_2; lon_3;];

% station name tag
for i = 1:3
name_tag_1{i} = ['������' num2str(i,'%02d')] 
end

name_tag_2{1} = ['������' num2str(1,'%02d')] 

name_tag_3{1} = ['��õ��' num2str(1,'%02d')] 

for i = 1:5
name_tag_4{i} = ['������' num2str(i,'%02d')] 
end

for i = 1:25
name_tag_5{i} = ['������' num2str(i,'%02d')] 
end

for i = 1:2
name_tag_6{i} = ['���ָ�' num2str(i,'%02d')] 
end

for i = 1:28
name_tag_7{i} = ['����' num2str(i,'%02d')] 
end

name_tag = name_tag_1'; 
name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
size_tag = length(name_tag);

for i = 1:length(name_tag_4)
name_tag{size_tag+i} = name_tag_4{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_5)
name_tag{size_tag+i} = name_tag_5{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_6)
name_tag{size_tag+i} = name_tag_6{i}; 
end
size_tag = length(name_tag);

for i = 1:length(name_tag_7)
name_tag{size_tag+i} = name_tag_7{i};
end

% plot masking
figure; hold on; pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; plot(lon_koem,lat_koem,'.','color','k')
for i=1:length(name_tag)
text(lon_koem(i),lat_koem(i), char(name_tag{i}));
end

%plot only in the domain
figure; hold on; pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; plot(lon_koem,lat_koem,'.','color','k')
for i=1:length(name_tag)
text(lon_koem(i),lat_koem(i), char(name_tag{i}));
end
xlim([min(min(lon)) max(max(lon))])
ylim([min(min(lat)) max(max(lat))])




