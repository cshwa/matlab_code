close all; clear; clc; 
lon=ncread('grid_sumjin_v1970.nc','lon_rho');
lat=ncread('grid_sumjin_v1970.nc','lat_rho');
mask=ncread('grid_sumjin_v1970.nc','mask_rho');
h=ncread('grid_sumjin_v1970.nc','h');

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
name_tag_1{i} = ['여수항' num2str(i,'%02d')] 
end

name_tag_2{1} = ['광양항' num2str(1,'%02d')] 

name_tag_3{1} = ['삼천포' num2str(1,'%02d')] 

for i = 1:5
name_tag_4{i} = ['가막만' num2str(i,'%02d')] 
end

for i = 1:25
name_tag_5{i} = ['섬진강' num2str(i,'%02d')] 
end

for i = 1:2
name_tag_6{i} = ['진주만' num2str(i,'%02d')] 
end

for i = 1:28
name_tag_7{i} = ['연안' num2str(i,'%02d')] 
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




