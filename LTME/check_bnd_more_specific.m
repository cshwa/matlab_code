close all; clear; clc;
cd D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary
load interp_bnd_fix_2001.mat
%south
for i = 1:12
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_temp_s_t(:,:,i)));h = colorbar;
set(get(h,'title'),'string','Temp.(^oC)'); shading flat
ylabel('depth (m)','fontsize',13);
xlabel('longitude (deg)','fontsize',13);
title(['south boundary (interped from yjtaks) = ',num2str(i),'month'],'fontsize',13);
saveas(gcf,['south_bnd_temp_',num2str(i),'mon.png'])
end

close all; clear; clc;
cd D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary
load interp_bnd_fix_2002.mat
%south
for i = 1:12
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_temp_s_t(:,:,i)));h = colorbar;
set(get(h,'title'),'string','Temp.(^oC)'); shading flat
ylabel('depth (m)','fontsize',13);
xlabel('longitude (deg)','fontsize',13);
title(['south boundary (interped from yjtaks) = ',num2str(i),'month'],'fontsize',13);
% saveas(gcf,['south_bnd_temp_',num2str(i),'mon.png'])
end



for i = 1:12
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_no3_s_t(:,:,i)));h = colorbar;
set(get(h,'title'),'string','NO3 (mmol N /m^3)'); shading flat
ylabel('depth (m)','fontsize',13);
xlabel('longitude (deg)','fontsize',13);
title(['south boundary (interped from yjtaks)= ',num2str(i),'mon'],'fontsize',11);
saveas(gcf,['south_bnd_NO3_',num2str(i),'mon.png'])
end

% nutri. east
no3_s=d3_no3_s_t;
nh4_s=d3_nh4_s_t;
no3_e=d3_no3_e_t;
nh4_e=d3_nh4_e_t;
chl_e=d3_chl_e_t;
chl_s=d3_chl_s_t; 


cd D:\장기생태\Dynamic\07_boundary_ts\pic\yjtak
for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(nh4_e(:,:,i))); shading flat
    caxis([0 1.5]); 
    h = colorbar; set(get(h,'title'),'string','NH4(uM)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from monthly mean'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_nh4_e'),'-dpng')
end

for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(no3_e(:,:,i)));colorbar; shading flat
    caxis([0 15])
    h = colorbar; set(get(h,'title'),'string','NO3(uM)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from monthly mean'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_no3_e'),'-dpng')
end


cd D:\장기생태\Dynamic\07_boundary_ts\pic\yjtak
for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(chl_e(:,:,i))); shading flat
    caxis([0 4]); 
    h = colorbar; set(get(h,'title'),'string','chl(mg/L)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from monthly mean'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_chl_e'),'-dpng')
end


% nutri. south

for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chl_s(:,:,i)));colorbar; shading flat
    caxis([0 4])
    h = colorbar; set(get(h,'title'),'string','chl(mg/L)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from monthly mean'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_chl_s'),'-dpng')
end


for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(nh4_s(:,:,i)));colorbar; shading flat
    caxis([0 1.5])
    h = colorbar; set(get(h,'title'),'string','NH4(uM)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from monthly mean'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_nh4_s'),'-dpng')
end


for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(no3_s(:,:,i)));colorbar; shading flat
    caxis([0 15])
    h = colorbar; set(get(h,'title'),'string','NO3(uM)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from monthly mean'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_no3_s'),'-dpng')
end


figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_salt_s_t(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_phy_s_t(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_chl_s_t(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_no3_s_t(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_nh4_s_t(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_zoo_s_t(:,:,1)));colorbar; shading flat
%east
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(d3_temp_e_t(:,:,12)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(d3_salt_e_t(:,:,12)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(d3_phy_e_t(:,:,12)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(d3_chl_e_t(:,:,12)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(d3_no3_e_t(:,:,12)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(d3_nh4_e_t(:,:,12)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(d3_zoo_e_t(:,:,12)));colorbar; shading flat


for i = 1:12
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(d3_no3_e_t(:,:,i)));h = colorbar;
set(get(h,'title'),'string','NO3 (mmol N /m^3)'); shading flat
ylabel('depth (m)','fontsize',13);
xlabel('longitude (deg)','fontsize',13);
title(['East boundary (interped from yjtaks)= ',num2str(i),'mon'],'fontsize',11);
saveas(gcf,['east_bnd_NO3_',num2str(i),'mon.png'])
end

% close all; clear; clc;
cd D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary
close all; clear; clc;
matf=load('interp_bnd_fix_2001.mat');
ref_lon_sm=matf.ref_lon_sm;
ref_dep_s=matf.ref_dep_s;
ref_lon_em=matf.ref_lon_em;
ref_dep_e=matf.ref_dep_e;
ref_lat_em=matf.ref_lat_em;

cd D:\장기생태\Dynamic\07_boundary_ts\Gwangyang_bry_ykang
mask=ncread('grid_gy_v11_s.nc','mask_rho');

mask_s=repmat(mask(:,1),1,20);  
mask_e=repmat(mask(end,:)',1,20);  

i=3
bnd_f = ['bndy_auto_NWP_1_10_test6_clim_',num2str(i),'regime.nc'];
clearvars temp_s temp_e salt_s salt_e chlo_s chlo_e
temp_s=ncread(bnd_f,'temp_south');
salt_s=ncread(bnd_f','salt_south');
chlo_s=ncread(bnd_f,'chlo_south');
no3_s=ncread(bnd_f,'NO3_south');
nh4_s=ncread(bnd_f','NH4_south');
po4_s=ncread(bnd_f','tPO4_south');
temp_e=ncread(bnd_f,'temp_east');
salt_e=ncread(bnd_f,'salt_east');
chlo_e=ncread(bnd_f,'chlo_east'); 
no3_e=ncread(bnd_f,'NO3_east');
nh4_e=ncread(bnd_f,'NH4_east');
chlo_e=ncread(bnd_f,'chlo_east'); 
po4_e=ncread(bnd_f,'tPO4_east'); 

%south
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(temp_s(:,:,1)).*(mask_s./mask_s));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,1)).*(mask_s./mask_s));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chlo_s(:,:,1)).*(mask_s./mask_s));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,1)).*(mask_s./mask_s));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(nh4_s(:,:,1)).*(mask_s./mask_s));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(po4_s(:,:,1)).*(mask_s./mask_s));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(no3_s(:,:,1)).*(mask_s./mask_s));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(temp_e(:,:,1)).*(mask_e./mask_e));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(salt_e(:,:,1)).*(mask_e./mask_e));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(chlo_e(:,:,1)).*(mask_e./mask_e));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(po4_e(:,:,1)).*(mask_e./mask_e));colorbar; shading flat


% close all; clear; clc;
close all; clear; clc;
cd D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary
matf=load('interp_bnd_fix_2001.mat');
ref_lon_sm=matf.ref_lon_sm;
ref_dep_s=matf.ref_dep_s;
ref_lon_em=matf.ref_lon_em;
ref_dep_e=matf.ref_dep_e;
ref_lat_em=matf.ref_lat_em;
% 
cd D:\장기생태\Dynamic\07_boundary_ts
bnd_f = 'bndy_2001_Gwangyang_P.nc'
clearvars temp_s temp_e salt_s salt_e chlo_s chlo_e
temp_s=ncread(bnd_f,'temp_south');
salt_s=ncread(bnd_f','salt_south');
chlo_s=ncread(bnd_f,'chlo_south');
no3_s=ncread(bnd_f,'NO3_south');
nh4_s=ncread(bnd_f','NH4_south');
po4_s=ncread(bnd_f','tPO4_south');
do_s=ncread(bnd_f','oxygen_south');
ldetp_s=ncread(bnd_f','LDeP_south');
sdetp_s=ncread(bnd_f','SDeP_south');
% ldetn_s=ncread(bnd_f','LDeN_south');
% sdetn_s=ncread(bnd_f','SDeN_south');

temp_e=ncread(bnd_f,'temp_east');
salt_e=ncread(bnd_f,'salt_east');
chlo_e=ncread(bnd_f,'chlo_east'); 
no3_e=ncread(bnd_f,'NO3_east');
nh4_e=ncread(bnd_f,'NH4_east');
chlo_e=ncread(bnd_f,'chlo_east'); 
po4_e=ncread(bnd_f','tPO4_east');
do_e=ncread(bnd_f','oxygen_east');
ldetp_e=ncread(bnd_f','LDeP_east');
sdetp_e=ncread(bnd_f','SDeP_east');
% ldetn_e=ncread(bnd_f','LDeN_east');
% sdetn_e=ncread(bnd_f','SDeN_east');

cd D:\장기생태\Dynamic\07_boundary_ts\pic

find(isnan(temp_s)==1)
find(isnan(salt_s)==1)
find(isnan(chlo_s)==1)
find(isnan(no3_s)==1)
find(isnan(nh4_s)==1)
find(isnan(po4_s)==1)
find(isnan(do_s)==1)
find(isnan(ldetp_s)==1)
find(isnan(sdetp_s)==1)
find(isnan(ldetn_s)==1)
find(isnan(sdetn_s)==1)

find(isnan(temp_e)==1)
find(isnan(salt_e)==1)
find(isnan(chlo_e)==1)
find(isnan(no3_e)==1)
find(isnan(nh4_e)==1)
find(isnan(chlo_e)==1)
find(isnan(po4_e)==1)
find(isnan(do_e)==1)
find(isnan(ldetp_e)==1)
find(isnan(sdetp_e)==1)
find(isnan(ldetn_e)==1)
find(isnan(sdetn_e)==1)

find(temp_s < 0)
find(salt_s< 0)
find(chlo_s<0)
find(no3_s<0)
find(nh4_s<0)
find(po4_s<0)
find(do_s<0)

find(temp_e <0)
find(salt_e<0)
find(chlo_e<0)
find(no3_e<0)
find(nh4_e<0)
find(chlo_e<0)
find(po4_e<0)
find(do_e<0)


%south
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(temp_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chlo_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(po4_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(nh4_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(no3_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(do_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(temp_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(salt_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(chlo_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(po4_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(nh4_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(no3_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(do_e(:,:,1)));colorbar; shading flat

cd D:\장기생태\Dynamic\07_boundary_ts\pic\originated_kodc

% east
for i = 1:12
    figure; pcolor(ref_lat_em,ref_dep_e,squeeze(po4_e(:,:,i)));colorbar; shading flat
    caxis([0 0.7])
end

for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(chlo_e(:,:,i))); shading flat
    caxis([0 4]); 
    h = colorbar; set(get(h,'title'),'string','Chl(mg/L)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from KODC'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_chl_e'),'-dpng')
end


for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(nh4_e(:,:,i))); shading flat
    caxis([0 1.5]); 
    h = colorbar; set(get(h,'title'),'string','NH4(uM)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from KODC'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_nh4_e'),'-dpng')
end

for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(no3_e(:,:,i)));colorbar; shading flat
    caxis([0 15])
    h = colorbar; set(get(h,'title'),'string','NO3(uM)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from KODC'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_no3_e'),'-dpng')
end

for i = 1:12
    figure; pcolor(ref_lat_em,ref_dep_e,squeeze(do_e(:,:,i)));colorbar; shading flat
%     caxis([0 200])
end


for i = 1:12
    figure; pcolor(ref_lat_em,ref_dep_e,squeeze(temp_e(:,:,i)));colorbar; shading flat
    caxis([0 30])
end

for i = 1:12
    figure; pcolor(ref_lat_em,ref_dep_e,squeeze(salt_e(:,:,i)));colorbar; shading flat
    caxis([30 34])
end

% south
for i = 1:12
    figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(po4_s(:,:,i)));colorbar; shading flat
%     caxis([0 0.7])
end

for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chlo_s(:,:,i)));colorbar; shading flat
    caxis([0 4])
    h = colorbar; set(get(h,'title'),'string','Chl(mg/L)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from KODC'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_chl_s'),'-dpng')
end

for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(nh4_s(:,:,i)));colorbar; shading flat
    caxis([0 1.5])
    h = colorbar; set(get(h,'title'),'string','NH4(uM)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from KODC'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_nh4_s'),'-dpng')
end


for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(no3_s(:,:,i)));colorbar; shading flat
    caxis([0 15])
    h = colorbar; set(get(h,'title'),'string','NO3(uM)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from KODC'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_no3_s'),'-dpng')
end


% figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(nh4_s(:,:,1)));colorbar; shading flat
% caxis([-inf inf])
% find(nh4_s < 0) 
% find(nh4_e < 0)


for i = 1:12
    figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(do_s(:,:,i)));colorbar; shading flat
%     caxis([0 200])
end


for i = 1:12
    figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(temp_s(:,:,i)));colorbar; shading flat
    caxis([0 30])
end

for i = 1:12
    figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,i)));colorbar; shading flat
    caxis([30 34])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clim_11_17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;
cd D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary
matf=load('interp_bnd_fix_2001.mat');
ref_lon_sm=matf.ref_lon_sm;
ref_dep_s=matf.ref_dep_s;
ref_lon_em=matf.ref_lon_em;
ref_dep_e=matf.ref_dep_e;
ref_lat_em=matf.ref_lat_em;
% 
cd D:\장기생태\Dynamic\07_boundary_ts\pic\yjtak_11to17
bnd_f = 'bndy_2001_Gwangyang_P_yjtak_11_17.nc'
clearvars temp_s temp_e salt_s salt_e chlo_s chlo_e
temp_s=ncread(bnd_f,'temp_south');
salt_s=ncread(bnd_f','salt_south');
chlo_s=ncread(bnd_f,'chlo_south');
no3_s=ncread(bnd_f,'NO3_south');
nh4_s=ncread(bnd_f','NH4_south');
po4_s=ncread(bnd_f','tPO4_south');
do_s=ncread(bnd_f','oxygen_south');
% ldetp_s=ncread(bnd_f','LDeP_south');
% sdetp_s=ncread(bnd_f','SDeP_south');
% ldetn_s=ncread(bnd_f','LDeN_south');
% sdetn_s=ncread(bnd_f','SDeN_south');

temp_e=ncread(bnd_f,'temp_east');
salt_e=ncread(bnd_f,'salt_east');
chlo_e=ncread(bnd_f,'chlo_east'); 
no3_e=ncread(bnd_f,'NO3_east');
nh4_e=ncread(bnd_f,'NH4_east');
chlo_e=ncread(bnd_f,'chlo_east'); 
po4_e=ncread(bnd_f','tPO4_east');
do_e=ncread(bnd_f','oxygen_east');
% ldetp_e=ncread(bnd_f','LDeP_east');
% sdetp_e=ncread(bnd_f','SDeP_east');
% ldetn_e=ncread(bnd_f','LDeN_east');
% sdetn_e=ncread(bnd_f','SDeN_east');

find(isnan(temp_s)==1)
find(isnan(salt_s)==1)
find(isnan(chlo_s)==1)
find(isnan(no3_s)==1)
find(isnan(nh4_s)==1)
find(isnan(po4_s)==1)
find(isnan(do_s)==1)
% find(isnan(ldetp_s)==1)
% find(isnan(sdetp_s)==1)
% find(isnan(ldetn_s)==1)
% find(isnan(sdetn_s)==1)

find(isnan(temp_e)==1)
find(isnan(salt_e)==1)
find(isnan(chlo_e)==1)
find(isnan(no3_e)==1)
find(isnan(nh4_e)==1)
find(isnan(chlo_e)==1)
find(isnan(po4_e)==1)
find(isnan(do_e)==1)
% find(isnan(ldetp_e)==1)
% find(isnan(sdetp_e)==1)
% find(isnan(ldetn_e)==1)
% find(isnan(sdetn_e)==1)

find(temp_s < 0)
find(salt_s< 0)
find(chlo_s<0)
find(no3_s<0)
find(nh4_s<0)
find(po4_s<0)
find(do_s<0)

find(temp_e <0)
find(salt_e<0)
find(chlo_e<0)
find(no3_e<0)
find(nh4_e<0)
find(chlo_e<0)
find(po4_e<0)
find(do_e<0)


% east
for i = 1:12
    figure; pcolor(ref_lat_em,ref_dep_e,squeeze(po4_e(:,:,i)));colorbar; shading flat
    caxis([0 0.7])
end

for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(chlo_e(:,:,i))); shading flat
    caxis([0 4]); 
    h = colorbar; set(get(h,'title'),'string','Chl(mg/L)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from 11~17 climate'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_chl_e'),'-dpng')
end

for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(nh4_e(:,:,i))); shading flat
    caxis([0 1.5]); 
    h = colorbar; set(get(h,'title'),'string','NH4(uM)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from 11~17 climate'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_nh4_e'),'-dpng')
end

for i = 1:12
    fig=figure; pcolor(ref_lat_em,ref_dep_e,squeeze(no3_e(:,:,i)));colorbar; shading flat
    caxis([0 15])
    h = colorbar; set(get(h,'title'),'string','NO3(uM)','fontsize',14,'fontweight','bold');
    title(['2001 east boundary from 11~17 climate'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_no3_e'),'-dpng')
end

for i = 1:12
    figure; pcolor(ref_lat_em,ref_dep_e,squeeze(do_e(:,:,i)));colorbar; shading flat
%     caxis([0 200])
end


for i = 1:12
    figure; pcolor(ref_lat_em,ref_dep_e,squeeze(temp_e(:,:,i)));colorbar; shading flat
    caxis([0 30])
end

for i = 1:12
    figure; pcolor(ref_lat_em,ref_dep_e,squeeze(salt_e(:,:,i)));colorbar; shading flat
    caxis([30 34])
end

% south
for i = 1:12
    figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(po4_s(:,:,i)));colorbar; shading flat
%     caxis([0 0.7])
end

for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chlo_s(:,:,i)));colorbar; shading flat
    caxis([0 4])
    h = colorbar; set(get(h,'title'),'string','Chl(mg/L)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from 11~17 climate'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_chl_s'),'-dpng')
end


for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(nh4_s(:,:,i)));colorbar; shading flat
    caxis([0 1.5])
    h = colorbar; set(get(h,'title'),'string','NH4(uM)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from 11~17 climate'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_nh4_s'),'-dpng')
end


% figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(nh4_s(:,:,1)));colorbar; shading flat
% caxis([-inf inf])
% find(nh4_s < 0) 
% find(nh4_e < 0)


for i = 1:12
    fig=figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(no3_s(:,:,i)));colorbar; shading flat
    caxis([0 15])
    h = colorbar; set(get(h,'title'),'string','NO3(uM)','fontsize',14,'fontweight','bold');
    title(['2001 south boundary from 11~17 climate'],'fontsize',14,'fontweight','bold');
    xlabel('lon','fontsize',16)
    ylabel('lat','fontsize',16)
    grid on
    set(gca,'fontsize',16,'fontweight','bold')
    print(fig,strcat('2001_', num2str(i),'_mon_no3_s'),'-dpng')
end

for i = 1:12
    figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(do_s(:,:,i)));colorbar; shading flat
%     caxis([0 200])
end


for i = 1:12
    figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(temp_s(:,:,i)));colorbar; shading flat
    caxis([0 30])
end

for i = 1:12
    figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,i)));colorbar; shading flat
    caxis([30 34])
end