close all; clear; clc;
load interp_bnd_fix_2002.mat
%south
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(d3_temp_s_t(:,:,1)));colorbar; shading flat
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

% close all; clear; clc;
bnd_f = 'bndy_2002_Gwangyang.nc'
clearvars temp_s temp_e salt_s salt_e chlo_s chlo_e
temp_s=ncread(bnd_f,'temp_south');
salt_s=ncread(bnd_f','salt_south');
chlo_s=ncread(bnd_f,'chlo_south');
no3_s=ncread(bnd_f,'NO3_south');
nh4_s=ncread(bnd_f','NH4_south');
temp_e=ncread(bnd_f,'temp_east');
salt_e=ncread(bnd_f,'salt_east');
chlo_e=ncread(bnd_f,'chlo_east'); 
no3_e=ncread(bnd_f,'NO3_east');
nh4_e=ncread(bnd_f,'NH4_east');
chlo_e=ncread(bnd_f,'chlo_east'); 


%south
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(temp_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chlo_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chlo_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(temp_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(salt_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(chlo_e(:,:,1)));colorbar; shading flat


% close all; clear; clc;
bnd_f = 'bndy_2001_Gwangyang_P.nc'
clearvars temp_s temp_e salt_s salt_e chlo_s chlo_e
temp_s=ncread(bnd_f,'temp_south');
salt_s=ncread(bnd_f','salt_south');
chlo_s=ncread(bnd_f,'chlo_south');
no3_s=ncread(bnd_f,'NO3_south');
nh4_s=ncread(bnd_f','NH4_south');
temp_e=ncread(bnd_f,'temp_east');
salt_e=ncread(bnd_f,'salt_east');
chlo_e=ncread(bnd_f,'chlo_east'); 
no3_e=ncread(bnd_f,'NO3_east');
nh4_e=ncread(bnd_f,'NH4_east');
chlo_e=ncread(bnd_f,'chlo_east'); 


%south
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(temp_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chlo_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(salt_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lon_sm,ref_dep_s,squeeze(chlo_s(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(temp_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(salt_e(:,:,1)));colorbar; shading flat
figure; pcolor(ref_lat_em,ref_dep_e,squeeze(chlo_e(:,:,1)));colorbar; shading flat

