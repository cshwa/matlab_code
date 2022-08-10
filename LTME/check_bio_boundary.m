close all; clear; clc;
no3_s=ncread('Gwangyang_bio_Y1980.nc','NO3_south');
no3_e=ncread('Gwangyang_bio_Y1980.nc','NO3_east');

nh4_s=ncread('Gwangyang_bio_Y1980.nc','NH4_south');
nh4_e=ncread('Gwangyang_bio_Y1980.nc','NH4_east');

do_s=ncread('Gwangyang_bio_Y1980.nc','oxygen_south');
do_e=ncread('Gwangyang_bio_Y1980.nc','oxygen_east');

chl_s=ncread('Gwangyang_bio_Y1980.nc','chlo_south');
chl_e=ncread('Gwangyang_bio_Y1980.nc','chlo_east');

load boundary_bio_input.mat

find((no3_s - no3_input) ~= 0)
find((no3_e - no3_input_we) ~= 0)

find((nh4_s - nh4_input_sn) ~= 0)
find((nh4_e - nh4_input_we) ~= 0)

find((chl_s - chl_input_sn) ~= 0)
find((chl_e - chl_input_we) ~= 0)

find((do_s - squeeze(do_input(:,:,1,:))) ~= 0)
find((do_e - squeeze(do_input_we(:,:,1,:))) ~= 0)



load interp_bnd_fix_2001.mat
i=1
pcolor(ref_lon_sm,ref_dep_s,squeeze(ref_temp(:,1,:))); colorbar; shading flat;




