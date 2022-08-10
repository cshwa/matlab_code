close all; clear; clc; 
do=ncread('monthly_2016_01_fix.nc','oxygen');
nh4=ncread('monthly_2016_01_fix.nc','NH4');
no3=ncread('monthly_2016_01_fix.nc','NO3');
chl=ncread('monthly_2016_01_fix.nc','chlorophyll');
phy=ncread('monthly_2016_01_fix.nc','phytoplankton');
zoo=ncread('monthly_2016_01_fix.nc','zooplankton');
detS=ncread('monthly_2016_01_fix.nc','SdetritusN');
detL=ncread('monthly_2016_01_fix.nc','LdetritusN');

find(do <= 0)
find(nh4 <= 0)
find(no3 <= 0)
find(chl <= 0)
find(zoo <= 0)
find(detS <= 0)
find(detL <= 0)


find(isnan(do) == 1)
find(isnan(nh4) == 1)
find(isnan(no3) == 1)
find(isnan(chl) == 1)
find(isnan(zoo) == 1)
find(isnan(detS) == 1)
find(isnan(detL) == 1)