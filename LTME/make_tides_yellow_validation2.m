close all; clear; clc; 
 

% cd /data1/cshwa/auto_fennel/input
cshwa='yellow_tide_1996.nc'
yjtak='..\input\roms_tide_1996.nc'


t_c_t=ncread(cshwa,'tide_Ephase');
t_y_t=ncread(yjtak,'tide_Ephase');

t_c_e=ncread(cshwa,'tide_Eamp');
t_y_e=ncread(yjtak,'tide_Eamp');

return