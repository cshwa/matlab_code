close all; clear; clc;
grd_file='D:\장기생태\Dynamic\02_grid_depth\grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');
h=ncread(grd_file, 'h');

% 1. right mouth (low)
 p1_1=[x(87,117), y(87,117)]
 p1_2=[x(83,98), y(83,98)]
 
% 2. right mouth (high)
 p2_1=[x(72,113), y(72,113)]
 p2_2=[x(83,121), y(83,121)] 
 
 % 3. left mouth 
 p3_1=[x(31,102), y(31,102)]
 p3_2=[x(42,104), y(42,104)] 
 
%  % 4. myodo (high)
%  p4_1=[x(49,77), y(49,77)]
%  p4_2=[x(49,67), y(49,67)]
%  
%   % 5. myodo (low)
%  p5_1=[x(49,85), y(49,85)]
%  p5_2=[x(49,97), y(49,97)] 
 
  % 6. left entrance (right on the myodo)
 p6_1=[x(78,90), y(78,90)]
 p6_2=[x(76,66), y(76,66)] 
 
  % 7. main channel
 p7_1=[x(77,66), y(77,66)]
 p7_2=[x(103,68), y(103,68)] 
 
  % 8. jinju channel
 p8_1=[x(131,122), y(131,122)]
 p8_2=[x(132,105), y(132,105)] 
 
%   % 9. small posco channel
%  p9_1=[x(73,111), y(73,111)]
%  p9_2=[x(73,105), y(73,105)] 
 
 save('picked_transport_cal_points_v2.mat')