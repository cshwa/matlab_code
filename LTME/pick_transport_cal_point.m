close all; clear; clc;
grd_file='D:\장기생태\Dynamic\02_grid_depth\grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');
h=ncread(grd_file, 'h');

% 1. right mouth (low)
 p1_1=[x(85,115), y(85,115)]
 p1_2=[x(82,98), y(82,98)]
 
% 2. right mouth (high)
 p2_1=[x(73,111), y(73,111)]
 p2_2=[x(81,119), y(81,119)] 
 
 % 3. left mouth 
 p3_1=[x(29,98), y(29,98)]
 p3_2=[x(40,98), y(40,98)] 
 
 % 4. myodo (high)
 p4_1=[x(49,77), y(49,77)]
 p4_2=[x(49,67), y(49,67)]
 
  % 5. myodo (low)
 p5_1=[x(49,85), y(49,85)]
 p5_2=[x(49,97), y(49,97)] 
 
  % 6. left entrance (right on the myodo)
 p6_1=[x(78,85), y(78,85)]
 p6_2=[x(77,66), y(77,66)] 
 
  % 7. main channel
 p7_1=[x(77,66), y(77,66)]
 p7_2=[x(101,66), y(101,66)] 
 
  % 8. jinju channel
 p8_1=[x(131,119), y(131,119)]
 p8_2=[x(132,109), y(132,109)] 
 
  % 9. small posco channel
 p9_1=[x(73,111), y(73,111)]
 p9_2=[x(73,105), y(73,105)] 
 
 save('picked_transport_cal_points.mat')