% old 2005
%    !copy roms_grid2_ADD_04_2.nc roms_grid2_ADD_05_2.nc
%--- 아래 과정을 반복하면서 원하는 정도로 smoothing
%--- for loop 를 사용하는 것도 반복하는 방법
f=netcdf('grid_gy_v11_s.nc','write')
% f=netcdf('grid_sumjin_v04_1st_masked.nc','write')
bathm=f{'h'}(:);
hmin=3.0; % 최소수심
rmax=0.1; % smoothing 정도
Ip=find(bathm>5000);
bathm(Ip)=5000;
newh=smoothgrid(bathm,hmin,rmax);
f{'h'}(:)=newh;
close(f)
