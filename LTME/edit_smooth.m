% old 2005
%    !copy roms_grid2_ADD_04_2.nc roms_grid2_ADD_05_2.nc
%--- �Ʒ� ������ �ݺ��ϸ鼭 ���ϴ� ������ smoothing
%--- for loop �� ����ϴ� �͵� �ݺ��ϴ� ���
f=netcdf('grid_gy_v11_s.nc','write')
% f=netcdf('grid_sumjin_v04_1st_masked.nc','write')
bathm=f{'h'}(:);
hmin=3.0; % �ּҼ���
rmax=0.1; % smoothing ����
Ip=find(bathm>5000);
bathm(Ip)=5000;
newh=smoothgrid(bathm,hmin,rmax);
f{'h'}(:)=newh;
close(f)
