clc; clear all; close all;
% filepath = 'D:\ROMS\data\Gwangyang\boundary\';
% grid_path = 'D:\ROMS\data\Gwangyang\grid\';
grid_name = 'grid_sumjin_v1970_fix_3m.nc';
% grdname = [grid_path, grid_name];
ROMS_title  = 'Gwangyang Model';
obc = [1 1 1 1];

g = grd('Gwangyang');
Lonr = g.lon_rho;
Latr = g.lat_rho;
grid_south_lon = Lonr(1,:);
grid_east_lat = Latr(:,end);
h = g.h;
theta_s = g.theta_s;
theta_b = g.theta_b;
hc = g.hc;
N = g.N;
zeta = g.zeta;

for year = 1980 : 1989
    y = num2str(year);
    ss = cumsum(eomday(year,1:12));
    cycle=ss(end);
    roms_time=15:30:cycle;
    
    oldfile = ['Gwangyang_Y',y,'.nc'];
    newfile = ['Gwangyang_bio_Y',y,'.nc'];
    
    copy_file_uv_bio(newfile,oldfile,grid_name,ROMS_title,obc,...
        theta_s,theta_b,hc,N,...
        roms_time,cycle,year,'clobber')
end

