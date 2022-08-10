clc;clear all;close all

[raw,txt]=xlsread('D:\포항산업\model\forcing\기상청150228-150507.xlsx','Sheet1','');
% 8 - air T, 9 - dir wind 1min, 11 - mag wind 1min, 12 - dir wind 10min,
% 14 - mag wind 10min, 15 - huminity, 16 - air P
dirwind=raw(:,9)*(pi/180);magwind=raw(:,11);
dirwind=flipud(dirwind);magwind=flipud(magwind);
[u,v] = pol2cart(dirwind,magwind);
for i=1:1:(length(magwind)-1)/12
data(i)=mean(magwind((i-1)*12+1:i*12));
end

time=[58.5:0.5:254];

grdfile='D:\포항산업\model\grid\grid_pohang.nc';

grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);
[n,m]=size(lat_rho);

fname=['pohang_frc_KMA_2014Uwind_12hourly.nc'];

nw = netcdf(fname, 'clobber');
result = redef(nw);

disp(['file is ',fname])

set_time=length(time);
%
%  Create dimensions
%
disp(['xi_rho is ',num2str(m)])
disp(['eta_rho is ',num2str(n)])
nw('eta_rho') = n; 
nw('xi_rho') = m; 
nw('eta_u') = n;  
nw('xi_u') = m-1;  
nw('eta_v') = n-1;
nw('xi_v') = m;   
nw('eta_psi') = n-1;
nw('xi_psi') = m-1; 
nw('Uwind_time') =set_time;
cycle=365;
%
%  Create variables and attributes
%

nw.type = ' ROMS forcing file ';
nw.title = ' Bulk Formular Forcing file ';
nw.source = 'KMA';
nw.author = 'Created by T';
nw.date = date ;

nw{'Uwind_time'} = ncfloat('Uwind_time');
nw{'Uwind_time'}.long_name = ncchar('Year 12hourly');
nw{'Uwind_time'}.long_name = 'Year 12hourly';
nw{'Uwind_time'}.units = ncchar('12hourly');
nw{'Uwind_time'}.units = '12hourly';
nw{'Uwind_time'}.cycle_length = ncdouble(cycle);
nw{'Uwind_time'}.cycle_length = cycle;
nw{'Uwind_time'}(:)=time;

 nw{'Uwind'} = ncfloat('Uwind_time', 'eta_rho', 'xi_rho');
 nw{'Uwind'}.long_name = ncchar('KMA Two meter U');
 nw{'Uwind'}.units = ncchar('m/s');
 nw{'Uwind'}.time = ncchar('Uwind_time');
  
for ii=1:1:length(time);
    nw{'Uwind'}(ii,:,:)=ones(n,m)*u(ii);
end
 close(nw);
clear U_daily grd Z lat_u lon_u Xi Yi Zi
