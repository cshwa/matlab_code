function create_roms_file_month2_NPZD(fname,temp,salt,zeta,u,v,NO3,phy,zoo,det)%,tframe,cycle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Create an empty netcdf frc file 
%       x: total number of rho points in x direction
%       y: total number of rho points in y direction
%       varname: name of field variable
%       fname: name of the ecmwf file
%       var: mean file
%
%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nw = netcdf(fname, 'clobber');
result = redef(nw);

disp(['file is ',fname])
[t,n,m]=size(temp);
%
%  Create dimensions
%
disp(['xi_rho is ',num2str(m)])
disp(['eta_rho is ',num2str(n)])
nw('eta_rho') = n ;
nw('xi_rho') = m ;
nw('eta_u') = n;  
nw('xi_u') = m-1;  
nw('eta_v') = n-1;
nw('xi_v') = m;   
nw('eta_psi') = n-1;
nw('xi_psi') = m-1; 
nw('level') = t;
nw('one')=1;
%
%  Create variables and attributes
%

nw.type = ' ROMS file ';
nw.title = ' averaging daily to monthly ';
nw.source = 'roms data ';
nw.author = 'Created by T';
nw.date = date;

% switch varname
% 
% case 'temp'
 nw{'temp'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'temp'}.long_name = ncchar('temperature');
 nw{'temp'}.units = ncchar('Celsius');
 nw{'temp'}(:)= temp;

% case 'salt'
 nw{'salt'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'salt'}.long_name = ncchar('salinity');
 nw{'salt'}.units = ncchar('psu');
 nw{'salt'}(:)= salt;

% case 'zeta'
 nw{'zeta'} = ncfloat(  'eta_rho', 'xi_rho');
 nw{'zeta'}.long_name = ncchar('time-averaged free-surface');
 nw{'zeta'}.units = ncchar('meter');
%  nw{'zeta'}.time = ncchar('level');
 nw{'zeta'}(:)= zeta;

% case 'u'
 nw{'u'} = ncfloat('level', 'eta_u', 'xi_u');
 nw{'u'}.long_name = ncchar('time-averaged u-momentum component');
 nw{'u'}.units = ncchar('meter second-1');
 nw{'u'}(:)= u;
 
% case 'v'
 nw{'v'} = ncfloat('level',  'eta_v', 'xi_v');
 nw{'v'}.long_name = ncchar('time-averaged v-momentum component');
 nw{'v'}.units = ncchar('meter second-1');
 nw{'v'}(:)= v;
 
 nw{'NO3'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'NO3'}.long_name = ncchar('time-averaged nitrate concentration');
 nw{'NO3'}.units = ncchar('millimole_N03 meter-3');
 nw{'NO3'}(:)= NO3;
 
 nw{'phytoplankton'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'phytoplankton'}.long_name = ncchar('time-averaged phytoplankton concentration');
 nw{'phytoplankton'}.units = ncchar('millimole_nitrogen meter-3');
 nw{'phytoplankton'}(:)= phy;
 
 nw{'zooplankton'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'zooplankton'}.long_name = ncchar('time-averaged zooplankton concentration');
 nw{'zooplankton'}.units = ncchar('millimole_nitrogen meter-3');
 nw{'zooplankton'}(:)= zoo;
 
 nw{'detritus'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'detritus'}.long_name = ncchar('time-averaged detritus concentration');
 nw{'detritus'}.units = ncchar('millimole_nitrogen meter-3');
 nw{'detritus'}(:)= det;
 
 close(nw);
% end
