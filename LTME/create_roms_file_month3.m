function create_roms_file_month3(fname,temp,salt,zeta,ubar,vbar,u,v,AKv,AKt,mm)%,tframe,cycle)
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
nw('s_w')=t+1;
nw('s_rho') = t;
nw('one')=1;
%
%  Create variables and attributes
%

nw.type = ' ROMS file ';
nw.title = ' averaging daily to monthly ';
nw.source = 'roms data ';
nw.author = 'Created by T';
nw.date = date;

nw{'time'} = ncdouble('one');
nw{'time'}.long_name = ncchar('Year Month (1: January, ... 12: December)');
nw{'time'}.long_name = 'Year Month (1: January, ... 12: December)';
nw{'time'}.units = ncchar('Month');
nw{'time'}.units = 'Month';
% nw{'time'}.cycle_length = ncdouble(cycle);
% nw{'time'}.cycle_length = cycle;
nw{'time'}(:)=mm;

% switch varname
% 
% case 'temp'
 nw{'temp'} = ncfloat('s_rho', 'eta_rho', 'xi_rho');
 nw{'temp'}.long_name = ncchar('temperature');
 nw{'temp'}.units = ncchar('Celsius');
%  nw{'temp'}.time = ncchar('s_rho');
 nw{'temp'}(:)= temp;

% case 'salt'
 nw{'salt'} = ncfloat('s_rho', 'eta_rho', 'xi_rho');
 nw{'salt'}.long_name = ncchar('salinity');
 nw{'salt'}.units = ncchar('psu');
%  nw{'salt'}.time = ncchar(var_time);
 nw{'salt'}(:)= salt;

% case 'zeta'
 nw{'zeta'} = ncfloat(  'eta_rho', 'xi_rho');
 nw{'zeta'}.long_name = ncchar('time-averaged free-surface');
 nw{'zeta'}.units = ncchar('meter');
%  nw{'zeta'}.time = ncchar(var_time);
 nw{'zeta'}(:)= zeta;

% case 'ubar'
 nw{'ubar'} = ncfloat(  'eta_u', 'xi_u');
 nw{'ubar'}.long_name = ncchar('time-averaged vertically integrated u-momentum component');
 nw{'ubar'}.units = ncchar('meter second-1');
%  nw{'ubar'}.time = ncchar(var_time);
 nw{'ubar'}(:)= ubar;
 
% case 'vbar'
 nw{'vbar'} = ncfloat(  'eta_v', 'xi_v');
 nw{'vbar'}.long_name = ncchar('time-averaged vertically integrated v-momentum component');
 nw{'vbar'}.units = ncchar('meter second-1');
%  nw{'vbar'}.time = ncchar(var_time);
 nw{'vbar'}(:)= vbar;
 
% case 'u'
 nw{'u'} = ncfloat('s_rho', 'eta_u', 'xi_u');
 nw{'u'}.long_name = ncchar('time-averaged u-momentum component');
 nw{'u'}.units = ncchar('meter second-1');
%  nw{'u'}.time = ncchar(var_time);
 nw{'u'}(:)= u;
 
% case 'v'
 nw{'v'} = ncfloat('s_rho',  'eta_v', 'xi_v');
 nw{'v'}.long_name = ncchar('time-averaged v-momentum component');
 nw{'v'}.units = ncchar('meter second-1');
%  nw{'v'}.time = ncchar(var_time);
 nw{'v'}(:)= v;
 
 % case 'AKv'
 nw{'AKv'} = ncfloat('s_w', 'eta_rho', 'xi_rho');
 nw{'AKv'}.long_name = ncchar('vertical viscosity coefficient');
 nw{'AKv'}.units = ncchar('meter2 second-1');
%  nw{'AKv'}.time = ncchar(var_time);
 nw{'AKv'}(:)= AKv;
 
  % case 'AKt'
 nw{'AKt'} = ncfloat('s_w', 'eta_rho', 'xi_rho');
 nw{'AKt'}.long_name = ncchar('vertical diffusion coefficient');
 nw{'AKt'}.units = ncchar('meter2 second-1');
%  nw{'AKt'}.time = ncchar(var_time);
 nw{'AKt'}(:)= AKt;
 
 % case 'shflux'
%  nw{'shflux'} = ncfloat(  'eta_rho', 'xi_rho');
%  nw{'shflux'}.long_name = ncchar('time-averaged surface net heat flux');
%  nw{'shflux'}.units = ncchar('watt meter-2');
% %  nw{'shflux'}.time = ncchar(var_time);
%  nw{'shflux'}(:)= shflux;
 
% % case 'omega'
%  nw{'omega'} = ncfloat (var_time_1, 'eta_rho', 'xi_rho');
%  nw{'omega'}.long_name = ncchar('time-averaged S-coordinate vertical momentum componen');
%  nw{'omega'}.units = ncchar('meter3 second-1');
% %  nw{'omega'}.time = ncchar(var_time);
%  nw{'omega'}(:)= omega;
 close(nw);
% end
