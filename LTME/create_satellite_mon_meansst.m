function create_satellite_mon_meansst(fname,sst,sstscfactor,lat,lon)%,tframe,cycle)
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
[t,n,m]=size(sst);
var_time=['time'];
%
%  Create dimensions
%
disp(['lon is ',num2str(m)])
disp(['lat is ',num2str(n)])
nw('lat') = n; 
nw('lon') = m; 
nw(var_time) = t;
%
%  Create variables and attributes
%

nw.type = ' Satellite monthly mean SST ';
nw.title = ' Satellite monthly mean SST ';
nw.source = 'NOAA/National Climatic Data Center';
nw.author = 'Created by T';
nw.date = date;

nw{var_time} = ncdouble(var_time);
nw{var_time}.long_name = ncchar('Year Month (1: January, ... 12: December)');
nw{var_time}.long_name = 'Year Month (1: January, ... 12: December)';
nw{var_time}.units = ncchar('DAYS');
nw{var_time}.units = 'DAYS';

 nw{'sst'} = ncfloat(var_time, 'lat', 'lon');
 nw{'sst'}.long_name = ncchar('Monthly sea surface temperature');
 nw{'sst'}.units = ncchar('degrees C');
 nw{'sst'}.time = ncchar(var_time);
 nw{'sst'}.scale_factor = ncdouble(sstscfactor);
% for i=1:1:yy
 nw{'sst'}(:)= sst;%(i,:,:);
% end

 nw{'lat'} = ncfloat('lat');
 nw{'lat'}.long_name = ncchar('Latitude');
 nw{'lat'}.units = ncchar('degrees_north');
 nw{'lat'}.grids= ncchar('Uniform grid from -89.875 to 89.875 by 0.25');
 nw{'lat'}(:)= lat;
 
 nw{'lon'} = ncfloat('lon');
 nw{'lon'}.long_name = ncchar('Longitude');
 nw{'lon'}.units = ncchar('degrees_east');
 nw{'lon'}.grids= ncchar('Uniform grid from 0.125 to 359.875 by 0.25');
 nw{'lon'}(:)= lon;
 close(nw);
