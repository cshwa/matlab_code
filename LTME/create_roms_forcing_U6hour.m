function create_roms_forcing_U(fname,Uwind,Uwind2,Uwind3,Uwind4,time,yy,cycle)%,tframe,cycle)
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
[t,n,m]=size(Uwind);

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
nw('Uwind_time') = t;
%
%  Create variables and attributes
%
t
nw.type = ' ROMS forcing file ';
nw.title = ' Bulk Formular Forcing file ';
nw.source = 'ECMWF-interim ';
nw.author = 'Created by T';
nw.date = date ;

nw{'Uwind_time'} = ncfloat('Uwind_time');
nw{'Uwind_time'}.long_name = ncchar('Year 6hourly');
nw{'Uwind_time'}.long_name = 'Year 6hourly';
nw{'Uwind_time'}.units = ncchar('6hourly');
nw{'Uwind_time'}.units = '6hourly';
nw{'Uwind_time'}.cycle_length = ncdouble(cycle);
nw{'Uwind_time'}.cycle_length = cycle;
nw{'Uwind_time'}(:)=time;

 nw{'Uwind'} = ncfloat('Uwind_time', 'eta_rho', 'xi_rho');
 nw{'Uwind'}.long_name = ncchar('ECMWF Ten meter U');
 nw{'Uwind'}.units = ncchar('m/s');
 nw{'Uwind'}.time = ncchar('Uwind_time');
for i=1:1:yy
 nw{'Uwind'}(i,:,:)= Uwind(i,:,:);
end
for i=yy+1:1:yy*2
 nw{'Uwind'}(i,:,:)= Uwind2(i-yy*1,:,:);
end
for i=yy*2+1:1:yy*3
 nw{'Uwind'}(i,:,:)= Uwind3(i-yy*2,:,:);
end
for i=yy*3+1:1:yy*4;
nw{'Uwind'}(i,:,:)= Uwind4(i-yy*3,:,:);
end
 close(nw);
