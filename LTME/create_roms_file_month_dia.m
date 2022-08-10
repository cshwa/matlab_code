function create_roms_file_month_dia(fname,ubar_cor,vbar_cor,ubar_hadv,vbar_hadv,ubar_xadv,vbar_xadv, ...
    ubar_yadv,vbar_yadv,ubar_hvisc,vbar_hvisc,ubar_xvisc,vbar_xvisc,ubar_yvisc,vbar_yvisc, ...
    ubar_prsgrd,vbar_prsgrd,ubar_sstr,vbar_sstr,ubar_bstr,vbar_bstr, ubar_accel,vbar_accel,  ... 
    u_prsgrd,v_prsgrd,u_vvisc,v_vvisc,u_hvisc,v_hvisc,temp_hdiff,temp_vdiff,salt_hdiff,salt_vdiff);
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
[t,n,m]=size(temp_hdiff);
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

 nw{'ubar_cor'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_cor'}.long_name = ncchar('2D u-momentum, Coriolis term');
 nw{'ubar_cor'}.units = ncchar('meter second-2');
 nw{'ubar_cor'}(:)= ubar_cor;
 
 nw{'vbar_cor'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_cor'}.long_name = ncchar('2D v-momentum, Coriolis term');
 nw{'vbar_cor'}.units = ncchar('meter second-2');
 nw{'vbar_cor'}(:)= vbar_cor;
 
 nw{'vbar_xadv'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_xadv'}.long_name = ncchar('2D v-momentum, horizontal XI-advection term');
 nw{'vbar_xadv'}.units = ncchar('meter second-2');
 nw{'vbar_xadv'}(:)= vbar_xadv;

 nw{'ubar_hadv'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_hadv'}.long_name = ncchar('2D u-momentum, horizontal advection term');
 nw{'ubar_hadv'}.units = ncchar('meter second-2');
 nw{'ubar_hadv'}(:)= ubar_hadv;
 
 nw{'vbar_hadv'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_hadv'}.long_name = ncchar('2D v-momentum, horizontal advection term');
 nw{'vbar_hadv'}.units = ncchar('meter second-2');
 nw{'vbar_hadv'}(:)= vbar_hadv;

 nw{'ubar_xadv'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_xadv'}.long_name = ncchar('2D u-momentum, horizontal XI-advection term');
 nw{'ubar_xadv'}.units = ncchar('meter second-2');
 nw{'ubar_xadv'}(:)= ubar_xadv;
 
 nw{'vbar_xadv'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_xadv'}.long_name = ncchar('2D v-momentum, horizontal XI-advection term');
 nw{'vbar_xadv'}.units = ncchar('meter second-2');
 nw{'vbar_xadv'}(:)= vbar_xadv;
 
 nw{'ubar_yadv'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_yadv'}.long_name = ncchar('2D u-momentum, horizontal ETA-advection term');
 nw{'ubar_yadv'}.units = ncchar('meter second-2');
 nw{'ubar_yadv'}(:)= ubar_yadv;
 
 nw{'vbar_yadv'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_yadv'}.long_name = ncchar('2D v-momentum, horizontal ETA-advection term');
 nw{'vbar_yadv'}.units = ncchar('meter second-2');
 nw{'vbar_yadv'}(:)= vbar_yadv;
 
 nw{'ubar_hvisc'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_hvisc'}.long_name = ncchar('2D u-momentum, horizontal viscosity term');
 nw{'ubar_hvisc'}.units = ncchar('meter second-2');
 nw{'ubar_hvisc'}(:)= ubar_hvisc;
 
 nw{'vbar_hvisc'} = ncfloat( 'eta_v', 'xi_v');
 nw{'vbar_hvisc'}.long_name = ncchar('2D v-momentum, horizontal viscosity term');
 nw{'vbar_hvisc'}.units = ncchar('meter second-2');
 nw{'vbar_hvisc'}(:)= vbar_hvisc;
 
 nw{'ubar_xvisc'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_xvisc'}.long_name = ncchar('2D u-momentum, XI-viscosity term');
 nw{'ubar_xvisc'}.units = ncchar('meter second-2');
 nw{'ubar_xvisc'}(:)= ubar_xvisc;
 
 nw{'vbar_xvisc'} = ncfloat( 'eta_v', 'xi_v');
 nw{'vbar_xvisc'}.long_name = ncchar('2D v-momentum, XI-viscosity term');
 nw{'vbar_xvisc'}.units = ncchar('meter second-2');
 nw{'vbar_xvisc'}(:)= vbar_xvisc;
 
 nw{'ubar_yvisc'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_yvisc'}.long_name = ncchar('2D u-momentum, ETA-viscosity term');
 nw{'ubar_yvisc'}.units = ncchar('meter second-2');
 nw{'ubar_yvisc'}(:)= ubar_yvisc;
 
 nw{'vbar_yvisc'} = ncfloat( 'eta_v', 'xi_v');
 nw{'vbar_yvisc'}.long_name = ncchar('2D v-momentum, ETA-viscosity term');
 nw{'vbar_yvisc'}.units = ncchar('meter second-2');
 nw{'vbar_yvisc'}(:)= vbar_yvisc;
 
 nw{'ubar_prsgrd'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_prsgrd'}.long_name = ncchar('2D u-momentum, pressure gradient term');
 nw{'ubar_prsgrd'}.units = ncchar('meter second-2');
 nw{'ubar_prsgrd'}(:)= ubar_prsgrd;
 
 nw{'vbar_prsgrd'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_prsgrd'}.long_name = ncchar('2D v-momentum, pressure gradient term');
 nw{'vbar_prsgrd'}.units = ncchar('meter second-2');
 nw{'vbar_prsgrd'}(:)= vbar_prsgrd;

 nw{'ubar_sstr'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_sstr'}.long_name = ncchar('2D u-momentum, surface stress term');
 nw{'ubar_sstr'}.units = ncchar('meter second-2');
 nw{'ubar_sstr'}(:)= ubar_sstr;
 
 nw{'vbar_sstr'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_sstr'}.long_name = ncchar('2D v-momentum, surface stress term');
 nw{'vbar_sstr'}.units = ncchar('meter second-2');
 nw{'vbar_sstr'}(:)= vbar_sstr;
 
 nw{'ubar_bstr'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_bstr'}.long_name = ncchar('2D u-momentum, bottom stress term');
 nw{'ubar_bstr'}.units = ncchar('meter second-2');
 nw{'ubar_bstr'}(:)= ubar_bstr;
 
 nw{'vbar_bstr'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_bstr'}.long_name = ncchar('2D v-momentum, bottom stress term');
 nw{'vbar_bstr'}.units = ncchar('meter second-2');
 nw{'vbar_bstr'}(:)= vbar_bstr;
 
 nw{'ubar_accel'} = ncfloat('eta_u', 'xi_u');
 nw{'ubar_accel'}.long_name = ncchar('2D u-momentum, acceleration term');
 nw{'ubar_accel'}.units = ncchar('meter second-2');
 nw{'ubar_accel'}(:)= ubar_accel;
 
 nw{'vbar_accel'} = ncfloat('eta_v', 'xi_v');
 nw{'vbar_accel'}.long_name = ncchar('2D v-momentum, acceleration term');
 nw{'vbar_accel'}.units = ncchar('meter second-2');
 nw{'vbar_accel'}(:)= vbar_accel;

 nw{'u_prsgrd'} = ncfloat('level', 'eta_u', 'xi_u');
 nw{'u_prsgrd'}.long_name = ncchar('3D u-momentum, pressure gradient term');
 nw{'u_prsgrd'}.units = ncchar('meter second-2');
 nw{'u_prsgrd'}(:)= u_prsgrd;
 
 nw{'v_prsgrd'} = ncfloat('level', 'eta_v', 'xi_v');
 nw{'v_prsgrd'}.long_name = ncchar('3D v-momentum, pressure gradient term');
 nw{'v_prsgrd'}.units = ncchar('meter second-2');
 nw{'v_prsgrd'}(:)= v_prsgrd;
 
 nw{'u_vvisc'} = ncfloat('level', 'eta_u', 'xi_u');
 nw{'u_vvisc'}.long_name = ncchar('3D u-momentum, vertical viscosity term');
 nw{'u_vvisc'}.units = ncchar('meter second-2');
 nw{'u_vvisc'}(:)= u_vvisc;
 
 nw{'v_vvisc'} = ncfloat('level', 'eta_v', 'xi_v');
 nw{'v_vvisc'}.long_name = ncchar('3D v-momentum, vertical viscosity term');
 nw{'v_vvisc'}.units = ncchar('meter second-2');
 nw{'v_vvisc'}(:)= v_vvisc;
 
 nw{'u_hvisc'} = ncfloat('level', 'eta_u', 'xi_u');
 nw{'u_hvisc'}.long_name = ncchar('3D u-momentum, horizontal viscosity term');
 nw{'u_hvisc'}.units = ncchar('meter second-2');
 nw{'u_hvisc'}(:)= u_hvisc;
 
 nw{'v_hvisc'} = ncfloat('level', 'eta_v', 'xi_v');
 nw{'v_hvisc'}.long_name = ncchar('3D v-momentum, horizontal viscosity term');
 nw{'v_hvisc'}.units = ncchar('meter second-2');
 nw{'v_hvisc'}(:)= v_hvisc;
 
 nw{'temp_hdiff'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'temp_hdiff'}.long_name = ncchar('potential temperature, horizontal diffusion term');
 nw{'temp_hdiff'}.units = ncchar('Celsius second-1');
 nw{'temp_hdiff'}(:)= temp_hdiff;
 
 nw{'temp_vdiff'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'temp_vdiff'}.long_name = ncchar('potential temperature, vertical diffusion term');
 nw{'temp_vdiff'}.units = ncchar('Celsius second-1');
 nw{'temp_vdiff'}(:)= temp_vdiff;
 
 nw{'salt_hdiff'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'salt_hdiff'}.long_name = ncchar('salinity, horizontal diffusion term');
 nw{'salt_hdiff'}.units = ncchar('nondimensional second-1');
 nw{'salt_hdiff'}(:)= salt_hdiff;
 
 nw{'salt_vdiff'} = ncfloat('level', 'eta_rho', 'xi_rho');
 nw{'salt_vdiff'}.long_name = ncchar('salinity, vertical diffusion term');
 nw{'salt_vdiff'}.units = ncchar('nondimensional second-1');
 nw{'salt_vdiff'}(:)= salt_vdiff;
 
 
 close(nw);
% end
