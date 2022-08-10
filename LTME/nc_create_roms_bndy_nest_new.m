% Create a netcdf file of boundary condition data for ROMS
%
% It is assumed the data to be written are on a ROMS 3-D grid, but only the
% boundary segments will be used. (OK, so it uses a lot of space but it
% makes it easy to recycle code).
%
% The ROMS 3-D data must be a structure named 'roms'
%   roms.time      = the times of the data
%   roms.base_date = COARDS format string describing base date
%   roms.temp      = temperature observations (such as from Ridgway's SBT
%                     method)
%   roms.salt
%   roms.{u,v,ubar,vbar} on their respective Arakawa-C grid locations
% 
%   add possibly also...
%   roms.temp_err  = normalized error variance of observations 
%   roms.salt_err
%   roms.Evel on the C-grid rho points <<< !!!!!!!!!! <<<<<<<<<<<<<<<<<
%
% Additional data required in the workspace are:
%
%   time_variable = name of the time variable (string)
%   out_file = the name of the netcdf file generated here
%   grd_file (for the record)
%
% Optional data in the workspace
%   titlestr
%   sourcestr = something aobut where the data came from
%   details = something about computations relevant to these data
%
% John Wilkin
% Updated for EAC 15-Aug-2003

disp(' ')
disp(['The output netcdf file (out_file) will be ' out_file])


switch noclobber
  case 1
    ncid = netcdf.create(out_file,'noclobber');
  case 0
    ncid = netcdf.create(out_file,'clobber');
  otherwise
end
if isempty(ncid)
  error(['Failed to open ' out_file])
end


% global attributes
varid = netcdf.getConstant('GLOBAL');
try
  nc.type = ncchar(roms.type);
catch
  nc.type = ncchar('ROMS BOUNDARY file');
end

if exist('titlestr')
  netcdf.putAtt(ncid,varid,'title',titlestr);
end
if exist('out_file')
  netcdf.putAtt(ncid,varid,'out_file',out_file);
end
if exist('grd_file')
  netcdf.putAtt(ncid,varid,'Grid',grd_file);
end
if exist('sourcestr')
  netcdf.putAtt(ncid,varid,'source',sourcestr);
end
if exist('details')
  netcdf.putAtt(ncid,varid,'details',details);
end
if exist('reference_string')
  netcdf.putAtt(ncid,varid,'reference',reference_string);
end
nc.history = ncchar(['Created by ' which(mfilename) ' - ' datestr(now)]);

% Dimensions

eta_rho_dimID = netcdf.defDim(ncid,'eta_rho',size(roms.grd.lon_rho,1));
xi_rho_dimID = netcdf.defDim(ncid,'xi_rho',size(roms.grd.lon_rho,2));
eta_u_dimID = netcdf.defDim(ncid,'eta_u',size(roms.grd.lon_u,1));
xi_u_dimID = netcdf.defDim(ncid,'xi_u',size(roms.grd.lon_u,2));
eta_v_dimID = netcdf.defDim(ncid,'eta_v',size(roms.grd.lon_v,1));
xi_v_dimID = netcdf.defDim(ncid,'xi_v',size(roms.grd.lon_v,2));
% time_dimID = netcdf.defDim(ncid, time_variable,netcdf.getConstant('NC_UNLIMITED'));
time_dimID = netcdf.defDim(ncid, time_variable,length(bndy_time));

s_rho_dimID = netcdf.defDim(ncid,'s_rho',length(roms.grd.sc_r));
s_w_dimID = netcdf.defDim(ncid,'s_w',length(roms.grd.sc_w));
one_dimID = netcdf.defDim(ncid,'one',1);

% The unlimited dimension
%time_dimension = 'time';
%nc(time_dimension) = 0; % UNLIMITED


% coordinates ------------------------------------------------------------

% THE TIME
theVarname = time_variable;
time_ID= netcdf.defVar(ncid,theVarname,'double',time_dimID);
netcdf.putAtt(ncid,time_ID,'long_name','subsurface temp/salt observations time');
netcdf.putAtt(ncid,time_ID,'units',time_variable_units);
netcdf.putAtt(ncid,time_ID,'cycle_length',cycle_length);


% THE VERTICAL COORDINATE VARIABLES

theVarname = 'theta_s';
theta_s_ID=netcdf.defVar(ncid,theVarname,'double',[one_dimID]);
netcdf.putAtt(ncid,theta_s_ID,'long_name','S-coordinate surface control parameter');
netcdf.putAtt(ncid,theta_s_ID,'units','nondimensional');

theVarname = 'theta_b';
theta_b_ID=netcdf.defVar(ncid,theVarname,'double',[one_dimID]);
netcdf.putAtt(ncid,theta_b_ID,'long_name','S-coordinate bottom control parameter');
netcdf.putAtt(ncid,theta_b_ID,'units','nondimensional');

theVarname = 'Tcline';
Tcline_ID=netcdf.defVar(ncid,theVarname,'double',[one_dimID]);
netcdf.putAtt(ncid,Tcline_ID,'long_name','S-coordinate surface/bottom layer width');
netcdf.putAtt(ncid,Tcline_ID,'units','meter');

theVarname = 'hc';
hc_ID=netcdf.defVar(ncid,theVarname,'double',[one_dimID]);
netcdf.putAtt(ncid,hc_ID,'long_name','S-coordinate parameter, critical depth');
netcdf.putAtt(ncid,hc_ID,'units','meters');


% THE CLIMATOLOGY VARIABLES

% NOTE: These variable names (and the associated time variable names) need
% to correspond the the entires in roms:init_scalars.F

theVarname = 'zeta_east';
zeta_e_ID=netcdf.defVar(ncid,theVarname,'double',[eta_rho_dimID time_dimID]);
netcdf.putAtt(ncid,zeta_e_ID,'long_name','free surface east boundary condition');
netcdf.putAtt(ncid,zeta_e_ID,'units','meter');
netcdf.putAtt(ncid,zeta_e_ID,'field','zeta_east, scalar, series');
netcdf.putAtt(ncid,zeta_e_ID, 'time', time_variable);

theVarname = 'zeta_west';
zeta_w_ID=netcdf.defVar(ncid,theVarname,'double',[eta_rho_dimID time_dimID]);
netcdf.putAtt(ncid,zeta_w_ID,'long_name','free surface west boundary condition');
netcdf.putAtt(ncid,zeta_w_ID,'units','meter');
netcdf.putAtt(ncid,zeta_w_ID,'field','zeta_west, scalar, series');
netcdf.putAtt(ncid,zeta_w_ID, 'time', time_variable);

theVarname = 'zeta_north';
zeta_n_ID=netcdf.defVar(ncid,theVarname,'double',[xi_rho_dimID time_dimID]);
netcdf.putAtt(ncid,zeta_n_ID,'long_name','free surface north boundary condition');
netcdf.putAtt(ncid,zeta_n_ID,'units','meter');
netcdf.putAtt(ncid,zeta_n_ID,'field','zeta_north, scalar, series');
netcdf.putAtt(ncid,zeta_n_ID, 'time', time_variable);

theVarname = 'zeta_south';
zeta_s_ID=netcdf.defVar(ncid,theVarname,'double',[xi_rho_dimID time_dimID ]);
netcdf.putAtt(ncid,zeta_s_ID,'long_name','free surface south boundary condition');
netcdf.putAtt(ncid,zeta_s_ID,'units','meter');
netcdf.putAtt(ncid,zeta_s_ID,'field','zeta_south, scalar, series');
netcdf.putAtt(ncid,zeta_s_ID, 'time', time_variable);

theVarname = 'ubar_east';
ubar_e_ID=netcdf.defVar(ncid,theVarname,'double',[eta_u_dimID time_dimID ]);
netcdf.putAtt(ncid,ubar_e_ID,'long_name','2D u-momentum east boundary condition');
netcdf.putAtt(ncid,ubar_e_ID,'units','meter secon-1');
netcdf.putAtt(ncid,ubar_e_ID,'field','ubar_east, scalar, series');
netcdf.putAtt(ncid,ubar_e_ID, 'time', time_variable);

theVarname = 'ubar_west';
ubar_w_ID=netcdf.defVar(ncid,theVarname,'double',[eta_u_dimID time_dimID ]);
netcdf.putAtt(ncid,ubar_w_ID,'long_name','2D u-momentum west boundary condition');
netcdf.putAtt(ncid,ubar_w_ID,'units','meter secon-1');
netcdf.putAtt(ncid,ubar_w_ID,'field','ubar_west, scalar, series');
netcdf.putAtt(ncid,ubar_w_ID, 'time', time_variable);

theVarname = 'ubar_north';
ubar_n_ID=netcdf.defVar(ncid,theVarname,'double',[xi_u_dimID time_dimID ]);
netcdf.putAtt(ncid,ubar_n_ID,'long_name','2D u-momentum north boundary condition');
netcdf.putAtt(ncid,ubar_n_ID,'units','meter secon-1');
netcdf.putAtt(ncid,ubar_n_ID,'field','ubar_north, scalar, series');
netcdf.putAtt(ncid,ubar_n_ID, 'time', time_variable);

theVarname = 'ubar_south';
ubar_s_ID=netcdf.defVar(ncid,theVarname,'double',[xi_u_dimID time_dimID ]);
netcdf.putAtt(ncid,ubar_s_ID,'long_name','2D u-momentum south boundary condition');
netcdf.putAtt(ncid,ubar_s_ID,'units','meter secon-1');
netcdf.putAtt(ncid,ubar_s_ID,'field','ubar_south, scalar, series');
netcdf.putAtt(ncid,ubar_s_ID, 'time', time_variable);

theVarname = 'vbar_east';
vbar_e_ID=netcdf.defVar(ncid,theVarname,'double',[eta_v_dimID time_dimID ]);
netcdf.putAtt(ncid,vbar_e_ID,'long_name','2D v-momentum east boundary condition');
netcdf.putAtt(ncid,vbar_e_ID,'units','meter secon-1');
netcdf.putAtt(ncid,vbar_e_ID,'field','vbar_east, scalar, series');
netcdf.putAtt(ncid,vbar_e_ID, 'time', time_variable);

theVarname = 'vbar_west';
vbar_w_ID=netcdf.defVar(ncid,theVarname,'double',[eta_v_dimID time_dimID]);
netcdf.putAtt(ncid,vbar_w_ID,'long_name','2D v-momentum west boundary condition');
netcdf.putAtt(ncid,vbar_w_ID,'units','meter secon-1');
netcdf.putAtt(ncid,vbar_w_ID,'field','vbar_west, scalar, series');
netcdf.putAtt(ncid,vbar_w_ID, 'time', time_variable);

theVarname = 'vbar_north';
vbar_n_ID=netcdf.defVar(ncid,theVarname,'double',[xi_v_dimID time_dimID]);
netcdf.putAtt(ncid,vbar_n_ID,'long_name','2D v-momentum north boundary condition');
netcdf.putAtt(ncid,vbar_n_ID,'units','meter secon-1');
netcdf.putAtt(ncid,vbar_n_ID,'field','vbar_north, scalar, series');
netcdf.putAtt(ncid,vbar_n_ID, 'time', time_variable);

theVarname = 'vbar_south';
vbar_s_ID=netcdf.defVar(ncid,theVarname,'double',[xi_v_dimID time_dimID]);
netcdf.putAtt(ncid,vbar_s_ID,'long_name','2D v-momentum south boundary condition');
netcdf.putAtt(ncid,vbar_s_ID,'units','meter secon-1');
netcdf.putAtt(ncid,vbar_s_ID,'field','vbar_south, scalar, series');
netcdf.putAtt(ncid,vbar_s_ID, 'time', time_variable);

theVarname = 'u_east';
u_e_ID=netcdf.defVar(ncid,theVarname,'double',[eta_u_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,u_e_ID,'long_name','3D u-momentum east boundary condition');
netcdf.putAtt(ncid,u_e_ID,'units','meter secon-1');
netcdf.putAtt(ncid,u_e_ID,'field','u_east, scalar, series');
netcdf.putAtt(ncid,u_e_ID, 'time', time_variable);

theVarname = 'u_west';
u_w_ID=netcdf.defVar(ncid,theVarname,'double',[eta_u_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,u_w_ID,'long_name','3D u-momentum west boundary condition');
netcdf.putAtt(ncid,u_w_ID,'units','meter secon-1');
netcdf.putAtt(ncid,u_w_ID,'field','u_west, scalar, series');
netcdf.putAtt(ncid,u_w_ID, 'time', time_variable);

theVarname = 'u_north';
u_n_ID=netcdf.defVar(ncid,theVarname,'double',[xi_u_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,u_n_ID,'long_name','3D u-momentum north boundary condition');
netcdf.putAtt(ncid,u_n_ID,'units','meter secon-1');
netcdf.putAtt(ncid,u_n_ID,'field','u_north, scalar, series');
netcdf.putAtt(ncid,u_n_ID, 'time', time_variable);

theVarname = 'u_south';
u_s_ID=netcdf.defVar(ncid,theVarname,'double',[xi_u_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,u_s_ID,'long_name','3D u-momentum south boundary condition');
netcdf.putAtt(ncid,u_s_ID,'units','meter secon-1');
netcdf.putAtt(ncid,u_s_ID,'field','u_south, scalar, series');
netcdf.putAtt(ncid,u_s_ID, 'time', time_variable);

theVarname = 'v_east';
v_e_ID=netcdf.defVar(ncid,theVarname,'double',[eta_v_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,v_e_ID,'long_name','3D v-momentum east boundary condition');
netcdf.putAtt(ncid,v_e_ID,'units','meter secon-1');
netcdf.putAtt(ncid,v_e_ID,'field','v_east, scalar, series');
netcdf.putAtt(ncid,v_e_ID, 'time', time_variable);

theVarname = 'v_west';
v_w_ID=netcdf.defVar(ncid,theVarname,'double',[eta_v_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,v_w_ID,'long_name','3D v-momentum west boundary condition');
netcdf.putAtt(ncid,v_w_ID,'units','meter secon-1');
netcdf.putAtt(ncid,v_w_ID,'field','v_west, scalar, series');
netcdf.putAtt(ncid,v_w_ID, 'time', time_variable);

theVarname = 'v_north';
v_n_ID=netcdf.defVar(ncid,theVarname,'double',[xi_v_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,v_n_ID,'long_name','3D v-momentum north boundary condition');
netcdf.putAtt(ncid,v_n_ID,'units','meter secon-1');
netcdf.putAtt(ncid,v_n_ID,'field','v_north, scalar, series');
netcdf.putAtt(ncid,v_n_ID, 'time', time_variable);

theVarname = 'v_south';
v_s_ID=netcdf.defVar(ncid,theVarname,'double',[xi_v_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,v_s_ID,'long_name','3D v-momentum south boundary condition');
netcdf.putAtt(ncid,v_s_ID,'units','meter secon-1');
netcdf.putAtt(ncid,v_s_ID,'field','v_south, scalar, series');
netcdf.putAtt(ncid,v_s_ID, 'time', time_variable);

theVarname = 'temp_east';
temp_e_ID=netcdf.defVar(ncid,theVarname,'double',[eta_rho_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,temp_e_ID,'long_name','potential temperature east boundary condition');
netcdf.putAtt(ncid,temp_e_ID,'units','Celcius');
netcdf.putAtt(ncid,temp_e_ID,'field','temp_east, scalar, series');
netcdf.putAtt(ncid,temp_e_ID, 'time', time_variable);

theVarname = 'temp_west';
temp_w_ID=netcdf.defVar(ncid,theVarname,'double',[eta_rho_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,temp_w_ID,'long_name','potential temperature west boundary condition');
netcdf.putAtt(ncid,temp_w_ID,'units','Celcius');
netcdf.putAtt(ncid,temp_w_ID,'field','temp_west, scalar, series');
netcdf.putAtt(ncid,temp_w_ID, 'time', time_variable);

theVarname = 'temp_north';
temp_n_ID=netcdf.defVar(ncid,theVarname,'double',[xi_rho_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,temp_n_ID,'long_name','potential temperature north boundary condition');
netcdf.putAtt(ncid,temp_n_ID,'units','Celcius');
netcdf.putAtt(ncid,temp_n_ID,'field','temp_north, scalar, series');
netcdf.putAtt(ncid,temp_n_ID, 'time', time_variable);

theVarname = 'temp_south';
temp_s_ID=netcdf.defVar(ncid,theVarname,'double',[xi_rho_dimID s_rho_dimID time_dimID]);
netcdf.putAtt(ncid,temp_s_ID,'long_name','potential temperature south boundary condition');
netcdf.putAtt(ncid,temp_s_ID,'units','Celcius');
netcdf.putAtt(ncid,temp_s_ID,'field','temp_south, scalar, series');
netcdf.putAtt(ncid,temp_s_ID, 'time', time_variable);

theVarname = 'salt_east';
salt_e_ID=netcdf.defVar(ncid,theVarname,'double',[eta_rho_dimID s_rho_dimID time_dimID ]);
netcdf.putAtt(ncid,salt_e_ID,'long_name','salinity east boundary condition');
netcdf.putAtt(ncid,salt_e_ID,'units','PSU');
netcdf.putAtt(ncid,salt_e_ID,'field','salt_east, scalar, series');
netcdf.putAtt(ncid,salt_e_ID, 'time', time_variable);

theVarname = 'salt_west';
salt_w_ID=netcdf.defVar(ncid,theVarname,'double',[eta_rho_dimID s_rho_dimID  time_dimID]);
netcdf.putAtt(ncid,salt_w_ID,'long_name','salinity west boundary condition');
netcdf.putAtt(ncid,salt_w_ID,'units','PSU');
netcdf.putAtt(ncid,salt_w_ID,'field','salt_west, scalar, series');
netcdf.putAtt(ncid,salt_w_ID, 'time', time_variable);

theVarname = 'salt_north';
salt_n_ID=netcdf.defVar(ncid,theVarname,'double',[xi_rho_dimID s_rho_dimID  time_dimID ]);
netcdf.putAtt(ncid,salt_n_ID,'long_name','salinity north boundary condition');
netcdf.putAtt(ncid,salt_n_ID,'units','PSU');
netcdf.putAtt(ncid,salt_n_ID,'field','salt_north, scalar, series');
netcdf.putAtt(ncid,salt_n_ID, 'time', time_variable);

theVarname = 'salt_south';
salt_s_ID=netcdf.defVar(ncid,theVarname,'double',[xi_rho_dimID s_rho_dimID  time_dimID ]);
netcdf.putAtt(ncid,salt_s_ID,'long_name','salinity south boundary condition');
netcdf.putAtt(ncid,salt_s_ID,'units','PSU');
netcdf.putAtt(ncid,salt_s_ID,'field','salt_south, scalar, series');
netcdf.putAtt(ncid,salt_s_ID, 'time', time_variable);

% specify additional tracers here

% force write-to-disk of the configuration and coordinates
% will need to re-open the file to write the data

netcdf.endDef(ncid);

%writing the values
netcdf.putVar(ncid, theta_s_ID, gn.theta_s);
netcdf.putVar(ncid, theta_b_ID, gn.theta_b);
netcdf.putVar(ncid, Tcline_ID, gn.Tcline);
netcdf.putVar(ncid, hc_ID, gn.hc);
netcdf.putVar(ncid, time_ID, bndy_time);