function grd = grd(location)
  
% if nargin == 0
%   location = 'eas'; % default
% end
    scoord = [5 0.4 50 20]; 
    
switch location
    case 'auto8'
    grd_file = '/data1/yjtak/auto_fennel/grid/roms_grd_auto_rdrg2_new8_smooth.nc'; % new 1/10
    scoord = [5 0.4 500 40];    
    case 'auto7'
    grd_file = '/data1/yjtak/auto_fennel/grid/roms_grd_auto_rdrg2_new7_smooth.nc'; % new 1/10
    scoord = [5 0.4 500 40];    
    case 'ideal9_sym_lin'
    grd_file = '/home/yjtak/roms/ideal9_7/input/ideal9_sym_linear_fplane_grid.nc'; % new 1/10
    scoord = [1 1 1 20];    
    case 'ideal8_sym_lin_c1'
    grd_file = '/home/yjtak/roms/ideal8_7/input/ideal8_sym_linear_connect1_grid.nc'; % new 1/10
    scoord = [1 1 1 20];
    case 'ideal8_sym_lin_c2'
    grd_file = '/home/yjtak/roms/ideal8_7/input/ideal8_sym_linear_connect2_grid.nc'; % new 1/10
    scoord = [1 1 1 20];
    case 'ideal3_sym_lin'
    grd_file = '/home/yjtak/roms/ideal3_7/input/ideal3_sym_linear_grid.nc'; % new 1/10
    scoord = [1 1 1 20];
    case 'NWP'
    grd_file = '/home/yjtak/make_frc_ecmwf/roms_add_grid2.nc'; % new 1/10
    scoord = [1 1 1 40];
    case 'auto'
    grd_file = '/data1/yjtak/auto_fennel/grid/roms_grd_auto_rdrg2_new5_smooth.nc'; % new 1/10
    scoord = [5 .4 500 40];
    case 'ADD2'
    grd_file = '/home/yjtak/roms/shift_run7/input/roms_add_grid2.nc'; % new 1/10
    scoord = [5 0.4 5 40];    
    case 'NWP_0p04'
    grd_file = '/data2/yjtak/roms_grid_combine2_0p04.nc'; % new 1/10
    scoord = [1 1 1 40];
    case 'NWP_jhjung'
    grd_file = '/data1/yjtak/make_obc_ROMS_NWP/roms_grid_NWP.nc'; 
    scoord = [10 0 250 40];
    case 'NWP_1_10_kimyy'
    grd_file = '/data1/RCM/CMIP5/CMIP5_input/nwp_1_20/input/test57/run/2006/roms_grid_nwp_1_20_test57.nc'; 
    scoord = [10 1 250 40];
    case 'GY_cshwa'
    grd_file = 'D:\장기생태\Dynamic\02_grid_depth\grid_gy_v11_s.nc'; 
    scoord = [1 1 1 20];
end
disp(' ')
disp([ 'Loading ROMS grd for application: ' location])
disp([ 'using grid file ' grd_file])
disp(' ')
grd = roms_get_grid(grd_file,scoord);

  
