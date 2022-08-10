function [lon_lim, lat_lim] = domain_J(location)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain_list = 'KODC_large',      lon_lim = [124 133]; lat_lim = [31 39]
%               'YECS_small',      lon_lim = [124 130]; lat_lim = [33 37]
%               'JWLee',           lon_lim = [124 130]; lat_lim = [33 38]
%               'YECS_large',      lon_lim = [117 130]; lat_lim = [27 41]
%               'YECS_flt',        lon_lim = [117.5 130.3]; lat_lim = [31 41]
%               'NWP',             lon_lim = [117 160]; lat_lim = [16 50]
%               'DA',              lon_lim = [123 135]; lat_lim = [30 41]
%               'TI',              lon_lim = [124 126.7]; lat_lim = [32 34]
%               'ECMWF_monthly'    lon_lim = [120 130]; lat_lim = [30 40]
%               'ECMWF_pressure'   lon_lim = [80 160]; lat_lim = [10 55]
%               'onshore'          lon_lim = [120 145]; lat_lim = [24 43]
%               'KS'               lon_lim = [125 132]; lat_lim = [30 36]
%               'Tsugaru'          lon_lim = [139.5 142]; lat_lim = [40.5 42.5]
%               'Taean'            lon_lim = [125 128]; lat_lim = [35 38]
%               'LT' 장기해양생태   lon_lim=[126.125 129.125]; lat_lim=[33.125 35.125]; 
%
% [lon_lim, lat_lim] = domain_J('Domain_list')
%
% J. Jung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch location
    case 'KODC_large'
        lon_lim = [124 133]; lat_lim = [31 39];
        
    case 'YECS_small'
        lon_lim = [124 130]; lat_lim = [33 37];
        
    case 'JWLee'
        lon_lim = [124 130]; lat_lim = [33 38];
        
    case 'YECS_large'
        lon_lim = [117 130]; lat_lim = [27 41];
        
    case 'YECS_flt'
        lon_lim = [117.5 130.3]; lat_lim = [31 41];
        
    case {'NWP', 'NWP_ET1K30'}
        lon_lim = [115 162]; lat_lim = [15 52];
        
    case 'DA'
        lon_lim = [123 135]; lat_lim = [30 41];
        
    case 'TI' % temperature inversion
        lon_lim = [124 126.7]; lat_lim = [32 34];
        
    case 'ECMWF_monthly'
        lon_lim = [119 134]; lat_lim = [30 42];
        
    case 'ECMWF_pressure'
        lon_lim = [80 160]; lat_lim = [10 55];
        
    case 'onshore'
        lon_lim = [120 145]; lat_lim = [24 43];
        
    case 'KS'
        lon_lim = [125 132]; lat_lim = [30 36];
        
    case 'Tsugaru'
        lon_lim = [139.5 142]; lat_lim = [40.5 42.5];
        
    case 'Taean'
        lon_lim = [125 128]; lat_lim = [35 38];
        
    case 'LT'     
        lon_lim=[126.125 129.125]; lat_lim=[33.125 35.125]; 
        
end
end