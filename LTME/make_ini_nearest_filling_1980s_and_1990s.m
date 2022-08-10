close all; clear all; clc; 

file_n = '1982to1986_monthly_2016_01_fix_po4.nc';

x=ncread(file_n,'lon_rho');
y=ncread(file_n,'lat_rho');
xv=ncread(file_n,'lon_v');
yv=ncread(file_n,'lat_v');
xu=ncread(file_n,'lon_u');
yu=ncread(file_n,'lat_u');
mask=ncread(file_n,'mask_rho');
maskv=ncread(file_n,'mask_v');
masku=ncread(file_n,'mask_u');

variable = {'temp','salt','NO3','NH4','tPO4','zeta',...
    'w','chlorophyll','phytoplankton','zooplankton','LdetritusN','SdetritusN',...
    'oxygen','SdetritusP','LdetritusP','ubar','vbar','u','v'};

for i= 1:numel(variable)
    
    var_name = variable{i};
    
    disp(var_name)
    
    switch var_name       
        case {'ubar','u'}
            xx=xu; yy=yu;
            mask_l = masku;
            
        case {'vbar','v'}
            xx=xv; yy=yv;
            mask_l = maskv;
            
        case {'temp','salt','NO3','NH4','tPO4','zeta',...
    'w','chlorophyll','phytoplankton','zooplankton','LdetritusN','SdetritusN',...
    'oxygen','SdetritusP','LdetritusP'}
            xx=x; yy=y;
            mask_l = mask;     
    end
%% load    
    pick_var=ncread(file_n,var_name);
%% interpolation
idx = (mask_l ~= 0); %logical mask to keep nonzeros in z
tempo_in = squeeze(NaN(size(pick_var,1),size(pick_var,2),size(pick_var,3)));
tempo = NaN(size(pick_var,1),size(pick_var,2));
for i = 1:size(pick_var,3)
    tempo=squeeze(pick_var(:,:,i));
    x_adj = xx(idx);
    y_adj = yy(idx);
    tempo_adj = tempo(idx);
    tempo_in(:,:,i)=griddata(x_adj,y_adj,tempo_adj,xx,yy,'nearest');
end
tempo_in= squeeze(tempo_in);

%% make form to put it in 
    switch var_name       
        case {'ubar','vbar','zeta'}
            tempo_in = permute(tempo_in,[2 1]);
            
        case {'temp','salt','NO3','NH4','tPO4','u','v'...
    'w','chlorophyll','phytoplankton','zooplankton','LdetritusN','SdetritusN',...
    'oxygen','SdetritusP','LdetritusP'}
            tempo_in = permute(tempo_in,[3 2 1]);
    end

nc=netcdf(file_n,'w');
nc{var_name}(:) = tempo_in;
close(nc);
end