close all; clear all; clc; 

file_n = 'monthly_2016_01_fix_po4.nc';

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

for 1:numel(variable)
    
    var_name = variable{i};
    
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
    
    pick_var=ncread(file_n,var_name);

idx = (mask_l ~= 0); %logical mask to keep nonzeros in z

tempo_in = squeeze(NaNs(size(pick_var,1),size(pick_var,2),size(pick_var,3)));
tempo = NaNs(size(pick_var,1),size(pick_var,2));
for i = 1:size(pick_var,3)
    tempo=squeeze(pick_var(:,:,i));
    x_adj = xx(idx);
    y_adj = yy(idx);
    tempo_adj = tempo(idx);
    tempo_in(:,:,i)=griddata(x_adj,y_adj,tempo_adj,xx,yy,'nearest');
end
tempo_in= squeeze(tempo_in);

nc=netcdf(file_n,'w');
nc{var_name}(:) = tempo_in;
close(nc);


close all; clear all; clc; 

nc=netcdf('gy2015_his_8785.nc','w');

z=nc{'zeta'}(:);
nc{'zeta'}(:) = zeros(size(z,1),size(z,2));

close(nc);


close all; clear all; clc; 
nc=netcdf('1997_tunn1_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('1997_tunn2_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('2000_mp_p_sewer_det_f_366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('1996_2nd_366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('1995_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('2002_tg_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('1996_tg_366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('1999_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('1995_tg_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('1994_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('1994_tg_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('1993_tg_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('1998_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('1993_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('1992_tg_366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('1997_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('1991_tg_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('1992_mp_p_sewer_det_f_366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('1996_mp_p_sewer_det_f_366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('1992_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('2009_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('2008_mp_p_sewer_det_f_366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('yellow_1999_ini_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);




close all; clear all; clc; 
nc=netcdf('yellow_1998_ini_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('2007_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


% nc{'tide_Eamp'}(:)=amp;
% nc{'tide_Ephase'}(:)=phase;

close all; clear all; clc; 
nc=netcdf('yellow_1997_ini_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('yellow_1996_ini_0366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('yellow_1995_ini_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('yellow_1994_ini_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('2006_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('yellow_1992_ini_0366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('yellow_1991_ini_1460.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('gy_1980_v3_bio_avg_0366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('2005_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc)


close all; clear all; clc; 
nc=netcdf('2004_mp_p_sewer_det_f_366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc)


close all; clear all; clc; 
nc=netcdf('2003_mp_p_sewer_det_f_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc)


close all; clear all; clc; 
nc=netcdf('2002_mp_p_sewer_det_365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc)


close all; clear all; clc; 
nc=netcdf('mp_p_sewer_det_2001_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc)


close all; clear all; clc; 
nc=netcdf('avg_365_2001.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc)


close all; clear all; clc; 
nc=netcdf('gy_1980_v3_bio_riv_avg_0366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('gy_1981_v3_bio_riv_avg_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('gy_1982_v3_bio_riv_avg_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);


close all; clear all; clc; 
nc=netcdf('gy_1983_v3_bio_riv_avg_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);



close all; clear all; clc; 
nc=netcdf('gy_1984_v3_bio_riv_avg_0366.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);

close all; clear all; clc; 
nc=netcdf('gy_1985_v3_bio_riv_avg_0365.nc','w');

time = nc{'ocean_time'}(:)

nc{'ocean_time'}(:) = 0;

close(nc);