clc;clear all;close all

nc=netcdf('ocean_his_9999.nc','w');

t=nc{'ocean_time'}(:);
% 
% h(40,71)=9;
% h(41,71)=9;

t1 =0
nc{'ocean_time'}(:)=t1;

close(nc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tide grid
clc;clear all;close all

nc=netcdf('grid_pohang_csh_fine_flat_tide.nc','w');

h=nc{'h'}(:);
mask=nc{'mask_rho'}(:);
% 
% h(40,71)=9;
% h(41,71)=9;

h(h<=10)=10;

nc{'h'}(:)=h;

close(nc);

pcolor(h.*(mask./mask))
colorbar
xlabel('i on the model','fontsize',20)
ylabel('j on the model','fontsize',20)
set(gca,'fontsize',20);
