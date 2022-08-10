close all; clear all; clc; 

nc=netcdf('tide_gy_v11_s_2021.nc','w');

amp=nc{'tide_Eamp'}(:);
phase=nc{'tide_Ephase'}(:);

% nc.components = ncchar(''M2 S2 N2 K2 K1 O1 P1 Q1'');

% tide correction
amp(1,:,:) = amp(1,:,:) + (-6.75132/100);
amp(2,:,:) = amp(2,:,:) + (-0.423672/100);
amp(3,:,:) = amp(3,:,:) + (-4.888992/100);
amp(4,:,:) = amp(4,:,:) + (0);
amp(5,:,:) = amp(5,:,:) + (2.877312/100);
amp(6,:,:) = amp(6,:,:) + (0.566928/100);
amp(7,:,:) = amp(7,:,:) + 0;
amp(8,:,:) = amp(8,:,:) + (0.67056/100);

phase(1,:,:) = phase(1,:,:) + (-124.36);
phase(2,:,:) = phase(2,:,:) + (14.91);
phase(3,:,:) = phase(3,:,:) + (-7.31);
phase(4,:,:) = phase(4,:,:) + 0;
phase(5,:,:) = phase(5,:,:) + (11.56);
phase(6,:,:) = phase(6,:,:) + 222.32;
phase(7,:,:) = phase(7,:,:) + 0;
phase(8,:,:) = phase(8,:,:) + (1.1);

nc{'tide_Eamp'}(:)=amp;
nc{'tide_Ephase'}(:)=phase;
close(nc);

%% comfirm correction
close all; clear all; clc; 

nc=netcdf('tide_gy_v11_s_2021.nc','r');
amp=nc{'tide_Eamp'}(:);
phase=nc{'tide_Ephase'}(:);
close(nc);

nc=netcdf('tide_gy_v11_s_0613.nc','r');
amp_ori=nc{'tide_Eamp'}(:);
phase_ori=nc{'tide_Ephase'}(:);
close(nc);

amp(:,74,24) - amp_ori(:,74,24) 
phase(:,74,24) - phase_ori(:,74,24)
