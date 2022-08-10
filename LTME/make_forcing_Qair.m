%   clc;close all;clear all;
% 
% for year=2017
%     clear t_time
%     P='Pair';T='Tair';D='Dair';Q='Qair';foot='.nc';
%     
%     yy=yeardays(year);stryy=num2str(yy);
%     
%     rawdata=['2017/frc_ecmwf_',num2str(year),P,'_',stryy,foot];
% %     rawdata=['easttest_',num2str(year),P,'_',stryy,foot];
%     Pnc=netcdf(rawdata);
%     rawdata=['2017/frc_ecmwf_',num2str(year),T,'_',stryy,foot];
% %     rawdata=['easttest_',num2str(year),T,'_',stryy,foot];
%     Tnc=netcdf(rawdata);
%     rawdata=['2017/frc_ecmwf_',num2str(year),D,'_',stryy,foot];
% %     rawdata=['easttest_',num2str(year),D,'_',stryy,foot];
%     Dnc=netcdf(rawdata);
% 
%     lon=Pnc{'lon_rho'}(:);
%     lat=Pnc{'lat_rho'}(:);
%     t_time = Pnc{'Pair_time'}(:);
%     ttt = length(t_time) ;
%         for dd=1:1:ttt;
%     
%             P_value = Pnc{P}(dd,:,:);
%             T_value = Tnc{T}(dd,:,:);
%             D_value = Dnc{D}(dd,:,:);
% 
%             Qair(dd,:,:) = (qsat_yg(D_value,P_value)./qsat_yg(T_value,P_value)).*100;
%             
%         end
%     close(Pnc);  close(Tnc);  close(Dnc);
% 
%     fname=['2017/frc_ecmwf_',num2str(year),'Qair_',stryy,foot];
% %     fname=['easttest_',num2str(year),'Qair_',stryy,foot];
%     create_roms_forcing_Qair(fname,Qair,t_time,yy,yy)
% end


clc;close all;clear all;
for year = 2001:2010
    clearvars -except year
grdfile='/home/cshwa/Long/02_grid_depth/smoothing/grid_gy_v11_s.nc';
tinterval=1;
ROMS_surface_forcing_Qair_1(grdfile, year, tinterval)
end