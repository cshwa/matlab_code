  clc;close all;clear all;

for year=2009
    clear t_time
    P='Pair';T='Tair';D='Dair';Q='Qair';foot='.nc';
    
    yy=yeardays(year);stryy=num2str(yy);
    
    rawdata=['frc_ecmwf_',num2str(year),P,'_',stryy,foot];
%     rawdata=['easttest_',num2str(year),P,'_',stryy,foot];
    Pnc=netcdf(rawdata);
    rawdata=['frc_ecmwf_',num2str(year),T,'_',stryy,foot];
%     rawdata=['easttest_',num2str(year),T,'_',stryy,foot];
    Tnc=netcdf(rawdata);
    rawdata=['frc_ecmwf_',num2str(year),D,'_',stryy,foot];
%     rawdata=['easttest_',num2str(year),D,'_',stryy,foot];
    Dnc=netcdf(rawdata);

    lon=Pnc{'lon_rho'}(:);
    lat=Pnc{'lat_rho'}(:);
    t_time = Pnc{'Pair_time'}(:);
    ttt = length(t_time) ;
        for dd=1:1:ttt;
    
            P_value = Pnc{P}(dd,:,:);
            T_value = Tnc{T}(dd,:,:);
            D_value = Dnc{D}(dd,:,:);

            Qair(dd,:,:) = (qsat_yg(D_value,P_value)./qsat_yg(T_value,P_value)).*100;
            
        end
    close(Pnc);  close(Tnc);  close(Dnc);

    fname=['frc_ecmwf_',num2str(year),'Qair_',stryy,foot];
%     fname=['easttest_',num2str(year),'Qair_',stryy,foot];
    create_roms_forcing_Qair(fname,Qair,t_time,yy,yy)
end