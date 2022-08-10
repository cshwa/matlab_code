clc;close all;clear all;
var = 'ssrd';

for year = 2012
    disp(year)

%     grdfile='D:\ROMS\eastsea\data\grid\2012\new_test9\new_etopo1_Eastsea_test9.nc';
%     inpath = ['D:\ROMS\eastsea\data\air_forcing\ECMWF\',var,'\'];
%     outpath = ['D:\ROMS\eastsea\data\air_forcing\',num2str(year),'\'];
%     rawdata=[inpath,'ECMWF_Interim_',var,'_',num2str(year),'.nc'];
    grdfile = 'D:\ROMS\eastsea\data\test\mix_test_grid_etopo1.nc';
    inpath = ['D:\ROMS\eastsea\data\air_forcing\2012\'];
    outpath = 'D:\ROMS\eastsea\data\air_forcing\2012\mix_test\';
    rawdata = '2012-3.nc';
    nc1=netcdf(rawdata);
    ly = (( mod(year,4) == 0 & mod(year,100) ~= 0 ) | mod(year,400) == 0);
    if ly,    yy = 366;
    else     yy = 365;
        
    end
    
    longitude=nc1{'longitude'}(:);
    latitude=nc1{'latitude'}(:);
    time=[0.5:1:yy-0.5];
    swrad=nc1{'ssrd'}(:);
    SWsfactor=nc1{'ssrd'}.scale_factor(:);
    SWaoffset=nc1{'ssrd'}.add_offset(:);
    size_SW=size(swrad);
    close(nc1);
    swrad=(swrad.*SWsfactor+SWaoffset);
    
    
    SW_daily=ones(yy,size_SW(2),size_SW(3))*NaN;
    i_c=1;
    for i=1:4:yy*4;
        SW_three(1+(i_c-1)*8,:,:)=swrad(i,:,:)/10800/4;  % 12시 부터 00시까지 누적 -00시
        SW_three(3+(i_c-1)*8,:,:)=swrad(i+1,:,:)/10800/2; % 불필요
        SW_three(5+(i_c-1)*8,:,:)=swrad(i+2,:,:)/10800/4; % 00시 부터 12시까지 누적 -12시
        SW_three(7+(i_c-1)*8,:,:)=swrad(i+3,:,:)/10800/2; % 불필요
        i_c=i_c+1;
    end
    i_c=1;
    for i=1+yy*4:4:yy*8;
        SW_three(2+(i_c-1)*8,:,:)=swrad(i,:,:)/10800/1; % 불필요
        SW_three(4+(i_c-1)*8,:,:)=swrad(i+1,:,:)/10800/3; % 불필요
        SW_three(6+(i_c-1)*8,:,:)=swrad(i+2,:,:)/10800/1; % 불필요
        SW_three(8+(i_c-1)*8,:,:)=swrad(i+3,:,:)/10800/3; % 불필요
        i_c=i_c+1;
    end
    
    for i=1:1:yy;
        %     aa=sum(swrad(1+4*(i-1):4*i,:,:));
        %     bb=sum(swrad((yy*4+1)+4*(i-1):(yy*4)+4*i,:,:));
        %     SW_daily(i,:,:)=(aa+bb)/8;
        SW_daily(i,:,:)=(SW_three(1+(i-1)*8,:,:)+SW_three(5+(i-1)*8,:,:))/2;
    end
    
    clear swrad nc1;
    
    grd=netcdf(grdfile);
    lat_rho=grd{'lat_rho'}(:);
    lon_rho=grd{'lon_rho'}(:);
    close(grd);
    
    size_inp_SW=size(lat_rho);
    inp_SW_daily=ones(yy,size_inp_SW(1),size_inp_SW(2))*NaN;
    [Xi,Yi]=meshgrid(longitude,latitude);
    
    for ii=1:1:yy;
        Zi=squeeze(SW_daily(ii,:,:));
        Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
        inp_SW_daily(ii,:,:)=Z;
    end
    
    clear SW_daily grd Z lat_u lon_u Xi Yi Zi
    fname=[outpath,'eastsea_',num2str(year),'_swrad.nc'];
    create_roms_forcing_SWD(fname,inp_SW_daily,time,yy,yy)
end


% i_c=1;
% for i=1:4:yy*4;
%     SW_three(1+(i_c-1)*8,:,:)=swrad(i,:,:)/10800;
%     SW_three(3+(i_c-1)*8,:,:)=(swrad(i+1,:,:)-swrad(i+yy*4,:,:))/10800;
%     SW_three(5+(i_c-1)*8,:,:)=(swrad(i+2,:,:)-swrad(i+1+yy*4,:,:))/10800;
%     SW_three(7+(i_c-1)*8,:,:)=(swrad(i+3,:,:)-swrad(i+2+yy*4,:,:))/10800;
%     i_c=i_c+1;
% end
% i_c=1;
% for i=1+yy*4:4:yy*8;
%     SW_three(2+(i_c-1)*8,:,:)=(swrad(i,:,:)-swrad(i-yy*4,:,:))/10800;
%     SW_three(4+(i_c-1)*8,:,:)=(swrad(i+1,:,:)-swrad(i+1-yy*4,:,:))/10800;
%     SW_three(6+(i_c-1)*8,:,:)=(swrad(i+2,:,:)-swrad(i+2-yy*4,:,:))/10800;
%     SW_three(8+(i_c-1)*8,:,:)=(swrad(i+3,:,:)-swrad(i+3-yy*4,:,:))/10800;
%     i_c=i_c+1;
% end
