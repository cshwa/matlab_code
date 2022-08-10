clc;clear all; close all

for year=2001:2012

tsuv_file=['G:\auto_fennel\lateral\NWP_1_10_test6\data\roms_bndy_auto_NWP_1_10_test6_',num2str(year),'.nc'];
% tsuv_file=['G:\auto_fennel\lateral\NWP_1_10_test6\data\roms_bndy_auto_NWP_1_10_test6_clim06_10.nc'];

bio_file=['G:\auto_fennel\lateral\NWP_1_10_test6\data\roms_bndy_auto_NWP_',num2str(year),'_fennel_woa13.nc'];
% bio_file=['G:\auto_fennel\lateral\NWP_1_10_test6\data\roms_bndy_auto_NWP_1_10_test6_clim06_10_fennel_woa13.nc'];

nc=netcdf(tsuv_file,'r');
nw=netcdf(bio_file,'write');

obc=[1 1 0 0];

for obcndx=1:4
    if obc(obcndx)==1
      if obcndx==1
        disp(' Processing southern boundary...')
    suffix='_south';
      elseif obcndx==2
        disp(' Processing eastern boundary...')
    suffix='_east';
      elseif obcndx==3
        disp(' Processing northern boundary...')
    suffix='_north';
      elseif obcndx==4
        disp(' Processing western boundary...')
    suffix='_west';
      end
    end
    tempname=['temp',suffix];
    saltname=['salt',suffix];
    uname=['u',suffix];
    vname=['v',suffix];
    ubarname=['ubar',suffix];
    vbarname=['vbar',suffix];
    zetaname=['zeta',suffix];
    temp=nc{tempname}(:);
    salt=nc{saltname}(:);
    u=nc{uname}(:);
    v=nc{vname}(:);
    ubar=nc{ubarname}(:);
    vbar=nc{vbarname}(:);
    zeta=nc{zetaname}(:);
    nw{tempname}(:)=temp;
    nw{saltname}(:)=salt;
    nw{uname}(:)=u;
    nw{vname}(:)=v;
    nw{ubarname}(:)=ubar;
    nw{vbarname}(:)=vbar;
    nw{zetaname}(:)=zeta;
end

close(nc);
close(nw);
end
