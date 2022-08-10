addpath(genpath(['/home/yjtak/roms/matlab/']));


for year=2001:2001;
    clearvars -except year
    ufile=strcat('/data1/yjtak/auto_fennel/input/ERA5/auto_ERA5_',num2str(year,'%04i'),'_Uwind.nc');
    vfile=strcat('/data1/yjtak/auto_fennel/input/ERA5/auto_ERA5_',num2str(year,'%04i'),'_Vwind.nc');
    g=grd('auto7');
    cosa=cos(g.angle);
    sina=sin(g.angle);

    unc=netcdf(ufile,'w');
    u=unc{'Uwind'}(:);

    vnc=netcdf(vfile,'w');
    v=vnc{'Vwind'}(:);

    cosa_mat = repmat(cosa,[1,1,length(u)]);
    cosa_mat = permute(cosa_mat,[3 1 2]);

    sina_mat = repmat(sina,[1,1,length(v)]);
    sina_mat = permute(sina_mat,[3 1 2]);

    urot =u.*cosa_mat + v.*sina_mat;
    vrot =v.*cosa_mat - u.*sina_mat;

    unc{'Uwind'}(:)=urot;
    vnc{'Vwind'}(:)=vrot;
    close(unc);
    close(vnc);
end
