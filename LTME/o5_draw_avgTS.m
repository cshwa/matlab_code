clear; clc; close all

fpath = 'E:\11.사업\장기생태_3단계\1차년도\Data\from_JJH_정선관측\KODC_mean\';
flist = dir(fullfile(fpath, '*.txt')); % filename with path

lon_lim = [124 130]; lat_lim = [33 37]; % small domain
%lon_lim = [124 133]; lat_lim = [31 39]; % large domain
lim = [lon_lim lat_lim];
fc = [.95 .95 .95 ];

for ts = 1:2
    if ts == 1
        contour_interval = [0:2:35];
        clim = [0 35];
        vari = 'TEMP';
        colorbarname = 'Temperature (deg C)';
    elseif ts == 2
        contour_interval = [30:1:35];
        clim = [30 35];
        vari = 'SALT';
        colorbarname = 'Salinity';
    end
    
    % load map data
    dum2=load('coast_sin.dat');
    coa_lon=dum2(:,1); coa_lat=dum2(:,2);
    
    % Bathymetry file
    zfile = load('E:\11.사업\장기생태_3단계\1차년도\Data\from_JJH_정선관측\KorBathy30s.mat');
    Zlon = zfile.xbathy; Zlat = zfile.ybathy; Zz = zfile.zbathy;
    
    for fi = 1:length(flist)
        fname = flist(fi).name;
                
        % load temperature data
        [ST, LONG, LAT, TEMP, SALT] = textread([fpath, fname],'%d %f %f %f %f');
        
        X1 = min(LONG); X2 = max(LONG);
        Y1 = min(LAT); Y2 = max(LAT);
        Xp = [X1:0.1:X2]; Yp = [Y1:0.1:Y2];
        [Xi, Yi] = meshgrid(Xp, Yp);
        
        zind = find(lon_lim(1) < Zlon & Zlon < lon_lim(2) & lat_lim(1) < Zlat & Zlat < lat_lim(2));
        Zlon2 = Zlon(zind); Zlat2 = Zlat(zind); Zz2 = Zz(zind);
        z_grid = griddata(Zlon2,Zlat2,Zz2, Xi, Yi);
        
        Zi = griddata(LONG, LAT, eval(vari), Xi, Yi);
        Zi(z_grid < str2num(fname(1:2))) = NaN; % Bathymetry
        
        for li = 1:length(Yp)
            lind = isnan(Zi(li,:));
            try
                lind2 = find(diff(lind) == 1);
                Zi(li,lind2(end):-1:lind2(end)-2) = nan;
            catch
            end
            
        end
        
        figure; hold on;
        set(gca,'ydir','nor');
        m_proj('miller','lon',[lim(1) lim(2)],'lat',[lim(3) lim(4)]);
        m_gshhs_i('color','k')
        m_pcolor(Xi, Yi, Zi); colormap('jet'); shading flat;
        [cs, h] = m_contour(Xi, Yi, Zi, contour_interval, 'w');
        clabel(cs, h);
        m_gshhs_i('patch',fc )
        c = colorbar; c.FontSize = 15;
        c.Label.String = colorbarname; c.Label.FontSize = 15;
        caxis(clim);
        m_grid('XTick',lim(1):2:lim(2),'YTick',lim(3):2:lim(4),'linewi',2,'linest','none','tickdir','out','fontsize',20, 'fontweight','bold','FontName', 'Times');
        
        titlename = [vari, ' Avg ', fname(5:6), ' ', fname(1:3),];
        title(titlename, 'fontsize', 25)
                
        % plot point
        for i = 1:length(LAT)
            m_plot(LONG(i), LAT(i), '.k', 'markersize', 5)
        end
        
        saveas(gcf,[fpath, fname(1:end-4), '_', vari],'png');
        
    end
end