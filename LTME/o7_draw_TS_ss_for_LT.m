clear; clc; close all

Year_all = [2011:2016];
fig_type = 'mean'; % Normal or Anomaly

casename = 'LT';
% lon_lim=[126.125 129.125]; lat_lim=[33.125 35.125]
[lon_lim, lat_lim] = domain_J(casename);

for ts = 1:1 % 수온과 염분중 선택
        if ts == 1
        contour_interval = [10:2:30];
        clim = [10 30];
        vari = 'TEMP';
        colorbarname = 'Temperature (deg C)';
        elseif ts == 2
        contour_interval = [30:1:35];
        clim = [30 35];
        vari = 'SALT';
        colorbarname = 'Salinity';
        end
    
        fpath1 = 'E:\11.사업\장기생태_3단계\1차년도\Data\from_JJH_정선관측\KODC_';
        titletype = [' '];
        colormap_style = 'jet';
   
        % Bathymetry file
        zfile = load('E:\11.사업\장기생태_3단계\1차년도\Data\from_JJH_정선관측\KorBathy30s.mat');
        Zlon = zfile.xbathy; Zlat = zfile.ybathy; Zz = zfile.zbathy;
    
    for yi = 1:length(Year_all); % 전체 수심

        year = Year_all(yi);
        fpath = [fpath1, num2char(year,4), '\'];        
        flist = dir(fullfile(fpath, '*.txt')); % filename with path
        
        for fi = 1:6       
            fname = flist(fi).name;
            
            % load temperature data
            [ST, LONG, LAT, TEMP, SALT] = textread([fpath, fname],'%d %f %f %f %f');
            
%             X1 = min(LONG); X2 = max(LONG);
%             Y1 = min(LAT); Y2 = max(LAT);
            X1 = 126.125; X2 = 129.125;
            Y1 = 33.125; Y2 = 35.125;
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
            LT_mean(fi,:,:) = Zi;
            ss = mean(LT_mean);
            ss2 = squeeze(ss);
    
        end
        

%             figure;
%%
Zi = ss2;
LT_mean2(yi,:,:) = Zi;
%             figure('Position', [100, 100, 500*1.5, 400*1.5])
%             map_J(casename)
%             m_pcolor(Xi, Yi, Zi); colormap(colormap_style); shading flat;
%             [cs, h] = m_contour(Xi, Yi, Zi, contour_interval, 'w');
%             clabel(cs, h);
%             c = colorbar; c.FontSize = 15;
%             c.Label.String = colorbarname; c.Label.FontSize = 15;
%             caxis(clim);
%             
%             titlename = [vari, titletype, num2str(year), ' ', fname(1:3)];
%             title(titlename, 'fontsize', 25)
%             
%             % plot point
%             for i = 1:length(LAT)
%                 m_plot(LONG(i), LAT(i), '.k', 'markersize', 5)
%             end
            
%             saveas(gcf,[fname(1:end-4), '_', vari],'png');
            %--- save figure --------------------------------------
%             outpath = '.\figure_';
%             out_vari = [fname(1:end-4), '_', vari];
% %             t = abs(xkm_sect);
%             out=[outpath,fig_type,'\',out_vari];
%             set(gcf,'renderer','painter');
%             set(gcf, 'PaperUnits', 'inches');
%             x_width=5*1.5;
%             y_width=4*1.5;
%             set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%             saveas(gcf,out,'tif');
            %------------------------------------------------------
    end
end


