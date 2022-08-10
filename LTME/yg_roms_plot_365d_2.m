%--------------------------------------------------------------------------
%
%
%                        Plotting ROMS OUT ( 365daily )
%
%                                             date : 2005. 11. 19
%                                             made by K.S. Moon
%
%                                             date : 2007. 10. 28
%                                             edited by Y.G.
%
%--------------------------------------------------------------------------
clear all;close all;
%==========================================================================
    variation = 1;   level  = 40;

current   = 1;   skip_v = 3;     size_v = 1;    color_v = 'k';

plot_pcolorjw = 1;    temp_lim = [0 35];    salt_lim = [20 34];

plot_contour  = 0;    color_c  = '-w';      level_c  = [20:2:34];

plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'jpg';

grdfile       = 'd:\eastsea model\1993\roms_grid_ADD_03.nc';

yy = 2003;
start_mm=1;
start_dd=1;

end_mm=8;
end_dd=15;

if yeardays(yy)==365
    days = [31 28 31 30 31 30 31 31 30 31 30 31];    
else
    days = [31 29 31 30 31 30 31 31 30 31 30 31];
end

aa=365;
start_day=sum(days(1:start_mm-1))+start_dd;
start_file_num=start_day+1;

end_day=sum(days(1:end_mm-1))+end_dd;
end_file_num=end_day+1;

mm=start_mm;
dd=start_dd;
%==========================================================================

gd       = netcdf(grdfile);
lon_rho  = gd{'lon_rho'}(:);
lat_rho  = gd{'lat_rho'}(:);
mask_rho = gd{'mask_rho'}(:);
angle    = gd{'angle'}(:);
close(gd)
warning off
mask_rho = mask_rho./mask_rho;
warning on
vname = {'temp','salt','zeta','ubar','vbar','u','v','omega'};
     
for i=start_file_num:end_file_num
    if dd>days(mm)
        mm=mm+1;
        dd=1;
    end
		mid=[num2char(i,4)];
        file = ['my_',mid,'.nc'];
        disp([file,' : ', num2char(mm,2),'/',num2char(dd,2)])
            nc=netcdf(file);

            switch variation
                case 1
                    value=nc{char(vname(variation))}(:);
                    val_name='temperature';
                    unit = '^oC';
                    if (level == 20)
                        val_title='Surface Temperature';
                    else
                        val_title='Bottom Temperature';
                    end

                    val_caxis=temp_lim;
                case 2
                    value=nc{char(vname(variation))}(:);
                    val_name='salinity';
                    unit = 'PSU';
                    if (level == 20)
                        val_title='Surface Salinity';
                    else
                        val_title='Bottom Salinity';
                    end
                    val_caxis=salt_lim;

            end
            close(nc)
            data=mask_rho.*squeeze(value(level,:,:));
            clear value;         

            if (current)
                val_title = [val_title,' & Current'];
                nc=netcdf(file);
                u=nc{'u'}(:); v=nc{'v'}(:);
                close(nc)

                u=squeeze(u(level,:,:));   v=squeeze(v(level,:,:));
                [uu,vv,lon,lat,mask]=uv_vec2rho(u,v,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu=mask_rho.*uu;  vv=mask_rho.*vv;
                w = (uu+sqrt(-1).*vv);
                clear u;clear v;
                psliceuv(lon_rho,lat_rho,w,skip_v,size_v,color_v);
                psliceuv(126.2,40.5,0.5,skip_v,size_v,color_v);
                text(127,40.5,'0.5 m/s','color',color_v,'FontSize',10,'fontweight','bold')
            end

            if (plot_pcolorjw)
                hold on
                pcolorjw(lon_rho,lat_rho,data);
            end


            if (plot_contour)
                hold on
                [C,h]=contour(lon_rho,lat_rho,data,level_c,color_c,'linewidth',1);
                clabel(C,h,'FontSize',10,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
            end

            if (plot_geoshow)
                bln=load('coast_yg.dat');
                lat=bln(:,2);
                lon=bln(:,1);
                geoshow(lat,lon,'DisplayType','polygon','facecolor',color_g)
            end

            axis equal;
            axis([min(min(lon_rho)) max(max(lon_rho)) min(min(lat_rho)) max(max(lat_rho))]);
            caxis(val_caxis);
            bar = colorbar('fontsize',12,'fontweight','bold');
            set(get(bar,'title'),'string',unit,'FontSize',15,'fontweight','bold')

            line([min(min(lon_rho)) min(min(lon_rho))],[min(min(lat_rho)) max(max(lat_rho))],'color','k','LineWidth',2.5)
            line([min(min(lon_rho)) max(max(lon_rho))],[max(max(lat_rho)) max(max(lat_rho))],'color','k','LineWidth',2.5)
            line([max(max(lon_rho)) max(max(lon_rho))],[min(min(lat_rho)) max(max(lat_rho))],'color','k','LineWidth',2.5)
            line([min(min(lon_rho)) max(max(lon_rho))],[min(min(lat_rho)) min(min(lat_rho))],'color','k','LineWidth',2.5)
            
            out_name=['SST-2002',num2char(mm,2),num2char(dd,2)];
            
            title([out_name],'FontSize',15,'fontweight','bold');


            if (switch_save)
                saveas(gcf,out_name,out_type);
            end
            clf
        dd=dd+1;
end
     close all
% % % end