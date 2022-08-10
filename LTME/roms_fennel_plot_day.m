% edit by G.H.SEO on JAN. 2008
%--------------------------------------------------------------------------
clc;clear all; close all;


addpath(genpath('/home/yjtak/roms/matlab'));
addpath(genpath('/home/yjtak/Dropbox/source/matlab/Common/m_map'));

syear = 2009; eyear = 2009;

spinup=2;

file_start=1+(365*(spinup-1));  
file_end =file_start;
file_end =120+(365*(spinup-1));


mid_name = 'ocean_avg_'; foot ='.nc';

%==========================================================================
variation =12;  %  1 = temperature , 2 = salinity , 3 = density , 4 = zeta , 5=rain 6=evapo  7=e-p 8 = nh4  9=no3  10=po4 11=do 12=chla  0 == none

level_T  =40;  % 20 = surface , 1 = bottom

flag_record_gif=1;

current =0; cur_interval=10;

skip_v = 1;     size_v = 2;     color_v = 'k';

plot_pcolorjw = 1;
temp_lim = [0 33];    salt_lim = [0 34];     zeta_lim = [-.3 .3]; dens_lim = [23 27]; rain_lim=[0 10];evap_lim=[0 10]; ep_lim=[-10 10];
nh4_lim=[0 2];no3_lim=[0 10];po4_lim=[0 0.5]; do_lim=[0 500]; chl_lim=[0 5];
level_a  = [-3:2:35]; level_b  = [20:1:35];    level_z  = [-.3 :0.01:.1 ]; level_d  = [23:0.5:27]; level_rain=[0:1:10];level_evap=[0:1:10];level_ep=[-10:1:10]; 
level_nh4=[-10:20:200];level_no3=[0:4:10];level_po4=[0:0.1:0.5];level_do=[0:50:500];level_chl=[0:1:10];
color_c  = '-k';  color_g = 'white';   out_type = 'tif';

out_vari='_';

%==========================================================================
grdfile       = '/data1/yjtak/auto_fennel/grid/roms_grd_auto_rdrg2_new3_smooth.nc';

gd       = netcdf(grdfile);
lon_rho  = gd{'lon_rho'}(:);
lat_rho  = gd{'lat_rho'}(:);
mask_rho = gd{'mask_rho'}(:);
angle    = gd{'angle'}(:);
close(gd)
warning off;
mask_rho = mask_rho./mask_rho;
warning on;

domaxis=[min(min(lon_rho)) max(max(lon_rho)) min(min(lat_rho)) max(max(lat_rho))+0.075];

%domaxis=[115 130.5 24.5 42];
%domaxis=[120 129 31 37];

Mon_T={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for yid  = syear:eyear
    in_path =['./spinup_',num2str(yid),'/daily_mean/'];
    out_path=['./spinup_',num2str(yid),'/daily_mean/'];
    y_num = num2str(yid);
    kkk=0;
    tc=1; % �ð� ���Կ� count
    for mid  = [file_start:file_end];
%     for mid = [22,154,225,319];
        kkk=kkk+1;
%         if kkk < 10
%             m_num = ['_',y_num,'_0',num2str(kkk)];
%         else
%             m_num = ['_',y_num,'_',num2str(kkk)];
%         end
%         filename = [in_path,mid_name,m_num,foot];
		midd=[num2char(mid,4)];
        filename = [in_path,'ocean_avg_',midd,'.nc'];
        nc=netcdf(filename);
        temp = nc{'temp'}(:);
        salt = nc{'salt'}(:);
        zeta = nc{'zeta'}(:);
                u=nc{'u'}(:);
                v=nc{'v'}(:);
                ubar=nc{'ubar'}(:);
                vbar=nc{'vbar'}(:);
        NH4 = nc{'NH4'}(:);
        NO3 = nc{'NO3'}(:);
%         tPO4 = nc{'tPO4'}(:);
        DO = nc{'oxygen'}(:);
        phyto=nc{'phytoplankton'}(:);
        chlo=nc{'chlorophyll'}(:);        
        close(nc)
        if variation == 1
                temp_T=mask_rho.*squeeze(temp(level_T,:,:));
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                u_T=ubar;v_T=vbar;
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='Temp.(^{o}C)';
                caxis_lim=temp_lim;
                cont_level=level_a;
                vari_name='SST';
                out_vari='SST';
        else if variation == 2
                temp_T=mask_rho.*squeeze(salt(level_T,:,:));
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='Sal.(psu)';
                caxis_lim=salt_lim;
                cont_level=level_b;
                vari_name='Sea Surface Salinity';
                out_vari='SSS';
                
            else if variation == 3
                dens=sw_dens(mask_rho.*squeeze(salt(level_T,:,:)),mask_rho.*squeeze(temp(level_T,:,:)),0)-1000;
                val_name='density';
                unit = '(kg/m^3-1000)';
                temp_T=abs(dens);
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='��';
                caxis_lim=dens_lim;
                cont_level=level_d;
                vari_name='Sea Surface Density';
                out_vari='SSD';
                
            else if variation == 4
                temp_T=mask_rho.*zeta;
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='zeta(cm)';
                caxis_lim=zeta_lim;
                cont_level=level_z;
                vari_name='Sea Surface Height';
                out_vari='SSH';
                color=cbrewer('div','RdBu',50);color=flipud(color);
                
            else if variation == 5
                temp_T=mask_rho.*rain;
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='rain rate(kg m^{-2}s^{-1})';
                caxis_lim=rain_lim;
                cont_level=level_rain;
                vari_name='Precipitation Rate';
                out_vari='rain';                
            else if variation == 6
                temp_T=mask_rho.*evaporation;
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='evap.rate(kg m^{-2}s^{-1})';
                caxis_lim=evap_lim;
                cont_level=level_evap;
                vari_name='Evaporation Rate';
                out_vari='Evap';
            else if variation == 7
                temp_T=mask_rho.*evaporation-mask_rho.*rain;
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='e-p.rate(kg m^{-2}s^{-1})';
                caxis_lim=ep_lim;
                cont_level=level_ep;
                vari_name='E-P Rate';
                out_vari='E-P';
                color=cbrewer('div','RdBu',50);color=flipud(color);
            else if variation == 8
                temp_T=mask_rho.*squeeze(NH4(level_T,:,:));
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='��M';
                caxis_lim=nh4_lim;
                cont_level=level_nh4;
                vari_name='ammonium';
                out_vari='NH4';
            else if variation == 9
                temp_T=mask_rho.*squeeze(NO3(level_T,:,:));
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='��M';
                caxis_lim=no3_lim;
                cont_level=level_no3;
                vari_name='nitrate';
                out_vari='NO3';
            else if variation == 10
                temp_T=mask_rho.*squeeze(tPO4(level_T,:,:));
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='��M';
                caxis_lim=po4_lim;
                cont_level=level_po4;
                vari_name='PO4';
                out_vari='PO4';               
            else if variation == 11
                temp_T=mask_rho.*squeeze(DO(level_T,:,:));
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='��M';
                caxis_lim=do_lim;
                cont_level=level_do;
                vari_name='DO';
                out_vari='DO';                    
            else if variation == 12
                temp_T=mask_rho.*squeeze(chlo(level_T,:,:));
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                bar_name='ug/L';
                caxis_lim=chl_lim;
                cont_level=level_chl;
                vari_name='chlorophyll';
                out_vari='chla';                      
           else
                u_T=squeeze(u(level_T,:,:));   v_T=squeeze(v(level_T,:,:));
                [uu_T,vv_T,lon,lat,mask]=uv_vec2rho(u_T,v_T,lon_rho,lat_rho,angle,mask_rho,1,[0 0 0 0]);
                uu_T=mask_rho.*uu_T;  vv_T=mask_rho.*vv_T;
                vari_name='Sea Surface Current';
                end
                end
                end
                end
                end
                end
                end
                end
                end
                end
            end
        end
        
        iinan = find(100<abs(uu_T) );
        uu_T(iinan)=NaN;
        vv_T(iinan)=NaN;
        
        ttime = [num2str(yid),'.',num2char(mid-(365*(spinup-1)),3)];
        figure('position',[400 100 600 450],'PaperUnits','inches','PaperPosition',[0 0 5 5.2]);
        set(gca,'Position',[0.15 0.08 0.75 0.85]);
%         subplot('position',[0. 0.14 width height])
        m_proj('mercator',...
        'lon',[domaxis(1) domaxis(2)],...
        'lat',[domaxis(3) domaxis(4)+0.01]);

%                         'lon',[115 142.5],...
%             'lat',[24 52]);


            hold on        
        if variation ~=0
            m_pcolor(lon_rho,lat_rho,temp_T);
            colormap('jet');
            if variation ==4 ||variation ==7 , colormap(color), end
            bar_h=colorbar;
            caxis(caxis_lim);
            title(bar_h,bar_name,'FontSize',15,'fontweight','bold');
            set(gca,'fontsize',15)
            shading flat
            [C1,h1]=m_contour(lon_rho,lat_rho,temp_T,cont_level,'color','k','linewidth',1);
            clabel(C1,h1,'Color','k','fontsize',13)
        end
%             hpv3 = m_vec(1, lon, lat, uu_T, vv_T, 'm', ...
%      'shaftwidth', 0.2, 'headlength', 2.5);
%           Q = m_quiver(lon,lat, uu_T,vv_T,'k');
%          qk = quiverkey(Q, 0.95, 1.05, 1, '1 m/s', labelpos='W')

ratio_x=(domaxis(4)-domaxis(3))/100;
ratio_y=(domaxis(2)-domaxis(1))/100;
%         m_gshhs_l('color','k');
%         m_gshhs_l('patch','w');
%         timename=[num2str(yid),'.',num2str(mid)];
        m_grid('linewi',2,'linest','none','tickdir','out','fontsize',15);
        m_text(domaxis(2)-40*ratio_x,domaxis(4)-(3*ratio_y),vari_name,'color','k','FontSize',16,'fontweight','bold')
        m_text(domaxis(2)-40*ratio_x,domaxis(4)-(10*ratio_y),{ttime},'fontsize',17,'color','k','fontweight','bold')
%         m_text(domaxis(1)+ratio_x,domaxis(4)-(9*ratio_y),Mon_T{kkk},'fontsize',16,'color','k','fontweight','bold')
        if (current)
    for j=max(find(domaxis(3)>lat_rho(:,1))):cur_interval:min(find(domaxis(4)<lat_rho(:,1)))
        for i=max(find(domaxis(1)>lon_rho(1,:))):cur_interval:min(find(domaxis(2)<lon_rho(1,:)))
            if isnan(uu_T(j,i)) == 0
                hh=m_quiver(lon(j,i), lat(j,i), uu_T(j,i), vv_T(j,i), 5, 'Color',[.6 .6 .6],'linewidth',1,'maxheadsize',100);
                %adjust_quiver_arrowhead_size(hh,8.1);

            else,continue
            end
        end
    end
%plot scale
m_quiver(domaxis(1)+(6*ratio_x),domaxis(4)-(18*ratio_y), 0.5,  0.0,2,'color',[.6 .6 .6],'linewidth',.5);
m_text(domaxis(1)+(6*ratio_x),domaxis(4)-(20*ratio_y), '0.5 m s^{-1}','fontsize',10,'fontweight','bold');

end

%         title('Oyashio Grid Test','fontsize',14,'fontweight','bold');
%         m_grid('box','fancy',...
%             'xtick',5,'ytick',5,'tickdir','out',...
%             'fontsize',7);'GSHHS\_I (intermediate)'
        hold off
%         colaxis=[0 35];
%         fixcolorbar([0.25 0.05 0.5 0.03],colaxis,...
%             'Topography',10)
        out=[out_path,'model_YSdaily',out_vari,num2str(yid),num2char(mid,4)];
        set(gcf,'renderer','painter')
        if (flag_record_gif)
            mov(kkk) =getframe(gcf);
            if kkk==1
                f=getframe(figure(1));
                [im,map]=rgb2ind(f.cdata,256,'nodither');
            end
            f = getframe(figure(1));
            im(:,:,1,kkk)=rgb2ind(f.cdata,map,'nodither');
        end
        drawnow;
        pause(3);
        saveas(gcf,out,out_type)
     clf;close all
    end
   if(flag_record_gif), imwrite(im,map,[out_path,'daily_',out_vari,'.gif'],'delaytime',0.4,'loopcount',inf); end
end