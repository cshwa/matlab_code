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
max_level=40;
variation =2;  %  1 = temperature , 2 = salinity, 3 =NO3',4 = NH4', 5= 'chlorophyll'
Vtransform=2;
Vstretching=4;

current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

plot_pcolorjw = 1;    
temp_lim = [6 15];    temp_c  =[6:2:28];
salt_lim = [31 34.5];   salt_c  =[26:.2:35];
NO3_lim = [0 10];    NO3_c  =[0:2:10];
NH4_lim = [0 1];    NH4_c  =[0:.1:1];
chlo_lim = [0 5];   chlo_c=[0:1:5];

plot_contour  = 1;    color_c  ='-k' ;      

% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'tif';

section = 6; % 0 => whole, 1=> yellowsea, 2=> Korean Strait 3=> tawain island 4 => southsea 5=> YellowSea trough
             % 6 => YSWC path (Liu et al.,2015)    

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';

yy = 2009;
start_mm=1;
end_mm=12;
time_step=1;

file_dir=['.\spinup_',num2str(yy),'\monthly_mean_rdrg2_hm6_8\'];

mm=start_mm;
%==========================================================================

gd = grd('auto');
lon_rho  = gd.lon_rho;
lat_rho  = gd.lat_rho; 
mask_rho = gd.mask_rho;
h=gd.h;
N = gd.N;
depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,Vtransform,Vstretching,'r');
% angle    = gd{'angle'}(:);
% mask_u = gd{'mask_u'}(:);
% mask_v = gd{'mask_v'}(:);
warning off
mask_rho = mask_rho./mask_rho;
% mask_u = mask_u./mask_u;
% mask_v = mask_v./mask_v;
warning on
vname = {'temp','salt','NO3','NH4','chlorophyll'};%,'zeta','ubar','vbar','u','v','omega'};
     
for im=start_mm:time_step:end_mm
		mid=[num2char(im,2)];
        file = [file_dir,'monthly1_',num2str(yy),'_',mid,'.nc'];
        file = [file_dir,'spinup3_monthly_',num2str(yy),'_',mid,'.nc'];
        disp([file,' : ', num2char(im,2)])
            nc=netcdf(file);
            date=[num2str(yy),'. ',num2str(im)];
            switch variation
                case 1
                    value=nc{char(vname(variation))}(:);
                    val_name='Temperature';
                    unit = '^oC';
                    out_name_1=['vert_diag_Temp-',num2str(yy),'-'];
                    val_caxis=temp_lim;
                    level_c=temp_c;
                    data=value;
                   clear value;         
                case 2
                    value=nc{char(vname(variation))}(:);
                    val_name='Salinity';
                    unit = 'psu';
                    out_name_1=['vertical_NS_Salt-',num2str(yy),'-'];
                    val_caxis=salt_lim;
                    level_c=salt_c;
                    data=value;
                   clear value;     
                case 3
                    value=nc{char(vname(variation))}(:);
                    val_name='Nitrate';
                    unit = 'mMole N m^{-3}';
                    out_name_1=['vert_diag_NO3',num2str(yy)];
                    val_caxis=NO3_lim;
                    level_c=NO3_c;
                    data=value;
                   clear value;    
                case 4
                    value=nc{char(vname(variation))}(:);
                    val_name='Ammonium';
                    unit = 'mMole N m^{-3}';
                    out_name_1=['vert_diag_NH4',num2str(yy)];
                    val_caxis=NH4_lim;
                    level_c=NH4_c;
                    data=value;
                   clear value;    
                case 5
                    value=nc{char(vname(variation))}(:);
                val_name='Chlorophyll';
                unit = '¥ìg/L';
                out_name_1=['vert_diag_chlo',num2str(yy)];
                    val_caxis=chlo_lim;
                    level_c=chlo_c;
                    data=value;
                   clear value;                       
            end
   

            if (plot_pcolorjw)
                for i=1:1:length(data(:,1))
                    for j=1:1:length(data(1,:))
                        if data(i,j) > 10000
                            data(i,j) = NaN;
                        end
                    end
                end
                
                switch section
                    case 1
                    domaxis=[36.0302 36.0302 129.422 129.431 -120 0]; 
                    case 2
                    domaxis=[ 35 33.7 128.78 130.1 -150 0];
%                     domaxis=[33.4263348468  35.2833167966 130.094555939 128.883481062 -150 0];
                    case 3
                    domaxis=[ 25.5 24.5 121. 124 -300 0];
                    case 4
                    domaxis=[ 35 33 127 128 -300 0];
                    case 5
                    domaxis=[ 34 32 124.5 127 -150 0];
                    case 6
                    domaxis=[ 32 35.5 128 120.5 -120 0];
                end
                dist=sqrt((lon_rho-domaxis(3)).^2+(lat_rho-domaxis(1)).^2);
                min_dist=min(min(dist));
                dist2=sqrt((lon_rho-domaxis(4)).^2+(lat_rho-domaxis(2)).^2);
                min_dist2=min(min(dist2));                
                [x1,y1]=find(dist==min_dist);
                [x2,y2]=find(dist2==min_dist2);                
                lat1=lat_rho(x1(1),y1(1));lon1=lon_rho(x1(1),y1(1));
                lat2=lat_rho(x2(1),y2(1));lon2=lon_rho(x2(1),y2(1));
                if (lon2-lon1) >= (lat2-lat1)
                    lon_line=[lon1:0.001:lon2];
                    lat_line=(lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1;
                    x=repmat(lon_line,gd.N,1);
                    x_label='Longitude(^oE)';
                    domaxis=[domaxis(3) domaxis(4) domaxis(5) domaxis(6) domaxis(5) domaxis(6)];
                else
                    lat_line=[min(lat1,lat2):0.001:max(lat1,lat2)];
                    lon_line=(lat_line-lat1)*((lon2-lon1)/(lat2-lat1))+lon1;
                    x=repmat(lat_line,gd.N,1);
                    x_label='Latitude(^oN)';
                end
                Temp=zeros(gd.N,length(lat_line));
                for k=1:1:gd.N
                    lon_range=lon_rho(min(x1,x2):max(x1,x2),min(y1,y2):max(y1,y2));
                    lat_range=lat_rho(min(x1,x2):max(x1,x2),min(y1,y2):max(y1,y2));
                    data_range=squeeze(data(k,min(x1,x2):max(x1,x2),min(y1,y2):max(y1,y2)));
                    depth_range=squeeze(depth(k,min(x1,x2):max(x1,x2),min(y1,y2):max(y1,y2)));
                    Temp(k,:)=griddata(lon_range,lat_range,data_range,lon_line,lat_line);
                    Yi(k,:)=griddata(lon_range,lat_range,depth_range,lon_line,lat_line);
                end
                data=Temp;
            end    
            
        figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
        set(gca,'Position',[0.2 0.15 0.73 0.75]);
        text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
        text_posi_y1=(domaxis(6)-domaxis(5))/20+domaxis(5);
        text_posi_y2=2*(domaxis(6)-domaxis(5))/20+domaxis(5);
        text_posi_y3=3*(domaxis(6)-domaxis(5))/20+domaxis(5);
            switch section
                case 1
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(5) domaxis(6)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel(x_label,'color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
%                 title('Data assimilation','fontsize',17);
                out_name_1=['YellowSea',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);                
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
                case 2
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(5) domaxis(6)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel(x_label,'color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
%                 title('Data assimilation','fontsize',17);
                out_name_1=['Korean_Strait',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);                
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
                case 3
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(5) domaxis(6)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel(x_label,'color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
%                 title('Data assimilation','fontsize',17);
                out_name_1=['Taiwan',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);                
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
                case 4
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(5) domaxis(6)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel(x_label,'color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
%                 title('Data assimilation','fontsize',17);
                out_name_1=['southsea',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);                
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
                case 5
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(5) domaxis(6)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel(x_label,'color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
%                 title('Data assimilation','fontsize',17);
                out_name_1=['YST',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);                
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');                  
                end
                case 6
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(5) domaxis(6)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel(x_label,'color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
%                 title('Data assimilation','fontsize',17);
                out_name_1=['YSWC',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);                
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');      
                end
            end
%             caxis(val_caxis);
            bar = colorbar('fontsize',17,'fontweight','bold');
            set(get(bar,'title'),'string',unit,'FontSize',17,'fontweight','bold')
            
            out_name=[file_dir,out_name_1,num2char(im,2)];

            if (switch_save)
                saveas(gcf,out_name,out_type);
            end
            close all
end
%      close all
% % % end