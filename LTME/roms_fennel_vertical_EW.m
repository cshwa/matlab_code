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
clear all;close all;clc
%==========================================================================

% 정선 307 - 36.925  // 308 - 36.33 // 309 - 35.855 // 310 - 35.335  // 311 -
% 34.716 // 312 - 33.975 ~ 34.0917 //313 - 33.4067 // 314 - 33
st_line=25;

section =1 ;  % 1 =yellow sea   //  2 =east sea

max_level= 40;

variation =6;  %  1 = temperature , 2 = salinity , 3 = V velocity , 4 = E velocity , 5 = density , 6 = NO3 , 7 = NH4 , 8 = chla

current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

plot_pcolorjw = 1; plot_contour  = 1;    color_c  ='-k' ;  

temp_lim = [0 30]; level_tem =[0:1:30]; 
salt_lim = [30 35]; level_sal=[30:0.2:35];
vel_lim=[-5 5];   level_vel=[-10:1:10];
vele_lim=[-10 10];  level_vele=[-10:2:10];
den_lim= [23 28];  level_den =[25.2:0.05:26];
NO3_lim= [0 10];   level_NO3=[0:1:10];
NH4_lim= [0 1];   level_NH4=[0:0.1:1];
chla_lim= [0 5.6];   level_chla=[0:1:6];
% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

temp_lim = [6 12];    salt_lim = [30.9 33.9]; den_lim = [19 27.4];
NO3_lim=[0 12]; NH4_lim=[0 1]; chla_lim=[0.1 2.5];
level_tem=[6:1:26]; level_sal=[32.4:0.5:32.9];
den_lim= [19 25.4];  level_den=[19:1:25];
level_NO3=[2:2:12];level_chla=[0.1:0.5:2.5];

switch_save   = 1;    out_type = 'tif';

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';
yy = 2009;
file_dir=['.\spinup_',num2str(yy),'\monthly_mean_rdrg2_hm6_3\'];
start_mm=1;
end_mm=12;
time_step=1;

mm=start_mm;
%==========================================================================

                
switch section
    case 1
    domaxis=[122 125.6 -100 0];
%     domaxis=[124 128 -100 0]; 
%     domaxis=[124.5 125.7 -100 0];
    domaxis=[121 123 -250 0];
    case 2
    domaxis=[129.5 131.2 -500 0];    
end

gd = grd('auto');
lon_rho  = gd.lon_rho;
lat_rho  = gd.lat_rho; 
mask_rho = gd.mask_rho;
h=gd.h;
N = gd.N;
depth=zlevs(gd.h,gd.zeta,gd.theta_s,gd.theta_b,gd.Tcline,N,2,4,'r');
% hz=gd.z_w(2:N+1,:,:)-gd.z_w(1:N,:,:); 
dist=abs(st_line-lat_rho);
[num_stline,temp1,temp2]=find(dist==min(min(dist)));   %%수정요망 (거리가 최소가 되게끔)
num_stline=num_stline(1);

% for level=1:1:N
%     depth(level,:,:)=gd.z_w(level+1,:,:)-(hz(level,:,:)/2);
% end

% angle    = gd{'angle'}(:);
% mask_u = gd{'mask_u'}(:);
% mask_v = gd{'mask_v'}(:);
warning off
mask_rho = mask_rho./mask_rho;
% mask_u = mask_u./mask_u;
% mask_v = mask_v./mask_v;
warning on
vname = {'temp','salt','v','u','NO3','NH4','chlorophyll'};%,'zeta','ubar','vbar','u','v','omega'};
str_month={'Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'};
for im=start_mm:time_step:end_mm
		mid=[num2char(im,2)];
        file = [file_dir,'monthly1_',num2str(yy),'_',mid,'.nc'];
        file = [file_dir,'spinup3_monthly_',num2str(yy),'_',mid,'.nc'];
        disp([file,' : ', num2char(im,2)])
        nc=netcdf(file);
        date=[str_month{im},' ',num2str(yy)];
        zeta=nc{'zeta'}(:);
        switch variation
            case 1
                value=nc{char(vname(variation))}(:);
                val_name='Temperature';
                unit = '^oC';
                out_name_1=['vertical_Temp-',num2str(yy),'-'];
                val_caxis=temp_lim;
                level_c=level_tem;
                data=squeeze(value(:,num_stline,:)); 
            case 2
                value=nc{char(vname(variation))}(:);
                val_name='Salinity';
                unit = 'psu';
                out_name_1=['Sal',num2str(yy)];
                val_caxis=salt_lim;
                level_c=level_sal;
                data=squeeze(value(:,num_stline,:)); 
             case 3
                value=nc{char(vname(variation))}(:);
                val_name= ' '; % 'N-ward Velocity';
                unit = 'cm/s';
                out_name_1=['Nvel',num2str(yy)];
                val_caxis=vel_lim;
                level_c=level_vel;
                data=squeeze(value(:,num_stline,:))*100; 
             case 4
                value=nc{char(vname(variation))}(:);
                val_name='E-ward Velocity';
                unit = 'cm/s';
                out_name_1=['Evel',num2str(yy)];
                val_caxis=vele_lim;
                level_c=level_vele;
                data=squeeze(value(:,num_stline,:))*100;                 
             case 5
                temp=nc{'temp'}(:);
                salt=nc{'salt'}(:);
                val_name='Density';
                unit = 'σ';
                out_name_1=['Den',num2str(yy)];
                val_caxis=den_lim;
                level_c=level_den;
                data=sw_dens(squeeze(salt(:,num_stline,:)),squeeze(temp(:,num_stline,:)),0)-1000;
            case 6
                value=nc{char(vname(variation-1))}(:);
                val_name='Nitrate';
                unit = 'μM';
                out_name_1=['NO3',num2str(yy)];
                val_caxis=NO3_lim;
                level_c=level_NO3;
                data=squeeze(value(:,num_stline,:)); 
            case 7
                value=nc{char(vname(variation-1))}(:);
                val_name='Ammonium';
                unit = 'μM';
                out_name_1=['NH4',num2str(yy)];
                val_caxis=NH4_lim;
                level_c=level_NH4;
                data=squeeze(value(:,num_stline,:)); 
            case 8
                value=nc{char(vname(variation-1))}(:);
                val_name='Chlorophyll';
                unit = 'μg/L';
                out_name_1=['Chla',num2str(yy)];
                val_caxis=chla_lim;
                level_c=level_chla;
                data=squeeze(value(:,num_stline,:));             
        end
        Yi=squeeze(depth(:,num_stline,:));
        clear value;         

        if (plot_pcolorjw)
            for i=1:1:length(data(:,1))
                for j=1:1:length(data(1,:))
                    if data(i,j) > 100
                        data(i,j) = NaN;
                    end
                end
            end                       
            x_1=find(lon_rho(1,:)>=domaxis(1));
            x_2=find(lon_rho(1,:)>=domaxis(2));
            x=lon_rho(num_stline,x_1(1):x_2(1));x=repmat(x,40,1);
            data=data(:,x_1(1):x_2(1));
            Yi=Yi(:,x_1(1):x_2(1));
        end               
        
        figure('position',[87 200 500 500],'PaperUnits','inches','PaperPosition',[1.22 0.66 8 6.02]);
        set(gca,'Position',[0.13 0.11 0.825 0.81]);
        text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
        text_posi_y1=(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y2=3*(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y3=6*(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y21=21*(domaxis(4)-domaxis(3))/20+domaxis(3);
            switch section
                case 0
                hold on
                pcolor(Xi,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)]);
                shading interp;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                lab_lat=['Latitude',num2str(lat_rho(num_stline,1)),'(^oN)'];
                xlabel('Longitude(^oE)','color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y1,lab_lat,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
%                 text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,-Yi,Zi,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(x,-Yi,Zi,level_c-1,color_c,'linewidth',1);
                  clabel(C,h,'FontSize',10,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
                case 1
                hold on
                pcolor(x,Yi,data);
                shading interp;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)]);
                lab_lat=['Latitude ',num2str(st_line),'^oN'];
                xlabel('Longitude (^oE)','color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth (m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y1,lab_lat,'color',color_v,'FontSize',22,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',20,'fontweight','bold')
%                 text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
%                 text(text_posi_x,text_posi_y21,'(c)','color',color_v,'FontSize',22,'fontweight','bold')
                out_name_1=['YellowSea',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);
%                   [C2,h2]=contour(x,Yi,data,level_c+1,color_c,'linewidth',1);
                  clabel(C,h,'Rotation',0,'Color','k','fontsize',24);

%                   clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
                case 2 
                hold on
                pcolor(x,Yi,data)
                shading interp;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)])
                lab_lat=['Latitude',num2str(st_line),'(^oN)'];
                xlabel('Longitude (^oE)','color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y1,lab_lat,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
                out_name_1=['EastSea',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);
%                   [C2,h2]=contour(x,Yi,data,level_c-1,color_c,'linewidth',1);
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
            end
            if variation ==3 ||variation ==4
                color=cbrewer('div','RdBu',50);color=flipud(color);
                colormap(color)
            end
            bar = colorbar('fontsize',22,'fontweight','bold');
            set(get(bar,'title'),'string',unit,'FontSize',22,'fontweight','bold')
%             title(['CONTROL  ',date],'fontsize',18)            
            out_name=[file_dir,num2str(floor(st_line)),out_name_1,num2char(im,2)];

            if (switch_save)
                saveas(gcf,out_name,out_type);
            end
            close all
end
%      close all
% % % end