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
addpath(genpath('/home/yjtak/roms/matlab'));
addpath(genpath('/home/yjtak/Dropbox/source/matlab/Common/m_map'));
%==========================================================================
max_level=40;
variation = 4;  %  1 = temperature , 2 = salinity, 3 = u 4 = v

Vtransform=2;
Vstretching=4;

current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

temp_lim = [6 28];    temp_c  =[6:2:28];
salt_lim = [32 35];   salt_c  =[26:1:35];
u_lim = [-.7 .7];    u_c  =[-.6:.2:.6];
v_lim = [-.7 .7];    v_c  =[-.6:.2:.6];

plot_contour  = 1;    color_c  ='-k' ;      

% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'tif';

section = 1; % 0 => whole, 1=> kuroshio, 2=> Korean Strait
bndy={'_south','_east'};

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';

yy = 2010;
start_mm=1;
end_mm=12;
time_step=1;

file_dir=['./'];

mm=start_mm;
%==========================================================================

gd = grd('auto');
lon  = gd.lon_rho;
lat  = gd.lat_rho; 
mask_rho = gd.mask_rho;
h=gd.h;
N = gd.N;
depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.Tcline,N,Vtransform,Vstretching,'r');
mask_rho = mask_rho./mask_rho;

warning on
vname = {'temp','salt','u','v'};%,'zeta','ubar','vbar','u','v','omega'};
file = [file_dir,'roms_bndy_auto_NWP_1_10_test6_clim.nc'];
for im=start_mm:time_step:end_mm
	switch section
    case 1
        maskusouth=repmat((gd.mask_u(1,:)./gd.mask_u(1,:)),[N,1]);
        %maskusouth=permute(maskusouth,[2,1]);
        maskvsouth=repmat((gd.mask_v(1,:)./gd.mask_v(1,:)),[N,1]);
        %maskvsouth=permute(maskvsouth,[2,1]);
        maskrsouth=repmat((gd.mask_rho(1,:)./gd.mask_rho(1,:)),[N,1]);
        %maskrsouth=permute(maskrsouth,[2,1]);
        cos_angle2D=repmat( cos(gd.angle(1,:)) ,[1 N ]);
        cos_angle2D=permute(cos_angle2D, [2 1]);
        sin_angle2D=repmat( sin(gd.angle(1,:)) ,[1 N ]);
        sin_angle2D=permute(sin_angle2D, [2 1]);
        Yi=squeeze(depth(:,1,:));
    case 2
        maskueast=repmat((gd.mask_u(:,end)./gd.mask_u(:,end)),[1,N]);
        maskueast=permute(maskueast,[2,1]);
        maskveast=repmat((gd.mask_v(:,end)./gd.mask_v(:,end)),[1,N]);
        maskveast=permute(maskveast,[2,1]);
        cos_angle2D=repmat( cos(gd.angle(:,1)) ,[1 N ]);
        cos_angle2D=permute(cos_angle2D, [2 1]);
        sin_angle2D=repmat( sin(gd.angle(:,1)) ,[1 N ]);
        sin_angle2D=permute(sin_angle2D, [2 1]);
        Yi=squeeze(depth(:,:,end));
    end
	
    mid=[num2char(im,2)];
    disp([file,' : ', num2char(im,2)])
    nc=netcdf(file);
    date=[num2str(yy),'. ',num2str(im)];
    switch variation
        case 1
            value=nc{char([vname{variation},bndy{section}])}(im,:,:);
            val_name='Temperature';
            unit = '^oC';
            out_name_1=['vert_bndy_Temp',bndy{section},num2char(im,2),'-'];
            val_caxis=temp_lim;
            level_c=temp_c;
            data=value.*maskrsouth;
           clear value;         
        case 2
            value=nc{char([vname{variation},bndy{section}])}(im,:,:);
            val_name='Salinity';
            unit = 'psu';
            out_name_1=['vert_bndy_Salt',bndy{section},num2char(im,2),'-'];
            val_caxis=salt_lim;
            level_c=salt_c;
            data=value.*maskrsouth;
           clear value;     
        case 3
            u=nc{char([vname{variation},bndy{section}])}(im,:,:);
            v=nc{char([vname{variation+1},bndy{section}])}(im,:,:);
            val_name='East.Vel.';
            unit = 'm s^{-1}';
            out_name_1=['vert_bndy_Uvel',bndy{section},num2char(im,2),'-'];
            val_caxis=u_lim;
            level_c=u_c;
            lat=gd.lat_u;
            lon=gd.lon_u;
            v=[v,v(:,end)];
            true_u=cos_angle2D.*u - sin_angle2D.*v;
            true_v=sin_angle2D.*u + cos_angle2D.*v;
            data=true_u.*maskueast;
           clear value;    
        case 4
            u=nc{char([vname{variation-1},bndy{section}])}(im,:,:);
            v=nc{char([vname{variation},bndy{section}])}(im,:,:);
            val_name='Nor.Vel';
            unit = 'm s^{-1}';
            out_name_1=['vert_bndy_Vvel',bndy{section},num2char(im,2),'-'];
            val_caxis=v_lim;
            level_c=v_c;
            lat=gd.lat_v;
            lon=gd.lon_v;
            u=u(:,1:length(v));
            cos_angle2D=cos_angle2D(:,1:length(v));
            sin_angle2D=sin_angle2D(:,1:length(v));
            true_u=cos_angle2D.*u - sin_angle2D.*v;
            true_v=sin_angle2D.*u + cos_angle2D.*v;
            data=true_v.*maskveast;
            clear value;    
    end


        for i=1:1:length(data(:,1))
            for j=1:1:length(data(1,:))
                if data(i,j) > 10000
                    data(i,j) = NaN;
                end
            end
        end
        switch section
            case 1
                domaxis=[117.6 130 -1000 0];       
                x_1=find(lon(1,:)>=domaxis(1));
                x_2=find(lon(1,:)>=domaxis(2));
                x=lon(1,x_1(1):x_2(1));x=repmat(x,40,1);
                data=data(:,x_1(1):x_2(1));
                Yi=Yi(:,x_1(1):x_2(1));
            case 2
                domaxis=[33 37 -150 0];       
                x_1=find(lat(:,end)>=domaxis(1));
                x_2=find(lat(:,end)>=domaxis(2));
                x=lat(x_1(1):x_2(1),end);x=repmat(x',40,1);
                data=data(:,x_1(1):x_2(1));
                Yi=Yi(:,x_1(1):x_2(1));
        end

            
        figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
        set(gca,'Position',[0.2 0.15 0.73 0.75]);
        text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
        text_posi_y1=(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y2=2*(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y3=3*(domaxis(4)-domaxis(3))/20+domaxis(3);
            switch section
                case 1
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel('longitudee(^{o}E','color',color_v,'FontSize',17,'fontweight','bold')
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
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel('latitude(^{o}N','color',color_v,'FontSize',17,'fontweight','bold')
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
            end
%             caxis(val_caxis);
            bar = colorbar('fontsize',17,'fontweight','bold');
            colormap('jet')
            set(get(bar,'title'),'string',unit,'FontSize',17,'fontweight','bold')
            
            out_name=[file_dir,'re8',out_name_1,num2char(im,2)];

            if (switch_save)
                saveas(gcf,out_name,out_type);
            end
            close all
end
%      close all
% % % end