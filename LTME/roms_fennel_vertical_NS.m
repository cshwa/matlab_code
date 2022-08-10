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
st_line  = 124.5;  %latitude of stations line
max_level= 40;
variation =6;  %  1 = temperature , 2 = salinity ,3 = density , 4 = NO3 , 5 = NH4 , 6 = chla

current   = 0;   skip_v = 3;     size_v = 2;    color_v = 'k';

plot_pcolorjw = 1;   

temp_lim = [5.8 9.8];    salt_lim = [32.2 34.2]; den_lim = [25.2 26.2];
NO3_lim=[0 15]; NH4_lim=[0 1]; chla_lim=[0 5.6];

level_temp =[6:.2:10]; level_c  =[32.1:0.5:34.5];level_dens  =[25.2:0.05:26];

level_NO3=[0:1:15];level_NH4=[0:.1:1];level_chla=[0:1:5.6];


temp_lim = [6 26];    salt_lim = [30.9 33.9]; den_lim = [19 27.4];
NO3_lim=[0 12]; NH4_lim=[0 1]; chla_lim=[0.1 2.5];
level_temp=[8,10,16,21]; level_sal=[32.4:0.5:32.9];
level_dens=[21.2 23.4 26.4];
level_NO3=[2:1:10];level_chla=[0.1:0.5:2.5];

plot_contour  = 1;    color_c  ='-k' ; 

% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'tif';

section = 1; % 0 => whole, 1=> yellowsea, 2=> eastsea

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';

yy = 2009;
start_mm=1;
end_mm=12;
time_step=1;

file_dir=['.\spinup_',num2str(yy),'\monthly_mean_rdrg2_hm6_7\'];

mm=start_mm;
%==========================================================================

gd = grd('auto');
lon_rho  = gd.lon_rho;
lat_rho  = gd.lat_rho; 
mask_rho = gd.mask_rho;
h=gd.h;
N = gd.N;
dist=abs(st_line-lon_rho);
[temp1,num_stline,temp2]=find(dist==min(min(dist)));%%%%%%%%%%수정요망 (거리가 최소가 되게끔)
num_stline=num_stline(1);
depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,2,4,'r');
% angle    = gd{'angle'}(:);
% mask_u = gd{'mask_u'}(:);
% mask_v = gd{'mask_v'}(:);
warning off
mask_rho = mask_rho./mask_rho;
% mask_u = mask_u./mask_u;
% mask_v = mask_v./mask_v;
warning on
vname = {'temp','salt','density','NO3','NH4','chlorophyll'};%,'zeta','ubar','vbar','u','v','omega'};
     
for im=start_mm:time_step:end_mm
		mid=[num2char(im,2)];
        file = [file_dir,'monthly1_',num2str(yy),'_',mid,'.nc'];
        file = [file_dir,'spinup2_monthly_',num2str(yy),'_',mid,'.nc'];
        disp([file,' : ', num2char(im,2)])
            nc=netcdf(file);
            date=[num2str(yy),'. ',num2str(im)];
            switch variation
                case 1
                    value=nc{char(vname(variation))}(:);
                    val_name='Temperature';
                    unit = '^oC';
                    out_name_1=['vertical_NS_Temp-',num2str(yy),'-'];
                    val_caxis=temp_lim;
                    data=squeeze(value(:,:,num_stline));
                    Yi=squeeze(depth(:,:,num_stline));
                    level_c=level_temp;
                   clear value;         

                case 2
                    value=nc{char(vname(variation))}(:);
                    val_name='Salinity';
                    unit = '';
                    out_name_1=['Sal',num2str(yy)];
                    val_caxis=salt_lim;
                    data=squeeze(value(:,:,num_stline));
                    Yi=squeeze(depth(:,:,num_stline));
                    clear value;         

                case 3
                    temp=nc{'temp'}(:);
                    salt=nc{'salt'}(:);
                    val_name='density';
                    unit = 'σ';
                    out_name_1=['Den',num2str(yy)];
                    val_caxis=den_lim;
                    level_c=level_dens;
                    tem=squeeze(temp(:,:,num_stline));
                    sal=squeeze(salt(:,:,num_stline));
data=-1000+999.842594+6.793952.*10.^(-2).*tem-9.095290*10^(-3).*tem.^(2)...
+1.001685*10^(-4).*tem.^(3)-1.120083*10^(-6).*tem.^(4)+6.536332*10^(-9).*tem.^(5)...
+8.24493*10^(-1).*sal-4.0899*10^(-3).*tem.*sal+7.6438*10^(-5).*tem.^2.*sal...
-8.2467*10^(-7).*tem.^3.*sal+5.3875*10^(-9).*tem.^(4).*sal-5.72466*10^(-3).*sal.^(3/2)...
+1.0227*10^(-4).*tem.*sal.^(3/2)-1.6546*10^(-6).*tem.^(2).*sal.^(3/2)+4.8314*10^(-4).*sal.^2; 
                   Yi=squeeze(depth(:,:,num_stline));
               case 4
                    value=nc{char(vname(variation))}(:);
                    val_name='Nitrate';
                    unit = 'μM';
                    out_name_1=['NO3',num2str(yy)];
                    val_caxis=NO3_lim;
                    level_c=level_NO3;
                    data=squeeze(value(:,:,num_stline));
                    Yi=squeeze(depth(:,:,num_stline));
                    clear value;  
                case 5
                    value=nc{char(vname(variation))}(:);
                    val_name='Ammonium';
                    unit = 'μM';
                    out_name_1=['NH4',num2str(yy)];
                    val_caxis=NH4_lim;
                    level_c=level_NH4;
                    data=squeeze(value(:,:,num_stline));
                    Yi=squeeze(depth(:,:,num_stline));
                    clear value;  
                case 6
                    value=nc{char(vname(variation))}(:);
                    val_name='Chlorophyll';
                    unit = 'μg/L';
                    out_name_1=['Chla',num2str(yy)];
                    val_caxis=chla_lim;
                    level_c=level_chla;
                    data=squeeze(value(:,:,num_stline));
                    Yi=squeeze(depth(:,:,num_stline));
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
                domaxis=[35 36 -100 0];       
                x_1=find(lat_rho(:,1)>=domaxis(1));
                x_2=find(lat_rho(:,1)>=domaxis(2));
                x=lat_rho(x_1(1):x_2(1),num_stline);x=repmat(x',40,1);
                data=data(:,x_1(1):x_2(1));
                Yi=Yi(:,x_1(1):x_2(1));
                end
               
            end               
        figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
        set(gca,'Position',[0.2 0.15 0.73 0.75]);
        text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
        text_posi_y1=(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y2=4*(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y3=8*(domaxis(4)-domaxis(3))/20+domaxis(3);
            switch section
                case 1
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                lab_lat=['Longitude',num2str(lon_rho(1,num_stline)),'(^oE)'];
                xlabel('Latitude(^oN)','color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
%                 text(text_posi_x,text_posi_y1,lab_lat,'color',color_v,'FontSize',17,'fontweight','bold')
%                 text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
%                 text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
%                 title(val_name,'fontsize',17);
                out_name_1=[out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
            end
%             caxis(val_caxis);
            bar = colorbar('fontsize',17,'fontweight','bold');
            set(get(bar,'title'),'string',unit,'FontSize',17,'fontweight','bold')
            
            out_name=[file_dir,num2str(floor(st_line)),out_name_1,num2char(im,2)];

            if (switch_save)
                saveas(gcf,out_name,out_type);
            end
            close all
end
%      close all
% % % end