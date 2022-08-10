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
max_level= 30;
out_filen = 'mp_p_sewer_det_10m_do.png';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variation =5; %  1 = temperature // 2 = salinity // 3 = NO3 // 4 = NH4 // 5 = chlorophyll
              %  6 = phytoplankton // 7 = zooplankton // 8 = total inorganic carbon  
              %  9 = oxygen // 10 = PO4 // 11 = Density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

% temp_lim = [6 25];
plot_pcolorjw = 1;    temp_lim = [0 35];    salt_lim = [28 33]; den_lim = [20 26]; no3_lim=[0 20]; nh4_lim=[0 15];
                      chl_lim=[0 10];   phy_lim=[0 3];    zoo_lim=[0 3];  oxy_lim=[0 14];
                      tic_lim=[-70 200]; po4_lim=[0 1];
% temp_c  =[6:1:25];
plot_contour  = 1;    color_c  ='-k' ;      temp_c  =[0:1:35];salt_c  =[20:.2:35]; no3_c=[0:2:200]; nh4_c=[0:.5:30];
                      chl_c=[0:.2:20];    phy_c=[0:.1:10];    zoo_c=[0:0.5:3];        oxy_c=[-200:.5:300];
                      tic_c=[-50:50:200];po4_c=[-0.1:0.1:1];  den_c=[20:1:26];

% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'tif';

section = 1; % 0 => whole, 1=> yellowsea, 2=> eastsea

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';

yy = 2001;  
start_mm=1;   %04.29 - 119 08.05 - 217  09.04-247   12.01-335    01.12-377
end_mm=12;
time_step=0;

%% 1.4.8. == 5/29 >? 1.4.9 .
%% 8.1. == 3/23 >? .8.2.
%% 7.2. == 3/14 >? 7.3.
%% 2.2.6. == 8/14
% file_dir=['\H:\E\포항모델\model\output\'];

%==========================================================================

% gd = grd('pohang_fine');
% gd = grd('pump_30_close');
fpath='D:\장기생태\Dynamic\result\2001\';
lon_rho  = ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'lon_rho');
lat_rho  = ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'lat_rho');
mask_rho = ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'mask_rho');
h=ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'h');
N = length(ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'s_rho'));
% depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,1,1,'r');
theta_s=ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'theta_s');
theta_b=ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'theta_b');
zeta=ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'zeta');
hc=ncread([fpath,'mp_p_sewer_det_monthly_2001_01.nc'],'hc');
depth=zlevs(h,zeta,theta_s,theta_b,hc,N,'r');
% angle    = gd{'angle'}(:);
% mask_u = gd{'mask_u'}(:);
% mask_v = gd{'mask_v'}(:);
warning off
mask_rho = mask_rho./mask_rho;
% mask_u = mask_u./mask_u;
% mask_v = mask_v./mask_v;
warning on
vname = {'temp','salt' ,'NO3','NH4','chlorophyll','phytoplankton','zooplankton','TIC','oxygen' 'tPO4'};%,'zeta','ubar','vbar','u','v','omega'};
     
% for im=start_mm:time_step:end_mm
for im=start_mm:end_mm %pohang
		mid=[num2str(im,'%02d')];
%         file = [file_dir,'ocean_avg_',mid,'.nc'];
file = [fpath,'mp_p_sewer_det_monthly_2001_',mid,'.nc'];
% file = ['pohang_avg_',mid,'.nc'];
%         file = ['pohang_wt8_fine_flat_3m_',mid,'.nc'];
        disp([file,' : ', num2str(im,'%02d')])
        nc=netcdf(file);
        date=[num2str(yy),'. ',num2str(im,'%02d')];
        switch variation
            case 1
                value=nc{char(vname(variation))}(:);
                val_name='Temperature';
                unit = '^oC';
                out_name_1=['vertical-GY-Temp-',num2str(yy),'-'];
                val_caxis=temp_lim;
                level_c=temp_c;
                data=value;
               clear value;         
            case 2
                value=nc{char(vname(variation))}(:);
                val_name='Salinity';
                unit = 'psu';
                out_name_1=['vertical-GY-Salt-',num2str(yy),'-'];
                val_caxis=salt_lim; 
                level_c=salt_c;
                data=value;
               clear value;     
           case 3
                value=nc{char(vname(variation))}(:);
                val_name='Nitrate';
                unit = 'μMol N L';
                out_name_1=['vertical-GY-NO3-',num2str(yy),'-'];
                val_caxis=no3_lim;
                level_c=no3_c;
                data=value; %.*14./1000; %% convert to mg/L
               clear value;
           case 4
                value=nc{char(vname(variation))}(:);
                val_name='Ammonium';
                unit = 'μMol N L';
                out_name_1=['vertical-GY-NH4-',num2str(yy),'-'];
                val_caxis=nh4_lim;
                level_c=nh4_c;
                data=value; %.*14./1000; %% convert to mg/L
               clear value;              
           case 5
                value=nc{char(vname(variation))}(:);
                val_name='chlorophyll';
                unit = 'μg C L';
                out_name_1=['vertical-GY-chl-',num2str(yy),'-'];
                val_caxis=chl_lim;
                level_c=chl_c;
                data=value;
               clear value;  
           case 6
                value=nc{char(vname(variation))}(:);
                val_name='phytoplankton';
                unit ='μMol N L';
                out_name_1=['vertical-GY-phy-',num2str(yy),'-'];
                val_caxis=phy_lim;
                level_c=phy_c;
                data=value;
               clear value;
           case 7
                value=nc{char(vname(variation))}(:);
                val_name='zootoplankton';
                unit ='μMol N L';
                out_name_1=['vertical-GY-zoo-',num2str(yy),'-'];
                val_caxis=zoo_lim;
                level_c=zoo_c;
                data=value;
               clear value;
           case 8
                value=nc{char(vname(variation))}(:);
                val_name='TIC';
                unit ='μMol C L';
                out_name_1=['vertical-GY-tic-',num2str(yy),'-'];
                val_caxis=tic_lim;
                level_c=tic_c;
                data=value;
               clear value;
           case 9
                value=nc{char(vname(variation))}(:);
                val_name='dissolved oxygen';
                unit ='mg/L';
                out_name_1=['vertical-GY-DO-',num2str(yy),'-'];
                val_caxis=oxy_lim;
                level_c=oxy_c;
                data=value; %.*(10/7).*(1/44.661); %% convert to mg/L 
               clear value;
           case 10
                value=nc{char(vname(variation))}(:);
                val_name='phosphate';
                unit ='μMol P L';
                out_name_1=['vertical-GY-PO4-',num2str(yy),'-'];
                val_caxis=po4_lim;
                level_c=po4_c;
                data=value;
               clear value;   
            case 11
                temp=nc{char(vname(1))}(:);salt=nc{char(vname(2))}(:);
                for N=1:1:20
                    value(N,:,:)=sw_dens(squeeze(salt(N,:,:)),squeeze(temp(N,:,:)),0)-1000;
                end
                val_name='density';
                unit = 'σ(kg/m^3)';
                out_name_1=['vertical-GY-Den-',num2str(yy),'-'];
                val_caxis=den_lim;
                level_c=den_c;
                data=value;
        end
   

            if (plot_pcolorjw)
                for i=1:1:length(data(:,1))
                    for j=1:1:length(data(1,:))
                        if data(i,j) > 10000
                            data(i,j) = NaN;
                        end
                    end
                end
            end    
   


dep=permute(depth,[2 3 1]);
data_3 = permute(data,[3 2 1]);
pick = 93;

dep(find(dep > 10000))=NaN;

%%% check line
% figure; hold on; 
% pcolor(lon_rho,lat_rho,squeeze(data_3(:,:,end)))
% shading flat; hold on; 
% plot(lon_rho(93,:),lat_rho(93,:),'r--','linew',2)
%%%%%%

for i= 1:N
    lat_rho_3(:,:,i) = lat_rho;
end
%         figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
%         set(gca,'Position',[0.2 0.15 0.73 0.75]);
%         figure('position',[400 100 550
%         550],'PaperUnits','inches','PaperPosition',[0 0 7.5 3.5]);
%         set(gca,'Position',[0.1 0.17 0.88 0.73]);
        set(gcf,'Position',[200 100 1200 700])
        set(gcf,'PaperPositionMode','auto')     
        pcolor(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(data_3(pick,:,:)))
        xlim([34.67 34.945]);
        shading flat;
        set(gca,'box','on','linewidth',1.5,'fontsize',13)
        xlabel('lat','color',color_v,'FontSize',30,'fontweight','bold')
        ylabel('Depth(m)','color',color_v,'FontSize',30,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',25,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',30,'fontweight','bold')
        out_name=[out_name_1,num2str(im,'%02d')];
        title(out_name)
        caxis(val_caxis)
        grid on;
        if (plot_contour)
                  hold on
                  [C,h]=contour(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(data_3(pick,:,:)),level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(data_3(pick,:,:)),[-1:2:1],'-w','linewidth',1);              
                  clabel(C,h,'FontSize',12,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
        end
            
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
                saveas(gca,[out_name,'.png']);
                close all;
            end
            
end
return

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
max_level= 30;
out_filen = 'mp_p_sewer_det_10m_do.png';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variation =5; %  1 = temperature // 2 = salinity // 3 = NO3 // 4 = NH4 // 5 = chlorophyll
              %  6 = phytoplankton // 7 = zooplankton // 8 = total inorganic carbon  
              %  9 = oxygen // 10 = PO4 // 11 = Density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

% temp_lim = [6 25];
plot_pcolorjw = 1;    temp_lim = [0 35];    salt_lim = [28 33]; den_lim = [20 26]; no3_lim=[0 20]; nh4_lim=[0 15];
                      chl_lim=[0 10];   phy_lim=[0 3];    zoo_lim=[0 3];  oxy_lim=[0 14];
                      tic_lim=[-70 200]; po4_lim=[0 1];
% temp_c  =[6:1:25];
plot_contour  = 1;    color_c  ='-k' ;      temp_c  =[0:1:35];salt_c  =[20:.2:35]; no3_c=[0:2:200]; nh4_c=[0:.5:30];
                      chl_c=[0:.2:20];    phy_c=[0:.1:10];    zoo_c=[0:0.5:3];        oxy_c=[-200:.5:300];
                      tic_c=[-50:50:200];po4_c=[-0.1:0.1:1];  den_c=[20:1:26];
                      
% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'tif';

section = 1; % 0 => whole, 1=> yellowsea, 2=> eastsea

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';

yy = 2001;  
start_mm=1;   %04.29 - 119 08.05 - 217  09.04-247   12.01-335    01.12-377
end_mm=12;
time_step=0;

%% 1.4.8. == 5/29 >? 1.4.9 .
%% 8.1. == 3/23 >? .8.2.
%% 7.2. == 3/14 >? 7.3.
%% 2.2.6. == 8/14
% file_dir=['\H:\E\포항모델\model\output\'];

%==========================================================================

% gd = grd('pohang_fine');
% gd = grd('pump_30_close');
fpath='D:\장기생태\Dynamic\result\2001\';
lon_rho  = ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'lon_rho');
lat_rho  = ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'lat_rho');
mask_rho = ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'mask_rho');
h=ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'h');
N = length(ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'s_rho'));
% depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,1,1,'r');
theta_s=ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'theta_s');
theta_b=ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'theta_b');
zeta=ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'zeta');
hc=ncread([fpath,'mp_p_sewer_det_wt8_monthly_2001_01.nc'],'hc');
depth=zlevs(h,zeta,theta_s,theta_b,hc,N,'r');
% angle    = gd{'angle'}(:);
% mask_u = gd{'mask_u'}(:);
% mask_v = gd{'mask_v'}(:);
warning off
mask_rho = mask_rho./mask_rho;
% mask_u = mask_u./mask_u;
% mask_v = mask_v./mask_v;
warning on
vname = {'temp','salt' ,'NO3','NH4','chlorophyll','phytoplankton','zooplankton','TIC','oxygen' 'tPO4'};%,'zeta','ubar','vbar','u','v','omega'};
     
% for im=start_mm:time_step:end_mm
for im=start_mm:end_mm %pohang
		mid=[num2str(im,'%02d')];
%         file = [file_dir,'ocean_avg_',mid,'.nc'];
file = [fpath,'mp_p_sewer_det_wt8_monthly_2001_',mid,'.nc'];
% file = ['pohang_avg_',mid,'.nc'];
%         file = ['pohang_wt8_fine_flat_3m_',mid,'.nc'];
        disp([file,' : ', num2str(im,'%02d')])
        nc=netcdf(file);
        date=[num2str(yy),'. ',num2str(im,'%02d')];
        switch variation
            case 1
                value=nc{char(vname(variation))}(:);
                val_name='Temperature';
                unit = '^oC';
                out_name_1=['vertical-GY-wt8-Temp-',num2str(yy),'-'];
                val_caxis=temp_lim;
                level_c=temp_c;
                data=value;
               clear value;         
            case 2
                value=nc{char(vname(variation))}(:);
                val_name='Salinity';
                unit = 'psu';
                out_name_1=['vertical-GY-wt8-Salt-',num2str(yy),'-'];
                val_caxis=salt_lim; 
                level_c=salt_c;
                data=value;
               clear value;     
           case 3
                value=nc{char(vname(variation))}(:);
                val_name='Nitrate';
                unit = 'μMol N L';
                out_name_1=['vertical-GY-NO3-wt8-',num2str(yy),'-'];
                val_caxis=no3_lim;
                level_c=no3_c;
                data=value; %.*14./1000; %% convert to mg/L
               clear value;
           case 4
                value=nc{char(vname(variation))}(:);
                val_name='Ammonium';
                unit = 'μMol N L';
                out_name_1=['vertical-GY-wt8-NH4-',num2str(yy),'-'];
                val_caxis=nh4_lim;
                level_c=nh4_c;
                data=value; %.*14./1000; %% convert to mg/L
               clear value;              
           case 5
                value=nc{char(vname(variation))}(:);
                val_name='chlorophyll';
                unit = 'μg C L';
                out_name_1=['vertical-GY-wt8-chl-',num2str(yy),'-'];
                val_caxis=chl_lim;
                level_c=chl_c;
                data=value;
               clear value;  
           case 6
                value=nc{char(vname(variation))}(:);
                val_name='phytoplankton';
                unit ='μMol N L';
                out_name_1=['vertical-GY-wt8-phy-',num2str(yy),'-'];
                val_caxis=phy_lim;
                level_c=phy_c;
                data=value;
               clear value;
           case 7
                value=nc{char(vname(variation))}(:);
                val_name='zootoplankton';
                unit ='μMol N L';
                out_name_1=['vertical-GY-wt8-zoo-',num2str(yy),'-'];
                val_caxis=zoo_lim;
                level_c=zoo_c;
                data=value;
               clear value;
           case 8
                value=nc{char(vname(variation))}(:);
                val_name='TIC';
                unit ='μMol C L';
                out_name_1=['vertical-GY-wt8-tic-',num2str(yy),'-'];
                val_caxis=tic_lim;
                level_c=tic_c;
                data=value;
               clear value;
           case 9
                value=nc{char(vname(variation))}(:);
                val_name='dissolved oxygen';
                unit ='mg/L';
                out_name_1=['vertical-GY-wt8-DO-',num2str(yy),'-'];
                val_caxis=oxy_lim;
                level_c=oxy_c;
                data=value; %.*(10/7).*(1/44.661); %% convert to mg/L 
               clear value;
           case 10
                value=nc{char(vname(variation))}(:);
                val_name='phosphate';
                unit ='μMol P L';
                out_name_1=['vertical-GY-wt8-PO4-',num2str(yy),'-'];
                val_caxis=po4_lim;
                level_c=po4_c;
                data=value;
               clear value;   
            case 11
                temp=nc{char(vname(1))}(:);salt=nc{char(vname(2))}(:);
                for N=1:1:20
                    value(N,:,:)=sw_dens(squeeze(salt(N,:,:)),squeeze(temp(N,:,:)),0)-1000;
                end
                val_name='density';
                unit = 'σ(kg/m^3)';
                out_name_1=['vertical-GY-wt8-Den-',num2str(yy),'-'];
                val_caxis=den_lim;
                level_c=den_c;
                data=value;
        end
   

            if (plot_pcolorjw)
                for i=1:1:length(data(:,1))
                    for j=1:1:length(data(1,:))
                        if data(i,j) > 10000
                            data(i,j) = NaN;
                        end
                    end
                end
            end    
   


dep=permute(depth,[2 3 1]);
data_3 = permute(data,[3 2 1]);
pick = 93;

dep(find(dep > 10000))=NaN;

%%% check line
% figure; hold on; 
% pcolor(lon_rho,lat_rho,squeeze(data_3(:,:,end)))
% shading flat; hold on; 
% plot(lon_rho(93,:),lat_rho(93,:),'r--','linew',2)
%%%%%%

for i= 1:N
    lat_rho_3(:,:,i) = lat_rho;
end
%         figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
%         set(gca,'Position',[0.2 0.15 0.73 0.75]);
%         figure('position',[400 100 550
%         550],'PaperUnits','inches','PaperPosition',[0 0 7.5 3.5]);
%         set(gca,'Position',[0.1 0.17 0.88 0.73]);
        set(gcf,'Position',[200 100 1200 700])
        set(gcf,'PaperPositionMode','auto')     
        pcolor(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(data_3(pick,:,:)))
        xlim([34.67 34.945]);
        shading flat;
        set(gca,'box','on','linewidth',1.5,'fontsize',13)
        xlabel('lat','color',color_v,'FontSize',30,'fontweight','bold')
        ylabel('Depth(m)','color',color_v,'FontSize',30,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',25,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',30,'fontweight','bold')
        out_name=[out_name_1,num2str(im,'%02d')];
        title(out_name)
        caxis(val_caxis)
        grid on;
        if (plot_contour)
                  hold on
                  [C,h]=contour(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(data_3(pick,:,:)),level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(data_3(pick,:,:)),[-1:2:1],'-w','linewidth',1);              
                  clabel(C,h,'FontSize',12,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
        end
            
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
                saveas(gca,[out_name,'.png']);
                close all;
            end
            
end
return

