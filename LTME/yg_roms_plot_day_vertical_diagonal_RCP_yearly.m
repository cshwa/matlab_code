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
cd F:\ROMS\roms_tools\Run
start
cd D:\장기생태\Dynamic\result
%==========================================================================
%-- appearance setting
% figPos = [200 100 1200 700];
figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

latLim = [34.67 34.945];
lonLim = [127.593 127.82];
% temp_lim = [6 25];
plot_pcolorjw = 1;    temp_lim = [14 23];    salt_lim = [28 35]; den_lim = [20 25]; no3_lim=[5 25]; nh4_lim=[0 2];
                      chl_lim=[0.5 2.5];   phy_lim=[0 3];    zoo_lim=[0 1];  oxy_lim=[240 260];
                      tic_lim=[-70 200];   po4_lim=[0 1];    v_lim =[-0.2 0];
% temp_c  =[6:1:25];
plot_contour  = 1;    color_c  ='-k' ;      temp_c  =[0:.5:35];salt_c  =[28:.5:35]; no3_c=[0:.2:7]; nh4_c=[0:.2:30];
                      chl_c=[0:.2:20];    phy_c=[0:.1:10];    zoo_c=[0:0.1:3];        oxy_c=[150:2:330];
                      tic_c=[-50:50:200]; po4_c=[-0.1:0.1:1];  den_c=[10:.5:30]; v_c=[-10:.05:10];

% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'tif';

section = 1; % 0 => whole, 1=> meridional, 2=> zonal

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';
  
% start_mm=1;   %04.29 - 119 08.05 - 217  09.04-247   12.01-335    01.12-377
% end_mm=12;
time_step=0;

%==========================================================================
% -- loop settings
% yy_range = 2021:10:2091;
yy_range = 2:3;
% yy_range = [1982, 1987, 1993];
% yy2_range = [1986, 1992, 2004];
% variation_range = [11,12]; %  1 = temperature // 2 = salinity // 3 = NO3 // 4 = NH4 // 5 = chlorophyll
%               %  6 = phytoplankton // 7 = zooplankton // 8 = total inorganic carbon  
%               %  9 = oxygen // 10 = PO4 // 11 = Density // 12 = DIN  //
%               %  13 = v-velocity
% variation_range = [1:5,9,11:12];
variation_range=12;
% no3_lim=[5 8];
no3_lim=[5 10];

if section == 1
% % meridional
% pick = 93;
pick = 93;
else
% % zonal
%     pick = 76;
pick = 74;
end
%==========================================================================
max_level= 30;
out_filen = 'mp_p_sewer_det_10m_do.png';
od_i = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for yy = yy_range
    od_i = od_i + 1;
    for variation = variation_range
        clearvars value
%     clearvars -except yy variation
%% 1.4.8. == 5/29 >? 1.4.9 .
%% 8.1. == 3/23 >? .8.2.
%% 7.2. == 3/14 >? 7.3.
%% 2.2.6. == 8/14
% file_dir=['\H:\E\포항모델\model\output\'];

%==========================================================================
% case_name='3reg_t4_cir_monthly_';
%% RCP
% case_name=['ncra_rcp85_',num2str(yy),'_v3_2yr'];
%% 2&3regime
case_name=['ncra_',num2str(yy),'reg_cir_2yr'];
%% 1982 to 2004
% case_name=['ncra_',num2str(yy),'to',num2str(yy2_range(od_i))];


% gd = grd('pohang_fine');
% gd = grd('pump_30_close');
fpath='J:\장기생태_2021\Dynamic\result\3regime\';
lon_rho  = ncread([fpath,case_name,'.nc'],'lon_rho');
lat_rho  = ncread([fpath,case_name,'.nc'],'lat_rho');
mask_rho = ncread([fpath,case_name,'.nc'],'mask_rho');
topo=ncread([fpath,case_name,'.nc'],'h');
N = length(ncread([fpath,case_name,'.nc'],'s_rho'));
% depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,1,1,'r');
theta_s=ncread([fpath,case_name,'.nc'],'theta_s');
theta_b=ncread([fpath,case_name,'.nc'],'theta_b');
hc=ncread([fpath,case_name,'.nc'],'hc');
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
% for im=start_mm:end_mm %pohang
% 		mid=[num2str(im,'%02d')];
%         file = [file_dir,'ocean_avg_',mid,'.nc'];
file = [fpath,case_name,'.nc'];
%% sigma2zeta
clearvars zeta depth
zeta=ncread(file,'zeta');
depth=zlevs(topo,zeta,theta_s,theta_b,hc,N,'r');
depth_w=zlevs(topo,zeta,theta_s,theta_b,hc,N,'w');
%% white color line range
white_range = [-1:2:1];
%%
% file = ['pohang_avg_',mid,'.nc'];
%         file = ['pohang_wt8_fine_flat_3m_',mid,'.nc'];
        disp([file,' : ','yearly'])
        nc=netcdf(file);
        date=[num2str(yy)];
        switch variation
            case 1
                value=nc{char(vname(variation))}(:);
                val_name='Temperature';
                unit = '^oC';
                out_name_1=['vertical-GY-Temp-',num2str(yy),'-'];
                val_caxis=temp_lim;
                level_c=temp_c;
%                 value=value;
%                clear value;         
            case 2
                value=nc{char(vname(variation))}(:);
                val_name='Salinity';
                unit = 'psu';
                out_name_1=['vertical-GY-Salt-',num2str(yy),'-'];
                val_caxis=salt_lim; 
                level_c=salt_c;
%                 value=value;
                white_range = [31.499999:0.000001:31.5000001];
%                clear value;     
           case 3
                value=nc{char(vname(variation))}(:);
                val_name='Nitrate';
                unit = 'μMol N L^-1';
                out_name_1=['vertical-GY-NO3-',num2str(yy),'-'];
                val_caxis=no3_lim;
                level_c=no3_c;
%                 value=value; %.*14./1000; %% convert to mg/L
%                clear value;
           case 4
                value=nc{char(vname(variation))}(:);
                val_name='Ammonium';
                unit = 'μMol N L^-1';
                out_name_1=['vertical-GY-NH4-',num2str(yy),'-'];
                val_caxis=nh4_lim;
                level_c=nh4_c;
%                 value=value; %.*14./1000; %% convert to mg/L
%                clear value;              
           case 5
                value=nc{char(vname(variation))}(:);
                val_name='chlorophyll';
                unit = 'μg C L^-1';
                out_name_1=['vertical-GY-chl-',num2str(yy),'-'];
                val_caxis=chl_lim;
                level_c=chl_c;
%                 value=value;
%                clear value;  
           case 6
                value=nc{char(vname(variation))}(:);
                val_name='phytoplankton';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-phy-',num2str(yy),'-'];
                val_caxis=phy_lim;
                level_c=phy_c;
%                 value=value;
%                clear value;
           case 7
                value=nc{char(vname(variation))}(:);
                val_name='zootoplankton';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-zoo-',num2str(yy),'-'];
                val_caxis=zoo_lim;
                level_c=zoo_c;
%                 value=value;
%                clear value;
           case 8
                value=nc{char(vname(variation))}(:);
                val_name='TIC';
                unit ='μMol C L^-1';
                out_name_1=['vertical-GY-tic-',num2str(yy),'-'];
                val_caxis=tic_lim;
                level_c=tic_c;
%                 value=value;
%                clear value;
           case 9
                value=nc{char(vname(variation))}(:);
                val_name='dissolved oxygen';
                unit ='μMol L^-1';
                out_name_1=['vertical-GY-DO-',num2str(yy),'-'];
                val_caxis=oxy_lim;
                level_c=oxy_c;
%                 value=value; %.*(10/7).*(1/44.661); %% convert to mg/L 
%                clear value;
           case 10
                value=nc{char(vname(variation))}(:);
                val_name='phosphate';
                unit ='μMol P L^-1';
                out_name_1=['vertical-GY-PO4-',num2str(yy),'-'];
                val_caxis=po4_lim;
                level_c=po4_c;
%                 value=value;
%                clear value;   
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
%                 value=value;
%                clear value;   
            case 12
                white_range = [5.2:0.2:5.4];
                red_range = [5.0:0.2:5.2];
                value1=nc{char(vname(3))}(:); %no3
                value2=nc{char(vname(4))}(:); %nh4
                value=value1 + value2;
                val_name='DIN';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-DIN-',num2str(yy),'-'];
                val_caxis=no3_lim;
                level_c=no3_c;
%                 value=value;
%                clear value;  
            case 13
                v=nc{char('v')}(:);
                val_name='v-velocity';
                unit ='m/s';
                out_name_1=['vertical-GY-v-vel-',num2str(yy),'-'];
                val_caxis=v_lim;
                level_c=v_c;
                [N M L]=size(v)
                 value(:,2:M,:) = mean(cat(4,v(:,1:M-1,:),v(:,2:M,:)),[4],'omitnan');
                 value(:,1,:)   = v(:,1,:);
                 value(:,M+1,:) = v(:,M,:);
%                 value=value;
%                clear value;    
        end
   

            if (plot_pcolorjw)
                for i=1:1:length(value(:,1))
                    for j=1:1:length(value(1,:))
                        if value(i,j) > 10000
                            value(i,j) = NaN;
                        end
                    end
                end
            end    
   


dep=permute(depth,[2 3 1]);
value = permute(value,[3 2 1]);

dep(find(dep > 10000))=NaN;

%%% check line
% figure; hold on; 
% pcolor(lon_rho,lat_rho,squeeze(data_3(:,:,end)))
% shading flat; hold on; 
% plot(lon_rho(93,:),lat_rho(93,:),'r--','linew',2)
%%%%%%



if section == 1
    xLim = latLim;

    lat_rho_3 = repmat(lat_rho,1,1,N);

    figure
%         figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
%         set(gca,'Position',[0.2 0.15 0.73 0.75]);
%         figure('position',[400 100 550
%         550],'PaperUnits','inches','PaperPosition',[0 0 7.5 3.5]);
%         set(gca,'Position',[0.1 0.17 0.88 0.73]);
        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);
%         set(gcf,'PaperPositionMode','auto')     
        pcolor(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(value(pick,:,:))); shading flat;
        xlim(xLim);
        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
        xlabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
        ylabel('Depth(m)','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
        caxis(val_caxis)
        out_name=[out_name_1,'yearly'];
        title(out_name)
        if (plot_contour)
                  hold on
                  [C,h]=contour(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(value(pick,:,:)),level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(value(pick,:,:)),white_range,'-w','linewidth',1);      
                  if exist('red_range','var') == 1
                      [C3,h3]=contour(squeeze(lat_rho_3(pick,:,:)),squeeze(dep(pick,:,:)),squeeze(value(pick,:,:)),red_range,'-r','linewidth',1);
                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                  end
                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
        end
            
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
%                 saveas(gca,[out_name,'.png']);
                print(gcf,['figures\',out_name,'.png'],'-dpng','-r200');
                close all;
            end
elseif section == 2
    xLim = lonLim;

    lon_rho_3 = repmat(lon_rho,1,1,N);

%         figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
%         set(gca,'Position',[0.2 0.15 0.73 0.75]);
%         figure('position',[400 100 550
%         550],'PaperUnits','inches','PaperPosition',[0 0 7.5 3.5]);
%         set(gca,'Position',[0.1 0.17 0.88 0.73]);
        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
%         set(gcf,'PaperPositionMode','auto')     
        pcolor(squeeze(lon_rho_3(:,pick,:)),squeeze(dep(:,pick,:)),squeeze(value(:,pick,:))); shading flat;
        xlim(xLim);
        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
        ylabel('Depth(m)','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
        caxis(val_caxis)
        out_name=[out_name_1,'yearly'];
        title(out_name)
        if (plot_contour)
                  hold on
                  [C,h]=contour(squeeze(lon_rho_3(:,pick,:)),squeeze(dep(:,pick,:)),squeeze(value(:,pick,:)),level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(squeeze(lon_rho_3(:,pick,:)),squeeze(dep(:,pick,:)),squeeze(value(:,pick,:)),white_range,'-w','linewidth',1);      
                  if exist('red_range','var') == 1
                      [C3,h3]=contour(squeeze(lon_rho_3(:,pick,:)),squeeze(dep(:,pick,:)),squeeze(value(:,pick,:)),red_range,'-r','linewidth',1);
                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                  end
                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
        end
            
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
%                 saveas(gca,[out_name,'_meridional','.png']);
                print(gcf,['figures\',out_name,'_meridional','.png'],'-dpng','-r200');
%                 close all;
            end   
    
            
end

% end


%% vertical integration
% switch variation
%             case 5
                
dep_w = permute(depth_w,[2 3 1]);

% for i = N:-1:1
%     if i == N
%         dep_l(:,:,i) =  dep_w(:,:,i+1) - dep(:,:,i);
%     else
%         dep_l(:,:,i) =  dep(:,:,i+1) - dep(:,:,i); 
%     end
% end

for i = N:-1:1
    if i == N
        dep_l(:,:,i) =  dep_w(:,:,i+1) - dep_w(:,:,i);
    else
        dep_l(:,:,i) =  dep_w(:,:,i+1) - dep_w(:,:,i); 
    end
end

depcchl= dep_l .* value;
int_chl =  sum(depcchl,3) ./ sum(dep_l,3);
% pcolor(int_chl'); shading flat; colorbar

figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
%         set(gcf,'PaperPositionMode','auto')     
        pcolor(lon_rho,lat_rho,int_chl); shading flat;
%         xlim(xLim);
        ylim([34.8 34.96])
        xlim([127.579 127.87])
        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
        caxis(val_caxis)
        out_name=[out_name_1,'yearly'];
        title(out_name)
        if (plot_contour)
                  hold on
                  [C,h]=contour(lon_rho,lat_rho,int_chl,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(lon_rho,lat_rho,int_chl,white_range,'-w','linewidth',1);      
                  if exist('red_range','var') == 1
                      [C3,h3]=contour(lon_rho,lat_rho,int_chl,red_range,'-r','linewidth',1);
                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                  end
                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
        end
            
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
%                 saveas(gca,[out_name,'_meridional','.png']);
                print(gcf,['figures\',out_name,'_vertical_integral','.png'],'-dpng','-r200');
%                 close all;
            end
            
if variation_range == 13            
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
%         set(gcf,'PaperPositionMode','auto')     
        pcolor(lon_rho,lat_rho,int_chl); shading flat;
%         xlim(xLim);
        ylim([34.9 35.02])
        xlim([127.68 127.87])
        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
        caxis(val_caxis)
        out_name=[out_name_1,'yearly'];
        title(out_name)
        if (plot_contour)
                  hold on
                  [C,h]=contour(lon_rho,lat_rho,int_chl,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(lon_rho,lat_rho,int_chl,white_range,'-w','linewidth',1);      
                  if exist('red_range','var') == 1
                      [C3,h3]=contour(lon_rho,lat_rho,int_chl,red_range,'-r','linewidth',1);
                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                  end
                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
        end
            
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
%                 saveas(gca,[out_name,'_meridional','.png']);
                print(gcf,['figures\',out_name,'_vertical_integral_upper','.png'],'-dpng','-r200');
%                 close all;
            end            
end

% end
    end
end
