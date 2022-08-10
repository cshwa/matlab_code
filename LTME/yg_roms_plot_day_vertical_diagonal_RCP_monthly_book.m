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
% temp_lim = [6 25]; chl_lim=[0 7]; den_c=[10:.5:100]; chl_c=[0:.5:30];
% no3_lim=[9 15]; % RCP
%% winter (feb)
plot_pcolorjw = 1;    temp_lim = [4 12];    salt_lim = [27 35]; den_lim = [10 25]; no3_lim=[8 30]; nh4_lim=[0 2];
                      chl_lim=[0 2];   phy_lim=[0 3];    zoo_lim=[0 1];  oxy_lim=[280 320];
                      tic_lim=[-70 200];   po4_lim=[0 1]; Ldet_lim=[0 1]; Sdet_lim=[0 4]; v_lim=[-1 1];
                      Sdet_up_lim=[0 5];
                      % temp_c  =[6:1:25];  no3_c=[0:2:30];
plot_contour  = 1;    color_c  ='-k' ;      temp_c =[0:.5:40];salt_c =[30:1:35]; no3_c=[0:5:30]; nh4_c=[0:.1:30];
                      chl_c=[0:.2:30];    phy_c=[0:.1:10];    zoo_c=[0:0.1:3];        oxy_c=[280:5:330];
                      tic_c=[-50:50:200]; po4_c=[-0.1:0.1:1];  den_c=[0:1:100]; Ldet_c=[0:.05:30]; Sdet_c=[0:.2:100];
                      v_c=[-10:0.05:1];                      
no3_lim=[9 15]; % RCP
                      
%% summer (aug)
% historical only
% salt_lim = [10 30]; oxy_lim=[205 280]; no3_lim=[4 60]; % historical only
% plot_pcolorjw = 1;    temp_lim = [27 33];    salt_lim = [20 31]; den_lim = [10 25]; no3_lim=[0 25]; nh4_lim=[0 2];
%                       chl_lim=[0 7];   phy_lim=[0 3];    zoo_lim=[0 1];  oxy_lim=[205 240];
%                       tic_lim=[-70 200];   po4_lim=[0 1]; Ldet_lim=[0 1]; Sdet_lim=[0 4]; v_lim=[-1 1];
%                       Sdet_up_lim=[0 5];
%                       
% plot_contour  = 1;    color_c  ='-k' ;      temp_c =[0:.5:40];salt_c =[10:2:35]; no3_c=[0:5:30]; nh4_c=[0:.1:30];
%                       chl_c=[0:1:30];    phy_c=[0:.1:10];    zoo_c=[0:0.1:3];        oxy_c=[205:5:240];
%                       tic_c=[-50:50:200]; po4_c=[-0.1:0.1:1];  den_c=[0:1:100]; Ldet_c=[0:.05:30]; Sdet_c=[0:.2:100];
%                       v_c=[-10:0.05:1];
                      
% salt_lim = [10 30]; oxy_lim=[205 280]; no3_lim=[4 60]; % historical only
% oxy_c=[205:5:280]; no3_c=[5:5:60];

                      
% no3_lim=[9 15];

% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'tif';

vertical =  0;  % vertical integration section 0 => no vertical integration (surf,  bot), 1 => plot vertical integration

section = 0; % 0 => no   plot, 1=> meridional, 2=> zonal

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';
  
% start_mm=1;   %04.29 - 119 08.05 - 217  09.04-247   12.01-335    01.12-377
% end_mm=12;
time_step=0;

switch_upper_river = 0; % 0 = off , 1 =on

%==========================================================================
%-- loop settings

yy_range = 3;
% yy_range = [1982];
% yy2 = [1986];
% yy_range = 2091;

% yy_range = 2021:10:2091;
% yy_range = 2:3;
% yy_range = [1982, 1987, 1993];
% yy2 = [1986, 1992, 2004];
% yy_range = 3;
% variation_range = [11,12]; %  1 = temperature // 2 = salinity // 3 = NO3 // 4 = NH4 // 5 = chlorophyll
%               %  6 = phytoplankton // 7 = zooplankton // 8 = total inorganic carbon  
%               %  9 = oxygen // 10 = PO4 // 11 = Density // 12 = DIN // 13 = LDet
%               %  14 = Sdet // 15 = Ldet + Sdet = total-det // 16 = vector(U&V)
%               %  17 = depth integration vector (Ubar & Vbar)
% variation_range = [1:5,9,11:12,15];
% variation_range = [1:2,5,9,12];
variation_range = [1];
% variation_range = [12];
% variation_range = [5];
% no3_lim=[5 8];
% no3_lim=[5 15];
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
% set_pic_path = 'figures\monthly\book\feb\';  %% historical
% set_pic_path = 'figures\monthly\book\RCP\feb\etc\';  %% RCP
% set_pic_path = 'figures\monthly\book\RCP\feb\';  %% RCP
% set_pic_path = 'figures\monthly\book\RCP\aug\etc\';  %% RCP 
% set_pic_path = 'figures\monthly\book\RCP\aug\';  %% RCP
% set_pic_path = 'figures\monthly\book\aug\';  %% RCP
set_pic_path = 'figures\monthly\book\';  %% RCP
od_i = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for yy = yy_range
    od_i = od_i +1;
%     for monthly_t = 1:12
%     for monthly_t = [2,8]
    for monthly_t = [8]
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
% case_name=['3reg_RCP_',num2str(yy),'_v3_2yr_monthly_', num2str(monthly_t,'%02d')];
%% 2&3regime
% case_name=['ncra_',num2str(yy),'reg_cir_2yr'];
case_name=[num2str(yy),'reg_cir_2yr_monthly_', num2str(monthly_t,'%02d')];
%% 1982to2004
% case_name=[num2str(yy),'to',num2str(yy2(od_i)),'_2yr_monthly_',num2str(monthly_t,'%02d')];
% 

% gd = grd('pohang_fine');
% gd = grd('pump_30_close');
fpath='J:\장기생태_2021\Dynamic\result\3regime\monthly\';
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
vname = {'temp','salt' ,'NO3','NH4','chlorophyll','phytoplankton','zooplankton','TIC','oxygen' 'tPO4','LdetritusN','SdetritusN'};%,'zeta','ubar','vbar','u','v','omega'};
     
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
        disp([file,' : ','monthly'])
        nc=netcdf(file);
        date=[num2str(yy)];
        switch variation
            case 1
                value=nc{char(vname(variation))}(:);
                val_name='수온';
                unit = '^oC';
                out_name_1=['vertical-GY-Temp-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=temp_lim;
                level_c=temp_c;
%                 value=value;
%                clear value;         
            case 2
                value=nc{char(vname(variation))}(:);
                val_name='염분';
                unit = 'g/kg';
                out_name_1=['vertical-GY-Salt-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=salt_lim; 
                level_c=salt_c;
%                 value=value;
%                 white_range = [31.499999:0.000001:31.5000001];
%                clear value;     
           case 3
                value=nc{char(vname(variation))}(:);
                val_name='질산염';
                unit = 'μM/L';
                out_name_1=['vertical-GY-NO3-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=no3_lim;
                level_c=no3_c;
%                 value=value; %.*14./1000; %% convert to mg/L
%                clear value;
           case 4
                value=nc{char(vname(variation))}(:);
                val_name='암모늄';
                unit = 'μM/L';
                out_name_1=['vertical-GY-NH4-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=nh4_lim;
                level_c=nh4_c;
%                 value=value; %.*14./1000; %% convert to mg/L
%                clear value;              
           case 5
                value=nc{char(vname(variation))}(:);
                val_name='클로로필';
                unit = 'μg/L';
                out_name_1=['vertical-GY-chl-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=chl_lim;
                level_c=chl_c;
%                 value=value;
%                clear value;  
           case 6
                value=nc{char(vname(variation))}(:);
                val_name='phytoplankton';
                unit ='μM/L';
                out_name_1=['vertical-GY-phy-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=phy_lim;
                level_c=phy_c;
%                 value=value;
%                clear value;
           case 7
                value=nc{char(vname(variation))}(:);
                val_name='zootoplankton';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-zoo-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=zoo_lim;
                level_c=zoo_c;
%                 value=value;
%                clear value;
           case 8
                value=nc{char(vname(variation))}(:);
                val_name='TIC';
                unit ='μMol C L^-1';
                out_name_1=['vertical-GY-tic-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=tic_lim;
                level_c=tic_c;
%                 value=value;
%                clear value;
           case 9
                value=nc{char(vname(variation))}(:);
                val_name='용존산소';
                unit ='μM/L';
                out_name_1=['vertical-GY-DO-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=oxy_lim;
                level_c=oxy_c;
%                 value=value; %.*(10/7).*(1/44.661); %% convert to mg/L 
%                clear value;
           case 10
                value=nc{char(vname(variation))}(:);
                val_name='phosphate';
                unit ='μMol P L^-1';
                out_name_1=['vertical-GY-PO4-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
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
                out_name_1=['vertical-GY-Den-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=den_lim;
                level_c=den_c;
%                 value=value;
%                clear value;   
            case 12
                white_range = [1:0.001:1.001];
                red_range = [-0.002:0:0.002];
                value1=nc{char(vname(3))}(:); %no3
                value2=nc{char(vname(4))}(:); %nh4
                value=value1 + value2;
                val_name='용존무기질소';
                unit ='μM/L';
                out_name_1=['vertical-GY-DIN-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=no3_lim;
                level_c=no3_c;
%                 value=value;
%                clear value;
           case 13
                value=nc{char(vname(11))}(:);
                val_name='LdetritusN';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-LdetN-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=Ldet_lim;
                level_c=Ldet_c;
%                 value=value;
%                clear value;  
           case 14
                value=nc{char(vname(12))}(:);
                val_name='SdetritusN';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-SdetN-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=Sdet_lim;
                level_c=Sdet_c;
%                 value=value;
%                clear value;   
          case 15
                value1=nc{char(vname(11))}(:);value2=nc{char(vname(12))}(:);
                value = value1 + value2;
                val_name='total-detN';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-total-detN-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=Sdet_lim;
                level_c=Sdet_c;
                val_up_caxis=Sdet_up_lim;
                level_up_c=Sdet_c;
%                 value=value;
%                clear value; 
          case 16
                value1_d=nc{char('u')}(:); value2_d=nc{char('v')}(:);
                %%
                [L M N]=size(value1_d);
                u =ones(L,M,N+1)*NaN;
                value1_d(value1_d>100)=0;
                u(:,:,2:N) = mean(cat(4,value1_d(:,:,1:N-1),value1_d(:,:,2:N)),[4],'omitnan');
                u(:,:,1)   = value1_d(:,:,1);
                u(:,:,N+1) = value1_d(:,:,N);
                [L M N]=size(value2_d);
                v =ones(L,M+1,N)*NaN;
                value2_d(value2_d>100)=0;
                v(:,2:M,:) = mean(cat(4,value2_d(:,1:M-1,:),value2_d(:,2:M,:)),[4],'omitnan');
                v(:,1,:)   = value2_d(:,1,:);
                v(:,M+1,:) = value2_d(:,M,:);
                u= permute(u,[3 2 1]); v= permute(v,[3 2 1]);
                value = v;
                interval = 3;
                %%
%                 value = value1 + value2;
                val_name='velocity';
                unit ='m/s';
                out_name_1=['vertical-GY-velocity-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=v_lim;
                level_c=v_c;
%                 value=value;
%                clear value; 
          case 17
                value1_d=nc{char('ubar')}(:); value2_d=nc{char('vbar')}(:);
                %%
                [M N]=size(value1_d);
                u =ones(M,N+1)*NaN;
                value1_d(value1_d>100)=0;
                u(:,2:N) =mean(cat(3,value1_d(:,1:N-1),value1_d(:,2:N)),[3],'omitnan');
                u(:,1)   = value1_d(:,1);
                u(:,N+1) = value1_d(:,N);
                [M N]=size(value2_d);
                v =ones(M+1,N)*NaN;
                value2_d(value2_d>100)=0;
                v(2:M,:) =mean(cat(3,value2_d(1:M-1,:),value2_d(2:M,:)),[3],'omitnan');
                v(1,:)   = value2_d(1,:);
                v(M+1,:) = value2_d(M,:);
                %%
%                 value = value1 + value2;
                val_name='velocity';
                unit ='m/s';
                out_name_1=['vertical-GY-velocity-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=v_lim;
                level_c=v_c;
                u=u'; v=v';
                value=v;
                interval = 3;
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
        out_name=[out_name_1,'monthly'];
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
                print(gcf,[set_pic_path,out_name,'_meridional.png'],'-dpng','-r200');
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
        out_name=[out_name_1,'monthly'];
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
                print(gcf,[set_pic_path,out_name,'_zonal','.png'],'-dpng','-r200');
                close all;
            end   
    
            
end

% end

if vertical == 1

%% vertical integration
switch variation
            case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
                
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
                if switch_upper_river == 0
                    ylim([34.8 34.96])
                    xlim([127.579 127.87])   
                elseif switch_upper_river == 1
                    ylim([34.95 35.1])
                    xlim([127.63 127.793])
                end
                        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
                        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                %         out_name_1=[out_name_1];
                         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
                         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')                       
                         if exist('val_up_caxis') == 1 & switch_upper_river == 1
                             caxis(val_up_caxis)
                         else
                             caxis(val_caxis)
                         end
                        out_name=[out_name_1];
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
                               if switch_upper_river == 0
                                       print(gcf,[set_pic_path,out_name,'_vertical_integral','.png'],'-dpng','-r200');
                                       close all; 
                                elseif switch_upper_river == 1
                                       print(gcf,[set_pic_path,out_name,'_vertical_integral_upper_river','.png'],'-dpng','-r200');
                                       close all;
                               end
                            end 
                            
                 case {16} %velocity (U and V) surf and bottom
     
                u_top = squeeze(u(:,:,20));
                v_top = squeeze(v(:,:,20));
                u_bot = squeeze(u(:,:,1));
                v_bot = squeeze(v(:,:,1));
                
                % referenc vector
                    scale_vec = 0.1;
                    clearvars *0_ref
                    u0_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    v0_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                if switch_upper_river == 0
                    x_pick = 11; y_pick = 120;  
                elseif switch_upper_river == 1
                    x_pick = 20; y_pick = 155;
                end
                    u0_ref(x_pick,y_pick) = 0.1 * scale_vec;

                        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
%                         quiver(lon_rho,lat_rho,u',v'); 
                        pcolor(lon_rho,lat_rho,v_top .*(mask_rho ./ mask_rho)); shading flat; hold on;
                        quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end), ...
    squeeze(u_top(1:interval:end,1:interval:end) .* scale_vec), ...
    squeeze(v_top(1:interval:end,1:interval:end) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
                        quiver(lon_rho,lat_rho,u0_ref, v0_ref, 'AutoScale','off','LineWidth',1,'Color','r');
                        text(lon_rho(x_pick,y_pick-5),lat_rho(x_pick,y_pick-5),'0.1m/s','fontsize',15,'Color','r');
                %         xlim(xLim);
                if switch_upper_river == 0
                    ylim([34.8 34.96])
                    xlim([127.579 127.87])   
                elseif switch_upper_river == 1
                    ylim([34.95 35.1])
                    xlim([127.63 127.793])
                end
                        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
                        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                %         out_name_1=[out_name_1];
                         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
                         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
                        caxis(val_caxis)
                        out_name=[out_name_1,'monthly'];
                        colormap('bluewhitered');
                        title(out_name)
%                         if (plot_contour)
%                                   hold on
%                                   [C,h]=contour(lon_rho,lat_rho,int_chl,level_c,color_c,'linewidth',1);
%                                   [C2,h2]=contour(lon_rho,lat_rho,int_chl,white_range,'-w','linewidth',1);      
%                                   if exist('red_range','var') == 1
%                                       [C3,h3]=contour(lon_rho,lat_rho,int_chl,red_range,'-r','linewidth',1);
%                                       clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
%                                   end
%                                   clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
%                                   clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
%                         end

                            if (switch_save)
                               if switch_upper_river == 0
                                       print(gcf,[set_pic_path,out_name,'_surf','.png'],'-dpng','-r200');
                                       close all; 
                                elseif switch_upper_river == 1
                                       print(gcf,[set_pic_path,out_name,'_surf_upper_river','.png'],'-dpng','-r200');
                                       close all;
                               end
                            end 
                            
                         figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
%                         quiver(lon_rho,lat_rho,u',v'); 
                        pcolor(lon_rho,lat_rho,v_bot .*(mask_rho ./ mask_rho)); shading flat; hold on;
                        quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end), ...
    squeeze(u_bot(1:interval:end,1:interval:end) .* scale_vec), ...
    squeeze(v_bot(1:interval:end,1:interval:end) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
                        quiver(lon_rho,lat_rho,u0_ref, v0_ref, 'AutoScale','off','LineWidth',1,'Color','r');
                        text(lon_rho(x_pick,y_pick-5),lat_rho(x_pick,y_pick-5),'0.1m/s','fontsize',15,'Color','r');
                %         xlim(xLim);
                if switch_upper_river == 0
                    ylim([34.8 34.96])
                    xlim([127.579 127.87])   
                elseif switch_upper_river == 1
                    ylim([34.95 35.1])
                    xlim([127.63 127.793])
                end
                        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
                        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                %         out_name_1=[out_name_1];
                         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
                         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
                        caxis(val_caxis)
                        out_name=[out_name_1,'monthly'];
                        colormap('bluewhitered');
                        title(out_name)
%                         if (plot_contour)
%                                   hold on
%                                   [C,h]=contour(lon_rho,lat_rho,int_chl,level_c,color_c,'linewidth',1);
%                                   [C2,h2]=contour(lon_rho,lat_rho,int_chl,white_range,'-w','linewidth',1);      
%                                   if exist('red_range','var') == 1
%                                       [C3,h3]=contour(lon_rho,lat_rho,int_chl,red_range,'-r','linewidth',1);
%                                       clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
%                                   end
%                                   clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
%                                   clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
%                         end

                            if (switch_save)
                               if switch_upper_river == 0
                                       print(gcf,[set_pic_path,out_name,'_bot','.png'],'-dpng','-r200');
                                       close all; 
                                elseif switch_upper_river == 1
                                       print(gcf,[set_pic_path,out_name,'_bot','.png'],'-dpng','-r200');
                                       close all;
                               end
                            end             
                            
            case {17} %velocity (Ubar and Vbar)
                
                % referenc vector
                    scale_vec = 0.2;
                    clearvars *_ref
                    u_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    v_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    u_ref(11,120) = 0.1 * scale_vec;

                        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
%                         quiver(lon_rho,lat_rho,u',v'); 
                        pcolor(lon_rho,lat_rho,v .*(mask_rho ./ mask_rho)); shading flat; hold on;
                        quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end), ...
    squeeze(u(1:interval:end,1:interval:end) .* scale_vec), ...
    squeeze(v(1:interval:end,1:interval:end) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
                        quiver(lon_rho,lat_rho,u_ref, v_ref, 'AutoScale','off','LineWidth',1,'Color','r');
                        text(lon_rho(11,115),lat_rho(11,115),'0.5m/s','fontsize',15,'Color','r');
                %         xlim(xLim);
                        ylim([34.8 34.96])
                        xlim([127.579 127.87])
                        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
                        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                %         out_name_1=[out_name_1];
                         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
                         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
%                         caxis(val_caxis)
                        out_name=[out_name_1,'monthly'];
                        colormap('bluewhitered');
                        title(out_name)
%                         if (plot_contour)
%                                   hold on
%                                   [C,h]=contour(lon_rho,lat_rho,int_chl,level_c,color_c,'linewidth',1);
%                                   [C2,h2]=contour(lon_rho,lat_rho,int_chl,white_range,'-w','linewidth',1);      
%                                   if exist('red_range','var') == 1
%                                       [C3,h3]=contour(lon_rho,lat_rho,int_chl,red_range,'-r','linewidth',1);
%                                       clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
%                                   end
%                                   clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
%                                   clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
%                         end

                            if (switch_save)
                %                 saveas(gcf,out_name,out_type);
                %                 saveas(gca,[out_name,'_meridional','.png']);
                                print(gcf,[set_pic_path,out_name,'_vertical_integral','.png'],'-dpng','-r200');
                                close all;
                            end    
                               
end %% variation
                                
elseif vertical == 0 %% surf and bottom
    switch variation
            case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
                clearvars value_plt
                value_plt = squeeze(value(:,:,20));
                 figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
                        pcolor(lon_rho,lat_rho,value_plt); shading flat;
                %         xlim(xLim);
                if switch_upper_river == 0
                    ylim([34.8 34.96])
                    xlim([127.579 127.87])   
                elseif switch_upper_river == 1
                    ylim([34.95 35.1])
                    xlim([127.63 127.793])
                end
                        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
%                         xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%                         ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                %         out_name_1=[out_name_1];
                         if od_i == length(yy_range) && yy_range(od_i) >= 2020
                         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
                         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')  
                         elseif yy_range(od_i) == 3 && od_i == 2
                         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
%                          bar = colorbar('fontsize',fontSizeLab,'fontweight','bold','Location','eastoutside');
                         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')      
                         end
                         
                             if exist('val_up_caxis') == 1 & switch_upper_river == 1
                                 caxis(val_up_caxis)
                             else
                                 caxis(val_caxis)
                             end
                         
                        out_name=[out_name_1];                      
%                         yy_range = [1982, 1987, 1993];
%                         yy2 = [1986, 1992, 2004];
%                         val_name
                        if exist('yy2','var') == 1
%                             title([val_name,' - ',num2str(yy_range(od_i)),'~',num2str(yy2(od_i))])
                            text(127.583, 34.82, [val_name,' ',num2str(yy_range(od_i)),'~',num2str(yy2(od_i))],'fontsize',13,'fontweight','bold')
                        else
                            if yy_range(od_i) <= 2020
                            ytt1=[2007, 2016];
                            ytt2=[2015, 2020];
                            if yy_range(od_i) == 3 && od_i == 1
                                ytt1=[2016];
                                ytt2=[2020];
                            end
%                             title([val_name,' - ',num2str(ytt1(od_i)),'~',num2str(ytt2(od_i))])
                            text(127.583, 34.82, [val_name,' ',num2str(ytt1(od_i)),'~',num2str(ytt2(od_i))],'fontsize',13,'fontweight','bold')
                            elseif yy_range(od_i) > 2020
%                             title([val_name,' - ',num2str(yy_range(od_i)),'~',num2str(yy_range(od_i)+9)]) 
                            text(127.583, 34.82, [val_name,' ',num2str(yy_range(od_i)),'~',num2str(yy_range(od_i)+9)],'fontsize',13,'fontweight','bold')
                            end
                        end

                        if (plot_contour)
                                  hold on
                                  [C,h]=contour(lon_rho,lat_rho,value_plt,level_c,color_c,'linewidth',1);
                                  [C2,h2]=contour(lon_rho,lat_rho,value_plt,white_range,'-w','linewidth',1);      
                                  if exist('red_range','var') == 1
                                      [C3,h3]=contour(lon_rho,lat_rho,value_plt,red_range,'-r','linewidth',1);
                                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                                  end
                                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
                        end
                        set(gca,'YTickLabel',[]);
                        set(gca,'XTickLabel',[]);
                        colormap('jet')
                            if (switch_save)
                               if switch_upper_river == 0
                                       print(gcf,[set_pic_path,out_name,'_surf','.png'],'-dpng','-r500');
                                       close all; 
                                elseif switch_upper_river == 1
                                       print(gcf,[set_pic_path,out_name,'_surf_upper_river','.png'],'-dpng','-r500');
                                       close all;
                               end
                            end 
                            
%                              clearvars value_plt
%                   value_plt = squeeze(value(:,:,1));
%                  figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
%                 %         set(gcf,'PaperPositionMode','auto')     
%                         pcolor(lon_rho,lat_rho,value_plt); shading flat;
%                 %         xlim(xLim);
%                 if switch_upper_river == 0
%                     ylim([34.8 34.96])
%                     xlim([127.579 127.87])   
%                 elseif switch_upper_river == 1
%                     ylim([34.95 35.1])
%                     xlim([127.63 127.793])
%                 end
%                         set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
% %                         xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
% %                         ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%                 %         out_name_1=[out_name_1];
%                          bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
%                          set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')                       
%                          if exist('val_up_caxis') == 1 & switch_upper_river == 1
%                              caxis(val_up_caxis)
%                          else
%                              caxis(val_caxis)
%                          end
%                         out_name=[out_name_1];
%                         title(out_name)
%                         if (plot_contour)
%                                   hold on
%                                   [C,h]=contour(lon_rho,lat_rho,value_plt,level_c,color_c,'linewidth',1);
%                                   [C2,h2]=contour(lon_rho,lat_rho,value_plt,white_range,'-w','linewidth',1);      
%                                   if exist('red_range','var') == 1
%                                       [C3,h3]=contour(lon_rho,lat_rho,value_plt,red_range,'-r','linewidth',1);
%                                       clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
%                                   end
%                                   clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
%                                   clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
%                         end
%                         set(gca,'YTickLabel',[]);
%                         set(gca,'XTickLabel',[]);
% 
%                             if (switch_save)
%                                if switch_upper_river == 0
%                                        print(gcf,[set_pic_path,out_name,'_bot','.png'],'-dpng','-r200');
%                                        close all; 
%                                 elseif switch_upper_river == 1
%                                        print(gcf,[set_pic_path,out_name,'_bot_upper_river','.png'],'-dpng','-r200');
%                                        close all;
%                                end
%                             end 

    end
end %% vertical end

    end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% final test for project figure
%% summer (aug)
% % % historical only
if monthly_t == 8
salt_lim = [10 30]; oxy_lim=[205 280]; no3_lim=[4 60]; % historical only
plot_pcolorjw = 1;    temp_lim = [27 33];    salt_lim = [20 31]; den_lim = [10 25]; no3_lim=[0 25]; nh4_lim=[0 2];
                      chl_lim=[0 7];   phy_lim=[0 3];    zoo_lim=[0 1];  oxy_lim=[205 240];
                      tic_lim=[-70 200];   po4_lim=[0 1]; Ldet_lim=[0 1]; Sdet_lim=[0 4]; v_lim=[-1 1];
                      Sdet_up_lim=[0 5];
                      
plot_contour  = 1;    color_c  ='-k' ;      temp_c =[0:.5:40];salt_c =[10:2:35]; no3_c=[0:5:30]; nh4_c=[0:.1:30];
                      chl_c=[0:1:30];    phy_c=[0:.1:10];    zoo_c=[0:0.1:3];        oxy_c=[205:5:240];
                      tic_c=[-50:50:200]; po4_c=[-0.1:0.1:1];  den_c=[0:1:100]; Ldet_c=[0:.05:30]; Sdet_c=[0:.2:100];
                      v_c=[-10:0.05:1];
val_caxis=temp_lim; 
end
    switch variation
            case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
                clearvars value_plt
                value_plt = squeeze(value(:,:,20));
                 figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto') 
                hold on;
                        pcolor(lon_rho,lat_rho,value_plt); shading flat;
                %         xlim(xLim);
                if switch_upper_river == 0
                    ylim([34.8 34.96])
                    xlim([127.579 127.87])   
                elseif switch_upper_river == 1
                    ylim([34.95 35.1])
                    xlim([127.63 127.793])
                end
                        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
%                         xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%                         ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                %         out_name_1=[out_name_1];
%                          if od_i == length(yy_range) && yy_range(od_i) >= 2020
% %                          bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
%                          set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')  
%                          elseif yy_range(od_i) == 3 && od_i == 2
% %                          bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
% %                          bar = colorbar('fontsize',fontSizeLab,'fontweight','bold','Location','eastoutside');
%                          set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')      
%                          end
                         
                             if exist('val_up_caxis') == 1 & switch_upper_river == 1
                                 caxis(val_up_caxis)
                             else
                                 caxis(val_caxis)
                             end
                         
                        out_name=[out_name_1];                      
%                         yy_range = [1982, 1987, 1993];
%                         yy2 = [1986, 1992, 2004];
%                         val_name
                        if exist('yy2','var') == 1
%                             title([val_name,' - ',num2str(yy_range(od_i)),'~',num2str(yy2(od_i))])
%                             text(127.583, 34.82, [val_name,' ',num2str(yy_range(od_i)),'~',num2str(yy2(od_i))],'fontsize',13,'fontweight','bold')
                              text(127.583, 34.82, [val_name,' ',num2str(yy_range(od_i)),' - ',num2str(monthly_t),'월'],'fontsize',13,'fontweight','bold')
%                               text(127.583, 34.82, [val_name,' ',num2str(yy_range(od_i)),' - 2월'],'fontsize',13,'fontweight','bold')
                        else
                            if yy_range(od_i) <= 2020
                            ytt1=[2007, 2016];
                            ytt2=[2015, 2020];
                            if yy_range(od_i) == 3 && od_i == 1
                                ytt1=[2016];
                                ytt2=[2020];
                            end
%                             title([val_name,' - ',num2str(ytt1(od_i)),'~',num2str(ytt2(od_i))])
%                             text(127.583, 34.82, [val_name,' ',num2str(ytt1(od_i)),'~',num2str(ytt2(od_i))],'fontsize',13,'fontweight','bold')
                              text(127.583, 34.82, [val_name,' ',num2str(ytt2(od_i)),' - ',num2str(monthly_t),'월'],'fontsize',13,'fontweight','bold')
%                               text(127.583, 34.82, [val_name,' ',num2str(ytt2(od_i)),' - 8월'],'fontsize',13,'fontweight','bold')
                            elseif yy_range(od_i) > 2020
%                             title([val_name,' - ',num2str(yy_range(od_i)),'~',num2str(yy_range(od_i)+9)]) 
                            text(127.583, 34.82, [val_name,' ',num2str(yy_range(od_i)+9),' - ',num2str(monthly_t),'월'],'fontsize',13,'fontweight','bold')
%                             text(127.583, 34.82, [val_name,' ',num2str(yy_range(od_i)+9),' - 8월'],'fontsize',13,'fontweight','bold')
                            end
                        end

                        if (plot_contour)
                                  hold on
                                  [C,h]=contour(lon_rho,lat_rho,value_plt,level_c,color_c,'linewidth',1);
                                  [C2,h2]=contour(lon_rho,lat_rho,value_plt,white_range,'-w','linewidth',1);      
                                  if exist('red_range','var') == 1
                                      [C3,h3]=contour(lon_rho,lat_rho,value_plt,red_range,'-r','linewidth',1);
                                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                                  end
                                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
                        end
                        set(gca,'YTickLabel',[]);
                        set(gca,'XTickLabel',[]);
                        colormap('jet')
                        
                            axin = gca;
%                             P = figPos*get(groot,'ScreenPixelsPerInch');
%                             P = [min(lon_rho,[],all) max(lon_rho,[],all) max(lon_rho,[],all) ;
%                             Px = [P(1), P(1)+P(3), P(1)+P(3), P(1), P(1)];
                            Px = [min(lon_rho,[],'all'), max(lon_rho,[],'all'), max(lon_rho,[],'all'), min(lon_rho,[],'all')];
%                             Py = [P(2), P(2), P(2)+P(4), P(2)+P(4), P(2)];
                            Py = [min(lat_rho,[],'all'), min(lat_rho,[],'all'), max(lat_rho,[],'all'), max(lat_rho,[],'all')];
                            hold(axin, 'on');
                            Pp = fill(Px, Py, [169/255 169/255 169/255]);
                            uistack(Pp, 'bottom');
                            ax = gca;
                            ax.XAxis.LineWidth = 5;
                            ax.YAxis.LineWidth = 5;

%                             hold(axin, 'off')
                        
                            if (switch_save)
                               if switch_upper_river == 0
                                       print(gcf,[set_pic_path,out_name,'_surf','.png'],'-dpng','-r2000');
%                                        close all; 
                                elseif switch_upper_river == 1
                                       print(gcf,[set_pic_path,out_name,'_surf_upper_river','.png'],'-dpng','-r2000');
%                                        close all;
                               end
                            end 
    end