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


latLim = [34.67 34.945]; %% meridional
lonLim = [127.593 127.82]; %% zonal


% temp_lim = [6 25];
plot_pcolorjw = 1;    temp_lim = [-8 8];    salt_lim = [-2 2]; den_lim = [-2 0]; no3_lim=[-5 5]; nh4_lim=[-.5 .5];
                      chl_lim=[-1.2 1.2];   phy_lim=[-3 3];    zoo_lim=[-1 1];  oxy_lim=[-30 30];  v_lim=([-0.03 0.03]);
                      tic_lim=[-70 200];   po4_lim=[-1 1]; tdet_lim = [-1 1]; ldet_lim =[-1 1]; sdet_lim =[-1 1];
                      DIN_up_lim=[-20 20];
% temp_c  =[6:1:25];
plot_contour  = 1;    color_c  ='-k' ;      temp_c  =[-35:.5:35]; salt_c  =[-35:.1:35]; no3_c=[-20:.2:20]; nh4_c=[-30:.1:30];
                      chl_c=[-20:.2:20];    phy_c=[-10:.1:10];    zoo_c=[-1.3:0.1:3];        oxy_c=[-330:1:330];
                      tic_c=[-50:50:200]; po4_c=[-0.1:0.1:1];  den_c=[-30:.1:30]; v_c=[-10:.01:10];
                      v_c=[-10:0.05:10];  tdet_c = [-60:.5:60]; ldet_c = [-60:.5:60]; sdet_c = [-60:.5:60];
                      DIN_up_c=[-20:1:20];
% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'tif';

section = 0; % 0 => no plot, 1=> meridional, 2=> zonal

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';
  
% start_mm=1;   %04.29 - 119 08.05 - 217  09.04-247   12.01-335    01.12-377
% end_mm=12;
time_step=0;

%==========================================================================
%-- loop settings
% yy_range = 2021:10:2091;
% yy_range = 2091;
yy_ref = 3;
% variation_range = [11,12]; %  1 = temperature // 2 = salinity // 3 = NO3 // 4 = NH4 // 5 = chlorophyll
%               %  6 = phytoplankton // 7 = zooplankton // 8 = total inorganic carbon  
%               %  9 = oxygen // 10 = PO4 // 11 = Density // 12 = DIN //
%               % 13 = v-vel  // 14 = ubar & vbar // 15 = u & v surf and
%               % bottom // 16 = Ldet // 17 = Sdet // 18 = total = Ldet + Sdet
variation_range = [1:5,9,11:12,18];
% variation_range = [15];
switch_upper_river = 1;  % 0 = off, 1 = on; 
% no3_lim=[-1.5 1.5]; % diff no3
% no3_lim=[5 8];
% no3_lim=[5 15];
if section == 1
% % meridional
% pick_m = 93;
pick_m = 87;
else
% % zonal
%     pick = 76;
pick = 74;
end
%==========================================================================
max_level= 30;
out_filen = 'mp_p_sewer_det_10m_do.png';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for yy = yy_range
    for monthly_t = 1:12
    for variation = variation_range
        clearvars value value_ref
%     clearvars -except yy variation
%% 1.4.8. == 5/29 >? 1.4.9 .
%% 8.1. == 3/23 >? .8.2.
%% 7.2. == 3/14 >? 7.3.
%% 2.2.6. == 8/14
% file_dir=['\H:\E\포항모델\model\output\'];

%==========================================================================
% case_name='3reg_t4_cir_monthly_';
%% RCP
case_name=['3reg_RCP_',num2str(yy),'_v3_2yr_monthly_', num2str(monthly_t,'%02d')];
%% 2&3regime
case_name_ref=[num2str(yy_ref),'reg_cir_2yr_monthly_', num2str(monthly_t,'%02d')];

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
vname = {'temp','salt' ,'NO3','NH4','chlorophyll','phytoplankton','zooplankton','TIC','oxygen' 'tPO4'};%,'zeta','ubar','vbar','u','v','omega'};
     
% for im=start_mm:time_step:end_mm
% for im=start_mm:end_mm %pohang
% 		mid=[num2str(im,'%02d')];
%         file = [file_dir,'ocean_avg_',mid,'.nc'];
file = [fpath,case_name,'.nc'];
file_ref = [fpath,case_name_ref,'.nc'];
%% sigma2zeta
clearvars zeta depth
zeta=ncread(file,'zeta');
depth=zlevs(topo,zeta,theta_s,theta_b,hc,N,'r');
depth_w=zlevs(topo,zeta,theta_s,theta_b,hc,N,'w');

zeta_ref=ncread(file_ref,'zeta');
depth_ref=zlevs(topo,zeta_ref,theta_s,theta_b,hc,N,'r');
depth_w_ref=zlevs(topo,zeta_ref,theta_s,theta_b,hc,N,'w');

%% white color line range
white_range = [-1:2:1];
%%
% file = ['pohang_avg_',mid,'.nc'];
%         file = ['pohang_wt8_fine_flat_3m_',mid,'.nc'];

%% control case


%% exp case
        disp([file,' : ','monthly'])
        nc_ref=netcdf(file_ref);
        nc=netcdf(file);
        date=[num2str(yy)];
        switch variation
            case 1
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='Temperature';
                unit = '^oC';
                out_name_1=['vertical-GY-Temp-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=temp_lim;
                level_c=temp_c;
%                 value=value;
%                clear value;         
            case 2
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='Salinity';
                unit = 'psu';
                out_name_1=['vertical-GY-Salt-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=salt_lim; 
                level_c=salt_c;
%                 value=value;
                white_range = [31.499999:0.000001:31.5000001];
%                clear value;     
           case 3
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='Nitrate';
                unit = 'μMol N L^-1';
                out_name_1=['vertical-GY-NO3-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=no3_lim;
                level_c=no3_c;
%                 value=value; %.*14./1000; %% convert to mg/L
%                clear value;
           case 4
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='Ammonium';
                unit = 'μMol N L^-1';
                out_name_1=['vertical-GY-NH4-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=nh4_lim;
                level_c=nh4_c;
%                 value=value; %.*14./1000; %% convert to mg/L
%                clear value;              
           case 5
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='chlorophyll';
                unit = 'μg C L^-1';
                out_name_1=['vertical-GY-chl-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=chl_lim;
                level_c=chl_c;
%                 value=value;
%                clear value;  
           case 6
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='phytoplankton';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-phy-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=phy_lim;
                level_c=phy_c;
%                 value=value;
%                clear value;
           case 7
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='zootoplankton';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-zoo-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=zoo_lim;
                level_c=zoo_c;
%                 value=value;
%                clear value;
           case 8
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='TIC';
                unit ='μMol C L^-1';
                out_name_1=['vertical-GY-tic-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=tic_lim;
                level_c=tic_c;
%                 value=value;
%                clear value;
           case 9
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='dissolved oxygen';
                unit ='μMol L^-1';
                out_name_1=['vertical-GY-DO-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=oxy_lim;
                level_c=oxy_c;
%                 value=value; %.*(10/7).*(1/44.661); %% convert to mg/L 
%                clear value;
           case 10
                value=nc{char(vname(variation))}(:);
                value_ref=nc_ref{char(vname(variation))}(:);
                val_name='phosphate';
                unit ='μMol P L^-1';
                out_name_1=['vertical-GY-PO4-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=po4_lim;
                level_c=po4_c;
%                 value=value;
%                clear value;   
            case 11
                temp=nc{char(vname(1))}(:);salt=nc{char(vname(2))}(:);
                temp_ref=nc_ref{char(vname(1))}(:);salt_ref=nc_ref{char(vname(2))}(:);
                for N=1:1:20
                    value(N,:,:)=sw_dens(squeeze(salt(N,:,:)),squeeze(temp(N,:,:)),0)-1000;
                    value_ref(N,:,:)=sw_dens(squeeze(salt_ref(N,:,:)),squeeze(temp_ref(N,:,:)),0)-1000;
                end
                val_name='density';
                unit = 'σ(kg/m^3)';
                out_name_1=['vertical-GY-Den-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=den_lim;
                level_c=den_c;
%                 value=value;
%                clear value;   
            case 12
                white_range = [5.2:0.2:5.4];
                red_range = [5.0:0.2:5.2];
                value1=nc{char(vname(3))}(:); %no3
                value2=nc{char(vname(4))}(:); %nh4
                value1_ref=nc_ref{char(vname(3))}(:); %no3
                value2_ref=nc_ref{char(vname(4))}(:); %nh4
                value=value1 + value2;
                value_ref=value1_ref + value2_ref;
                val_name='DIN';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-DIN-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=no3_lim;
                level_c=no3_c;
                val_up_caxis=DIN_up_lim;
                level_up_c=DIN_up_c;
%                 value=value;
%                clear value;

            case 13
                v=nc{char('v')}(:);
                v_ref=nc_ref{char('v')}(:);
                val_name='v-velocity';
                unit ='m/s';
                out_name_1=['vertical-GY-v-vel-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=v_lim;
                level_c=v_c;
                [L M NN]=size(v);
                 value(:,2:M,:) = mean(cat(4,v(:,1:M-1,:),v(:,2:M,:)),[4],'omitnan');
                 value(:,1,:)   = v(:,1,:);
                 value(:,M+1,:) = v(:,M,:);
                [L M NN]=size(v_ref);
                 value_ref(:,2:M,:) = mean(cat(4,v_ref(:,1:M-1,:),v_ref(:,2:M,:)),[4],'omitnan');
                 value_ref(:,1,:)   = v_ref(:,1,:);
                 value_ref(:,M+1,:) = v_ref(:,M,:);
%                 value=value;
%                clear value;
            case 14
                value1_d=nc{char('ubar')}(:); value2_d=nc{char('vbar')}(:);
                value1_d_ref=nc_ref{char('ubar')}(:); value2_d_ref=nc_ref{char('vbar')}(:);
                %%
                [M N]=size(value1_d);
                u =ones(M,N+1)*NaN;
                value1_d(value1_d>100)=0;
                u(:,2:N) =mean(cat(3,value1_d(:,1:N-1),value1_d(:,2:N)),[3],'omitnan');
                u(:,1)   = value1_d(:,1);
                u(:,N+1) = value1_d(:,N);
                
                u_ref =ones(M,N+1)*NaN;
                value1_d_ref(value1_d_ref>100)=0;
                u_ref(:,2:N) =mean(cat(3,value1_d_ref(:,1:N-1),value1_d_ref(:,2:N)),[3],'omitnan');
                u_ref(:,1)   = value1_d_ref(:,1);
                u_ref(:,N+1) = value1_d_ref(:,N);
                
                [M N]=size(value2_d);
                v =ones(M+1,N)*NaN;
                value2_d(value2_d>100)=0;
                v(2:M,:) =mean(cat(3,value2_d(1:M-1,:),value2_d(2:M,:)),[3],'omitnan');
                v(1,:)   = value2_d(1,:);
                v(M+1,:) = value2_d(M,:);
                
                v_ref =ones(M+1,N)*NaN;
                value2_d_ref(value2_d_ref>100)=0;
                v_ref(2:M,:) =mean(cat(3,value2_d_ref(1:M-1,:),value2_d_ref(2:M,:)),[3],'omitnan');
                v_ref(1,:)   = value2_d_ref(1,:);
                v_ref(M+1,:) = value2_d_ref(M,:);               
                %%
%                 value = value1 + value2;
                val_name='velocity';
                unit ='m/s';
                out_name_1=['vertical-GY-velocity-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=v_lim;
                level_c=v_c;
                u=u'; v=v';
                u_ref=u_ref'; v_ref=v_ref';
                value=v - v_ref;
                interval = 3;
%                clear value; 

            case 15
                value1_d=nc{char('u')}(:); value2_d=nc{char('v')}(:);
                value1_d_ref=nc_ref{char('u')}(:); value2_d_ref=nc_ref{char('v')}(:);
                %%
                [L M N]=size(value1_d);
                u =ones(L,M,N+1)*NaN;
                value1_d(value1_d>100)=0;
            %     u(:,:,2:N) = ( u_dump(:,:,1:N-1)+u_dump(:,:,2:N) ) * 0.5;
                u(:,:,2:N) = mean(cat(4,value1_d(:,:,1:N-1),value1_d(:,:,2:N)),[4],'omitnan');
                u(:,:,1)   = value1_d(:,:,1);
                u(:,:,N+1) = value1_d(:,:,N);
                               
                [L M N]=size(value1_d_ref);
                u_ref =ones(L,M,N+1)*NaN;
                value1_d_ref(value1_d_ref>100)=0;
            %     u(:,:,2:N) = ( u_dump(:,:,1:N-1)+u_dump(:,:,2:N) ) * 0.5;
                u_ref(:,:,2:N) = mean(cat(4,value1_d_ref(:,:,1:N-1),value1_d_ref(:,:,2:N)),[4],'omitnan');
                u_ref(:,:,1)   = value1_d_ref(:,:,1);
                u_ref(:,:,N+1) = value1_d_ref(:,:,N);
                
                
                [L M N]=size(value2_d);
                v =ones(L,M+1,N)*NaN;
                value2_d(value2_d>100)=0;
            %     v(:,2:M,:) = ( v_dump(:,1:M-1,:)+v_dump(:,2:M,:) ) * 0.5;
                v(:,2:M,:) = mean(cat(4,value2_d(:,1:M-1,:),value2_d(:,2:M,:)),[4],'omitnan');
                v(:,1,:)   = value2_d(:,1,:);
                v(:,M+1,:) = value2_d(:,M,:);
                
                [L M N]=size(value2_d_ref);
                v_ref =ones(L,M+1,N)*NaN;
                value2_d_ref(value2_d_ref>100)=0;
            %     v(:,2:M,:) = ( v_dump(:,1:M-1,:)+v_dump(:,2:M,:) ) * 0.5;
                v_ref(:,2:M,:) = mean(cat(4,value2_d_ref(:,1:M-1,:),value2_d_ref(:,2:M,:)),[4],'omitnan');
                v_ref(:,1,:)   = value2_d_ref(:,1,:); 
                v_ref(:,M+1,:) = value2_d_ref(:,M,:);
                                          
                %%
%                 value = value1 + value2;
                val_name='velocity';
                unit ='m/s';
                out_name_1=['vertical-GY-velocity-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=v_lim;
                level_c=v_c;
                u= permute(u,[3 2 1]); v= permute(v,[3 2 1]);
                u_ref= permute(u_ref,[3 2 1]); v_ref=permute(v_ref,[3 2 1]);
                value = v - v_ref;
                interval = 3;
%                clear value; 

            case 16
%                 white_range = [5.2:0.2:5.4];
%                 red_range = [5.0:0.2:5.2];
                value=nc{char('LdetritusN')}(:); %no3
                value_ref=nc_ref{char('LdetritusN')}(:); %no3
                val_name='Ldet';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-Ldet-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=ldet_lim;
                level_c=ldet_c;
%                 value=value;
%                clear value;

            case 17
%                 white_range = [5.2:0.2:5.4];
%                 red_range = [5.0:0.2:5.2];
                value=nc{char('SdetritusN')}(:); %no3
                value_ref=nc_ref{char('SdetritusN')}(:); %no3
                val_name='Sdet';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-Sdet-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=sdet_lim;
                level_c=sdet_c;
%                 value=value;
%                clear value;

            case 18
%                 white_range = [5.2:0.2:5.4];
%                 red_range = [5.0:0.2:5.2];
                value1=nc{char('LdetritusN')}(:); %no3
                value2=nc{char('SdetritusN')}(:); %nh4
                value1_ref=nc_ref{char('LdetritusN')}(:); %no3
                value2_ref=nc_ref{char('SdetritusN')}(:); %nh4
                value=value1 + value2;
                value_ref=value1_ref + value2_ref;
                val_name='total-det';
                unit ='μMol N L^-1';
                out_name_1=['vertical-GY-total-det-',num2str(yy),'-',num2str(monthly_t,'%02d'),'mon-'];
                val_caxis=tdet_lim;
                level_c=tdet_c;
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
dep(find(dep > 10000))=NaN;

switch variation
            case {1,2,3,4,5,6,7,8,9,10,11,12,13,16,17,18}
                value = permute(value,[3 2 1]);
                value_ref = permute(value_ref,[3 2 1]);
end

%%% check line
% figure; hold on; 
% pcolor(lon_rho,lat_rho,squeeze(data_3(:,:,end)))
% shading flat; hold on; 
% plot(lon_rho(93,:),lat_rho(93,:),'r--','linew',2)
%%%%%%



if section == 1
    yLim = latLim;

    lat_rho_3 = repmat(lat_rho,1,1,N); %make it 3d
    
    pick_diff = squeeze( value(pick_m,:,:) - value_ref(pick_m,:,:) );
    
%         figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
%         set(gca,'Position',[0.2 0.15 0.73 0.75]);
%         figure('position',[400 100 550
%         550],'PaperUnits','inches','PaperPosition',[0 0 7.5 3.5]);
%         set(gca,'Position',[0.1 0.17 0.88 0.73]);
        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);
%         set(gcf,'PaperPositionMode','auto')     
        pcolor(squeeze(lat_rho_3(pick_m,:,:)),squeeze(dep(pick_m,:,:)),pick_diff); shading flat;
        xlim(yLim);
%         ylim([-5 0]);
        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
        xlabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
        ylabel('Depth(m)','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
        caxis(val_caxis)
        out_name=[out_name_1];
        title(out_name)
        if (plot_contour)
                  hold on
                  [C,h]=contour(squeeze(lat_rho_3(pick_m,:,:)),squeeze(dep(pick_m,:,:)),pick_diff,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(squeeze(lat_rho_3(pick_m,:,:)),squeeze(dep(pick_m,:,:)),pick_diff,white_range,'-w','linewidth',1);      
                  if exist('red_range','var') == 1
                      [C3,h3]=contour(squeeze(lat_rho_3(pick_m,:,:)),squeeze(dep(pick_m,:,:)),pick_diff,red_range,'-r','linewidth',1);
                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                  end
                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
        end
        colormap(bluewhitered);
        set(gca,'Color',[210/255 180/255 140/255]); %background color
        set(gcf,'InvertHardCopy','Off'); %remove default "white background"
            
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
%                 saveas(gca,[out_name,'.png']);
                print(gcf,['figures\diff\monthly\',out_name,'.png'],'-dpng','-r200');
%                 close all;
            end
elseif section == 2
    xLim = lonLim;

    lon_rho_3 = repmat(lon_rho,1,1,N); %make it 3d
    
    pick_diff = squeeze( value(:,pick,:) - value_ref(:,pick,:) );

%         figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
%         set(gca,'Position',[0.2 0.15 0.73 0.75]);
%         figure('position',[400 100 550
%         550],'PaperUnits','inches','PaperPosition',[0 0 7.5 3.5]);
%         set(gca,'Position',[0.1 0.17 0.88 0.73]);
        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
%         set(gcf,'PaperPositionMode','auto')     
        pcolor(squeeze(lon_rho_3(:,pick,:)),squeeze(dep(:,pick,:)),pick_diff); shading flat;
        xlim(xLim);
        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
        ylabel('Depth(m)','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
        caxis(val_caxis)
        out_name=[out_name_1];
        title(out_name)
        if (plot_contour)
                  hold on
                  [C,h]=contour(squeeze(lon_rho_3(:,pick,:)),squeeze(dep(:,pick,:)),pick_diff,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(squeeze(lon_rho_3(:,pick,:)),squeeze(dep(:,pick,:)),pick_diff,white_range,'-w','linewidth',1);      
                  if exist('red_range','var') == 1
                      [C3,h3]=contour(squeeze(lon_rho_3(:,pick,:)),squeeze(dep(:,pick,:)),pick_diff,red_range,'-r','linewidth',1);
                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                  end
                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
        end
        colormap(bluewhitered);
        set(gca,'Color',[210/255 180/255 140/255]); %background color
        set(gcf,'InvertHardCopy','Off'); %remove default "white background"
        
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
%                 saveas(gca,[out_name,'_meridional','.png']);
                print(gcf,['figures\diff\monthly\',out_name,'_meridional','.png'],'-dpng','-r200');
%                 close all;
            end   
    
            
end

% end


%% vertical integration
switch variation
            case {1,2,3,4,5,6,7,8,9,10,11,12,13,16,17,18}
                
dep_w = permute(depth_w,[2 3 1]);
dep_w_ref = permute(depth_w_ref,[2 3 1]);


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

for i = N:-1:1
    if i == N
        dep_l_ref(:,:,i) =  dep_w_ref(:,:,i+1) - dep_w_ref(:,:,i);
    else
        dep_l_ref(:,:,i) =  dep_w_ref(:,:,i+1) - dep_w_ref(:,:,i); 
    end
end

depcval= dep_l .* value;
int_value =  sum(depcval,3) ./ sum(dep_l,3);

depcval_ref= dep_l_ref .* value_ref;
int_value_ref =  sum(depcval_ref,3) ./ sum(dep_l_ref,3);

int_diff = int_value - int_value_ref;

% pcolor(int_chl'); shading flat; colorbar

figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
%         set(gcf,'PaperPositionMode','auto')     
        pcolor(lon_rho,lat_rho,int_diff); shading flat;
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
        out_name=[out_name_1];
        title(out_name)
        if (plot_contour)
                  hold on
                  [C,h]=contour(lon_rho,lat_rho,int_diff,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(lon_rho,lat_rho,int_diff,white_range,'-w','linewidth',1);      
                  if exist('red_range','var') == 1
                      [C3,h3]=contour(lon_rho,lat_rho,int_diff,red_range,'-r','linewidth',1);
                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                  end
                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
        end
        colormap(bluewhitered);
        set(gca,'Color',[210/255 180/255 140/255]); %background color
        set(gcf,'InvertHardCopy','Off'); %remove default "white background"
            
            if (switch_save)
%                 saveas(gcf,out_name,out_type);
%                 saveas(gca,[out_name,'_meridional','.png']);
                print(gcf,['figures\diff\monthly\',out_name,'_vertical_integral','.png'],'-dpng','-r200');
                close all;
            end
            


        case {14} %velocity (Ubar and Vbar)
                 
                %diff
                diff_u = u - u_ref;
                diff_v = v - v_ref;
                
                % referenc vector
                    scale_vec = 10;
                    clearvars *0_ref
                    u0_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    v0_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    u0_ref(11,120) = 0.01 * scale_vec;

                        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
%                         quiver(lon_rho,lat_rho,u',v'); 
                        pcolor(lon_rho,lat_rho,diff_v .*(mask_rho ./ mask_rho)); shading flat; hold on;
                        quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end), ...
    squeeze(diff_u(1:interval:end,1:interval:end) .* scale_vec), ...
    squeeze(diff_v(1:interval:end,1:interval:end) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
                        quiver(lon_rho,lat_rho,u0_ref, v0_ref, 'AutoScale','off','LineWidth',1,'Color','r');
                        text(lon_rho(11,115),lat_rho(11,115),'0.01m/s','fontsize',15,'Color','r');
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
                        out_name=[out_name_1];
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
                                print(gcf,['figures\diff\monthly\',out_name,'_vertical_integral','.png'],'-dpng','-r200');
                                close all;
                            end 
                            
                case {15} %velocity (U and V) surf and bottom
                 
                %diff
                diff_u_top = squeeze(u(:,:,20) - u_ref(:,:,20));
                diff_v_top = squeeze(v(:,:,20) - v_ref(:,:,20));
                diff_u_bot = squeeze(u(:,:,1) - u_ref(:,:,1));
                diff_v_bot = squeeze(v(:,:,1) - v_ref(:,:,1));
                
                % referenc vector
                    scale_vec = 1;
                    clearvars *0_ref
                    u0_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    v0_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    u0_ref(11,120) = 0.01 * scale_vec;

                        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
%                         quiver(lon_rho,lat_rho,u',v'); 
                        pcolor(lon_rho,lat_rho,diff_v_top .*(mask_rho ./ mask_rho)); shading flat; hold on;
                        quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end), ...
    squeeze(diff_u_top(1:interval:end,1:interval:end) .* scale_vec), ...
    squeeze(diff_v_top(1:interval:end,1:interval:end) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
                        quiver(lon_rho,lat_rho,u0_ref, v0_ref, 'AutoScale','off','LineWidth',1,'Color','r');
                        text(lon_rho(11,115),lat_rho(11,115),'0.01m/s','fontsize',15,'Color','r');
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
                        out_name=[out_name_1];
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
                                print(gcf,['figures\diff\monthly\',out_name,'_surf','.png'],'-dpng','-r200');
                                close all;
                            end 
                            
                         figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
%                         quiver(lon_rho,lat_rho,u',v'); 
                        pcolor(lon_rho,lat_rho,diff_v_bot .*(mask_rho ./ mask_rho)); shading flat; hold on;
                        quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end), ...
    squeeze(diff_u_bot(1:interval:end,1:interval:end) .* scale_vec), ...
    squeeze(diff_v_bot(1:interval:end,1:interval:end) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
                        quiver(lon_rho,lat_rho,u0_ref, v0_ref, 'AutoScale','off','LineWidth',1,'Color','r');
                        text(lon_rho(11,115),lat_rho(11,115),'0.01m/s','fontsize',15,'Color','r');
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
                        out_name=[out_name_1];
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
                                print(gcf,['figures\diff\monthly\',out_name,'_bot','.png'],'-dpng','-r200');
                                close all;
                            end 
end    

   %% switch to turn on for drawing upper river plot
if switch_upper_river == 1
    switch variation
        case {1,2,3,4,5,6,7,8,9,10,11,12,13,16,17,18}
                
dep_w = permute(depth_w,[2 3 1]);
dep_w_ref = permute(depth_w_ref,[2 3 1]);

for i = N:-1:1
    if i == N
        dep_l(:,:,i) =  dep_w(:,:,i+1) - dep_w(:,:,i);
    else
        dep_l(:,:,i) =  dep_w(:,:,i+1) - dep_w(:,:,i); 
    end
end

for i = N:-1:1
    if i == N
        dep_l_ref(:,:,i) =  dep_w_ref(:,:,i+1) - dep_w_ref(:,:,i);
    else
        dep_l_ref(:,:,i) =  dep_w_ref(:,:,i+1) - dep_w_ref(:,:,i); 
    end
end

depcval= dep_l .* value;
int_value =  sum(depcval,3) ./ sum(dep_l,3);

depcval_ref= dep_l_ref .* value_ref;
int_value_ref =  sum(depcval_ref,3) ./ sum(dep_l_ref,3);

int_diff = int_value - int_value_ref;

% pcolor(int_chl'); shading flat; colorbar

figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
%         set(gcf,'PaperPositionMode','auto')     
        pcolor(lon_rho,lat_rho,int_diff); shading flat;
%         xlim(xLim);
        ylim([34.95 35.1])
        xlim([127.63 127.793])
        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
%         out_name_1=[out_name_1];
         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')         
         if exist('val_up_caxis') == 1
             caxis(val_up_caxis)
             level_c = level_up_c;
         else
             caxis(val_caxis)
         end
        out_name=[out_name_1];
        title(out_name)
        if (plot_contour)
                  hold on
                  [C,h]=contour(lon_rho,lat_rho,int_diff,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(lon_rho,lat_rho,int_diff,white_range,'-w','linewidth',1);      
                  if exist('red_range','var') == 1
                      [C3,h3]=contour(lon_rho,lat_rho,int_diff,red_range,'-r','linewidth',1);
                      clabel(C3,h3,'FontSize',fontSizeCLab,'Color','r','labelspacing',100000,'fontweight','bold');
                  end
                  clabel(C,h,'FontSize',fontSizeCLab,'Color','k','labelspacing',100000,'fontweight','bold');
                  clabel(C2,h2,'FontSize',fontSizeCLab,'Color','w','labelspacing',100000,'fontweight','bold');
        end
        colormap(bluewhitered);
        set(gca,'Color',[210/255 180/255 140/255]); %background color
        set(gcf,'InvertHardCopy','Off'); %remove default "white background"
            
            if (switch_save)
                print(gcf,['figures\diff\monthly\',out_name,'_vertical_integral_upper_river','.png'],'-dpng','-r200');
                close all;
            end
        
        
    case {15} %velocity (U and V) surf and bottom
     %diff
                diff_u_top = squeeze(u(:,:,20) - u_ref(:,:,20));
                diff_v_top = squeeze(v(:,:,20) - v_ref(:,:,20));
                diff_u_bot = squeeze(u(:,:,1) - u_ref(:,:,1));
                diff_v_bot = squeeze(v(:,:,1) - v_ref(:,:,1));
                
                % referenc vector
                    scale_vec = 5;
                    clearvars *0_ref
                    u0_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    v0_ref = zeros(size(lon_rho,1),size(lat_rho,2));
                    u0_ref(20,155) = 0.01 * scale_vec;

                        figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
%                         quiver(lon_rho,lat_rho,u',v'); 
                        pcolor(lon_rho,lat_rho,diff_v_top .*(mask_rho ./ mask_rho)); shading flat; hold on;
                        quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end), ...
    squeeze(diff_u_top(1:interval:end,1:interval:end) .* scale_vec), ...
    squeeze(diff_v_top(1:interval:end,1:interval:end) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
                        quiver(lon_rho,lat_rho,u0_ref, v0_ref, 'AutoScale','off','LineWidth',1,'Color','r');
                        text(lon_rho(20,150),lat_rho(20,150),'0.01m/s','fontsize',15,'Color','r');
                %         xlim(xLim);
        ylim([34.95 35.1])
        xlim([127.63 127.793])
                        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
                        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                %         out_name_1=[out_name_1];
                         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
                         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
                        caxis(val_caxis)
                        out_name=[out_name_1];
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
                                print(gcf,['figures\diff\monthly\',out_name,'_surf_upper_river','.png'],'-dpng','-r200');
                                close all;
                            end 
                            
                         figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)
                %         set(gcf,'PaperPositionMode','auto')     
%                         quiver(lon_rho,lat_rho,u',v'); 
                        pcolor(lon_rho,lat_rho,diff_v_bot .*(mask_rho ./ mask_rho)); shading flat; hold on;
                        quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end), ...
    squeeze(diff_u_bot(1:interval:end,1:interval:end) .* scale_vec), ...
    squeeze(diff_v_bot(1:interval:end,1:interval:end) .* scale_vec), ...
    'AutoScale','off','LineWidth',0.5,'Color','k');
                        quiver(lon_rho,lat_rho,u0_ref, v0_ref, 'AutoScale','off','LineWidth',1,'Color','r');
                        text(lon_rho(20,150),lat_rho(20,150),'0.01m/s','fontsize',15,'Color','r');
                %         xlim(xLim);
        ylim([34.95 35.1])
        xlim([127.63 127.793])
                        set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
                        xlabel('lon','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                        ylabel('lat','color',color_v,'FontSize',fontSizeLab,'fontweight','bold')
                %         out_name_1=[out_name_1];
                         bar = colorbar('fontsize',fontSizeLab,'fontweight','bold');
                         set(get(bar,'title'),'string',unit,'FontSize',fontSizeLab,'fontweight','bold')
                        caxis(val_caxis)
                        out_name=[out_name_1];
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
                                print(gcf,['figures\diff\monthly\',out_name,'_bot_upper_river','.png'],'-dpng','-r200');
                                close all;
                            end 
    end
end
                            
                            
                            
end



    end
end
