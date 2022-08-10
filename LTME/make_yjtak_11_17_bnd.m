close all; clear; clc; 

% for i = 2001:2010
for i = 2001:2001
    for j = 1:12
        clearvars -except i j d3_*_s_t d3_*_e_t
        path = ['D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary\clim_11_17\']
        cd(path)
    %load raw file
%     rawfile = ['stdep1_avg_mon_', num2str(i),'_', num2str(j,'%02d'),'.nc'];
    rawfile = ['clim_monthly' num2str(j,'%02d') '.nc'];
    raw_temp=ncread(rawfile,'temp');
    raw_salt=ncread(rawfile,'salt');
    raw_u=ncread(rawfile,'u');
    raw_v=ncread(rawfile,'v');
    raw_chl=ncread(rawfile,'chlorophyll');
    raw_phy=ncread(rawfile,'phytoplankton');
    raw_zoo=ncread(rawfile,'zooplankton');
    raw_no3=ncread(rawfile,'NO3');
    raw_nh4=ncread(rawfile,'NH4');
    raw_do=ncread(rawfile,'oxygen');
%srho to depth
    %   raw_dep=ncread(rawfile,'depth');
    raw_h=ncread(rawfile,'h');
    raw_zeta=ncread(rawfile,'zeta');
    raw_srho=ncread(rawfile,'s_rho');
    raw_thetas=ncread(rawfile,'theta_s');
    raw_thetab=ncread(rawfile,'theta_b');
    raw_hc = ncread(rawfile,'hc');
    raw_dep = zlevs(raw_h, raw_zeta, raw_thetas, raw_thetab, raw_hc, length(raw_srho), 'r'); 
    raw_dep = permute(raw_dep,[2 3 1]);
    
    raw_dep_w = zlevs(raw_h, raw_zeta, raw_thetas, raw_thetab, raw_hc, length(raw_srho), 'u'); 
    raw_dep_w = permute(raw_dep_w,[2 3 1]);
    %%
    
    raw_lon=ncread(rawfile,'lon_rho');
    raw_lat=ncread(rawfile,'lat_rho');
    raw_lonu=ncread(rawfile,'lon_u');
    raw_latu=ncread(rawfile,'lat_u');
    raw_lonv=ncread(rawfile,'lon_v');
    raw_latv=ncread(rawfile,'lat_v');
    
    for i = 1:length(raw_srho)
        raw_lon_m(:,:,i)=raw_lon;
        raw_lat_m(:,:,i)=raw_lat;
        raw_lonu_m(:,:,i)=raw_lonu;
        raw_latu_m(:,:,i)=raw_latu;
        raw_lonv_m(:,:,i)=raw_lonv;
        raw_latv_m(:,:,i)=raw_latv;
    end
    
%     [raw_lat_m raw_lon_m]=meshgrid(raw_lat,raw_lon);
%     [raw_lat_m raw_lon_m]=meshgrid(raw_latu,raw_lonu);

    %load ref bnd file
    cd D:\장기생태\Dynamic\02_grid_depth
    bnd_reffile = 'grid_gy_v11_s.nc';
    ref_lon=ncread(bnd_reffile,'lon_rho');
    ref_lat=ncread(bnd_reffile,'lat_rho');
    ref_lonu=ncread(bnd_reffile,'lon_u');
    ref_latu=ncread(bnd_reffile,'lat_u');
    ref_lonv=ncread(bnd_reffile,'lon_v');
    ref_latv=ncread(bnd_reffile,'lat_v');
    
    ref_lon(:,1)
    
    cd D:\장기생태\Dynamic\result\2013\Input
    reffile = 'ocean_his_8785.nc';
    ref_srho=ncread(reffile,'s_rho');
    ref_hc=ncread(reffile,'hc');
    ref_theta_s=ncread(reffile,'theta_s');
    ref_theta_b=ncread(reffile,'theta_b');
    ref_h=ncread(reffile,'h');
    ref_temp=ncread(reffile,'temp');
    
    ref_dep = zlevs(ref_h, zeros(size(ref_h,1),size(ref_h,2)),ref_theta_s,ref_theta_b,ref_hc,length(ref_srho),'r');
    ref_dep = permute(ref_dep,[2 3 1]);
    
    for i = 1:length(ref_srho)
        ref_lon_m(:,:,i)=ref_lon;
        ref_lat_m(:,:,i)=ref_lat;
        ref_lonu_m(:,:,i)=ref_lonu;
        ref_latu_m(:,:,i)=ref_latu;
        ref_lonv_m(:,:,i)=ref_lonv;
        ref_latv_m(:,:,i)=ref_latv;
    end
    
%     figure; pcolor(squeeze(sum(ref_dep,1))'); colorbar
    figure; pcolor(ref_lon,ref_lat,squeeze(ref_dep(:,:,1))); colorbar; hold on;
    shading flat;
    plot(ref_lon(:,1),ref_lat(:,1),'ro');
    plot(ref_lon(end,:),ref_lat(end,:),'ko');
    
    ref_lon_s = ref_lon(:,1); ref_lat_s = ref_lat(:,1); % south boundary
    ref_lonu_s = ref_lonu(:,1); ref_latv_s = ref_latv(:,1); % south boundary vel.
    ref_lon_e = ref_lon(end,:); ref_lat_e = ref_lat(end,:); % east boundary
    ref_lonu_e = ref_lonu(end,:); ref_latv_e = ref_latv(end,:); % east boundary vel.
    
    figure; pcolor(ref_lon,ref_lat,squeeze(ref_dep(:,:,1))); colorbar; hold on;
    shading flat;
%     plot(ref_lon_s,ref_lat_s,'ro'); plot(ref_lon_e,ref_lat_e,'ko');
%       plot(ref_lonu_s,ref_latv_s,'ro'); plot(ref_lonu_e,ref_latv_e,'ko');
      
    ref_lon_s = ref_lon(:,1); ref_dep_s = squeeze(ref_dep(:,1,:)); ref_lat_s = ref_lat(:,1); % south boundary
    ref_lat_e = ref_lat(end,:); ref_dep_e = squeeze(ref_dep(end,:,:)); ref_lon_e = ref_lon(end,:); % east boundary
    
%     [ref_lat_sm ref_lon_sm]= meshgrid(ref_lat_s,ref_lon_s);
    
    % confirm plot boundary side
%     ref_temp(ref_temp > 10000) =NaN;
%     ref_lon_sm=repmat(ref_lon_s,1,20);
%     pcolor(ref_lon_sm,ref_dep_s,squeeze(ref_temp(:,1,:))); colorbar; shading flat;
%     
%     ref_lat_em=repmat(ref_lat_e',1,20);
%     pcolor(ref_lat_em,ref_dep_e,squeeze(ref_temp(end,:,:))); colorbar; shading flat;    
    
    %% 
    
    %interp 2d
    raw_temp(raw_temp>(10^20))=NaN;
    raw_salt(raw_salt>(10^20))=NaN;
    raw_u(raw_u>(10^20))=NaN;
    raw_v(raw_v>(10^20))=NaN;
    raw_no3(raw_no3>(10^20))=NaN;
    raw_nh4(raw_nh4>(10^20))=NaN;
    raw_zoo(raw_zoo>(10^20))=NaN;
    raw_phy(raw_phy>(10^20))=NaN;
    raw_chl(raw_chl>(10^20))=NaN;
    raw_do(raw_do>(10^20))=NaN;
    
%     for k = 1:length(raw_dep)
%         %south
%         d2_temp_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_temp(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_salt_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_salt(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_u_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_u(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_v_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_v(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_no3_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_no3(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_nh4_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_nh4(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_zoo_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_zoo(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_phy_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_phy(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_chl_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_chl(:,:,k)),ref_lon_s,ref_lat_s);
%         %east
%         d2_temp_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_temp(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_salt_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_salt(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_u_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_u(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_v_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_v(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_no3_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_no3(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_nh4_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_nh4(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_zoo_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_zoo(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_phy_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_phy(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_chl_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_chl(:,:,k)),ref_lon_e,ref_lat_e);
%     end

    
%     %interp 3d
% 
%      % south
%      ref_lon_sm=repmat(ref_lon_s,1,20); ref_lat_sm=repmat(ref_lat_s,1,20);
%      raw_lon_sm=repmat(ref_lon_s,1,length(raw_dep)); raw_dep_sm = repmat(raw_dep',length(ref_lon_s),1);
%      
%      d3_temp_s = griddata(raw_lon_sm,raw_dep_sm,d2_temp_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_salt_s = griddata(raw_lon_sm,raw_dep_sm,d2_salt_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_u_s = griddata(raw_lon_sm,raw_dep_sm,d2_u_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_v_s = griddata(raw_lon_sm,raw_dep_sm,d2_v_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_no3_s = griddata(raw_lon_sm,raw_dep_sm,d2_no3_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_nh4_s = griddata(raw_lon_sm,raw_dep_sm,d2_nh4_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_zoo_s = griddata(raw_lon_sm,raw_dep_sm,d2_zoo_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_phy_s = griddata(raw_lon_sm,raw_dep_sm,d2_phy_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_chl_s = griddata(raw_lon_sm,raw_dep_sm,d2_chl_s,ref_lon_sm,ref_dep_s,'nearest');
%      
%      % east
%      ref_lon_em=repmat(ref_lon_e',1,20); ref_lat_em=repmat(ref_lat_e',1,20);
%      raw_lat_em=repmat(ref_lat_e',1,length(raw_dep)); raw_dep_em = repmat(raw_dep',length(ref_lon_e),1);
%      
%      d3_temp_e = griddata(raw_lat_em,raw_dep_em,d2_temp_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_salt_e = griddata(raw_lat_em,raw_dep_em,d2_salt_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_u_e = griddata(raw_lat_em,raw_dep_em,d2_u_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_v_e = griddata(raw_lat_em,raw_dep_em,d2_v_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_no3_e = griddata(raw_lat_em,raw_dep_em,d2_no3_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_nh4_e = griddata(raw_lat_em,raw_dep_em,d2_nh4_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_zoo_e = griddata(raw_lat_em,raw_dep_em,d2_zoo_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_phy_e = griddata(raw_lat_em,raw_dep_em,d2_phy_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_chl_e = griddata(raw_lat_em,raw_dep_em,d2_chl_e,ref_lat_em,ref_dep_e,'nearest');

%     %interp 3d
% 
     % south
     ref_lon_sm=repmat(ref_lon_s,1,20); ref_lat_sm=repmat(ref_lat_s,1,20);
%      raw_lon_sm=repmat(ref_lon_s,1,length(raw_dep)); raw_dep_sm = repmat(raw_dep',length(ref_lon_s),1);
find(isnan(raw_lon_m)==1)
find(isnan(raw_lat_m)==1)
find(isnan(raw_lonu_m)==1)
find(isnan(raw_latu_m)==1)
find(isnan(raw_lonv_m)==1)
find(isnan(raw_latv_m)==1)
length(find(isnan(raw_dep)==1))
pcolor(raw_lon,raw_lat,squeeze(raw_dep(:,:,1))); colorbar; shading flat;  %there is NaN on the land
raw_dep(find(isnan(raw_dep)==1))=0;

pcolor(raw_lon,raw_lat,squeeze(raw_temp(:,:,1))); colorbar; shading flat;


find(isnan(ref_lon_m)==1)
find(isnan(ref_lat_m)==1)
find(isnan(ref_lonu_m)==1)
find(isnan(ref_latu_m)==1)
find(isnan(ref_lonv_m)==1)
find(isnan(ref_latv_m)==1)
find(isnan(ref_dep)==1)

for i = 1:length(raw_srho)
    raw_depu(:,:,i)=griddata(raw_lon,raw_lat,raw_dep(:,:,i),raw_lonu,raw_latu);
    raw_depv(:,:,i)=griddata(raw_lon,raw_lat,raw_dep(:,:,i),raw_lonv,raw_latv);
end

for i = 1:length(ref_srho)
    ref_depu(:,:,i)=griddata(ref_lon,ref_lat,ref_dep(:,:,i),ref_lonu,ref_latu);
    ref_depv(:,:,i)=griddata(ref_lon,ref_lat,ref_dep(:,:,i),ref_lonv,ref_latv);
end

degtometer=111000; % 1deg = 111km = 111000m

     d3_temp_s = griddata(raw_lon_m.*degtometer ,raw_lat_m.*degtometer,raw_dep,raw_temp,ref_lon_m.*degtometer,ref_lat_m.*degtometer,ref_dep);
     d3_salt_s = griddata(raw_lon_m.*degtometer,raw_lat_m.*degtometer,raw_dep,raw_salt,ref_lon_m.*degtometer,ref_lat_m.*degtometer,ref_dep);
     d3_u_s = griddata(raw_lonu_m.*degtometer,raw_latu_m.*degtometer,raw_depu,raw_u,ref_lonu_m.*degtometer,ref_latu_m.*degtometer,ref_depu);
     d3_v_s = griddata(raw_lonv_m.*degtometer,raw_latv_m.*degtometer,raw_depv,raw_v,ref_lonv_m.*degtometer,ref_latv_m.*degtometer,ref_depv);
     d3_no3_s = griddata(raw_lon_m.*degtometer,raw_lat_m.*degtometer,raw_dep,raw_no3,ref_lon_m.*degtometer,ref_lat_m.*degtometer,ref_dep);
     d3_nh4_s = griddata(raw_lon_m.*degtometer,raw_lat_m.*degtometer,raw_dep,raw_nh4,ref_lon_m.*degtometer,ref_lat_m.*degtometer,ref_dep);
     d3_zoo_s = griddata(raw_lon_m.*degtometer,raw_lat_m.*degtometer,raw_dep,raw_zoo,ref_lon_m.*degtometer,ref_lat_m.*degtometer,ref_dep);
     d3_phy_s = griddata(raw_lon_m.*degtometer,raw_lat_m.*degtometer,raw_dep,raw_phy,ref_lon_m.*degtometer,ref_lat_m.*degtometer,ref_dep);
     d3_chl_s = griddata(raw_lon_m.*degtometer,raw_lat_m.*degtometer,raw_dep,raw_chl,ref_lon_m.*degtometer,ref_lat_m.*degtometer,ref_dep);
     d3_do_s = griddata(raw_lon_m.*degtometer,raw_lat_m.*degtometer,raw_dep,raw_do,ref_lon_m.*degtometer,ref_lat_m.*degtometer,ref_dep);
    
     
     %fill missing value

     % south
     d3_temp_s(isnan(d3_temp_s)) = griddata(ref_lon_m(~isnan(d3_temp_s)).*degtometer ,ref_lat_m(~isnan(d3_temp_s)).*degtometer,ref_dep(~isnan(d3_temp_s)),d3_temp_s(~isnan(d3_temp_s)),ref_lon_m(isnan(d3_temp_s)).*degtometer,ref_lat_m(isnan(d3_temp_s)).*degtometer,ref_dep(isnan(d3_temp_s)),'nearest');
     d3_salt_s(isnan(d3_salt_s)) = griddata(ref_lon_m(~isnan(d3_salt_s)).*degtometer,ref_lat_m(~isnan(d3_salt_s)).*degtometer,ref_dep(~isnan(d3_salt_s)),d3_salt_s(~isnan(d3_salt_s)),ref_lon_m(isnan(d3_salt_s)).*degtometer,ref_lat_m(isnan(d3_salt_s)).*degtometer,ref_dep(isnan(d3_salt_s)),'nearest');
     d3_u_s(isnan(d3_u_s)) = griddata(ref_lonu_m(~isnan(d3_u_s)).*degtometer,ref_latu_m(~isnan(d3_u_s)).*degtometer,ref_depu(~isnan(d3_u_s)),d3_u_s(~isnan(d3_u_s)),ref_lonu_m(isnan(d3_u_s)).*degtometer,ref_latu_m(isnan(d3_u_s)).*degtometer,ref_depu(isnan(d3_u_s)),'nearest');
     d3_v_s(isnan(d3_v_s)) = griddata(ref_lonv_m(~isnan(d3_v_s)).*degtometer,ref_latv_m(~isnan(d3_v_s)).*degtometer,ref_depv(~isnan(d3_v_s)),d3_v_s(~isnan(d3_v_s)),ref_lonv_m(isnan(d3_v_s)).*degtometer,ref_latv_m(isnan(d3_v_s)).*degtometer,ref_depv(isnan(d3_v_s)),'nearest');
     d3_no3_s(isnan(d3_no3_s)) = griddata(ref_lon_m(~isnan(d3_no3_s)).*degtometer,ref_lat_m(~isnan(d3_no3_s)).*degtometer,ref_dep(~isnan(d3_no3_s)),d3_no3_s(~isnan(d3_no3_s)),ref_lon_m(isnan(d3_no3_s)).*degtometer,ref_lat_m(isnan(d3_no3_s)).*degtometer,ref_dep(isnan(d3_no3_s)),'nearest');
     d3_nh4_s(isnan(d3_nh4_s)) = griddata(ref_lon_m(~isnan(d3_nh4_s)).*degtometer,ref_lat_m(~isnan(d3_nh4_s)).*degtometer,ref_dep(~isnan(d3_nh4_s)),d3_nh4_s(~isnan(d3_nh4_s)),ref_lon_m(isnan(d3_nh4_s)).*degtometer,ref_lat_m(isnan(d3_nh4_s)).*degtometer,ref_dep(isnan(d3_nh4_s)),'nearest');
     d3_zoo_s(isnan(d3_zoo_s)) = griddata(ref_lon_m(~isnan(d3_zoo_s)).*degtometer,ref_lat_m(~isnan(d3_zoo_s)).*degtometer,ref_dep(~isnan(d3_zoo_s)),d3_zoo_s(~isnan(d3_zoo_s)),ref_lon_m(isnan(d3_zoo_s)).*degtometer,ref_lat_m(isnan(d3_zoo_s)).*degtometer,ref_dep(isnan(d3_zoo_s)),'nearest');
     d3_phy_s(isnan(d3_phy_s)) = griddata(ref_lon_m(~isnan(d3_phy_s)).*degtometer,ref_lat_m(~isnan(d3_phy_s)).*degtometer,ref_dep(~isnan(d3_phy_s)),d3_phy_s(~isnan(d3_phy_s)),ref_lon_m(isnan(d3_phy_s)).*degtometer,ref_lat_m(isnan(d3_phy_s)).*degtometer,ref_dep(isnan(d3_phy_s)),'nearest');
     d3_chl_s(isnan(d3_chl_s)) = griddata(ref_lon_m(~isnan(d3_chl_s)).*degtometer,ref_lat_m(~isnan(d3_chl_s)).*degtometer,ref_dep(~isnan(d3_chl_s)),d3_chl_s(~isnan(d3_chl_s)),ref_lon_m(isnan(d3_chl_s)).*degtometer,ref_lat_m(isnan(d3_chl_s)).*degtometer,ref_dep(isnan(d3_chl_s)),'nearest');
     d3_do_s(isnan(d3_do_s)) = griddata(ref_lon_m(~isnan(d3_do_s)).*degtometer,ref_lat_m(~isnan(d3_do_s)).*degtometer,ref_dep(~isnan(d3_do_s)),d3_do_s(~isnan(d3_do_s)),ref_lon_m(isnan(d3_do_s)).*degtometer,ref_lat_m(isnan(d3_do_s)).*degtometer,ref_dep(isnan(d3_do_s)),'nearest');    
    
%      d3_temp_s(isnan(d3_temp_s)) = griddata(ref_lon_m(~isnan(d3_temp_s)).*degtometer,ref_lat_m(~isnan(d3_temp_s)).*degtometer,ref_dep_s(~isnan(d3_temp_s)),d3_temp_s(~isnan(d3_temp_s)),ref_lon_m(isnan(d3_temp_s)),ref_dep_s(isnan(d3_temp_s)),'nearest');
%      d3_salt_s(isnan(d3_salt_s)) = griddata(ref_lon_m(~isnan(d3_salt_s)).*degtometer,ref_lat_m(~isnan(d3_temp_s)).*degtometer,ref_dep_s(~isnan(d3_salt_s)),d3_salt_s(~isnan(d3_salt_s)),ref_lon_m(isnan(d3_salt_s)),ref_dep_s(isnan(d3_salt_s)),'nearest');
%      d3_u_s(isnan(d3_u_s)) = griddata(ref_lon_sm(~isnan(d3_u_s)).*degtometer,ref_dep_s(~isnan(d3_u_s)),d3_u_s(~isnan(d3_u_s)),ref_lon_sm(isnan(d3_u_s)),ref_dep_s(isnan(d3_u_s)),'nearest');
%      d3_v_s(isnan(d3_v_s)) = griddata(ref_lon_sm(~isnan(d3_v_s)).*degtometer,ref_dep_s(~isnan(d3_v_s)),d3_v_s(~isnan(d3_v_s)),ref_lon_sm(isnan(d3_v_s)),ref_dep_s(isnan(d3_v_s)),'nearest');
%      d3_no3_s(isnan(d3_no3_s)) = griddata(ref_lon_m(~isnan(d3_no3_s)).*degtometer,ref_dep_s(~isnan(d3_no3_s)),d3_no3_s(~isnan(d3_no3_s)),ref_lon_m(isnan(d3_no3_s)),ref_dep_s(isnan(d3_no3_s)),'nearest');
%      d3_nh4_s(isnan(d3_nh4_s)) = griddata(ref_lon_m(~isnan(d3_nh4_s)).*degtometer,ref_dep_s(~isnan(d3_nh4_s)),d3_nh4_s(~isnan(d3_nh4_s)),ref_lon_m(isnan(d3_nh4_s)),ref_dep_s(isnan(d3_nh4_s)),'nearest');
%      d3_zoo_s(isnan(d3_zoo_s)) = griddata(ref_lon_m(~isnan(d3_zoo_s)).*degtometer,ref_dep_s(~isnan(d3_zoo_s)),d3_zoo_s(~isnan(d3_zoo_s)),ref_lon_m(isnan(d3_zoo_s)),ref_dep_s(isnan(d3_zoo_s)),'nearest');
%      d3_phy_s(isnan(d3_phy_s)) = griddata(ref_lon_m(~isnan(d3_phy_s)).*degtometer,ref_dep_s(~isnan(d3_phy_s)),d3_phy_s(~isnan(d3_phy_s)),ref_lon_m(isnan(d3_phy_s)),ref_dep_s(isnan(d3_phy_s)),'nearest');
%      d3_chl_s(isnan(d3_chl_s)) = griddata(ref_lon_m(~isnan(d3_chl_s)).*degtometer,ref_dep_s(~isnan(d3_chl_s)),d3_chl_s(~isnan(d3_chl_s)),ref_lon_m(isnan(d3_chl_s)),ref_dep_s(isnan(d3_chl_s)),'nearest');
%      d3_do_s(isnan(d3_do_s)) = griddata(ref_lon_m(~isnan(d3_do_s)).*degtometer,ref_dep_s(~isnan(d3_do_s)),d3_do_s(~isnan(d3_do_s)),ref_lon_m(isnan(d3_do_s)),ref_dep_s(isnan(d3_do_s)),'nearest');
     
     % add time dim.
     % south
     d3_temp_s_t(:,:,:,j) = d3_temp_s;
     d3_salt_s_t(:,:,:,j) = d3_salt_s;
     d3_u_s_t(:,:,:,j) = d3_u_s;
     d3_v_s_t(:,:,:,j) = d3_v_s;
     d3_no3_s_t(:,:,:,j) = d3_no3_s;
     d3_nh4_s_t(:,:,:,j) = d3_nh4_s;
     d3_zoo_s_t(:,:,:,j) = d3_zoo_s;
     d3_phy_s_t(:,:,:,j) = d3_phy_s;
     d3_chl_s_t(:,:,:,j) = d3_chl_s;
     d3_do_s_t(:,:,:,j) = d3_do_s;
     
     close all; 
     disp(j)
        if j == 12
           savefile = ['interp_bnd_fix_11to17_', num2str(i),'.mat'];
           save(savefile,'d3_*_s_t','-v7.3');
        end
    end
end

% load('KODC_line_cruse.mat','total','txt1');