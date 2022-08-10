close all; clear; clc; 

% for i = 2001:2010
for i = 2006:2006
    for j = 1:12
        clearvars -except i j d3_*_s_t d3_*_e_t
        path = ['D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary\',num2str(i)]
        cd(path)
    %load raw file
    rawfile = ['stdep1_avg_mon_', num2str(i),'_', num2str(j,'%02d'),'.nc'];
    raw_temp=ncread(rawfile,'temp');
    raw_salt=ncread(rawfile,'salt');
    raw_u=ncread(rawfile,'u');
    raw_v=ncread(rawfile,'v');
    raw_chl=ncread(rawfile,'chlorophyll');
    raw_phy=ncread(rawfile,'phytoplankton');
    raw_zoo=ncread(rawfile,'zooplankton');
    raw_no3=ncread(rawfile,'NO3');
    raw_nh4=ncread(rawfile,'NH4');
    raw_dep=ncread(rawfile,'depth');
    raw_lon=ncread(rawfile,'lon');
    raw_lat=ncread(rawfile,'lat');
    [raw_lat_m raw_lon_m]=meshgrid(raw_lat,raw_lon);

    %load ref bnd file
    cd D:\장기생태\Dynamic\02_grid_depth
    bnd_reffile = 'grid_gy_v11_s.nc';
    ref_lon=ncread(bnd_reffile,'lon_rho');
    ref_lat=ncread(bnd_reffile,'lat_rho');
    ref_lonu=ncread(bnd_reffile,'lon_u');
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
%     raw_no3(raw_no3>(10^20))=NaN;
%     raw_nh4(raw_nh4>(10^20))=NaN;
%     raw_zoo(raw_zoo>(10^20))=NaN;
%     raw_phy(raw_phy>(10^20))=NaN;
%     raw_chl(raw_chl>(10^20))=NaN;
    for k = 1:length(raw_dep)
        %south
        d2_temp_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_temp(:,:,k)),ref_lon_s,ref_lat_s);
        d2_salt_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_salt(:,:,k)),ref_lon_s,ref_lat_s);
        d2_u_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_u(:,:,k)),ref_lon_s,ref_lat_s);
        d2_v_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_v(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_no3_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_no3(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_nh4_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_nh4(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_zoo_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_zoo(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_phy_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_phy(:,:,k)),ref_lon_s,ref_lat_s);
%         d2_chl_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_chl(:,:,k)),ref_lon_s,ref_lat_s);
        %east
        d2_temp_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_temp(:,:,k)),ref_lon_e,ref_lat_e);
        d2_salt_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_salt(:,:,k)),ref_lon_e,ref_lat_e);
        d2_u_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_u(:,:,k)),ref_lon_e,ref_lat_e);
        d2_v_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_v(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_no3_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_no3(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_nh4_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_nh4(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_zoo_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_zoo(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_phy_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_phy(:,:,k)),ref_lon_e,ref_lat_e);
%         d2_chl_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_chl(:,:,k)),ref_lon_e,ref_lat_e);
    end
    
    %interp 3d

     % south
     ref_lon_sm=repmat(ref_lon_s,1,20); ref_lat_sm=repmat(ref_lat_s,1,20);
     raw_lon_sm=repmat(ref_lon_s,1,length(raw_dep)); raw_dep_sm = repmat(raw_dep',length(ref_lon_s),1);
     
     d3_temp_s = griddata(raw_lon_sm,raw_dep_sm,d2_temp_s,ref_lon_sm,ref_dep_s,'nearest');
     d3_salt_s = griddata(raw_lon_sm,raw_dep_sm,d2_salt_s,ref_lon_sm,ref_dep_s,'nearest');
     d3_u_s = griddata(raw_lon_sm,raw_dep_sm,d2_u_s,ref_lon_sm,ref_dep_s,'nearest');
     d3_v_s = griddata(raw_lon_sm,raw_dep_sm,d2_v_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_no3_s = griddata(raw_lon_sm,raw_dep_sm,d2_no3_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_nh4_s = griddata(raw_lon_sm,raw_dep_sm,d2_nh4_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_zoo_s = griddata(raw_lon_sm,raw_dep_sm,d2_zoo_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_phy_s = griddata(raw_lon_sm,raw_dep_sm,d2_phy_s,ref_lon_sm,ref_dep_s,'nearest');
%      d3_chl_s = griddata(raw_lon_sm,raw_dep_sm,d2_chl_s,ref_lon_sm,ref_dep_s,'nearest');
     
     % east
     ref_lon_em=repmat(ref_lon_e',1,20); ref_lat_em=repmat(ref_lat_e',1,20);
     raw_lat_em=repmat(ref_lat_e',1,length(raw_dep)); raw_dep_em = repmat(raw_dep',length(ref_lon_e),1);
     
     d3_temp_e = griddata(raw_lat_em,raw_dep_em,d2_temp_e,ref_lat_em,ref_dep_e,'nearest');
     d3_salt_e = griddata(raw_lat_em,raw_dep_em,d2_salt_e,ref_lat_em,ref_dep_e,'nearest');
     d3_u_e = griddata(raw_lat_em,raw_dep_em,d2_u_e,ref_lat_em,ref_dep_e,'nearest');
     d3_v_e = griddata(raw_lat_em,raw_dep_em,d2_v_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_no3_e = griddata(raw_lat_em,raw_dep_em,d2_no3_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_nh4_e = griddata(raw_lat_em,raw_dep_em,d2_nh4_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_zoo_e = griddata(raw_lat_em,raw_dep_em,d2_zoo_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_phy_e = griddata(raw_lat_em,raw_dep_em,d2_phy_e,ref_lat_em,ref_dep_e,'nearest');
%      d3_chl_e = griddata(raw_lat_em,raw_dep_em,d2_chl_e,ref_lat_em,ref_dep_e,'nearest');
    
     
     %fill missing value

     % south
    
     d3_temp_s(isnan(d3_temp_s)) = griddata(ref_lon_sm(~isnan(d3_temp_s)),ref_dep_s(~isnan(d3_temp_s)),d3_temp_s(~isnan(d3_temp_s)),ref_lon_sm(isnan(d3_temp_s)),ref_dep_s(isnan(d3_temp_s)),'nearest');
     d3_salt_s(isnan(d3_salt_s)) = griddata(ref_lon_sm(~isnan(d3_salt_s)),ref_dep_s(~isnan(d3_salt_s)),d3_salt_s(~isnan(d3_salt_s)),ref_lon_sm(isnan(d3_salt_s)),ref_dep_s(isnan(d3_salt_s)),'nearest');
     d3_u_s(isnan(d3_u_s)) = griddata(ref_lon_sm(~isnan(d3_u_s)),ref_dep_s(~isnan(d3_u_s)),d3_u_s(~isnan(d3_u_s)),ref_lon_sm(isnan(d3_u_s)),ref_dep_s(isnan(d3_u_s)),'nearest');
     d3_v_s(isnan(d3_v_s)) = griddata(ref_lon_sm(~isnan(d3_v_s)),ref_dep_s(~isnan(d3_v_s)),d3_v_s(~isnan(d3_v_s)),ref_lon_sm(isnan(d3_v_s)),ref_dep_s(isnan(d3_v_s)),'nearest');
%      d3_no3_s(isnan(d3_no3_s)) = griddata(ref_lon_sm(~isnan(d3_no3_s)),ref_dep_s(~isnan(d3_no3_s)),d3_no3_s(~isnan(d3_no3_s)),ref_lon_sm(isnan(d3_no3_s)),ref_dep_s(isnan(d3_no3_s)),'nearest');
%      d3_nh4_s(isnan(d3_nh4_s)) = griddata(ref_lon_sm(~isnan(d3_nh4_s)),ref_dep_s(~isnan(d3_nh4_s)),d3_nh4_s(~isnan(d3_nh4_s)),ref_lon_sm(isnan(d3_nh4_s)),ref_dep_s(isnan(d3_nh4_s)),'nearest');
%      d3_zoo_s(isnan(d3_zoo_s)) = griddata(ref_lon_sm(~isnan(d3_zoo_s)),ref_dep_s(~isnan(d3_zoo_s)),d3_zoo_s(~isnan(d3_zoo_s)),ref_lon_sm(isnan(d3_zoo_s)),ref_dep_s(isnan(d3_zoo_s)),'nearest');
%      d3_phy_s(isnan(d3_phy_s)) = griddata(ref_lon_sm(~isnan(d3_phy_s)),ref_dep_s(~isnan(d3_phy_s)),d3_phy_s(~isnan(d3_phy_s)),ref_lon_sm(isnan(d3_phy_s)),ref_dep_s(isnan(d3_phy_s)),'nearest');
%      d3_chl_s(isnan(d3_chl_s)) = griddata(ref_lon_sm(~isnan(d3_chl_s)),ref_dep_s(~isnan(d3_chl_s)),d3_chl_s(~isnan(d3_chl_s)),ref_lon_sm(isnan(d3_chl_s)),ref_dep_s(isnan(d3_chl_s)),'nearest');
     
     % east
     
     d3_temp_e(isnan(d3_temp_e)) = griddata(ref_lat_em(~isnan(d3_temp_e)),ref_dep_e(~isnan(d3_temp_e)),d3_temp_e(~isnan(d3_temp_e)),ref_lat_em(isnan(d3_temp_e)),ref_dep_e(isnan(d3_temp_e)),'nearest');
     d3_salt_e(isnan(d3_salt_e)) = griddata(ref_lat_em(~isnan(d3_salt_e)),ref_dep_e(~isnan(d3_salt_e)),d3_salt_e(~isnan(d3_salt_e)),ref_lat_em(isnan(d3_salt_e)),ref_dep_e(isnan(d3_salt_e)),'nearest');
     d3_u_e(isnan(d3_u_e)) = griddata(ref_lat_em(~isnan(d3_u_e)),ref_dep_e(~isnan(d3_u_e)),d3_u_e(~isnan(d3_u_e)),ref_lat_em(isnan(d3_u_e)),ref_dep_e(isnan(d3_u_e)),'nearest');
     d3_v_e(isnan(d3_v_e)) = griddata(ref_lat_em(~isnan(d3_v_e)),ref_dep_e(~isnan(d3_v_e)),d3_v_e(~isnan(d3_v_e)),ref_lat_em(isnan(d3_v_e)),ref_dep_e(isnan(d3_v_e)),'nearest');
%      d3_no3_e(isnan(d3_no3_e)) = griddata(ref_lat_em(~isnan(d3_no3_e)),ref_dep_e(~isnan(d3_no3_e)),d3_no3_e(~isnan(d3_no3_e)),ref_lat_em(isnan(d3_no3_e)),ref_dep_e(isnan(d3_no3_e)),'nearest');
%      d3_nh4_e(isnan(d3_nh4_e)) = griddata(ref_lat_em(~isnan(d3_nh4_e)),ref_dep_e(~isnan(d3_nh4_e)),d3_nh4_e(~isnan(d3_nh4_e)),ref_lat_em(isnan(d3_nh4_e)),ref_dep_e(isnan(d3_nh4_e)),'nearest');
%      d3_zoo_e(isnan(d3_zoo_e)) = griddata(ref_lat_em(~isnan(d3_zoo_e)),ref_dep_e(~isnan(d3_zoo_e)),d3_zoo_e(~isnan(d3_zoo_e)),ref_lat_em(isnan(d3_zoo_e)),ref_dep_e(isnan(d3_zoo_e)),'nearest');
%      d3_phy_e(isnan(d3_phy_e)) = griddata(ref_lat_em(~isnan(d3_phy_e)),ref_dep_e(~isnan(d3_phy_e)),d3_phy_e(~isnan(d3_phy_e)),ref_lat_em(isnan(d3_phy_e)),ref_dep_e(isnan(d3_phy_e)),'nearest');
%      d3_chl_e(isnan(d3_chl_e)) = griddata(ref_lat_em(~isnan(d3_chl_e)),ref_dep_e(~isnan(d3_chl_e)),d3_chl_e(~isnan(d3_chl_e)),ref_lat_em(isnan(d3_chl_e)),ref_dep_e(isnan(d3_chl_e)),'nearest');
     
     % add time dim.
     % south
     d3_temp_s_t(:,:,j) = d3_temp_s;
     d3_salt_s_t(:,:,j) = d3_salt_s;
     d3_u_s_t(:,:,j) = d3_u_s;
     d3_v_s_t(:,:,j) = d3_v_s;
%      d3_no3_s_t(:,:,j) = d3_no3_s;
%      d3_nh4_s_t(:,:,j) = d3_nh4_s;
%      d3_zoo_s_t(:,:,j) = d3_zoo_s;
%      d3_phy_s_t(:,:,j) = d3_phy_s;
%      d3_chl_s_t(:,:,j) = d3_chl_s;
     %east
     d3_temp_e_t(:,:,j) = d3_temp_e;
     d3_salt_e_t(:,:,j) = d3_salt_e;
     d3_u_e_t(:,:,j) = d3_u_e;
     d3_v_e_t(:,:,j) = d3_v_e;
%      d3_no3_e_t(:,:,j) = d3_no3_e;
%      d3_nh4_e_t(:,:,j) = d3_nh4_e;
%      d3_zoo_e_t(:,:,j) = d3_zoo_e;
%      d3_phy_e_t(:,:,j) = d3_phy_e;
%      d3_chl_e_t(:,:,j) = d3_chl_e;
     close all;
     disp(j)
        if j == 12
           savefile = ['interp_bnd_fix_', num2str(i),'.mat'];
           save(savefile,'-v7.3');
        end
    end
end

% load('KODC_line_cruse.mat','total','txt1');