close all; clear; clc; 

for i = 2001:2001
% for i = 2001:2001
    for j = 1:12
        clearvars -except i j d3_*_s_t d3_*_e_t
        path = ['D:\장기생태\Dynamic\07_boundary_ts\yjtak_boundary\',num2str(i)]
        cd(path)
    %load raw file
    rawfile = ['spinup1_monthly', num2str(j,'%02d'),'.nc'];
    
    raw_do=ncread(rawfile,'oxygen');
    raw_h=ncread(rawfile,'h');
    raw_lon_m=ncread(rawfile,'lon_rho');
    raw_lat_m=ncread(rawfile,'lat_rho');
    raw_srho=ncread(rawfile,'s_rho');
    raw_hc=ncread(rawfile,'hc');
    raw_theta_s=ncread(rawfile,'theta_s');
    raw_theta_b=ncread(rawfile,'theta_b');
%     [raw_lat_m raw_lon_m]=meshgrid(raw_lat,raw_lon);

    raw_dep = zlevs(raw_h, zeros(size(raw_h,1),size(raw_h,2)),raw_theta_s,raw_theta_b,raw_hc,length(raw_srho),'r');

    %load ref bnd file
    cd D:\장기생태\Dynamic\02_grid_depth
    bnd_reffile = 'grid_gy_v11_s.nc';
    ref_lon=ncread(bnd_reffile,'lon_rho');
    ref_lat=ncread(bnd_reffile,'lat_rho');
    ref_lonu=ncread(bnd_reffile,'lon_u');
    ref_latv=ncread(bnd_reffile,'lat_v');
    
    
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
    raw_do(raw_do>(10^20))=NaN;
        
    for k = 1:length(raw_srho)
        %south
        d2_do_s(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_do(:,:,k)),ref_lon_s,ref_lat_s);
        
        %east
        d2_do_e(:,k)= griddata(raw_lon_m,raw_lat_m,squeeze(raw_do(:,:,k)),ref_lon_e,ref_lat_e);
        
    end
      
     %fill missing value
     % south
     d2_do_s(isnan(d2_do_s))= griddata(ref_lon_s(~isnan(d2_do_s)),ref_lat_s(~isnan(d2_do_s)),d2_do_s(~isnan(d2_do_s)),ref_lon_s(isnan(d2_do_s)),ref_lat_s(isnan(d2_do_s)));
     
     % east
     d2_do_e(isnan(d2_do_e))= griddata(ref_lon_e(~isnan(d2_do_e)),ref_lat_e(~isnan(d2_do_e)),d2_do_e(~isnan(d2_do_e)),ref_lon_e(isnan(d2_do_e)),ref_lat_e(isnan(d2_do_e)));
     
     d3_do_s=d2_do_s(:,1:2:end);
     d3_do_e=d2_do_e(:,1:2:end);
     
     % add time dim.
     % south
     d3_do_s_t(:,:,j) = d3_do_s;

     %east
     d3_do_e_t(:,:,j) = d3_do_e;
     disp(j)
        if j == 12
           savefile = ['interp_bnd_do_', num2str(i),'.mat'];
           save(savefile,'-v7.3');
        end
    end
end

% load('KODC_line_cruse.mat','total','txt1');