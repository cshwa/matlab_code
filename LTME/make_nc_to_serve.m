
close all;clear all; clc;

cd 'D:\장기생태\Dynamic\result\1991'

minyear=1991;
maxyear=2010;

list = dir('mp_p_sewer_det_f_monthly_*');
in_file = 'mp_p_sewer_det_f_monthly_';

NO3_1 = ncread(list(1).name,'NO3');
u_1 = ncread(list(1).name,'u');
v_1 = ncread(list(1).name,'v');
[size11 size12 size13]= size(NO3_1);
[sizeu1 sizeu2 sizeu3]= size(u_1);
[sizev1 sizev2 sizev3]= size(v_1);
clearvars NO3_1 u_1 v_1
%%% concatenate array along specified dimension
% hur_sd = cat(3, hur_sd_1, hur_sd_2, hur_sd_3, hur_sd_4); %combine matrix into one
% 
% t_sd=size(hur_sd,3)

tstep=0;

 for j = minyear:maxyear
    fpath = ['D:\장기생태\Dynamic\result\',num2str(j)];
    cd(fpath)
    if j ~= 2001 && j~= 1996 && j~= 2002
        list = dir('mp_p_sewer_det_f_monthly_*');
    elseif j == 2001 | j== 1996 
        list = dir('mp_p_sewer_det_f_2nd_monthly_*');
    elseif j == 2002
        list = dir('mp_p_sewer_det_monthly_*');
    end

    for k = 1:12
        if k == 1
            clearvars mod_*
            mod_lon=ncread(list(k).name,'lon_rho');
            mod_lat=ncread(list(k).name,'lat_rho');
            mod_lonu=ncread(list(k).name,'lon_u');
            mod_latu=ncread(list(k).name,'lat_u');
            mod_lonv=ncread(list(k).name,'lon_v');
            mod_latv=ncread(list(k).name,'lat_v');
            mod_mask=ncread(list(k).name,'mask_rho');
            mod_masku=ncread(list(k).name,'mask_u');
            mod_maskv=ncread(list(k).name,'mask_v');
            mod_time=ncread(list(k).name,'ocean_time');
        end
    mod_temp=ncread(list(k).name,'temp');
    mod_salt=ncread(list(k).name,'salt');
    mod_no3=ncread(list(k).name,'NO3');
    mod_nh4=ncread(list(k).name,'NH4');
    mod_po4=ncread(list(k).name,'tPO4');
    mod_chl=ncread(list(k).name,'chlorophyll');
    mod_u=ncread(list(k).name,'u');
    mod_v=ncread(list(k).name,'v');
    mod_phy=ncread(list(k).name,'phytoplankton');
    mod_zoo=ncread(list(k).name,'zooplankton');
   
    
    mod_int_temp=mean(mod_temp,3) .* (mod_mask./mod_mask);
    mod_int_salt=mean(mod_salt,3) .* (mod_mask./mod_mask);
    mod_int_no3=mean(mod_no3,3) .* (mod_mask./mod_mask);
    mod_int_nh4=mean(mod_nh4,3) .* (mod_mask./mod_mask);
    mod_int_po4=mean(mod_po4,3) .* (mod_mask./mod_mask);
    mod_int_chl=mean(mod_chl,3) .* (mod_mask./mod_mask);
    mod_int_u=mean(mod_u,3) .* (mod_masku./mod_masku);
    mod_int_v=mean(mod_v,3).* (mod_maskv./mod_maskv);
    mod_int_phy=mean(mod_phy,3) .* (mod_mask./mod_mask);
    mod_int_zoo=mean(mod_zoo,3) .* (mod_mask./mod_mask);
   
    ncid = netcdf.create(strcat('LTME_depth_averaged_',num2str(j),'_',num2str(k,'%02d'),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    lonu_dimid = netcdf.defDim(ncid,'lonu',sizeu1);
    latu_dimid = netcdf.defDim(ncid, 'latu', sizeu2);
    lonv_dimid = netcdf.defDim(ncid,'lonv',sizev1);
    latv_dimid = netcdf.defDim(ncid, 'latv', sizev2);
    time_dimid = netcdf.defDim(ncid, 'time',1);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', 'LTME modeling results');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'monthly mean value');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'Seoul Nation Univ. MEPL');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude of RHO-points');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','latitude of RHO-points');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');
    
    lonuvarid = netcdf.defVar(ncid, 'lonu', 'NC_float', lonu_dimid);
    netcdf.putAtt(ncid,lonuvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonuvarid,'long_name','longitude of U-points');
    netcdf.putAtt(ncid,lonuvarid,'units','degrees_east');

    latuvarid = netcdf.defVar(ncid, 'latu', 'NC_float', latu_dimid);
    netcdf.putAtt(ncid,latuvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latuvarid,'long_name','latitude of U-points');
    netcdf.putAtt(ncid,latuvarid,'units','degrees_north');
    
    lonvvarid = netcdf.defVar(ncid, 'lonv', 'NC_float', lonv_dimid);
    netcdf.putAtt(ncid,lonvvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvvarid,'long_name','longitude of V-points');
    netcdf.putAtt(ncid,lonvvarid,'units','degrees_east');

    latvvarid = netcdf.defVar(ncid, 'latv', 'NC_float', latv_dimid);
    netcdf.putAtt(ncid,latvvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvvarid,'long_name','latitude of V-points');
    netcdf.putAtt(ncid,latvvarid,'units','degrees_north');
    
    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    tempvarid = netcdf.defVar(ncid, 'temp', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,tempvarid,'standard_name','temp');
    netcdf.putAtt(ncid,tempvarid,'long_name','temperature');
    netcdf.putAtt(ncid,tempvarid,'units','Celsius');
    
    saltvarid = netcdf.defVar(ncid, 'salt', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,saltvarid,'standard_name','salt');
    netcdf.putAtt(ncid,saltvarid,'long_name','salinity');
    netcdf.putAtt(ncid,saltvarid,'units','psu');
    
    no3varid = netcdf.defVar(ncid, 'no3', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,no3varid,'standard_name','no3');
    netcdf.putAtt(ncid,no3varid,'long_name','nitrate concentration');
    netcdf.putAtt(ncid,no3varid,'units','millimole_nitrogen meter-3');

    nh4varid = netcdf.defVar(ncid, 'nh4', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,nh4varid,'standard_name','nh4');
    netcdf.putAtt(ncid,nh4varid,'long_name','nitrate concentration');
    netcdf.putAtt(ncid,nh4varid,'units','millimole_nitrogen meter-3');

    po4varid = netcdf.defVar(ncid, 'po4', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,po4varid,'standard_name','po4');
    netcdf.putAtt(ncid,po4varid,'long_name','phosphate concentration');
    netcdf.putAtt(ncid,po4varid,'units','millimole_P meter-3');

    chlvarid = netcdf.defVar(ncid, 'chl', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,chlvarid,'standard_name','chl');
    netcdf.putAtt(ncid,chlvarid,'long_name','chlorophyll concentration');
    netcdf.putAtt(ncid,chlvarid,'units','milligrams_chlorophyll meter-3');
    
    phyvarid = netcdf.defVar(ncid, 'phy', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,phyvarid,'standard_name','phy');
    netcdf.putAtt(ncid,phyvarid,'long_name','phytoplankton concentration');
    netcdf.putAtt(ncid,phyvarid,'units','millimole_nitrogen meter-3');
    
    zoovarid = netcdf.defVar(ncid, 'zoo', 'NC_float', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,zoovarid,'standard_name','zoo');
    netcdf.putAtt(ncid,zoovarid,'long_name','zooplankton concentration');
    netcdf.putAtt(ncid,zoovarid,'units','millimole_nitrogen meter-3');
    
    uvarid = netcdf.defVar(ncid, 'u', 'NC_float', [lonu_dimid latu_dimid time_dimid]);
    netcdf.putAtt(ncid,uvarid,'standard_name','u');
    netcdf.putAtt(ncid,uvarid,'long_name','u-velocity');
    netcdf.putAtt(ncid,uvarid,'units','meter per sec');
    
    vvarid = netcdf.defVar(ncid, 'v', 'NC_float', [lonv_dimid latv_dimid time_dimid]);
    netcdf.putAtt(ncid,vvarid,'standard_name','v');
    netcdf.putAtt(ncid,vvarid,'long_name','v-velocity');
    netcdf.putAtt(ncid,vvarid,'units','meter per sec');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(mod_lon(:,1))), squeeze(mod_lon(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(mod_lat(1,:))), squeeze(mod_lat(1,:)));
    netcdf.putVar(ncid, lonuvarid, 0, length(squeeze(mod_lonu(:,1))), squeeze(mod_lonu(:,1)));
    netcdf.putVar(ncid, latuvarid, 0, length(squeeze(mod_latu(1,:))), squeeze(mod_latu(1,:)));
    netcdf.putVar(ncid, lonvvarid, 0, length(squeeze(mod_lonv(:,1))), squeeze(mod_lonv(:,1)));
    netcdf.putVar(ncid, latvvarid, 0, length(squeeze(mod_latv(1,:))), squeeze(mod_latv(1,:)));
    netcdf.putVar(ncid, timevarid, 0, length(mod_time), mod_time);
    netcdf.putVar(ncid, tempvarid, [0 0 0], [size11 size12 1], mod_int_temp);
    netcdf.putVar(ncid, saltvarid, [0 0 0], [size11 size12 1], mod_int_salt);
    netcdf.putVar(ncid, no3varid, [0 0 0], [size11 size12 1], mod_int_no3);
    netcdf.putVar(ncid, nh4varid, [0 0 0], [size11 size12 1], mod_int_nh4);
    netcdf.putVar(ncid, po4varid, [0 0 0], [size11 size12 1], mod_int_po4);
    netcdf.putVar(ncid, chlvarid, [0 0 0], [size11 size12 1], mod_int_chl);
    netcdf.putVar(ncid, phyvarid, [0 0 0], [size11 size12 1], mod_int_phy);
    netcdf.putVar(ncid, zoovarid, [0 0 0], [size11 size12 1], mod_int_zoo);
    netcdf.putVar(ncid, uvarid, [0 0 0], [sizeu1 sizeu2 1], mod_int_u);
    netcdf.putVar(ncid, vvarid, [0 0 0], [sizev1 sizev2 1], mod_int_v);
    netcdf.close(ncid);
 end
 end
 
 return
 

close all;clear all; clc;

cd 'D:\장기생태\Dynamic\result\1991'

minyear=1991;
maxyear=2010;

list = dir('mp_p_sewer_det_f_monthly_*');
in_file = 'mp_p_sewer_det_f_monthly_';

NO3_1 = ncread(list(1).name,'NO3');
u_1 = ncread(list(1).name,'u');
v_1 = ncread(list(1).name,'v');
ref_srho=ncread(list(1).name,'s_rho');
ref_hc=ncread(list(1).name,'hc');
ref_theta_s=ncread(list(1).name,'theta_s');
ref_theta_b=ncread(list(1).name,'theta_b');
% mask=ncread(list(1).name,'mask_rho');
ref_h=ncread(list(1).name,'h');
ref_dep = zlevs(ref_h, zeros(size(ref_h,1),size(ref_h,2)),ref_theta_s,ref_theta_b,ref_hc,length(ref_srho),'r');
ref_dep = permute(ref_dep,[2 3 1]);
[size11 size12 size13]= size(NO3_1);
[sizeu1 sizeu2 sizeu3]= size(u_1);
[sizev1 sizev2 sizev3]= size(v_1);
clearvars NO3_1 u_1 v_1

std_dep_pre=[-1 -10 -20 -30 -35];
std_dep = zeros(size11,size12,length(std_dep_pre));

for i = 1:length(std_dep_pre)
std_dep(:,:,i) = std_dep_pre(i);
end

tstep=0;

 for j = minyear:maxyear
%  for j = minyear:minyear
    fpath = ['D:\장기생태\Dynamic\result\',num2str(j)];
    cd(fpath)
    if j ~= 2001 && j~= 1996 && j~= 2002
        list = dir('mp_p_sewer_det_f_monthly_*');
    elseif j == 2001 | j== 1996 
        list = dir('mp_p_sewer_det_f_2nd_monthly_*');
    elseif j == 2002
        list = dir('mp_p_sewer_det_monthly_*');
    end

    for k = 1:12
        if k == 1
            clearvars mod_*
            mod_lon=ncread(list(k).name,'lon_rho');
            mod_lat=ncread(list(k).name,'lat_rho');
            mod_lonu=ncread(list(k).name,'lon_u');
            mod_latu=ncread(list(k).name,'lat_u');
            mod_lonv=ncread(list(k).name,'lon_v');
            mod_latv=ncread(list(k).name,'lat_v');
            mod_mask=ncread(list(k).name,'mask_rho');
            mod_masku=ncread(list(k).name,'mask_u');
            mod_maskv=ncread(list(k).name,'mask_v');
            mod_time=ncread(list(k).name,'ocean_time');
        end
    mod_temp=ncread(list(k).name,'temp');
    mod_salt=ncread(list(k).name,'salt');
    mod_no3=ncread(list(k).name,'NO3');
    mod_nh4=ncread(list(k).name,'NH4');
    mod_po4=ncread(list(k).name,'tPO4');
    mod_chl=ncread(list(k).name,'chlorophyll');
    mod_u=ncread(list(k).name,'u');
    mod_v=ncread(list(k).name,'v');
    mod_phy=ncread(list(k).name,'phytoplankton');
    mod_zoo=ncread(list(k).name,'zooplankton');

for isrho = 1:size13    
    mod_temp(:,:,isrho)=squeeze(mod_temp(:,:,isrho)) .* (mod_mask./mod_mask);
    mod_salt(:,:,isrho)=squeeze(mod_salt(:,:,isrho)) .* (mod_mask./mod_mask);
    mod_no3(:,:,isrho)=squeeze(mod_no3(:,:,isrho)) .* (mod_mask./mod_mask);
    mod_nh4(:,:,isrho)=squeeze(mod_nh4(:,:,isrho)) .* (mod_mask./mod_mask);
    mod_po4(:,:,isrho)=squeeze(mod_po4(:,:,isrho)) .* (mod_mask./mod_mask);
    mod_chl(:,:,isrho)=squeeze(mod_chl(:,:,isrho)) .* (mod_mask./mod_mask);
    mod_u(:,:,isrho)=squeeze(mod_u(:,:,isrho)) .* (mod_masku./mod_masku);
    mod_v(:,:,isrho)=squeeze(mod_v(:,:,isrho)).* (mod_maskv./mod_maskv);
    mod_phy(:,:,isrho)=squeeze(mod_phy(:,:,isrho)) .* (mod_mask./mod_mask);
    mod_zoo(:,:,isrho)=squeeze(mod_zoo(:,:,isrho)) .* (mod_mask./mod_mask);
end

% u, v to rho points
for isho2 = 1:size13
 temp_u=griddata(mod_lonu,mod_latu,mod_u(:,:,isho2),mod_lon,mod_lat,'nearest');
 temp_v=griddata(mod_lonv,mod_latv,mod_v(:,:,isho2),mod_lon,mod_lat,'nearest');
 mod_u_rho(:,:,isho2) = temp_u .* (mod_mask./mod_mask);
 mod_v_rho(:,:,isho2) = temp_v .* (mod_mask./mod_mask);
end
    
    mod_int_temp=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_temp,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_salt=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_salt,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_no3=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_no3,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_nh4=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_nh4,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_po4=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_po4,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_chl=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_chl,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_phy=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_phy,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_zoo=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_zoo,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_u=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_u_rho,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');
    mod_int_v=griddata(repmat(mod_lon,1,1,size13),repmat(mod_lat,1,1,size13),ref_dep,mod_v_rho,repmat(mod_lon,1,1,length(std_dep_pre)),repmat(mod_lat,1,1,length(std_dep_pre)),std_dep,'nearest');

for iisrho = 1:length(std_dep_pre) 
    if iisrho == 1
        mod_int_temp(:,:,iisrho)=squeeze(mod_temp(:,:,end)) .* (mod_mask./mod_mask);
        mod_int_salt(:,:,iisrho)=squeeze(mod_salt(:,:,end)) .* (mod_mask./mod_mask);
        mod_int_no3(:,:,iisrho)=squeeze(mod_no3(:,:,end)) .* (mod_mask./mod_mask);
        mod_int_nh4(:,:,iisrho)=squeeze(mod_nh4(:,:,end)) .* (mod_mask./mod_mask);
        mod_int_po4(:,:,iisrho)=squeeze(mod_po4(:,:,end)) .* (mod_mask./mod_mask);
        mod_int_chl(:,:,iisrho)=squeeze(mod_chl(:,:,end)) .* (mod_mask./mod_mask);
        mod_int_u(:,:,iisrho)=squeeze(mod_u_rho(:,:,end)) .* (mod_mask./mod_mask);
        mod_int_v(:,:,iisrho)=squeeze(mod_v_rho(:,:,end)).* (mod_mask./mod_mask);
        mod_int_phy(:,:,iisrho)=squeeze(mod_phy(:,:,end)) .* (mod_mask./mod_mask);
        mod_int_zoo(:,:,iisrho)=squeeze(mod_zoo(:,:,end)) .* (mod_mask./mod_mask);    
    else
    mod_int_temp(:,:,iisrho)=squeeze(mod_int_temp(:,:,iisrho)) .* (mod_mask./mod_mask);
    mod_int_salt(:,:,iisrho)=squeeze(mod_int_salt(:,:,iisrho)) .* (mod_mask./mod_mask);
    mod_int_no3(:,:,iisrho)=squeeze(mod_int_no3(:,:,iisrho)) .* (mod_mask./mod_mask);
    mod_int_nh4(:,:,iisrho)=squeeze(mod_int_nh4(:,:,iisrho)) .* (mod_mask./mod_mask);
    mod_int_po4(:,:,iisrho)=squeeze(mod_int_po4(:,:,iisrho)) .* (mod_mask./mod_mask);
    mod_int_chl(:,:,iisrho)=squeeze(mod_int_chl(:,:,iisrho)) .* (mod_mask./mod_mask);
    mod_int_u(:,:,iisrho)=squeeze(mod_int_u(:,:,iisrho)) .* (mod_mask./mod_mask);
    mod_int_v(:,:,iisrho)=squeeze(mod_int_v(:,:,iisrho)).* (mod_mask./mod_mask);
    mod_int_phy(:,:,iisrho)=squeeze(mod_int_phy(:,:,iisrho)) .* (mod_mask./mod_mask);
    mod_int_zoo(:,:,iisrho)=squeeze(mod_int_zoo(:,:,iisrho)) .* (mod_mask./mod_mask);
    end 
end
    

    ncid = netcdf.create(strcat('LTME_10m_depth_interval_',num2str(j),'_',num2str(k,'%02d'),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);
    depth_dimid = netcdf.defDim(ncid, 'depth', length(std_dep_pre));
    time_dimid = netcdf.defDim(ncid, 'time',1);  % does not exist leap year 


    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', 'LTME modeling results');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'monthly mean value');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'Seoul Nation Univ. MEPL');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude of RHO-points');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','latitude of RHO-points');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');
    
    depthvarid = netcdf.defVar(ncid, 'depth', 'NC_float', depth_dimid);
    netcdf.putAtt(ncid,depthvarid,'standard_name','depth');
    netcdf.putAtt(ncid,depthvarid,'long_name','ocean depth');
    netcdf.putAtt(ncid,depthvarid,'units','meter');
    
    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days from 1976-01-01');

    tempvarid = netcdf.defVar(ncid, 'temp', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,tempvarid,'standard_name','temp');
    netcdf.putAtt(ncid,tempvarid,'long_name','temperature');
    netcdf.putAtt(ncid,tempvarid,'units','Celsius');
    
    saltvarid = netcdf.defVar(ncid, 'salt', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,saltvarid,'standard_name','salt');
    netcdf.putAtt(ncid,saltvarid,'long_name','salinity');
    netcdf.putAtt(ncid,saltvarid,'units','psu');
    
    no3varid = netcdf.defVar(ncid, 'no3', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,no3varid,'standard_name','no3');
    netcdf.putAtt(ncid,no3varid,'long_name','nitrate concentration');
    netcdf.putAtt(ncid,no3varid,'units','millimole_nitrogen meter-3');

    nh4varid = netcdf.defVar(ncid, 'nh4', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,nh4varid,'standard_name','nh4');
    netcdf.putAtt(ncid,nh4varid,'long_name','nitrate concentration');
    netcdf.putAtt(ncid,nh4varid,'units','millimole_nitrogen meter-3');

    po4varid = netcdf.defVar(ncid, 'po4', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,po4varid,'standard_name','po4');
    netcdf.putAtt(ncid,po4varid,'long_name','phosphate concentration');
    netcdf.putAtt(ncid,po4varid,'units','millimole_P meter-3');

    chlvarid = netcdf.defVar(ncid, 'chl', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,chlvarid,'standard_name','chl');
    netcdf.putAtt(ncid,chlvarid,'long_name','chlorophyll concentration');
    netcdf.putAtt(ncid,chlvarid,'units','milligrams_chlorophyll meter-3');
    
    phyvarid = netcdf.defVar(ncid, 'phy', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,phyvarid,'standard_name','phy');
    netcdf.putAtt(ncid,phyvarid,'long_name','phytoplankton concentration');
    netcdf.putAtt(ncid,phyvarid,'units','millimole_nitrogen meter-3');
    
    zoovarid = netcdf.defVar(ncid, 'zoo', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,zoovarid,'standard_name','zoo');
    netcdf.putAtt(ncid,zoovarid,'long_name','zooplankton concentration');
    netcdf.putAtt(ncid,zoovarid,'units','millimole_nitrogen meter-3');
    
    uvarid = netcdf.defVar(ncid, 'u', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,uvarid,'standard_name','u');
    netcdf.putAtt(ncid,uvarid,'long_name','u-velocity');
    netcdf.putAtt(ncid,uvarid,'units','meter per sec');
    
    vvarid = netcdf.defVar(ncid, 'v', 'NC_float', [lon_dimid lat_dimid depth_dimid time_dimid]);
    netcdf.putAtt(ncid,vvarid,'standard_name','v');
    netcdf.putAtt(ncid,vvarid,'long_name','v-velocity');
    netcdf.putAtt(ncid,vvarid,'units','meter per sec');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(mod_lon(:,1))), squeeze(mod_lon(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(mod_lat(1,:))), squeeze(mod_lat(1,:)));
    netcdf.putVar(ncid, depthvarid, 0, length(std_dep_pre),std_dep_pre);
    netcdf.putVar(ncid, timevarid, 0, length(mod_time), mod_time);
    netcdf.putVar(ncid, tempvarid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_temp);
    netcdf.putVar(ncid, saltvarid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_salt);
    netcdf.putVar(ncid, no3varid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_no3);
    netcdf.putVar(ncid, nh4varid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_nh4);
    netcdf.putVar(ncid, po4varid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_po4);
    netcdf.putVar(ncid, chlvarid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_chl);
    netcdf.putVar(ncid, phyvarid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_phy);
    netcdf.putVar(ncid, zoovarid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_zoo);
    netcdf.putVar(ncid, uvarid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_u);
    netcdf.putVar(ncid, vvarid, [0 0 0 0], [size11 size12 length(std_dep_pre) 1], mod_int_v);
    netcdf.close(ncid);
 end
 end
    return
 
    
    
    %% make depth file 
    close all;clear all; clc;

cd 'D:\장기생태\Dynamic\result\1991'

list = dir('mp_p_sewer_det_f_monthly_*');
in_file = 'mp_p_sewer_det_f_monthly_';

NO3_1 = ncread(list(1).name,'NO3');
mod_lon = ncread(list(1).name,'lon_rho');
mod_lat = ncread(list(1).name,'lat_rho');
u_1 = ncread(list(1).name,'u');
v_1 = ncread(list(1).name,'v');
ref_srho=ncread(list(1).name,'s_rho');
ref_hc=ncread(list(1).name,'hc');
ref_theta_s=ncread(list(1).name,'theta_s');
ref_theta_b=ncread(list(1).name,'theta_b');
mask=ncread(list(1).name,'mask_rho');
ref_h=ncread(list(1).name,'h');
ref_dep = zlevs(ref_h, zeros(size(ref_h,1),size(ref_h,2)),ref_theta_s,ref_theta_b,ref_hc,length(ref_srho),'r');
ref_dep = permute(ref_dep,[2 3 1]);
[size11 size12 size13]= size(NO3_1);
[sizeu1 sizeu2 sizeu3]= size(u_1);
[sizev1 sizev2 sizev3]= size(v_1);
clearvars NO3_1 u_1 v_1

    end_dep=ref_dep(:,:,1).* (mask./mask);

    ncid = netcdf.create(strcat('LTME_model_depth_file.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',size11);
    lat_dimid = netcdf.defDim(ncid, 'lat', size12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', 'LTME modeling results');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'monthly mean value');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'Seoul Nation Univ. MEPL');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude of RHO-points');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','latitude of RHO-points');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');
    
    depthvarid = netcdf.defVar(ncid, 'depth', 'NC_float', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,depthvarid,'standard_name','depth');
    netcdf.putAtt(ncid,depthvarid,'long_name','ocean depth');
    netcdf.putAtt(ncid,depthvarid,'units','meter');

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, length(squeeze(mod_lon(:,1))), squeeze(mod_lon(:,1)));
    netcdf.putVar(ncid, latvarid, 0, length(squeeze(mod_lat(1,:))), squeeze(mod_lat(1,:)));
    netcdf.putVar(ncid, depthvarid, [0 0], [size11 size12],end_dep);
    netcdf.close(ncid);
 
