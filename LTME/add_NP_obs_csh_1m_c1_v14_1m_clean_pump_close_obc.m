clear all;close all;clc;
load('KODC_line_cruse.mat','total','txt1');
load KOEM_clim.mat
linear_inx=find(total(:,1) > 15);total(linear_inx,:)=[]; %OBC only has 10m dep.

pre_temp_sur=NaN(1,12); pre_temp_bot=NaN(1,12); pre_salt_sur=NaN(1,12); pre_salt_bot=NaN(1,12); 
pre_oxy_sur=NaN(1,12); pre_oxy_bot=NaN(1,12); pre_NO3_sur=NaN(1,12);

pre_temp_sur(2:2:12)= total(1:2:12,2); pre_temp_bot(2:2:12)= total(2:2:12,2);
pre_salt_sur(2:2:12)= total(1:2:12,4); pre_salt_bot(2:2:12)= total(2:2:12,4);
pre_oxy_sur(2:2:12)= total(1:2:12,6); pre_oxy_bot(2:2:12)= total(2:2:12,6); %ml/L
pre_NO3_sur(2:2:12)= total(1:2:12,11); % ug-at/L
t=1:12;
pre_temp_sur(isnan(pre_temp_sur))=interp1(t(~isnan(pre_temp_sur)), pre_temp_sur(~isnan(pre_temp_sur)), t(isnan(pre_temp_sur))); pre_temp_sur(1) = (pre_temp_sur(2)+pre_temp_sur(end))/2;
pre_temp_bot(isnan(pre_temp_bot))=interp1(t(~isnan(pre_temp_bot)), pre_temp_bot(~isnan(pre_temp_bot)), t(isnan(pre_temp_bot))); pre_temp_bot(1) = (pre_temp_bot(2)+pre_temp_bot(end))/2;
pre_salt_sur(isnan(pre_salt_sur))=interp1(t(~isnan(pre_salt_sur)), pre_salt_sur(~isnan(pre_salt_sur)), t(isnan(pre_salt_sur))); pre_salt_sur(1) = (pre_salt_sur(2)+pre_salt_sur(end))/2;
pre_salt_bot(isnan(pre_salt_bot))=interp1(t(~isnan(pre_salt_bot)), pre_salt_bot(~isnan(pre_salt_bot)), t(isnan(pre_salt_bot))); pre_salt_bot(1) = (pre_salt_bot(2)+pre_salt_bot(end))/2;
pre_oxy_sur(isnan(pre_oxy_sur))=interp1(t(~isnan(pre_oxy_sur)), pre_oxy_sur(~isnan(pre_oxy_sur)), t(isnan(pre_oxy_sur))); pre_oxy_sur(1) = (pre_oxy_sur(2)+pre_oxy_sur(end))/2;
pre_oxy_bot(isnan(pre_oxy_bot))=interp1(t(~isnan(pre_oxy_bot)), pre_oxy_bot(~isnan(pre_oxy_bot)), t(isnan(pre_oxy_bot))); pre_oxy_bot(1) = (pre_oxy_bot(2)+pre_oxy_bot(end))/2;
pre_NO3_sur(isnan(pre_NO3_sur))=interp1(t(~isnan(pre_NO3_sur)), pre_NO3_sur(~isnan(pre_NO3_sur)), t(isnan(pre_NO3_sur))); pre_NO3_sur(1) = (pre_NO3_sur(2)+pre_NO3_sur(end))/2;
pre_oxy_sur=pre_oxy_sur*0.7*44.661;
pre_oxy_bot=pre_oxy_bot*0.7*44.661; % mg/L -> millimole_oxygen meter-3;  
%¥ìg-at/l = mmol/m3 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% plot(pre_temp_sur,'linew',2); grid on; set(gca,'linew',2);xlabel('month','fontsize',20);
% ylabel('Temp.(^oC)','fontsize',20); xlim([1 12]);set(gca,'fontsize',20);
% 
% plot(pre_salt_sur,'linew',2); grid on; set(gca,'linew',2);xlabel('month','fontsize',20);
% ylabel('salt.(PSU)','fontsize',20); xlim([1 12]);set(gca,'fontsize',20);
% 
% plot(pre_oxy_sur,'linew',2); grid on; set(gca,'linew',2);xlabel('month','fontsize',20); 
% ylabel('DO.(ml/L)','fontsize',20); xlim([1 12]); set(gca,'fontsize',20);
% 
% plot(pre_NO3_sur,'linew',2); grid on; set(gca,'linew',2);xlabel('month','fontsize',20);
% ylabel('NO3.(mmol/m^3)','fontsize',20); xlim([1 12]);set(gca,'fontsize',20);

filename='pohang_closed_obs_1m_c1_v11_runoff_clean_pump_obc.nc';

nc=netcdf(filename,'write');
time=nc{'temp_time'}(:);
cycle=365;

% load new_lateral_BC_mean_koem&OBS.mat
% salt_OBS_bt=load('new_lateral_BC_mean_OBS.mat','pre_salt_bot');
% temp_OBS_bt=load('new_lateral_BC_mean_OBS.mat','pre_temp_bot');

%% raw maximum of observation  
% temp_obs_bt_max1=[NaN NaN 5.17 NaN 7.21 NaN NaN 9.87 NaN 15.49 NaN NaN];
% salt_obs_bt_max1=[NaN NaN 32.67 NaN 32.93 NaN NaN 32.84 NaN 31.81 NaN NaN];

temp_obs_bt_max1=[15.49 NaN NaN NaN NaN 5.17 NaN 7.21 NaN NaN 9.87 NaN 15.49 NaN NaN NaN NaN 5.17];
salt_obs_bt_max1=[31.81 NaN NaN NaN NaN 32.67 NaN 32.93 NaN NaN 32.84 NaN 31.81 NaN NaN NaN NaN 32.67];

t=1:18;
temp_obs_bt_max1(isnan(temp_obs_bt_max1))=interp1(t(~isnan(temp_obs_bt_max1)), temp_obs_bt_max1(~isnan(temp_obs_bt_max1)), t(isnan(temp_obs_bt_max1))); 
salt_obs_bt_max1(isnan(salt_obs_bt_max1))=interp1(t(~isnan(salt_obs_bt_max1)), salt_obs_bt_max1(~isnan(salt_obs_bt_max1)), t(isnan(salt_obs_bt_max1)));

temp_obs_bt_max=temp_obs_bt_max1(4:15); temp_obs_bt_max([1,2]) = [2 4];
salt_obs_bt_max=salt_obs_bt_max1(4:15);

gd=grd('v11_pump_closed_obc');
lon_rho=gd.lon_rho;
N=gd.N;
[y,x]=size(lon_rho);

temp=ones(length(time),N,x);
salt=ones(length(time),N,x);
oxygen=ones(length(time),N,x);
temp_w=ones(length(time),N,y);
salt_w=ones(length(time),N,y);
oxygen_w=ones(length(time),N,y);
no3=ones(length(time),N,x);
nh4=ones(length(time),N,x);
po4=ones(length(time),N,x);
chl=ones(length(time),N,x);

temp_ref=temp_obs_bt_max;
salt_ref=salt_obs_bt_max;

oxygen_ref=pre_oxy_sur;
no3_ref=pre_NO3_sur;

nh4_ref=[0.9 0.9 0.9 0.5 0.3 0.1 0.01 0.01 0.01 0.05 0.1 0.7];
po4_ref=[0.4 0.3 0.2 0.1 0.1 0. 0. 0. 0. 0.1 0.2 0.4];
chl_ref=[1.5 0.9 2 3 2 1.5 1 5 5 5.5 6.0 7.5 ];

% for i=1:1:12
%     temp(i,:,:)=temp(i,:,:)*temp_ref(i);
%     salt(i,:,:)=salt(i,:,:)*salt_ref(i);
%     oxygen(i,:,:)=oxygen(i,:,:)*oxygen_ref(i);
%     no3(i,:,:)=no3(i,:,:)*no3_ref(i);
%     nh4(i,:,:)=nh4(i,:,:)*nh4_ref(i);
%     po4(i,:,:)=po4(i,:,:)*po4_ref(i);
%     chl(i,:,:)=chl(i,:,:)*chl_ref(i);
% end

for i=1:1:12
    temp(i,:,:)=temp(i,:,:)*temp_ref(i);
    salt(i,:,:)=salt(i,:,:)*salt_ref(i);
    oxygen(i,:,:)=oxygen(i,:,:).*0;
    temp_w1(i,:,:)=temp_w(i,:,:)*temp_ref(i);
    salt_w1(i,:,:)=salt_w(i,:,:)*salt_ref(i);
    oxygen_w1(i,:,:)=oxygen_w(i,:,:).*0;
    temp_e1(i,:,:)=temp_w(i,:,:)*temp_ref(i);
    salt_e1(i,:,:)=salt_w(i,:,:)*salt_ref(i);
    oxygen_e1(i,:,:)=oxygen_w(i,:,:).*300;
    no3(i,:,:)=no3(i,:,:).*0;
    nh4(i,:,:)=nh4(i,:,:).*0;
    po4(i,:,:)=po4(i,:,:).*0;
    chl(i,:,:)=chl(i,:,:).*0;
end

load 'vel_pump_OBC_fix.mat'
ob_u = zeros(12,30,92); %west obc eta_rho : 50 (x:1,y:50)
ob_u1 = zeros(12,30,92);

for i = 1:12
    ob_u(i,1:10,35) = vel_bt;
    ob_u(i,11:30,35) = vel_surf;
    ob_u1(i,1:10,57) = vel_bt;
    ob_u1(i,11:30,57) = vel_surf;
end 
%  plot(squeeze(ob_u(50,:,1)),1:30)
ubar_w=zeros(92,12);
ubar_w(35,:)=mean(squeeze(ob_u(:,:,35)),2);
nc{'temp_east'}(:)=temp_e1;
nc{'salt_east'}(:)=salt_e1;
nc{'u_east'}(:)=-1.*ob_u1;
nc{'ubar_east'}(:)=-1.*ubar_w';

nc{'temp_west'}(:)=temp_w1;
nc{'salt_west'}(:)=salt_w1;
nc{'u_west'}(:)=-1.*ob_u1;
nc{'ubar_west'}(:)=-1.*ubar_w';

% nc('NO3_time') = length(time);
% nc{'NO3_time'} = ncdouble('NO3_time') ;
% nc{'NO3_time'}.long_name = ncchar('time for Nitrate');
% nc{'NO3_time'}.long_name = 'time for Nitrate';
% nc{'NO3_time'}.units = ncchar('day');
% nc{'NO3_time'}.units = 'day';
% nc{'NO3_time'}.calendar = ncchar('360.0 days in every year');
% nc{'NO3_time'}.calendar = '360.0 days in every year';
% nc{'NO3_time'}.cycle_length = cycle;
% nc{'NO3_time'}(:)=time;
% 
% nc{'NO3_north'} = ncdouble('NO3_time','s_rho','xi_rho') ;
% nc{'NO3_north'}.long_name = ncchar('northern boundary Nitrate');
% nc{'NO3_north'}.long_name = 'northern boundary Nitrate';
% nc{'NO3_north'}.units = ncchar('mMol N m-3');
% nc{'NO3_north'}.units = 'mMol N m-3';
% nc{'NO3_north'}.coordinates = ncchar('lon_rho s_rho NO3_time');
% nc{'NO3_north'}.coordinates = 'lon_rho s_rho NO3_time';
% nc{'NO3_north'}(:)=no3;
% 
% nc('NH4_time') = length(time);
% nc{'NH4_time'} = ncdouble('NH4_time') ;
% nc{'NH4_time'}.long_name = ncchar('time for Ammonium');
% nc{'NH4_time'}.long_name = 'time for Ammonium';
% nc{'NH4_time'}.units = ncchar('day');
% nc{'NH4_time'}.units = 'day';
% nc{'NH4_time'}.calendar = ncchar('360.0 days in every year');
% nc{'NH4_time'}.calendar = '360.0 days in every year';
% nc{'NH4_time'}.cycle_length = cycle;
% nc{'NH4_time'}(:)=time;
% 
% nc{'NH4_north'} = ncdouble('NH4_time','s_rho','xi_rho') ;
% nc{'NH4_north'}.long_name = ncchar('northern boundary Ammonium');
% nc{'NH4_north'}.long_name = 'northern boundary Ammonium';
% nc{'NH4_north'}.units = ncchar('mMol N m-3');
% nc{'NH4_north'}.units = 'mMol N m-3';
% nc{'NH4_north'}.coordinates = ncchar('lon_rho s_rho NH4_time');
% nc{'NH4_north'}.coordinates = 'lon_rho s_rho NH4_time';
% nc{'NH4_north'}(:)=nh4;

% nc('tPO4_time') = length(time);
% 
% nc{'tPO4_time'} = ncdouble('PO4_time') ;
% nc{'tPO4_time'}.long_name = ncchar('time for Phosphate');
% nc{'tPO4_time'}.long_name = 'time for Phosphate';
% nc{'tPO4_time'}.units = ncchar('day');
% nc{'tPO4_time'}.units = 'day';
% nc{'tPO4_time'}.calendar = ncchar('360.0 days in every year');
% nc{'tPO4_time'}.calendar = '360.0 days in every year';
% nc{'tPO4_time'}.cycle_length = cycle;
% nc{'tPO4_time'}(:)=time;
% 
% nc{'tPO4_north'} = ncdouble('tPO4_time','s_rho','xi_rho') ;
% nc{'tPO4_north'}.long_name = ncchar('northern boundary Phosphate');
% nc{'tPO4_north'}.long_name = 'northern boundary Phosphate';
% nc{'tPO4_north'}.units = ncchar('milimole PO4 m-3');
% nc{'tPO4_north'}.units = 'milimole PO4 m-3';
% nc{'tPO4_north'}.coordinates = ncchar('lon_rho s_rho tPO4_time');
% nc{'tPO4_north'}.coordinates = 'lon_rho s_rho tPO4_time';
% nc{'tPO4_north'}(:)=po4;

nc('oxygen_time') = length(time);
nc{'oxygen_time'} = ncdouble('oxygen_time') ;
nc{'oxygen_time'}.long_name = ncchar('time for oxygen');
nc{'oxygen_time'}.long_name = 'time for oxygen';
nc{'oxygen_time'}.units = ncchar('day');
nc{'oxygen_time'}.units = 'day';
nc{'oxygen_time'}.calendar = ncchar('360.0 days in every year');
nc{'oxygen_time'}.calendar = '360.0 days in every year';
nc{'oxygen_time'}.cycle_length = cycle;
nc{'oxygen_time'}(:)=time;
% 
nc{'oxygen_west'} = ncdouble('oxygen_time','s_rho','eta_rho') ;
nc{'oxygen_west'}.long_name = ncchar('western boundary oxygen');
nc{'oxygen_west'}.long_name = 'western boundary oxygen';
nc{'oxygen_west'}.units = ncchar('mg/L');
nc{'oxygen_west'}.units = 'mg/L';
nc{'oxygen_west'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
nc{'oxygen_west'}.coordinates = 'lon_rho s_rho oxygen_time';
nc{'oxygen_west'}(:)=oxygen_w1;

nc{'oxygen_east'} = ncdouble('oxygen_time','s_rho','eta_rho') ;
nc{'oxygen_east'}.long_name = ncchar('western boundary oxygen');
nc{'oxygen_east'}.long_name = 'western boundary oxygen';
nc{'oxygen_east'}.units = ncchar('mg/L');
nc{'oxygen_east'}.units = 'mg/L';
nc{'oxygen_east'}.coordinates = ncchar('lon_rho s_rho oxygen_time');
nc{'oxygen_east'}.coordinates = 'lon_rho s_rho oxygen_time';
nc{'oxygen_east'}(:)=oxygen_e1;
% 
% nc('chlo_time') = length(time);
% nc{'chlo_time'} = ncdouble('chlo_time') ;
% nc{'chlo_time'}.long_name = ncchar('time for chlo');
% nc{'chlo_time'}.long_name = 'time for chlo';
% nc{'chlo_time'}.units = ncchar('day');
% nc{'chlo_time'}.units = 'day';
% nc{'chlo_time'}.calendar = ncchar('360.0 days in every year');
% nc{'chlo_time'}.calendar = '360.0 days in every year';
% nc{'chlo_time'}.cycle_length = cycle;
% nc{'chlo_time'}(:)=time;
% 
% nc{'chlo_north'} = ncdouble('chlo_time','s_rho','xi_rho') ;
% nc{'chlo_north'}.long_name = ncchar('northern boundary chlo');
% nc{'chlo_north'}.long_name = 'northern boundary chlo';
% nc{'chlo_north'}.units = ncchar('mg/L');
% nc{'chlo_north'}.units = 'mg/L';
% nc{'chlo_north'}.coordinates = ncchar('lon_rho s_rho chlo_time');
% nc{'chlo_north'}.coordinates = 'lon_rho s_rho chlo_time';
% nc{'chlo_north'}(:)=chl;
% 
% nc('zoop_time') = length(time);
% nc{'zoop_time'} = ncdouble('zoop_time') ;
% nc{'zoop_time'}.long_name = ncchar('time for zoop');
% nc{'zoop_time'}.long_name = 'time for zoop';
% nc{'zoop_time'}.units = ncchar('day');
% nc{'zoop_time'}.units = 'day';
% nc{'zoop_time'}.calendar = ncchar('360.0 days in every year');
% nc{'zoop_time'}.calendar = '360.0 days in every year';
% nc{'zoop_time'}.cycle_length = cycle;
% nc{'zoop_time'}(:)=time;
% 
% nc{'zoop_north'} = ncdouble('zoop_time','s_rho','xi_rho') ;
% nc{'zoop_north'}.long_name = ncchar('northern boundary zoop');
% nc{'zoop_north'}.long_name = 'northern boundary zoop';
% nc{'zoop_north'}.units = ncchar('mg/L');
% nc{'zoop_north'}.units = 'mg/L';
% nc{'zoop_north'}.coordinates = ncchar('lon_rho s_rho zoop_time');
% nc{'zoop_north'}.coordinates = 'lon_rho s_rho zoop_time';
% nc{'zoop_north'}(:)=0;
% 
% nc('phyt_time') = length(time);
% nc{'phyt_time'} = ncdouble('phyt_time') ;
% nc{'phyt_time'}.long_name = ncchar('time for phyt');
% nc{'phyt_time'}.long_name = 'time for phyt';
% nc{'phyt_time'}.units = ncchar('day');
% nc{'phyt_time'}.units = 'day';
% nc{'phyt_time'}.calendar = ncchar('360.0 days in every year');
% nc{'phyt_time'}.calendar = '360.0 days in every year';
% nc{'phyt_time'}.cycle_length = cycle;
% nc{'phyt_time'}(:)=time;
% 
% nc{'phyt_north'} = ncdouble('phyt_time','s_rho','xi_rho') ;
% nc{'phyt_north'}.long_name = ncchar('northern boundary phyt');
% nc{'phyt_north'}.long_name = 'northern boundary phyt';
% nc{'phyt_north'}.units = ncchar('mg/L');
% nc{'phyt_north'}.units = 'mg/L';
% nc{'phyt_north'}.coordinates = ncchar('lon_rho s_rho phyt_time');
% nc{'phyt_north'}.coordinates = 'lon_rho s_rho phyt_time';
% nc{'phyt_north'}(:)=0;
% 
% nc('LDeN_time') = length(time);
% nc{'LDeN_time'} = ncdouble('LDeN_time') ;
% nc{'LDeN_time'}.long_name = ncchar('time for LDeN');
% nc{'LDeN_time'}.long_name = 'time for LDeN';
% nc{'LDeN_time'}.units = ncchar('day');
% nc{'LDeN_time'}.units = 'day';
% nc{'LDeN_time'}.calendar = ncchar('360.0 days in every year');
% nc{'LDeN_time'}.calendar = '360.0 days in every year';
% nc{'LDeN_time'}.cycle_length = cycle;
% nc{'LDeN_time'}(:)=time;
% 
% nc{'LDeN_north'} = ncdouble('LDeN_time','s_rho','xi_rho') ;
% nc{'LDeN_north'}.long_name = ncchar('northern boundary LDeN');
% nc{'LDeN_north'}.long_name = 'northern boundary LDeN';
% nc{'LDeN_north'}.units = ncchar('mg/L');
% nc{'LDeN_north'}.units = 'mg/L';
% nc{'LDeN_north'}.coordinates = ncchar('lon_rho s_rho LDeN_time');
% nc{'LDeN_north'}.coordinates = 'lon_rho s_rho LDeN_time';
% nc{'LDeN_north'}(:)=0;
% 
% nc('SDeN_time') = length(time);
% nc{'SDeN_time'} = ncdouble('SDeN_time') ;
% nc{'SDeN_time'}.long_name = ncchar('time for SDeN');
% nc{'SDeN_time'}.long_name = 'time for SDeN';
% nc{'SDeN_time'}.units = ncchar('day');
% nc{'SDeN_time'}.units = 'day';
% nc{'SDeN_time'}.calendar = ncchar('360.0 days in every year');
% nc{'SDeN_time'}.calendar = '360.0 days in every year';
% nc{'SDeN_time'}.cycle_length = cycle;
% nc{'SDeN_time'}(:)=time;
% 
% nc{'SDeN_north'} = ncdouble('SDeN_time','s_rho','xi_rho') ;
% nc{'SDeN_north'}.long_name = ncchar('northern boundary SDeN');
% nc{'SDeN_north'}.long_name = 'northern boundary SDeN';
% nc{'SDeN_north'}.units = ncchar('mg/L');
% nc{'SDeN_north'}.units = 'mg/L';
% nc{'SDeN_north'}.coordinates = ncchar('lon_rho s_rho SDeN_time');
% nc{'SDeN_north'}.coordinates = 'lon_rho s_rho SDeN_time';
% nc{'SDeN_north'}(:)=0;
close(nc);