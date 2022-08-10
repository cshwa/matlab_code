close all; clear; clc;
cd D:\장기생태\Dynamic\06_river
rt=ncread('river_1993to2004_reg_bio_det_nosew_mixed_daily.nc','river_transport');
cd J:\장기생태_2021\Dynamic\result\3regime\daily

plot(rt(1,:))

u=ncread('1993to2004_2yr.nc','u');
v=ncread('1993to2004_2yr.nc','v');
z=ncread('1993to2004_2yr.nc','zeta');

grd_file='D:\장기생태\Dynamic\KOEM\grid_gy_v11_s.nc';
lon_psi=ncread(grd_file, 'lon_psi');
lat_psi=ncread(grd_file, 'lat_psi');
lon=ncread(grd_file, 'lon_rho');
lat=ncread(grd_file, 'lat_rho');
h=ncread(grd_file, 'h');
mask=ncread(grd_file, 'mask_rho');
lon_u=ncread(grd_file, 'lon_u');
lat_u=ncread(grd_file, 'lat_u');
mask_u=ncread(grd_file, 'mask_u');
lon_v=ncread(grd_file, 'lon_v');
lat_v=ncread(grd_file, 'lat_v');
mask_v=ncread(grd_file, 'mask_v');
dx=1./ncread(grd_file, 'pm');
dy=1./ncread(grd_file, 'pn');
N = length(ncread('J:\장기생태_2021\Dynamic\result\3regime\monthly\2reg_cir_2yr_monthly_01.nc','s_rho'));
% depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,1,1,'r');
theta_s=ncread('J:\장기생태_2021\Dynamic\result\3regime\monthly\2reg_cir_2yr_monthly_01.nc','theta_s');
theta_b=ncread('J:\장기생태_2021\Dynamic\result\3regime\monthly\2reg_cir_2yr_monthly_01.nc','theta_b');
hc=ncread('J:\장기생태_2021\Dynamic\result\3regime\monthly\2reg_cir_2yr_monthly_01.nc','hc');

% River.lon = 12; % 열
% River.lat = 174;  % 향
figure;
pcolor(lon,lat,h.*(mask./mask)); shading flat; hold on;
% hold on; plot(lon_u(12,174),lat_u(12,174),'+')
hold on; plot(lon_u(13,175),lat_u(13,175),'r+')
hold on; plot(lon_u(13,175),lat_u(13,176),'r+')
% hold on; plot(lon_u(12,175),lat_u(12,175),'r+')
% hold on; plot(lon_u(14,175),lat_u(14,175),'r+')
text(lon_u(13,175),lat_u(13,175),'+13')
text(lon_u(13,176),lat_u(13,176),'+13')
plot(lon(13,175:176),lat(13,175:176),'c+')
plot(lon(14,175:176),lat(14,175:176),'c+')
plot(lon_v(13,175),lat_v(13,175),'g+')
plot(lon_v(14,175),lat_v(14,175),'g+')
plot(lon(13,175),lat(13,175),'k+')
plot(lon(14,175),lat(14,175),'k+')
plot(lon(14,176),lat(14,176),'k+')
plot(lon(13,176),lat(13,176),'k+')

plot(lon_u(41,175),lat_u(41,175),'r+')
plot(lon_u(42,175),lat_u(42,175),'r+')
plot(lon_v(42,174),lat_v(42,174),'r+')

    plot(lon(:,end),lat(:,end),'r');

% for i = 1:252
%     plot(lon(i,:),lat(i,:),'k');
% end
% 
% for i = 1:176
%     plot(lon_u(:,i),lat_u(:,i),'k');
% end
% for i = 1:251
%     plot(lon_u(i,:),lat_u(i,:),'k');
% end
% 
% for i = 1:175
%     plot(lon_v(:,i),lat_v(:,i),'k');
% end
% for i = 1:252
%     plot(lon_v(i,:),lat_v(i,:),'k');
% end

for i = 1:251
    plot(lon_psi(i,:),lat_psi(i,:),'k');
end

for i = 1:175
    plot(lon_psi(:,i),lat_psi(:,i),'k');
end


% lon_psi(13,:)-lon_psi(12,:)

for i = 1:251
    plot(lon_psi(i,:),lat_psi(i,:),'k');
end

for i = 1:175
    plot(lon_psi(:,i),lat_psi(:,i),'k');
end

ix=13; jx=175;
figure;
plot(squeeze(u(ix,jx,:,:))');

figure;
plot(squeeze(v(ix,jx,:,:))');

figure;
plot(squeeze(u(ix,jx+1,:,:))');

ix=13:41; jx=175;
clearvars depth depth_w dep_l
for i = 1:365
depth(:,:,i)=zlevs(squeeze(h(ix,jx)),squeeze(z(ix,jx,i)),theta_s,theta_b,hc,N,'r');
depth_w(:,:,i)=zlevs(squeeze(h(ix,jx)),squeeze(z(ix,jx,i)),theta_s,theta_b,hc,N,'w');
end

for i = N:-1:1
        dep_l(i,:,:) =  depth_w(i+1,:,:) - depth_w(i,:,:);
end

cal_trans = (dy(ix,jx)) .* dep_l .* squeeze(u(ix,jx,:,:));
cal_trans_sum = squeeze(sum(cal_trans,1));

figure; plot(cal_trans_sum,'r'); hold on; plot(rt(1,:));

figure; plot(cal_trans_sum - rt(1,:),'r'); 
yline(mean(cal_trans_sum - rt(1,:)),'c'); 

lon(ix,jx) - lon(ix+1,jx)                 
lat(ix,jx) - lat(ix,jx+1)  

u(13,175,20,:)
v(13,175,20,:)
h(ix,jx)

