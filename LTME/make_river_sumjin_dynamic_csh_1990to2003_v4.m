close all; clear all; clc
cd D:\장기생태\Dynamic\06_river

for t_year = 1990:2003 %before regime change
   
grd_file='grid_gy_v11_s.nc';
% Fname='pohang_tide_inout.nc';
Fname =['river_' num2str(t_year) '_realts_biofennel_Gwangyang.nc'];
ss = cumsum(eomday(t_year,1:12));
cycle=ss(end);

grd=roms_get_grid(grd_file);
time=[1:1:cycle];
ntime=length(time);
s_rho = 20;
river = 2;
create_river_NH4(grd,Fname,cycle,grd_file,s_rho,river,ntime);
ncout=netcdf(Fname,'write');
% ncout=netcdf(Fname,'clobber');


%% 송정 유량, 금강댐 방류량 load
cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
load songjung_discharge_1980to2018.mat  % pre_merg_dis is discharge

cd D:\장기생태\Dynamic\06_river
% approximation : Sumjin(songjung) & Namgang have same discharge during 1980~1989.  

load('gawha_recons_water_temp_present.mat','merg_recon_w_c');
gahwa_re_w_c = merg_recon_w_c;
clearvars merg_recon_w_c
load('sumjin_recons_water_temp_present.mat','merg_recon_w_c');
sumjin_re_w_c = merg_recon_w_c;

load('sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig.mat','yp_w_*');   % polynomial_fitting_on_ecological_variable_to2004_new_v5_3sigma.m 
% ~2003 : advanced 전체기간
% 2004~ : 04~ data's climate
% load('sumjin(songjung)_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime

sumjin_do = yp_w_do_04'*0.7*44.661;
sumjin_chl = yp_w_chl_04';
sumjin_no3 = yp_w_no3_04'*1000/14;
sumjin_nh4 = yp_w_nh4_04'*1000/14;

load D:\장기생태\Dynamic\06_river\환경과학원\sumjin(songjung)_polynomial_climate_to2011_3sig_po4.mat % ~2010 : 1st regime(yp_w_po4_04), 2011~2018 : 2nd(yp_w_po4_af)
sumjin_po4 = yp_w_po4_04'.*1000./30.973762;

clearvars yp_w_*
load('Namgang_polynomial_climate_to2004_advanced(v4)_3sig.mat','yp_w_*');
% load('Namgang_polynomial_climate_to2004_3sig.mat','yp_w_*'); % ~2004 : 1st regime, 2004~2018 : 2nd regime
nam_do = yp_w_do_04'*0.7*44.661;  %% mg/L to mM/m^3;
nam_chl = yp_w_chl_04'; 
% nam_chl(364:366) = nam_chl(363); % negative concen. to be closest sample.
nam_no3 = yp_w_no3_04'*1000/14; %% mg/L to mM/m^3;
nam_nh4 = yp_w_nh4_04'*1000/14; %% mg/L to mM/m^3;

clearvars yp_w_po4_04
load D:\장기생태\Dynamic\06_river\환경과학원\Namgang_polynomial_climate_to2011_3sig_po4.mat % ~2010 : 1st regime(yp_w_po4_04), 2011~2018 : 2nd(yp_w_po4_af)
nam_po4 = yp_w_po4_04'.*1000./30.973762;

%if it's not leap yr.. delete 2/29
if cycle == 365
    sumjin_do(60)= []; 
    sumjin_chl(60)= []; 
    sumjin_no3(60)= []; 
    sumjin_nh4(60)= []; 
    sumjin_po4(60)= []; 
    nam_do(60)= []; 
    nam_chl(60)= []; 
    nam_no3(60)= []; 
    nam_nh4(60)= []; 
    nam_po4(60)= []; 
end

temp_out_sumjin=sumjin_re_w_c{t_year-1989};
temp_out_gahwa=gahwa_re_w_c{t_year-1989};
salt_out = zeros(cycle,1); 
% salt_out = 1*ones(1,365);  %silver
sj_trans_out = dis_pre_total{t_year-1979}; %dis_pre_total

if cycle == 365 & length(sj_trans_out) == 366
    sj_trans_out(60) = [];
end

N=20; % vertical layer

r=1;  % 송정 river flow
Name='outflow_songjung';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 12; % 열
River.lat(r) = 174;  % 향
River.dir(r) = 0; % 0:u-face, 1:v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=sj_trans_out;
River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
%--- NPZD -----------------------------
River.NH4(:,:,r)=repmat(sumjin_nh4,1,20);
River.NO3(:,:,r)=repmat(sumjin_no3,1,20);
River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
River.chlorophyll(:,:,r)=repmat(sumjin_chl,1,20);
River.tPO4(:,:,r)=repmat(sumjin_po4,1,20);
%---------------------------------------
River.vshape(r,1:N)=1/N;

r=2; % 남강댐 방류량
Name='outflow_sacheon_bay';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 177; % 열
River.lat(r) = 174;  % 행
River.dir(r) = 0; % v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=sj_trans_out.*0.7;
% River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out_gahwa,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(nam_nh4,1,20);
River.NO3(:,:,r)=repmat(nam_no3,1,20);
River.Oxyg(:,:,r)=repmat(nam_do,1,20);
River.chlorophyll(:,:,r)=repmat(nam_chl,1,20);
River.tPO4(:,:,r)=repmat(nam_po4,1,20);
River.vshape(r,1:N)=1/N;

theVarname = 'river';
ncout{theVarname}(:) = [1:r];

theVarname = 'river_Xposition';
ncout{theVarname}(:) = River.lon;

theVarname = 'river_Eposition';
ncout{theVarname}(:) = River.lat;

theVarname = 'river_direction';
ncout{theVarname}(:) = River.dir;

% theVarname = 'river_flag';
% ncout{theVarname}(:) = River.flag;

theVarname = 'river_Vshape';
ncout{theVarname}(:,:) = River.vshape;

theVarname = 'river_time';
ncout{theVarname}(:) = time';

theVarname = 'river_transport';
ncout{theVarname}(:,:) = River.trans;

theVarname = 'river_temp';
ncout{theVarname}(:,:,:) = River.temp;

theVarname = 'river_salt';
ncout{theVarname}(:,:,:) = River.salt;

theVarname = 'river_NH4';
ncout{theVarname}(:,:,:) = River.NH4;

theVarname = 'river_NO3';
ncout{theVarname}(:,:,:) = River.NO3;

theVarname = 'river_Oxyg';
ncout{theVarname}(:,:,:) = River.Oxyg;
 
theVarname = 'river_chlorophyll';
ncout{theVarname}(:,:,:) = River.chlorophyll;

theVarname = 'river_tPO4';
ncout{theVarname}(:,:,:) = River.tPO4;
close(ncout)

clearvars -except t_year
cd D:\장기생태\Dynamic\06_river
end


%% NH4 input

% 'river_NH4'                                        ! Input
%   'river runoff NH4'
%   'millimole_nitrogen meter-3'                     ! [millimole/m3]
%   'river_NH4, scalar, series'
%   'river_time'
%   'idRtrc(iNH4_)'
%   'nulvar'
%   1.0d0


%% check

temp=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_temp');
dis=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_temp');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');



temp=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_Oxyg');
dis=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_Oxyg');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');


temp=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_NH4');
dis=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_NH4');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');


temp=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_NO3');
dis=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_NO3');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');


temp=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_chlorophyll');
dis=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_chlorophyll');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

temp=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_tPO4');
dis=ncread('river_1990_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_tPO4');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');
