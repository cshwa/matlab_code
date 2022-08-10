% dynamic model 용 유량 만들기 

close all; clear all; clc

% addpath(genpath('G:\ROMS\mathlab_setpath/silver/romsplot'));
% addpath(genpath('G:\ROMS\mathlab_setpath/silver/netcdf_old/'));

%--- set target year ------------------------------------------------------
t_year = 2001;
t_mm = 11;
ss = cumsum(eomday(t_year,1:12));
cycle = ss(end);
cycle= 365;

%--- get grid file --------------------------------------------------------
grd_file='../domain/grid_gy_v11_s_9_pollution.nc'; %0704는 수심만 깊어진것
grd = roms_get_grid(grd_file);
s_rho = 20;
%--- file name for new river nc file --------------------------------------
% Fname =['result_nc/river_gy_v11_thermal_' num2str(t_year) '_1204.nc'];

time=[1:1:cycle]';
ntime=length(time);

%--- total number of river sources ----------------------------------------
river = 11;

load('gy_pollution.mat');

%% load river discharge rate

ncout=netcdf('./river_2001_nemuro_con.nc');
river_transport = ncout{'river_transport'}(:);
river_temp = ncout{'river_temp'}(:);
river_salt = ncout{'river_salt'}(:);
sj_trans_out = abs(river_transport(:,1));
nam_trans_out = abs(river_transport(:,2));
temp_out = squeeze(river_temp(:,20,1));
salt_out = squeeze(river_salt(:,20,1));

%%
N=20; % vertical layer
%--- 송정 -----------------------------------------------------------------
r = 1;  
Name='outflow_songjung';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 57; % 열 % org:
River.lat(r) = 174;  % 행
River.dir(r) = 0; % 0:u-face, 1:v-face
River.trans(:,r)=-riv_s_dis_f;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv_s_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv_s_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv_s_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv_s_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) = (1/N);

%--- 남강댐 ---------------------------------------------------------------
r = 2; 
Name='outflow_namgang-dam';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 177; % 열
River.lat(r) = 174;  % 행
River.dir(r) = 0; % 0:u-face, 1:v-face
River.trans(:,r)=nam_trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
river_Vshape(1:N,r) = (1/N);

%--- 광양산업단지폐수종말처리장----------------------------------------------
r = 3; 
Name='gy_ind';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 68;
River.lat(r) = 112;
River.dir(r) = 1; % v-face
River.trans(:,r)=-riv1_dis_f/86400;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv1_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv1_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv1_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv1_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) = 1/N;

%---광양중앙하수종말처리장---------------------------------------------------
r = 4; 
Name='gy_ind2';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 29;
River.lat(r) = 99;  
River.dir(r) = 0; % u-face
River.trans(:,r)=riv2_dis_f/86400;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv2_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv2_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv2_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv2_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) =  1/N;

%--- 광양하수종말처리장 ----------------------------------------------------
r = 5; 
Name='gy_ind3';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 3; % 열
River.lat(r) = 106;  % 행
River.dir(r) = 1; % v-face
River.trans(:,r)=-riv3_dis_f/86400;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv3_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv3_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv3_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv3_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) =  1/N;

%--- 광영하수종말처리장-----------------------------------------------------
r = 6; 
Name='gy_ind4';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 49; % 열
River.lat(r) = 121;  % 행
River.dir(r) = 0; % v-face
River.trans(:,r)=riv4_dis_f/86400;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv4_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv4_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv4_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv4_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) =  1/N;
%--- 여수월내 -----------------------------------------
r = 7; 
Name='ys_ind1';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 56; % 열
River.lat(r) = 63;  % 행
River.dir(r) = 1; % 0:u-face, 1:v-face
River.trans(:,r) = riv5_dis_f/86400; % positive: northward flow 
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv5_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv5_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv5_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv5_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) =  1/N;

%--- 여수율촌 -------------------------------------------
r = 8; 
Name='ys_ind2';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 5; % 열
River.lat(r) = 84;  % 행
River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
River.trans(:,r) = -riv6_dis_f/86400; % negative: sourthward flow 
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv6_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv6_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv6_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv6_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) =  1/N;

%--- 여수중흥 -------------------------------------------
r = 9; 
Name='ys_ind3';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 24; % 열
River.lat(r) = 56;  % 행
River.dir(r) = 0; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
River.trans(:,r) = riv7_dis_f/86400; % positive: northward flow 
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv7_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv7_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv7_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv7_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) =  1/N;

%--- 여수하수 -------------------------------------------
r = 10;
Name='ys_in4';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 28; % 열
River.lat(r) = 23;  % 행
River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
River.trans(:,r) = -riv8_dis_f/86400; % positive: northward flow 
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv8_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv8_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv8_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv8_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) =  1/N;

%--- 진월하수 -------------------------------------------
r = 11; 
Name='ys_ind5';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 61; % 열
River.lat(r) = 126;  % 행
River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
River.trans(:,r) = -riv9_dis_f/86400; % positive: northward flow 
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(riv9_TN_f*(140/1000000)*(1000/14),1,20);
River.NO3(:,:,r)=repmat(riv9_TN_f*(136/1000000)*(1000/14),1,20);
River.DON(:,:,r)=repmat(riv9_TN_f*(999724/2/1000000)*(1000/14),1,20);
River.PON(:,:,r)=repmat(riv9_TN_f*(999724/2/1000000)*(1000/14),1,20);
river_Vshape(1:N,r) =  1/N;

%%
%--- make empty river source nc file---------------------------------------
% Fname =['result_nc/river_dye_' num2str(2019) '_posco_1212.nc'];
Fname =['gy_river_9pollution_lowNH4.nc'];
create_river_NH4(grd,Fname,cycle,grd_file,s_rho,river,ntime);
ncout=netcdf(Fname,'write');

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
% ncout{theVarname}(:,:) = River.vshape;
ncout{'river_Vshape'}(:,:) = river_Vshape;

theVarname = 'river_time';
ncout{theVarname}(:) = time;

theVarname = 'river_transport';
ncout{theVarname}(:,:) = River.trans;

theVarname = 'river_temp';
ncout{theVarname}(:,:,:) = River.temp;

theVarname = 'river_salt';
ncout{theVarname}(:,:,:) = River.salt;

theVarname = 'river_NO3';
ncout{theVarname}(:,:,:) = River.NO3;

theVarname = 'river_NH4';
ncout{theVarname}(:,:,:) = River.NH4*(1/5);

theVarname = 'river_DON';
ncout{theVarname}(:,:,:) = River.DON*(1/5);

theVarname = 'river_PON';
ncout{theVarname}(:,:,:) = River.PON*(1/5);

close(ncout)
