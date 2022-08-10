close all
clear all
clc

file=2013;

grd_file='../02_grid_depth\smoothing\grid_gy_v11_s.nc';
% Fname='pohang_tide_inout.nc';
Fname = 'ocean_rivers2013_ts_real.nc';

cycle=365;
grd=roms_get_grid(grd_file);
time=[1:1:365];
ntime=length(time);
s_rho = 20;
river = 2;
create_river_mks(grd,Fname,cycle,grd_file);
ncout=netcdf(Fname,'write');
% ncout=netcdf(Fname,'clobber');
%% 송정 유량, 금강댐 방류량 load
[num,txt,raw] = xlsread('data\sj_2013.xls');
num = flipud(num(1:365,:));
sj_trans_out =  num(:,2)';

[num,txt,raw] = xlsread('data\nam_2013.xls');
num = flipud(num(1:365,:));
nam_trans_out =  num(:,3)'; % 2:유입량 3: 방류량

% temp_out = 10*ones(1,366);
temp_out= load('data/sj_Tair2012.dat');
salt_out = 1*ones(1,365);

N=20; % vertical layer

r=1;  % 송정 river flow
Name='outflow_songjung';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 12; % 열
River.lat(r) = 174;  % 향
River.dir(r) = 0; % 0:u-face, 1:v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=sj_trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.vshape(r,1:N)=1/N;

r=2; % 남강댐 방류량
Name='outflow_namgang-dam';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 177; % 열
River.lat(r) = 174;  % 행
River.dir(r) = 0; % v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=nam_trans_out;
% River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.vshape(r,1:N)=1/N;

%%
% River.vshape(:,r)=ones(20,1)*(0);
% River.vshape(1:2,r)=0.5;

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
ncout{theVarname}(:) = time;

theVarname = 'river_transport';
ncout{theVarname}(:,:) = River.trans;

theVarname = 'river_temp';
ncout{theVarname}(:,:,:) = River.temp;

theVarname = 'river_salt';
ncout{theVarname}(:,:,:) = River.salt;

close(ncout)


