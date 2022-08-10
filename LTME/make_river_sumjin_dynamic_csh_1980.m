close all; clear all; clc
cd E:\장기생태\Dynamic\06_river

for t_year = 1980:1989

grd_file='grid_sumjin_v1970_fix_3m.nc';
% Fname='pohang_tide_inout.nc';
Fname =['river_' num2str(t_year) '_realts_Gwangyang.nc'];
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
cd E:\장기생태\Dynamic\06_river\data\sj_1980s
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
load songjung_estimate_data_ykang.mat  % pre_merg_dis is discharge

cd E:\장기생태\Dynamic\06_river
% approximation : Sumjin(songjung) & Namgang have same discharge during 1980~1989.  

cd E:\장기생태\Dynamic\06_river\data\기상청기온\gy_airtemp\
fn1 = ['남해_1980s.xlsx']
[num1,txt1,raw1] = xlsread(fn1,num2str(t_year),'');
temp_out_namhe=num1(:,3);
fn2 = ['진주_1980s.xlsx']
[num2,txt2,raw2] = xlsread(fn2,num2str(t_year),'');
temp_out_jinju=num2(:,3);

find(isnan(temp_out_namhe)==1)
find(isnan(temp_out_jinju)==1)

temp_out_namhe(isnan(temp_out_namhe))=[];
temp_out_jinju(isnan(temp_out_jinju))=[];
salt_out = zeros(cycle,1); 
% salt_out = 1*ones(1,365);  %silver
sj_trans_out = pre_merg_dis{t_year-1979};


N=20; % vertical layer

r=1;  % 송정 river flow
Name='outflow_songjung';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 12; % 열
River.lat(r) = 174;  % 향
River.dir(r) = 0; % 0:u-face, 1:v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=sj_trans_out;
River.temp(:,:,r)=repmat(temp_out_namhe,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
%--- NPZD -----------------------------
% River.NH4(:,:,r)=repmat(NH4_out,1,20);
% River.NO3(:,:,r)=repmat(NO3_out,1,20);
% River.Oxyg(:,:,r)=ones(ntime,20)*0;
% River.tPO4(:,:,r)=repmat(PO4_out,1,20);
%---------------------------------------
River.vshape(r,1:N)=1/N;

r=2; % 남강댐 방류량
Name='outflow_sacheon_bay';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 177; % 열
River.lat(r) = 174;  % 행
River.dir(r) = 0; % v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=sj_trans_out;
% River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out_jinju,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
% River.NH4(:,:,r)=repmat(NH4_out,1,20);
% River.NO3(:,:,r)=repmat(NO3_out,1,20);
% River.Oxyg(:,:,r)=ones(ntime,20)*0;
% River.tPO4(:,:,r)=rnnepmat(PO4_out,1,20);
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

% theVarname = 'river_NH4';
% ncout{theVarname}(:,:,:) = River.NH4;
% 
% theVarname = 'river_NO3';
% ncout{theVarname}(:,:,:) = River.NO3;
% 
% theVarname = 'river_Oxyg';
% ncout{theVarname}(:,:,:) = River.Oxyg;
% 
% theVarname = 'river_tPO4';
% ncout{theVarname}(:,:,:) = River.tPO4;
close(ncout)

clearvars -except t_year
cd E:\장기생태\Dynamic\06_river
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
temp=ncread('river_2009_realts_Gwangyang.nc','river_temp');
dis=ncread('river_2009_realts_Gwangyang.nc','river_transport');

temp=ncread('river_1980_realts_Gwangyang.nc','river_temp');
dis=ncread('river_1980_realts_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1989_realts_Gwangyang.nc','river_temp');
dis=ncread('river_1989_realts_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');