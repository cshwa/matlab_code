close all
clear all
clc

t_year=2009;

grd_file='F:\ROMS\Sumjin\02_grid_depth\smoothing\grid_sumjin_v11_smoo.nc';
% Fname='pohang_tide_inout.nc';
Fname =['river_' num2str(t_year) '_realts_Gwangyang.nc'];

ss = cumsum(eomday(t_year,1:12));
cycle=ss(end);
cycle= 365;
grd=roms_get_grid(grd_file);
time=[1:1:cycle];
ntime=length(time);
s_rho = 20;
river = 2;
create_river_NH4(grd,Fname,cycle,grd_file,s_rho,river,ntime);
ncout=netcdf(Fname,'write');
% ncout=netcdf(Fname,'clobber');
%% ¼ÛÁ¤ À¯·®, ±Ý°­´ï ¹æ·ù·® load
load D:\mepl\data\03_Sumjin_river\songjung\sj_00_16_year\to_rd.mat
sj_trans_out = [to_rd(:,t_year-1996)]; % 2010:14, 2011:15,2012:16,2013:17,2014:18

%--- 2016 ÀÌÈÄ
% [num,txt,raw] = xlsread('F:\ROMS\Sumjin\07_river\data\nam_2015.xls');
% nam_trans_out =  str2double(txt(5:4+cycle,5)); % 2:À¯ÀÔ·® 3: ¹æ·ù·®

%--- 2015 ±îÁö 
fn = ['F:\ROMS\Sumjin\07_river\data\nam_',num2str(t_year),'.xls'];
[num,txt,raw] = xlsread(fn);
if t_year == 2010 
    num = flipud(num(1:365,:));
    nam_trans_out =  num(:,1);
elseif t_year == 2011
    num = str2double(txt(5:4+cycle,5));
    nam_trans_out =  num(:,1);
elseif t_year >=2012 && t_year <2016
    num = flipud(num(1:365,3));
    nam_trans_out = num(:,1);
elseif t_year >= 2016
    nam_trans_out =  str2double(txt(5:4+cycle,5));
end
 % 2:À¯ÀÔ·® 3: ¹æ·ù·®

% temp_out = 10*ones(1,366);
% f = ['D:\mepl\data\03_OtherData\data4lt\database\±â»óÃ»±â¿Â/atemp' num2str(t_year) '.dat'];
temp_out= load(['D:\mepl\data\03_OtherData\data4lt\database\±â»óÃ»±â¿Â/atemp' num2str(t_year) '.dat']);
temp_out = temp_out(1:31*12);
temp_out(isnan(temp_out))=[];
temp_out = temp_out(1:365);
salt_out = 1*ones(1,365);

N=20; % vertical layer

r=1;  % ¼ÛÁ¤ river flow
Name='outflow_songjung';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 12; % ¿­
River.lat(r) = 174;  % Çâ
River.dir(r) = 0; % 0:u-face, 1:v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=sj_trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
%--- NPZD -----------------------------
% River.NH4(:,:,r)=repmat(NH4_out,1,20);
% River.NO3(:,:,r)=repmat(NO3_out,1,20);
% River.Oxyg(:,:,r)=ones(ntime,20)*0;
% River.tPO4(:,:,r)=repmat(PO4_out,1,20);
%---------------------------------------
River.vshape(r,1:N)=1/N;

r=2; % ³²°­´ï ¹æ·ù·®
Name='outflow_namgang-dam';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 177; % ¿­
River.lat(r) = 174;  % Çà
River.dir(r) = 0; % v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=nam_trans_out;
% River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
% River.NH4(:,:,r)=repmat(NH4_out,1,20);
% River.NO3(:,:,r)=repmat(NO3_out,1,20);
% River.Oxyg(:,:,r)=ones(ntime,20)*0;
% River.tPO4(:,:,r)=rnnepmat(PO4_out,1,20);
River.vshape(r,1:N)=1/N;

%%
%{
% - ÃÑ ¹èÃâ·® : ÀÏÁÖÀÏ¿¡ 5000Åæ, ´ë·« 0.01m^3/s
% - ´ã¼ö À¯ÀÔ·® : 0.001m^3/s ·Î ÇØº¸ÀÚ. ´ã¼ö ÇØ¼ö ºñÀ² :1:9
% r=1;
% Name='Inflow_fresh';
% River.Name(r,1:length(Name)) = Name;
% River.lon(r) = 111; % 112 ¹øÂ° Ä­  
% River.lat(r) = 3;  % 4 ¹øÂ° Ä­ 
% River.dir(r)=0; % u-face
% River.flag(r) = 3; %temp, salt
% River.trans(:,r)=ones(ntime,1)*0.001;
% River.temp(:,:,r)=ones(ntime,20)*10;
% River.salt(:,:,r)=ones(ntime,20)*0;
% River.NH4(:,:,r)=ones(ntime,20)*100;
% River.NO3(:,:,r)=ones(ntime,20)*1000;
% % River.vshape(:,r)=ones(20,1)*(0.4/16);
% % River.vshape(1:4,r)=0.15;
% River.vshape(:,r)=ones(20,1)*(1/20);
% theVarname = 'river';
% ncout{theVarname}(:) = r;
%}

%{
trans_out=0.03;
temp_out=[10 10 10 10 10 10 10 10 10 10 10 10];
salt_out=[33 33 33 33 33 33 33 33 33 33 33 33];
NH4_out=[0 0 0 0 0 0 0 0 0 0 0 0];
NH4_out=[0 0 30 30 30 30 30 30 30 30 30 0];
NO3_out=[7.5 7.5 7.5 7.5 7.5 7.5 7.5  5 0.1 0.1 0.4 7.5];
PO4_out=[0.09 0.09 0.08 0.07 0.06 0.05 0.04 0.04 0.04 0.05 0.07 0.08];


trans_in=0.06;
temp_in=10*ones(365,1);
salt_in=[33 33 33 33 33 33 33 33 33 33 33 33];
NH4_in=[0 0 0 0 0 0 0 0 0 0 0 0];
NO3_in=[7.5 7.5 7.5 5 0.1 0.1 0.1 0.1 0.1 0.1 0.1 7.5];
PO4_in=[0.05 0.05 0.05 0.05 0.04 0.03 0.02 0.02 0.02 0.03 0.04 0.05];
    
r=1;
Name='outflow_pond1';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 54; % 72 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % u-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=2;
Name='outflow_pond2';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 56; % 72 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=3;
Name='outflow_pond3';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 58; % 72 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;
r=4;
Name='outflow_pond4';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 60; % 72 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=5;
Name='outflow_pond5';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 62; % 72 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=6;
Name='outflow_pond6';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 64; % 72 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=7;
Name='outflow_sea';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 66; % 66 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 0; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=8;
Name='outflow_sea';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 68; % 66 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 0; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=9;
Name='outflow_sea';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 70; % 66 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 0; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=10;
Name='outflow_sea';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 72; % 73 ¹øÂ° Ä­
River.lat(r) = 47;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 0; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_out;
River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(NH4_out,1,20);
River.NO3(:,:,r)=repmat(NO3_out,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_out,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=11;
Name='Inflow_sea1';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 33; % 
River.lat(r) = 40;  %
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_in;
River.temp(:,:,r)=repmat(temp_in,1,20);
River.salt(:,:,r)=repmat(salt_in,1,20);
River.NH4(:,:,r)=repmat(NH4_in,1,20);
River.NO3(:,:,r)=repmat(NO3_in,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_in,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=12;
Name='Inflow_sea2';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 39; % 70 ¹øÂ° Ä­
River.lat(r) = 39;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_in;
River.temp(:,:,r)=repmat(temp_in,1,20);
River.salt(:,:,r)=repmat(salt_in,1,20);
River.NH4(:,:,r)=repmat(NH4_in,1,20);
River.NO3(:,:,r)=repmat(NO3_in,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_in,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=13;
Name='Inflow_sea3';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 47; % 70 ¹øÂ° Ä­
River.lat(r) = 38;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_in;
River.temp(:,:,r)=repmat(temp_in,1,20);
River.salt(:,:,r)=repmat(salt_in,1,20);
River.NH4(:,:,r)=repmat(NH4_in,1,20);
River.NO3(:,:,r)=repmat(NO3_in,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_in,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=14;
Name='Inflow_sea4';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 55; % 70 ¹øÂ° Ä­
River.lat(r) = 38;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_in;
River.temp(:,:,r)=repmat(temp_in,1,20);
River.salt(:,:,r)=repmat(salt_in,1,20);
River.NH4(:,:,r)=repmat(NH4_in,1,20);
River.NO3(:,:,r)=repmat(NO3_in,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_in,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=15;
Name='Inflow_sea5';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 63; % 70 ¹øÂ° Ä­
River.lat(r) = 37;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_in;
River.temp(:,:,r)=repmat(temp_in,1,20);
River.salt(:,:,r)=repmat(salt_in,1,20);
River.NH4(:,:,r)=repmat(NH4_in,1,20);
River.NO3(:,:,r)=repmat(NO3_in,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_in,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;

r=16;
Name='Inflow_sea6';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 63; % 70 ¹øÂ° Ä­
River.lat(r) = 36;  % 39 ¹øÂ° Ä­
River.dir(r)=1; % v-face
River.flag(r) = 3; %temp, salt
River.trans(:,r)=ones(ntime,1)*trans_in;
River.temp(:,:,r)=repmat(temp_in,1,20);
River.salt(:,:,r)=repmat(salt_in,1,20);
River.NH4(:,:,r)=repmat(NH4_in,1,20);
River.NO3(:,:,r)=repmat(NO3_in,1,20);
River.Oxyg(:,:,r)=ones(ntime,20)*0;
River.tPO4(:,:,r)=repmat(PO4_in,1,20);
River.vshape(:,r)=zeros(20,1);
River.vshape(1:10,r)=0.1;
% 
% r=17;
% Name='Inflow_sea7';
% River.Name(r,1:length(Name)) = Name;
% River.lon(r) = 64; % 70 ¹øÂ° Ä­
% River.lat(r) = 39;  % 39 ¹øÂ° Ä­
% River.dir(r)=1; % v-face
% River.flag(r) = 3; %temp, salt
% River.trans(:,r)=ones(ntime,1)*trans_in;
% River.temp(:,:,r)=repmat(temp_in,1,20);
% River.salt(:,:,r)=repmat(salt_in,1,20);
% River.NH4(:,:,r)=repmat(NH4_in,1,20);
% River.NO3(:,:,r)=repmat(NO3_in,1,20);
% River.Oxyg(:,:,r)=ones(ntime,20)*0;
% River.tPO4(:,:,r)=repmat(PO4_in,1,20);
% River.vshape(:,r)=zeros(20,1);
% River.vshape(1:10,r)=0.1;
% 
% r=18;
% Name='Inflow_sea8';
% River.Name(r,1:length(Name)) = Name;
% River.lon(r) = 66; % 70 ¹øÂ° Ä­
% River.lat(r) = 39;  % 39 ¹øÂ° Ä­
% River.dir(r)=1; % v-face
% River.flag(r) = 3; %temp, salt
% River.trans(:,r)=ones(ntime,1)*trans_in;
% River.temp(:,:,r)=repmat(temp_in,1,20);
% River.salt(:,:,r)=repmat(salt_in,1,20);
% River.NH4(:,:,r)=repmat(NH4_in,1,20);
% River.NO3(:,:,r)=repmat(NO3_in,1,20);
% River.Oxyg(:,:,r)=ones(ntime,20)*0;
% River.tPO4(:,:,r)=repmat(PO4_in,1,20);
% River.vshape(:,r)=zeros(20,1);
% River.vshape(1:10,r)=0.1;
% 
% r=19;
% Name='Inflow_sea9';
% River.Name(r,1:length(Name)) = Name;
% River.lon(r) = 68; % 70 ¹øÂ° Ä­
% River.lat(r) = 39;  % 39 ¹øÂ° Ä­
% River.dir(r)=1; % v-face
% River.flag(r) = 3; %temp, salt
% River.trans(:,r)=ones(ntime,1)*trans_in;
% River.temp(:,:,r)=repmat(temp_in,1,20);
% River.salt(:,:,r)=repmat(salt_in,1,20);
% River.NH4(:,:,r)=repmat(NH4_in,1,20);
% River.NO3(:,:,r)=repmat(NO3_in,1,20);
% River.Oxyg(:,:,r)=ones(ntime,20)*0;
% River.tPO4(:,:,r)=repmat(PO4_in,1,20);
% River.vshape(:,r)=zeros(20,1);
% River.vshape(1:10,r)=0.1;
% 
% r=20;
% Name='Inflow_sea10';
% River.Name(r,1:length(Name)) = Name;
% River.lon(r) = 70; % 70 ¹øÂ° Ä­
% River.lat(r) = 39;  % 39 ¹øÂ° Ä­
% River.dir(r)=1; % v-face
% River.flag(r) = 3; %temp, salt
% River.trans(:,r)=ones(ntime,1)*trans_in;
% River.temp(:,:,r)=repmat(temp_in,1,20);
% River.salt(:,:,r)=repmat(salt_in,1,20);
% River.NH4(:,:,r)=repmat(NH4_in,1,20);
% River.NO3(:,:,r)=repmat(NO3_in,1,20);
% River.Oxyg(:,:,r)=ones(ntime,20)*0;
% River.tPO4(:,:,r)=repmat(PO4_in,1,20);
% River.vshape(:,r)=zeros(20,1);
% River.vshape(1:10,r)=0.1;

%}

% r=2;
% Name='outflow_sea';
% River.Name(r,1:length(Name)) = Name;
% River.lon(r) = 65; % 66 ¹øÂ° Ä­
% River.lat(r) = 38;  % 39 ¹øÂ° Ä­
% River.dir(r)=1; % v-face
% River.flag(r) = 0; %temp, salt
% % temp1_trans=ones(ntime/2,1)*(0);
% temp2_trans=ones(ntime,1)*(-0.01);
% % temp_trans=[temp1_trans';temp2_trans'];
% % f_trans=reshape(temp_trans,1,[]);
% River.trans(:,r)=temp2_trans';
% 
% temp=ones(ntime,20)*(10);
% River.temp(:,:,r)=temp;
% salt=ones(ntime,20)*(33);
% River.salt(:,:,r)=salt;
% NH4=ones(ntime,20)*(120);
% River.NH4(:,:,r)=NH4;
% NO3=ones(ntime,20)*(10);
% River.NO3(:,:,r)=NO3;
% DO=ones(ntime,20)*0;
% River.DO(:,:,r)=DO;
% PO4=ones(ntime,20)*0;
% River.PO4(:,:,r)=PO4;



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




%% NH4 input

% 'river_NH4'                                        ! Input
%   'river runoff NH4'
%   'millimole_nitrogen meter-3'                     ! [millimole/m3]
%   'river_NH4, scalar, series'
%   'river_time'
%   'idRtrc(iNH4_)'
%   'nulvar'
%   1.0d0
