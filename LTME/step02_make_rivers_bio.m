%%% river input - river_NO3, river_NH4, river_SDeN,
%%% river_DIC, river_alkalinity, river_Oxyg, 

close all
clear all
clc

NO3_han=load('NO3_han_1996-2017(kimpo)_input.txt');NH4_han=load('NH4_han_1996-2017(kimpo)_input.txt');DO_han=load('DO_han_1996-2017(kimpo)_input.txt');PO4_han=load('PO4_han_1996-2017(kimpo)_input.txt');
NO3_nak=load('NO3_nak_1996-2017(koo)_input.txt');NH4_nak=load('NH4_nak_1996-2017(koo)_input.txt');DO_nak=load('DO_nak_1996-2017(koo)_input.txt');PO4_nak=load('PO4_nak_1996-2017(koo)_input.txt');
NO3_keum=load('NO3_keum_1996-2017(kang)_input.txt');NH4_keum=load('NH4_keum_1996-2017(kang)_input.txt');DO_keum=load('DO_keum_1996-2017(kang)_input.txt');PO4_keum=load('PO4_keum_1996-2017(kang)_input.txt');
NO3_sum=load('NO3_sum_1996-2017(hadong)_input.txt');NH4_sum=load('NH4_sum_1996-2017(hadong)_input.txt');DO_sum=load('DO_sum_1996-2017(hadong)_input.txt');PO4_sum=load('PO4_sum_1996-2017(hadong)_input.txt');
NO3_young=load('NO3_young_1996-2017(moo1)_input.txt');NH4_young=load('NH4_young_1996-2017(moo1)_input.txt');DO_young=load('DO_young_1996-2017(moo1)_input.txt');PO4_young=load('PO4_young_1996-2017(moo1)_input.txt');
startyy=1996;
file=2017;
grd_file='g:\auto_fennel\grid\roms_grd_auto_rdrg2_new8_smooth.nc';
% Fname='pohang_tide_inout.nc';
% Fname = ['roms_river_auto_new8_',num2str(file),'_fennel_po4_8rivers.nc'];
Fname = ['roms_river_auto_new8_',num2str(file),'_yangtze80s_fennel_po4_8rivers.nc'];
% Fname = ['roms_river_auto_new8_',num2str(file),'_inYS_fennel_po4_8rivers.nc'];

cycle=yeardays(file);
grd=roms_get_grid(grd_file);
time=[15:30:cycle];
ntime=length(time);
s_rho = 40;
river = 12;
create_river_NH4(grd,Fname,cycle,grd_file,s_rho,river,ntime);
ncout=netcdf(Fname,'write');
%%% DIN concentration [Dai et al. 2011]
%%% 1960s 1970s 1980s 1990s 2000s
%%% 18.4  21.6  57.1  105.1 123.6 [μmol/L]
% NO3_out_yangze= [50; 60; 70; 80; 90; 75; 75; 70; 60; 60; 65; 60];  % duan et al., 2008
% NO3_out_yangze= ones(1,ntime)'*123.6*(0.95); % 2000s
NO3_out_yangze= ones(1,ntime)'*57.1*(0.95); %1980s
NO3_out_yellow= ones(1,ntime)'*313.7; %Gong et al. 2015 2002~2004년자료
NO3_out_yalu= ones(1,ntime)'*293.03;   % Liu and zhang et al.2004 JCR (Nutrient dynamics in the macro tidal yalujiang estuary) 1992,94,96년 자료
% 단위 환산 mg/L -> millimole N m-3 (*1000/14) 
NO3_out_han= NO3_han((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NO3_out_keum= NO3_keum((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NO3_out_young= NO3_young((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NO3_out_sum= NO3_sum((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NO3_out_nak= NO3_nak((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NO3_out_0 = zeros(1,ntime)';

% NH4_out_yangze= ones(1,ntime)'*123.6*(0.05); % Zhou et al., 2017 에 wang et al. 2011 관측치 기반.(dai et al. 2011 DIN의 1/20)
NH4_out_yangze= ones(1,ntime)'*57.1*(0.05); % Zhou et al., 2017 에 wang et al. 2011 관측치 기반.(dai et al. 2011 DIN의 1/20)
NH4_out_yellow= ones(1,ntime)'*21.6; %Gong et al. 2015 2002~2004년자료
NH4_out_yalu= ones(1,ntime)'*4.8; % Liu and zhang et al.2004 JCR
% 단위 환산 mg/L -> millimole N m-3 (*1000/14) 
NH4_out_han= NH4_han((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NH4_out_keum= NH4_keum((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NH4_out_young= NH4_young((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NH4_out_sum= NH4_sum((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NH4_out_nak= NH4_nak((file-startyy)*12+1:(file-startyy+1)*12)*1000/14;
NH4_out_0 = zeros(1,ntime)';

Oxyg_out_yangze= [300; 300; 300; 300; 300; 300; 300; 300; 300; 300; 300; 300];
% 단위 환산 mg/L -> millimole N m-3  (*0.7*44.661);
Oxyg_out_han= DO_han((file-startyy)*12+1:(file-startyy+1)*12)*0.7*44.661;
Oxyg_out_keum= DO_keum((file-startyy)*12+1:(file-startyy+1)*12)*0.7*44.661;
Oxyg_out_young= DO_young((file-startyy)*12+1:(file-startyy+1)*12)*0.7*44.661;
Oxyg_out_sum= DO_sum((file-startyy)*12+1:(file-startyy+1)*12)*0.7*44.661;
Oxyg_out_nak= DO_nak((file-startyy)*12+1:(file-startyy+1)*12)*0.7*44.661;
Oxyg_out_0 = ones(1,ntime)'*300;

tPO4_out_yangze= ones(1,ntime)'*1.42; % Dai et al., 2011 DIP 농도  
tPO4_out_yellow= ones(1,ntime)'*0.34;  %Gong et al. 2015 2002~2004년자료
tPO4_out_yalu= ones(1,ntime)'*0.11;  % Liu and zhang et al.2004 JCR
% 단위 환산 mg/L -> millimole P m-3 (*1000./30.9)
tPO4_out_han= PO4_han((file-startyy)*12+1:(file-startyy+1)*12)*1000./30.9;
tPO4_out_keum= PO4_keum((file-startyy)*12+1:(file-startyy+1)*12)*1000./30.9;
tPO4_out_young= PO4_young((file-startyy)*12+1:(file-startyy+1)*12)*1000./30.9;
tPO4_out_sum= PO4_sum((file-startyy)*12+1:(file-startyy+1)*12)*1000./30.9;
tPO4_out_nak= PO4_nak((file-startyy)*12+1:(file-startyy+1)*12)*1000./30.9;
tPO4_out_0 = zeros(1,ntime)';
%%%
detritus_out= 0*ones(1,ntime)';
phytoplankton_out= 0*ones(1,ntime)';
zooplankton_out= 0*ones(1,ntime)';

N=s_rho; % vertical layer

r=1;  %양자강
River.detritus(:,:,r)=repmat(detritus_out,1,N);
River.NO3(:,:,r)=repmat(NO3_out_yangze,1,N);
River.NH4(:,:,r)=repmat(NH4_out_yangze,1,N);
River.Oxyg(:,:,r)=repmat(Oxyg_out_yangze,1,N);
River.tPO4(:,:,r)=repmat(tPO4_out_yangze,1,N);
% River.vshape(r,1:N)=1/N;

r=2; %황하
River.detritus(:,:,r)=repmat(detritus_out,1,N);
River.NO3(:,:,r)=repmat(NO3_out_yellow,1,N);
River.NH4(:,:,r)=repmat(NH4_out_yellow,1,N);
River.Oxyg(:,:,r)=repmat(Oxyg_out_0,1,N);
River.tPO4(:,:,r)=repmat(tPO4_out_yellow,1,N);

for r=3:7
    River.detritus(:,:,r)=repmat(detritus_out,1,N);
    River.NO3(:,:,r)=repmat(NO3_out_0,1,N);
    River.NH4(:,:,r)=repmat(NH4_out_0,1,N);
    River.Oxyg(:,:,r)=repmat(Oxyg_out_0,1,N);
    River.tPO4(:,:,r)=repmat(tPO4_out_0,1,N);
%     River.vshape(r,1:N)=1/N;
end

r=7; %압록강
River.detritus(:,:,r)=repmat(detritus_out,1,N);
River.NO3(:,:,r)=repmat(NO3_out_yalu,1,N);
River.NH4(:,:,r)=repmat(NH4_out_yalu,1,N);
River.Oxyg(:,:,r)=repmat(Oxyg_out_0,1,N);
River.tPO4(:,:,r)=repmat(tPO4_out_yalu,1,N);

r=8; %한강
River.detritus(:,:,r)=repmat(detritus_out,1,N);
River.NO3(:,:,r)=repmat(NO3_out_han,1,N);
River.NH4(:,:,r)=repmat(NH4_out_han,1,N);
River.Oxyg(:,:,r)=repmat(Oxyg_out_han,1,N);
River.tPO4(:,:,r)=repmat(tPO4_out_han,1,N);

r=9; %금강
River.detritus(:,:,r)=repmat(detritus_out,1,N);
River.NO3(:,:,r)=repmat(NO3_out_keum,1,N);
River.NH4(:,:,r)=repmat(NH4_out_keum,1,N);
River.Oxyg(:,:,r)=repmat(Oxyg_out_keum,1,N);
River.tPO4(:,:,r)=repmat(tPO4_out_keum,1,N);

r=10; %영산
River.detritus(:,:,r)=repmat(detritus_out,1,N);
River.NO3(:,:,r)=repmat(NO3_out_young,1,N);
River.NH4(:,:,r)=repmat(NH4_out_young,1,N);
River.Oxyg(:,:,r)=repmat(Oxyg_out_young,1,N);
River.tPO4(:,:,r)=repmat(tPO4_out_young,1,N);

r=11; %섬진
River.detritus(:,:,r)=repmat(detritus_out,1,N);
River.NO3(:,:,r)=repmat(NO3_out_sum,1,N);
River.NH4(:,:,r)=repmat(NH4_out_sum,1,N);
River.Oxyg(:,:,r)=repmat(Oxyg_out_sum,1,N);
River.tPO4(:,:,r)=repmat(tPO4_out_sum,1,N);

r=12; %낙동
River.detritus(:,:,r)=repmat(detritus_out,1,N);
River.NO3(:,:,r)=repmat(NO3_out_nak,1,N);
River.NH4(:,:,r)=repmat(NH4_out_nak,1,N);
River.Oxyg(:,:,r)=repmat(Oxyg_out_nak,1,N);
River.tPO4(:,:,r)=repmat(tPO4_out_nak,1,N);

theVarname = 'river';
ncout{theVarname}(:) = [1:r];

% theVarname = 'river_Xposition';
% ncout{theVarname}(:) = River.lon;
% 
% theVarname = 'river_Eposition';
% ncout{theVarname}(:) = River.lat;
% 
% theVarname = 'river_direction';
% ncout{theVarname}(:) = River.dir;

% theVarname = 'river_flag';
% ncout{theVarname}(:) = River.flag;

% theVarname = 'river_Vshape';
% ncout{theVarname}(:,:) = River.vshape;

% theVarname = 'river_time';
% ncout{theVarname}(:) = time';

% theVarname = 'river_transport';
% ncout{theVarname}(:,:) = River.trans;

% theVarname = 'river_temp';
% ncout{theVarname}(:,:,:) = River.temp;
% 
% theVarname = 'river_salt';
% ncout{theVarname}(:,:,:) = River.salt;

theVarname = 'river_SDeN';
ncout{theVarname}(:,:,:) = River.detritus;

theVarname = 'river_LDeN';
ncout{theVarname}(:,:,:) = River.detritus;

theVarname = 'river_NO3';
ncout{theVarname}(:,:,:) = River.NO3;

theVarname = 'river_NH4';
ncout{theVarname}(:,:,:) = River.NH4;

theVarname = 'river_Oxyg';
ncout{theVarname}(:,:,:) = River.Oxyg;

theVarname = 'river_tPO4';
ncout{theVarname}(:,:,:) = River.tPO4;


close(ncout)


