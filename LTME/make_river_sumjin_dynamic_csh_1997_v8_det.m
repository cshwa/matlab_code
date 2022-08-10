close all; clear all; clc
cd F:\ROMS\roms_tools\Run
start

cd D:\장기생태\Dynamic\06_river


for t_year = 1997:1997 %before regime change
   
grd_file='grid_gy_v11_s.nc';
% Fname='pohang_tide_inout.nc';
Fname =['river_' num2str(t_year) '_realts_biofennel_GY_sewer_det_v8.nc'];
ss = cumsum(eomday(t_year,1:12));
cycle=ss(end);

grd=roms_get_grid(grd_file);
time=[1:1:cycle];
ntime=length(time);
s_rho = 20;
% river = 2;
river = 11;
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

%TN
sumjin_tn=load('D:\장기생태\Dynamic\06_river\환경과학원\sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig_TN.mat','yp_w_tn_04');
nam_tn=load('D:\장기생태\Dynamic\06_river\환경과학원\Namgang_polynomial_climate_to2004_advanced(v4)_3sig_TN.mat','yp_w_tn_04');


%PO4
clearvars yp_w_po4_04
load D:\장기생태\Dynamic\06_river\환경과학원\Namgang_polynomial_climate_to2011_3sig_po4.mat % ~2010 : 1st regime(yp_w_po4_04), 2011~2018 : 2nd(yp_w_po4_af)
nam_po4 = yp_w_po4_04'.*1000./30.973762;

%% 하수종말처리장
cd D:\장기생태\Dynamic\06_river\
load('sewer_monthly_climate_fix.mat');

cd D:\장기생태\Dynamic\06_river
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
    sumjin_tn.yp_w_tn_04(60) =[];
    nam_tn.yp_w_tn_04(60)=[];
            %dis
        clim_dis_g_d(60)=[];
        clim_dis_gc_d(60)=[];
        clim_dis_gh_d(60)=[];
        clim_dis_gh2_d(60)=[];
        clim_dis_d(60)=[];
        clim_dis_y_d(60)=[];
        clim_dis_j_d(60)=[];
        clim_dis_ye_d(60)=[];
        clim_dis_jw_d(60)=[];
        %TN
        clim_tn_g_ns_d(60)=[];
        clim_tn_gc_ns_d(60)=[];
        clim_tn_gh_ns_d(60)=[];
        clim_tn_gh2_ns_d(60)=[];
        clim_tn_w_1st_ns_d(60)=[];
        clim_tn_w_2nd_ns_d(60)=[];
        clim_tn_w_3rd_ns_d(60)=[];
        clim_tn_y_ns_d(60)=[];
        clim_tn_j_1st_ns_d(60)=[];
        clim_tn_j_2nd_ns_d(60)=[];
        clim_tn_j_3rd_ns_d(60)=[];
        clim_tn_ye_ns_d(60)=[];
        clim_tn_jw_ns_d(60)=[];
        %TP
        clim_tp_g_ns_d(60)=[];
        clim_tp_gc_ns_d(60)=[];
        clim_tp_gh_ns_d(60)=[];
        clim_tp_gh2_ns_d(60)=[];
        clim_tp_w_1st_ns_d(60)=[];
        clim_tp_w_2nd_ns_d(60)=[];
        clim_tp_w_3rd_ns_d(60)=[];
        clim_tp_y_ns_d(60)=[];
        clim_tp_j_1st_ns_d(60)=[];
        clim_tp_j_2nd_ns_d(60)=[];
%         clim_tp_j_3rd_ns_d(60)=[];
        clim_tp_ye_ns_d(60)=[];
        clim_tp_jw_ns_d(60)=[];
end

temp_out_sumjin=sumjin_re_w_c{t_year-1989};
temp_out_gahwa=gahwa_re_w_c{t_year-1989};
salt_out = zeros(cycle,1); 
% salt_out = 1*ones(1,365);  %silver
sj_trans_out = dis_pre_total{t_year-1979}; %dis_pre_total

if cycle == 365 & length(sj_trans_out) == 366
    sj_trans_out(60) = [];
end

sumjin_det_n=(sumjin_tn.yp_w_tn_04.*1000/14)' - (sumjin_no3 + sumjin_nh4);
nam_det_n= (nam_tn.yp_w_tn_04.*1000/14)' - (nam_no3 + nam_nh4);

N=20; % vertical layer

r=1;  % 송정 river flow
sumjin_po4(sumjin_po4 < 0 ) = 0.01;
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
River.LDeN(:,:,r)=repmat(sumjin_det_n.*0.7,1,20);
River.SDeN(:,:,r)=repmat(sumjin_det_n.*0.3,1,20);
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
River.temp(:,:,r)=repmat(temp_out_gahwa,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(nam_nh4,1,20);
River.NO3(:,:,r)=repmat(nam_no3,1,20);
River.Oxyg(:,:,r)=repmat(nam_do,1,20);
River.chlorophyll(:,:,r)=repmat(nam_chl,1,20);
River.tPO4(:,:,r)=repmat(nam_po4,1,20);
River.LDeN(:,:,r)=repmat(nam_det_n.*0.7,1,20);
River.SDeN(:,:,r)=repmat(nam_det_n.*0.3,1,20);
River.vshape(r,1:N)=1/N;

if cycle == 365

        %--- 광양산업단지폐수종말처리장----------------------------------------------
        % 1988-05-30
        r = 3; 
        Name='gy_ind';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 68;
        River.lat(r) = 112;
        River.dir(r) = 1; % v-face
        River.trans(:,r)=-clim_dis_g_d/86400;
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(clim_tn_g_ns_d'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(clim_tn_g_ns_d'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(clim_tp_g_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        % River.DON(:,:,r)=repmat(riv1_TN_f'*(999724/1000000)*(1000/14),1,20);
        River.LDeN(:,:,r)=repmat(clim_tn_g_ns_d'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(clim_tn_g_ns_d'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;

        % %---광양중앙하수종말처리장---------------------------------------------------
        % % 2004-07-30
%         % r = 4; 
        r=4;
        Name='gy_ind2';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 29;
        River.lat(r) = 99;  
        River.dir(r) = 0; % u-face
        River.trans(:,r)=clim_dis_gc_d/86400;
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(clim_tn_gc_ns_d'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(clim_tn_gc_ns_d'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(clim_tp_gc_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        % River.DON(:,:,r)=repmat(riv1_TN_f'*(999724/1000000)*(1000/14),1,20);
        River.LDeN(:,:,r)=repmat(clim_tn_gc_ns_d'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(clim_tn_gc_ns_d'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;

        % %--- 광양하수종말처리장 ----------------------------------------------------
        % % 2002-06-11
%         % r = 5; 
        r=5;
        Name='gy_ind3';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 3; % 열
        River.lat(r) = 106;  % 행
        River.dir(r) = 1; % v-face
        River.trans(:,r)=-clim_dis_gh_d/86400;
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(clim_tn_gh_ns_d'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(clim_tn_gh_ns_d'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(clim_tp_gh_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        River.LDeN(:,:,r)=repmat(clim_tn_gh_ns_d'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(clim_tn_gh_ns_d'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;

        %--- 광영하수종말처리장-----------------------------------------------------
        % 1993-02-04
        r = 6; 
%         r = 4;
        Name='gy_ind4';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 49; % 열
        River.lat(r) = 121;  % 행
        River.dir(r) = 0; % v-face
        River.trans(:,r)=clim_dis_gh2_d/86400;
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(clim_tn_gh2_ns_d'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(clim_tn_gh2_ns_d'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(clim_tp_gh2_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        % River.DON(:,:,r)=repmat(riv4_TN_f'*(999724/1000000)*(1000/14),1,20);
        River.LDeN(:,:,r)=repmat(clim_tn_gh2_ns_d'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(clim_tn_gh2_ns_d'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;
        %--- 여수월내 -----------------------------------------
        % 1990-12-31
        r = 7; 
        if t_year <= 1999
            temp_clim_tp=clim_tp_w_1st_ns_d;
        elseif t_year >= 2000 & t_year <= 2008
            temp_clim_tp=clim_tp_w_2nd_ns_d;
        elseif  t_year >= 2009
            temp_clim_tp=clim_tp_w_3rd_ns_d;
        end     
        
        if t_year <= 1999
            temp_clim_tn=clim_tn_w_1st_ns_d;
        elseif t_year >= 2000 & t_year <= 2003
            temp_clim_tn=clim_tn_w_2nd_ns_d;
        elseif t_year >= 2004
            temp_clim_tn=clim_tn_w_3rd_ns_d;
        end 
        
%         r = 5; 
        Name='ys_ind1';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 56; % 열
        River.lat(r) = 63;  % 행
        River.dir(r) = 1; % 0:u-face, 1:v-face
        River.trans(:,r) = clim_dis_d/86400; % positive: northward flow 
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(temp_clim_tn'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(temp_clim_tn'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(temp_clim_tp'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        % River.DON(:,:,r)=repmat(riv5_TN_f'*(999724/1000000)*(1000/14),1,20);
        River.LDeN(:,:,r)=repmat(temp_clim_tn'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(temp_clim_tn'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;

        % %--- 여수율촌 -------------------------------------------
        % % 2009-11-07
        r = 8; 
        Name='ys_ind2';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 5; % 열
        River.lat(r) = 84;  % 행
        River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
        River.trans(:,r) = -clim_dis_y_d/86400; % negative: sourthward flow 
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(clim_tn_y_ns_d'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(clim_tn_y_ns_d'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(clim_tp_y_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        % River.DON(:,:,r)=repmat(riv5_TN_f'*(999724/1000000)*(1000/14),1,20);
        River.LDeN(:,:,r)=repmat(clim_tn_y_ns_d'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(clim_tn_y_ns_d'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;

        %--- 여수중흥 -------------------------------------------
        % 1997-07-12
        r = 9; 
%         r =6;
        if t_year <= 2007
            temp_clim_tp=clim_tp_w_1st_ns_d;
        elseif t_year >= 2008
            temp_clim_tp=clim_tp_w_2nd_ns_d;
        end
        
        if t_year <= 2003
            temp_clim_tn=clim_tn_w_1st_ns_d;
        elseif t_year >= 2004 & t_year <= 2008
            temp_clim_tn=clim_tn_w_2nd_ns_d;
        elseif t_year >= 2009
            temp_clim_tn=clim_tn_w_3rd_ns_d;
        end 

        Name='ys_ind3';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 24; % 열
        River.lat(r) = 56;  % 행
        River.dir(r) = 0; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
        River.trans(:,r) = clim_dis_j_d/86400; % positive: northward flow 
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(temp_clim_tn'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(temp_clim_tn'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(temp_clim_tp'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        % River.DON(:,:,r)=repmat(temp_clim_tn'*(999724/1000000)*(1000/14),1,20);
        % River.PON(:,:,r)=repmat(riv7_TN_f'*(999724/1000000)*(1000/14),1,20);
        River.LDeN(:,:,r)=repmat(temp_clim_tn'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(temp_clim_tn'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;

        % %--- 여수하수 -------------------------------------------
        % % 2005-01-01
        r = 10;
        Name='ys_in4';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 28; % 열
        River.lat(r) = 23;  % 행
        River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
        River.trans(:,r) = -clim_dis_ye_d/86400; % positive: northward flow 
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(clim_tn_ye_ns_d'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(clim_tn_ye_ns_d'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(clim_tp_ye_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        % River.DON(:,:,r)=repmat(temp_clim_tn'*(999724/1000000)*(1000/14),1,20);
        % River.PON(:,:,r)=repmat(riv7_TN_f'*(999724/1000000)*(1000/14),1,20);
        River.LDeN(:,:,r)=repmat(clim_tn_ye_ns_d'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(clim_tn_ye_ns_d'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;

        % %--- 진월하수 -------------------------------------------
        % % 2011-06-01
        r = 11; 
        Name='ys_ind5';
        River.Name(r,1:length(Name)) = Name;
        River.lon(r) = 61; % 열
        River.lat(r) = 126;  % 행
        River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
        River.trans(:,r) = -clim_dis_jw_d/86400; % positive: northward flow 
        River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
        River.salt(:,:,r)=repmat(salt_out,1,20);
        River.NH4(:,:,r)=repmat(clim_tn_jw_ns_d'*(140/1000000)*(1000/14),1,20);
        River.NO3(:,:,r)=repmat(clim_tn_jw_ns_d'*(136/1000000)*(1000/14),1,20);
        River.Oxyg(:,:,r)=repmat(sumjin_do,1,20);
        River.tPO4(:,:,r)=repmat(clim_tp_jw_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
        River.chlorophyll(:,:,r)=repmat(sumjin_chl.* 0,1,20);
        % River.DON(:,:,r)=repmat(temp_clim_tn'*(999724/1000000)*(1000/14),1,20);
        % River.PON(:,:,r)=repmat(riv7_TN_f'*(999724/1000000)*(1000/14),1,20);
        River.LDeN(:,:,r)=repmat(clim_tn_jw_ns_d'*(999724/1000000)*(1000/14).*0.7,1,20);
        River.SDeN(:,:,r)=repmat(clim_tn_jw_ns_d'*(999724/1000000)*(1000/14).*0.3,1,20);
        River.vshape(r,1:N)=1/N;
           
end

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

theVarname = 'River.vshape';
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

theVarname = 'river_LDeN';
ncout{theVarname}(:,:,:) = River.LDeN;

theVarname = 'river_SDeN';
ncout{theVarname}(:,:,:) = River.SDeN;
close(ncout)


clearvars -except t_year
cd D:\장기생태\Dynamic\06_river
end

return
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

temp=ncread('river_1997_realts_biofennel_GY_sewer_det_v8.nc','river_temp');
dis=ncread('river_1997_realts_biofennel_GY_sewer_det_v8.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1997_realts_biofennel_GY_sewer_det_v8.nc','river_temp');
dis=ncread('river_1997_realts_biofennel_GY_sewer_det_v8.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');



temp=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_Oxyg');
dis=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_Oxyg');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');


temp=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_NH4');
dis=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');plot(squeeze(temp(3,20,:)),'g');plot(squeeze(temp(4,20,:)),'k');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_NH4');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');


temp=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_NO3');
dis=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');plot(squeeze(temp(3,20,:)),'g');plot(squeeze(temp(4,20,:)),'k');


temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_NO3');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');


temp=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_chlorophyll');
dis=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_chlorophyll');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

temp=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_tPO4');
dis=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_tPO4');
dis=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp2(1,20,:)));hold on; plot(squeeze(temp2(2,20,:)),'r');

figure; hold on;plot(squeeze(temp2(1,20,:)),'k'); plot(squeeze(temp2(2,20,:)),'c');
plot(squeeze(temp(1,20,:))); plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_2001_realts_biofennel_GY_sewer.nc','river_tPO4');
dis=ncread('river_2001_realts_biofennel_GY_sewer.nc','river_transport');
figure;hold on; plot(dis(3,:),'r');
plot(dis(4,:),'g');
plot(dis(5,:),'b');
plot(dis(6,:),'m');

figure;hold on; plot(temp2(3,:),'r');
plot(temp2(4,:),'g');
plot(temp2(5,:),'b');
plot(temp2(6,:),'m');

temp2 = temp2 ./ 30.973762 %ug/L
figure;hold on; plot(temp2(3,:),'r');
plot(temp2(4,:),'g');
plot(temp2(5,:),'b');
plot(temp2(6,:),'m');

temp2=ncread('river_2001_realts_biofennel_GY_sewer.nc','river_NO3');
figure;hold on; plot(temp2(3,:),'r');
plot(temp2(4,:),'g');
plot(temp2(5,:),'b');
plot(temp2(6,:),'m');

temp2=ncread('river_2001_realts_biofennel_GY_sewer.nc','river_NH4');
figure;hold on; plot(temp2(3,:),'r');
plot(temp2(4,:),'g');
plot(temp2(5,:),'b');
plot(temp2(6,:),'m');

