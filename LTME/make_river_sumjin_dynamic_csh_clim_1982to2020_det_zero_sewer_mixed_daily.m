close all; clear all; clc
cd F:\ROMS\roms_tools\Run
start

cd D:\장기생태\Dynamic\06_river

% year_txt = {'1982to1986','1987to1996','1997to2006'};
year_txt = {'1982to1986','1987to1992','1993to2004'};

for t_year = 1:3 %before regime change
    clearvars -except t_year year_txt
grd_file='grid_gy_v11_s.nc';
% Fname='pohang_tide_inout.nc';
Fname =['river_', char(year_txt(t_year)), '_reg_bio_det_nosew_mixed_daily.nc'];
ss = 1:365;
cycle=ss(end);

grd=roms_get_grid(grd_file);
time=[1:1:cycle];
ntime=length(time);
s_rho = 20;
% river = 11;
river = 2;
create_river_NH4(grd,Fname,cycle,grd_file,s_rho,river,ntime);
ncout=netcdf(Fname,'write');
% ncout=netcdf(Fname,'clobber');

%% 송정 유량, 금강댐 방류량 load
cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
%river discharge 1980~1989, calculated from water lv. & discharge relation
%from ykang
sj = load('songjung_climate_days_bio_1982to2020_high_airT_mixed_daily_v2.mat'); % 3regime climate (1997~2006, 2007~2015, 2016~2019);
ng = load('namgang_climate_days_bio_1982to2020_high_airT_mixed_daily_v2.mat');  % 3regime climate (1997~2006, 2007~2015, 2016~2019);

% ng_dis=load('D:\장기생태\Dynamic\06_river\data\Namgang_dam\gawha_gate_discharge_and_compare_with_namgang.mat');
% namgang_climate_days_bio_to07to16_high_airT.mat
%TN
sumjin_tn=load('D:\장기생태\Dynamic\06_river\환경과학원\sumjin(songjung)_polynomial_climate_to2004_advanced(v4)_3sig_TN.mat','yp_w_tn_04');
nam_tn=load('D:\장기생태\Dynamic\06_river\환경과학원\Namgang_polynomial_climate_to2004_advanced(v4)_3sig_TN.mat','yp_w_tn_04');

%% 하수종말처리장
cd D:\장기생태\Dynamic\06_river\하수종말처리장\여수_중흥\
load('sewer_monthly_climate_fix_to06to15.mat');

cd D:\장기생태\Dynamic\06_river
%if it's not leap yr.. delete 2/29
if cycle == 365
%     sumjin_do(60)= []; 
%     sumjin_chl(60)= []; 
%     sumjin_no3(60)= []; 
%     sumjin_nh4(60)= []; 
%     sumjin_po4(60)= []; 
%     nam_do(60)= []; 
%     nam_chl(60)= []; 
%     nam_no3(60)= []; 
%     nam_nh4(60)= []; 
%     nam_po4(60)= []; 
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
        clim_tp_j_3rd_ns_d(60)=[];
        clim_tp_ye_ns_d(60)=[];
        clim_tp_jw_ns_d(60)=[];
end

if t_year == 1 
    trans_out_sumjin=sj.clim_day_trans_1;
    temp_out_sumjin=sj.clim_day_riv_temp_1;
%     trans_out_gawha=ng_dis.climate_out_gawha_1st;
    temp_out_gahwa=ng.clim_day_riv_temp_1;
        do_out_sumjin=sj.day_clim_do_1;
        do_out_gahwa=ng.day_clim_do_1;
        chl_out_sumjin=sj.day_clim_chl_1;
        chl_out_gahwa=ng.day_clim_chl_1;
        nh4_out_sumjin=sj.day_clim_nh4_1;
        nh4_out_gahwa=ng.day_clim_nh4_1;
        no3_out_sumjin=sj.day_clim_no3_1;
        no3_out_gahwa=ng.day_clim_no3_1;
        po4_out_sumjin=sj.day_clim_po4_1;
        po4_out_gahwa=ng.day_clim_po4_1;
        wout_clim_tp=clim_tp_w_1st_ns_d;
        wout_clim_tn=clim_tn_w_1st_ns_d;
        jout_clim_tp=clim_tp_j_1st_ns_d;
        jout_clim_tn=clim_tn_j_1st_ns_d;
elseif t_year == 2
    trans_out_sumjin=sj.clim_day_trans_2;
%     trans_out_gawha=ng_dis.climate_out_gawha_1st;
    temp_out_sumjin=sj.clim_day_riv_temp_2;
    temp_out_gahwa=ng.clim_day_riv_temp_2;
        do_out_sumjin=sj.day_clim_do_1;
        do_out_gahwa=ng.day_clim_do_1;
        chl_out_sumjin=sj.day_clim_chl_1;
        chl_out_gahwa=ng.day_clim_chl_1;
        nh4_out_sumjin=sj.day_clim_nh4_1;
        nh4_out_gahwa=ng.day_clim_nh4_1;
        no3_out_sumjin=sj.day_clim_no3_1;
        no3_out_gahwa=ng.day_clim_no3_1;
        po4_out_sumjin=sj.day_clim_po4_1;
        po4_out_gahwa=ng.day_clim_po4_1;
        wout_clim_tp=clim_tp_w_1st_ns_d;
        wout_clim_tn=clim_tn_w_1st_ns_d;
        jout_clim_tp=clim_tp_j_1st_ns_d;
        jout_clim_tn=clim_tn_j_1st_ns_d;
elseif t_year == 3
    trans_out_sumjin=sj.clim_day_trans_3;
%     trans_out_gawha=ng_dis.climate_out_gawha_2nd;
    temp_out_sumjin=sj.clim_day_riv_temp_3;
    temp_out_gahwa=ng.clim_day_riv_temp_3;
       do_out_sumjin=sj.day_clim_do_1;
        do_out_gahwa=ng.day_clim_do_1;
        chl_out_sumjin=sj.day_clim_chl_1;
        chl_out_gahwa=ng.day_clim_chl_1;
        nh4_out_sumjin=sj.day_clim_nh4_1;
        nh4_out_gahwa=ng.day_clim_nh4_1;
        no3_out_sumjin=sj.day_clim_no3_1;
        no3_out_gahwa=ng.day_clim_no3_1;
        po4_out_sumjin=sj.day_clim_po4_1;
        po4_out_gahwa=ng.day_clim_po4_1;
        wout_clim_tp=clim_tp_w_1st_ns_d;
        wout_clim_tn=clim_tn_w_1st_ns_d;
        jout_clim_tp=clim_tp_j_1st_ns_d;
        jout_clim_tn=clim_tn_j_1st_ns_d;
end

salt_out = zeros(cycle,1); 
% salt_out = 1*ones(1,365);  %silver

sumjin_det_n=(sumjin_tn.yp_w_tn_04.*1000/14)' - (no3_out_sumjin + nh4_out_sumjin);
nam_det_n= (nam_tn.yp_w_tn_04.*1000/14)' - (no3_out_gahwa + nh4_out_gahwa);

N=20; % vertical layer

%% gawha discharge factor ( songjung dis. * factor = gawha dis. )
switch t_year
    case {1,2}
        factor = 0.7;
    case {3}
        factor = 1.0;
end
    

r=1;  % 송정 river flow
Name='outflow_songjung';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 12; % 열
River.lat(r) = 174;  % 향
River.dir(r) = 0; % 0:u-face, 1:v-face
% River.flag(r) = 3; %temp, salt
River.trans(:,r)=trans_out_sumjin;
River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
%--- NPZD -----------------------------
River.NH4(:,:,r)=repmat(nh4_out_sumjin,1,20);
River.NO3(:,:,r)=repmat(no3_out_sumjin,1,20);
River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
River.Chlo(:,:,r)=repmat(chl_out_sumjin,1,20);
River.tPO4(:,:,r)=repmat(po4_out_sumjin,1,20);
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
River.trans(:,r)=(trans_out_sumjin .* factor);
% River.trans(:,r)=ones(ntime,1)*trans_out;
River.temp(:,:,r)=repmat(temp_out_gahwa,1,20);
River.salt(:,:,r)=repmat(salt_out,1,20);
River.NH4(:,:,r)=repmat(nh4_out_gahwa,1,20);
River.NO3(:,:,r)=repmat(no3_out_gahwa,1,20);
River.Oxyg(:,:,r)=repmat(do_out_gahwa,1,20);
River.Chlo(:,:,r)=repmat(chl_out_gahwa,1,20);
River.tPO4(:,:,r)=repmat(po4_out_gahwa,1,20);
River.LDeN(:,:,r)=repmat(nam_det_n.*0.7,1,20);
River.SDeN(:,:,r)=repmat(nam_det_n.*0.3,1,20);
River.vshape(r,1:N)=1/N;

% if cycle == 365
% 
%         %--- 광양산업단지폐수종말처리장----------------------------------------------
%         % 1988-05-30
%         r = 3; 
%         Name='gy_ind';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 68;
%         River.lat(r) = 112;
%         River.dir(r) = 1; % v-face
%         River.trans(:,r)=-clim_dis_g_d/86400;
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(clim_tn_g_ns_d'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(clim_tn_g_ns_d'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(clim_tp_g_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         % River.DON(:,:,r)=repmat(riv1_TN_f'*(999724/2/1000000)*(1000/14),1,20);
%         River.LDeN(:,:,r)=repmat(clim_tn_g_ns_d'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(clim_tn_g_ns_d'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
% 
%         % %---광양중앙하수종말처리장---------------------------------------------------
%         % % 2004-07-30
% %         % r = 4; 
%         r=4;
%         Name='gy_ind2';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 29;
%         River.lat(r) = 99;  
%         River.dir(r) = 1; % v-face
%         River.trans(:,r)=-clim_dis_gc_d/86400;
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(clim_tn_gc_ns_d'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(clim_tn_gc_ns_d'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(clim_tp_gc_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         % River.DON(:,:,r)=repmat(riv1_TN_f'*(999724/2/1000000)*(1000/14),1,20);
%         River.LDeN(:,:,r)=repmat(clim_tn_gc_ns_d'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(clim_tn_gc_ns_d'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
% 
%         % %--- 광양하수종말처리장 ----------------------------------------------------
%         % % 2002-06-11
% %         % r = 5; 
%         r=5;
%         Name='gy_ind3';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 3; % 열
%         River.lat(r) = 106;  % 행
%         River.dir(r) = 1; % v-face
%         River.trans(:,r)=-clim_dis_gh_d/86400;
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(clim_tn_gh_ns_d'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(clim_tn_gh_ns_d'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(clim_tp_gh_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         River.LDeN(:,:,r)=repmat(clim_tn_gh_ns_d'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(clim_tn_gh_ns_d'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
% 
%         %--- 광영하수종말처리장-----------------------------------------------------
%         % 1993-02-04
%         r = 6; 
% %         r = 4;
%         Name='gy_ind4';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 49; % 열
%         River.lat(r) = 121;  % 행
%         River.dir(r) = 1; % v-face
%         River.trans(:,r)=-clim_dis_gh2_d/86400;
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(clim_tn_gh2_ns_d'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(clim_tn_gh2_ns_d'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(clim_tp_gh2_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         % River.DON(:,:,r)=repmat(riv4_TN_f'*(999724/2/1000000)*(1000/14),1,20);
%         River.LDeN(:,:,r)=repmat(clim_tn_gh2_ns_d'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(clim_tn_gh2_ns_d'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
%         %--- 여수월내 -----------------------------------------
%         % 1990-12-31
%         r = 7;     
% %         r = 5; 
%         Name='ys_ind1';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 56; % 열
%         River.lat(r) = 63;  % 행
%         River.dir(r) = 1; % 0:u-face, 1:v-face
%         River.trans(:,r) = clim_dis_d/86400; % positive: northward flow 
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(wout_clim_tn'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(wout_clim_tn'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(wout_clim_tp'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         % River.DON(:,:,r)=repmat(riv5_TN_f'*(999724/2/1000000)*(1000/14),1,20);
%         River.LDeN(:,:,r)=repmat(wout_clim_tn'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(wout_clim_tn'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
% 
%         % %--- 여수율촌 -------------------------------------------
%         % % 2009-11-07
%         r = 8; 
%         Name='ys_ind2';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 5; % 열
%         River.lat(r) = 84;  % 행
%         River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
%         River.trans(:,r) = -clim_dis_y_d/86400; % negative: sourthward flow 
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(clim_tn_y_ns_d'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(clim_tn_y_ns_d'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(clim_tp_y_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         % River.DON(:,:,r)=repmat(riv5_TN_f'*(999724/2/1000000)*(1000/14),1,20);
%         River.LDeN(:,:,r)=repmat(clim_tn_y_ns_d'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(clim_tn_y_ns_d'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
% 
%         %--- 여수중흥 -------------------------------------------
%         % 1997-07-12
%         r = 9; 
% %         r =6;
%         Name='ys_ind3';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 24; % 열
%         River.lat(r) = 56;  % 행
%         River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
%         River.trans(:,r) = clim_dis_j_d/86400; % positive: northward flow 
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(jout_clim_tn'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(jout_clim_tn'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(jout_clim_tp'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         % River.DON(:,:,r)=repmat(temp_clim_tn'*(999724/2/1000000)*(1000/14),1,20);
%         % River.PON(:,:,r)=repmat(riv7_TN_f'*(999724/2/1000000)*(1000/14),1,20);
%         River.LDeN(:,:,r)=repmat(jout_clim_tn'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(jout_clim_tn'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
% 
%         % %--- 여수하수 -------------------------------------------
%         % % 2005-01-01
%         r = 10;
%         Name='ys_in4';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 28; % 열
%         River.lat(r) = 23;  % 행
%         River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
%         River.trans(:,r) = -clim_dis_ye_d/86400; % positive: northward flow 
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(clim_tn_ye_ns_d'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(clim_tn_ye_ns_d'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(clim_tp_ye_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         % River.DON(:,:,r)=repmat(temp_clim_tn'*(999724/2/1000000)*(1000/14),1,20);
%         % River.PON(:,:,r)=repmat(riv7_TN_f'*(999724/2/1000000)*(1000/14),1,20);
%         River.LDeN(:,:,r)=repmat(clim_tn_ye_ns_d'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(clim_tn_ye_ns_d'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
% 
%         % %--- 진월하수 -------------------------------------------
%         % % 2011-06-01
%         r = 11; 
%         Name='ys_ind5';
%         River.Name(r,1:length(Name)) = Name;
%         River.lon(r) = 61; % 열
%         River.lat(r) = 126;  % 행
%         River.dir(r) = 1; % 남으로 흐름 -> River.lon(r) 값을 edtimark 값에 +1
%         River.trans(:,r) = -clim_dis_jw_d/86400; % positive: northward flow 
%         River.temp(:,:,r)=repmat(temp_out_sumjin,1,20);
%         River.salt(:,:,r)=repmat(salt_out,1,20);
%         River.NH4(:,:,r)=repmat(clim_tn_jw_ns_d'*(140/1000000)*(1000/14),1,20);
%         River.NO3(:,:,r)=repmat(clim_tn_jw_ns_d'*(136/1000000)*(1000/14),1,20);
%         River.Oxyg(:,:,r)=repmat(do_out_sumjin,1,20);
%         River.tPO4(:,:,r)=repmat(clim_tp_jw_ns_d'.*0.8 .*(1000/30.973762)  ,1,20);
%         River.Chlo(:,:,r)=repmat(chl_out_sumjin.* 0,1,20);
%         % River.DON(:,:,r)=repmat(temp_clim_tn'*(999724/2/1000000)*(1000/14),1,20);
%         % River.PON(:,:,r)=repmat(riv7_TN_f'*(999724/2/1000000)*(1000/14),1,20);
%         River.LDeN(:,:,r)=repmat(clim_tn_jw_ns_d'*(999724/2/1000000)*(1000/14).*0.7,1,20);
%         River.SDeN(:,:,r)=repmat(clim_tn_jw_ns_d'*(999724/2/1000000)*(1000/14).*0.3,1,20);
%         River.vshape(r,1:N)=1/N;
%            
% end

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
 
theVarname = 'river_Chlo';
ncout{theVarname}(:,:,:) = River.Chlo;

theVarname = 'river_tPO4';
ncout{theVarname}(:,:,:) = River.tPO4;

theVarname = 'river_LDeN';
ncout{theVarname}(:,:,:) = River.LDeN;

theVarname = 'river_SDeN';
ncout{theVarname}(:,:,:) = River.SDeN;
close(ncout)


clearvars -except t_year year_txt
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

clearvars -except Fname
Fname =['river_', '1982to1986', '_reg_bio_det_nosew_mixed_daily.nc'];
ncload(Fname);

figure('position',[100 100 800 900]);
subplot(611)
plot(river_transport); xlim([1 365])
subplot(612)
plot(squeeze(river_temp(:,20,:))); xlim([1 365])
subplot(613)
plot(squeeze(river_salt(:,20,:))); xlim([1 365])
subplot(614)
plot(squeeze(river_NO3(:,20,:))); xlim([1 365])
subplot(615)
plot(squeeze(river_NH4(:,20,:))); xlim([1 365])
subplot(616)
plot(squeeze(river_Chlo(:,20,:))); xlim([1 365])

clearvars -except Fname
Fname =['river_', '1987to1992', '_reg_bio_det_nosew_mixed_daily.nc'];
ncload(Fname);

figure('position',[100 100 800 900]);
subplot(611)
plot(river_transport); xlim([1 365])
subplot(612)
plot(squeeze(river_temp(:,20,:))); xlim([1 365])
subplot(613)
plot(squeeze(river_salt(:,20,:))); xlim([1 365])
subplot(614)
plot(squeeze(river_NO3(:,20,:))); xlim([1 365])
subplot(615)
plot(squeeze(river_NH4(:,20,:))); xlim([1 365])
subplot(616)
plot(squeeze(river_Chlo(:,20,:))); xlim([1 365])

clearvars -except Fname
Fname =['river_', '1993to2004', '_reg_bio_det_nosew_mixed_daily.nc'];
ncload(Fname);
figure('position',[100 100 800 900]);
subplot(611)
plot(river_transport); xlim([1 365])
subplot(612)
plot(squeeze(river_temp(:,20,:))); xlim([1 365])
subplot(613)
plot(squeeze(river_salt(:,20,:))); xlim([1 365])
subplot(614)
plot(squeeze(river_NO3(:,20,:))); xlim([1 365])
subplot(615)
plot(squeeze(river_NH4(:,20,:))); xlim([1 365])
subplot(616)
plot(squeeze(river_Chlo(:,20,:))); xlim([1 365])


%% check
close all; clear; clc;
cd D:\장기생태\Dynamic\06_river
po4_v7=ncread('river_1997_realts_biofennel_GY_sewer_det_v7.nc','river_tPO4');
dis_v7=ncread('river_1997_realts_biofennel_GY_sewer_det_v7.nc','river_transport');
po4_v8=ncread('river_1997_realts_biofennel_GY_sewer_det_v8.nc','river_tPO4');
dis_v8=ncread('river_1997_realts_biofennel_GY_sewer_det_v8.nc','river_transport');
po4_v7_07=ncread('river_2007_realts_biofennel_GY_sewer_det_v7.nc','river_tPO4');
dis_v7_07=ncread('river_2007_realts_biofennel_GY_sewer_det_v7.nc','river_transport');
po4_v8_07=ncread('river_2007_realts_biofennel_GY_sewer_det_v8.nc','river_tPO4');
dis_v8_07=ncread('river_2007_realts_biofennel_GY_sewer_det_v8.nc','river_transport');
cd C:\Users\user\Desktop\장기생태_하수_여수월내
po4_v9=ncread('river_1997_realts_biofennel_GY_sewer_det_v9.nc','river_tPO4');
dis_v9=ncread('river_1997_realts_biofennel_GY_sewer_det_v9.nc','river_transport');
po4_v9_07=ncread('river_2007_realts_biofennel_GY_sewer_det_v9.nc','river_tPO4');
dis_v9_07=ncread('river_2007_realts_biofennel_GY_sewer_det_v9.nc','river_transport');
po4_v9_2=ncread('river_1997_realts_biofennel_GY_sewer_det_v9_2.nc','river_tPO4');
dis_v9_2=ncread('river_1997_realts_biofennel_GY_sewer_det_v9_2.nc','river_transport');

figure; hold on; plot(squeeze(po4_v7(5,20,:)),'b'); plot(squeeze(po4_v8(7,20,:)),'r'); plot(squeeze(po4_v9(7,20,:)),'k'); % wallne
figure; hold on; plot(squeeze(po4_v7(6,20,:)),'b'); plot(squeeze(po4_v8(9,20,:)),'r'); plot(squeeze(po4_v9(9,20,:)),'k'); % jungheung

figure; hold on; plot(squeeze(po4_v7(5,20,:)),'b'); plot(squeeze(po4_v8(7,20,:)),'r'); plot(squeeze(po4_v9(7,20,:)),'k'); % wallne
plot(squeeze(po4_v7_07(5,20,:)),'m'); plot(squeeze(po4_v8_07(7,20,:)),'g'); plot(squeeze(po4_v9_07(7,20,:)),'c'); % wallne

figure; hold on; plot(squeeze(po4_v7(6,20,:)),'b'); plot(squeeze(po4_v8(9,20,:)),'r'); plot(squeeze(po4_v9(9,20,:)),'k'); % jungheung
plot(squeeze(po4_v7_07(6,20,:)),'m'); plot(squeeze(po4_v8_07(9,20,:)),'g'); plot(squeeze(po4_v9_07(9,20,:)),'c'); % jungheung

figure; hold on; plot(dis_v7(5,:),'b'); plot(dis_v8(7,:),'r'); plot(dis_v9(7,:),'k'); % wallne
figure; hold on; plot(dis_v7(6,:),'b'); plot(dis_v8(9,:),'r'); plot(dis_v9(9,:),'k'); % jungheung

figure; hold on; plot(squeeze(dis_v7(5,:)),'b'); plot(squeeze(dis_v8(7,:)),'r'); plot(squeeze(dis_v9(7,:)),'k'); % wallne
plot(squeeze(dis_v7_07(5,:)),'m'); plot(squeeze(dis_v8_07(7,:)),'g'); plot(squeeze(dis_v9_07(7,:)),'c'); % wallne

figure; hold on; plot(squeeze(dis_v7(6,:)),'b'); plot(squeeze(dis_v8(9,:)),'r'); plot(squeeze(dis_v9(9,:)),'k'); % jungheung
plot(squeeze(dis_v7_07(6,:)),'m'); plot(squeeze(dis_v8_07(9,:)),'g'); plot(squeeze(dis_v9_07(9,:)),'c'); % jungheung



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


temp=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_Chlo');
dis=ncread('river_1996_realts_biofennel_GY_sewer_det.nc','river_transport');
figure;plot(dis(1,:));hold on; plot(dis(2,:),'r');
figure;plot(squeeze(temp(1,20,:)));hold on; plot(squeeze(temp(2,20,:)),'r');

temp2=ncread('river_1999_realts_biofennel_Gwangyang.nc','river_Chlo');
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

