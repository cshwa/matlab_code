close all; clc; clear;
% clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_9points.mat'); %default
clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_9points_no_regime_cut.mat'); %no cut regime
% clim1=load('D:\장기생태\Dynamic\06_river\koem_climate_to06to15_v5_sig_gy_8points_no_regime_cut.mat'); %no cut regime 8 points

yoonja=load('yoonjakangs_koem_data_monthly_v3_16points_2021.mat','mon_d','time');
sp_9p_st = [1:6, 14:16];
mon_d=yoonja.mon_d(sp_9p_st,:);
time_hr=yoonja.time(sp_9p_st,:);

plt_time = zeros(1,12);
plt_time(1,[2,5,8,11]) = [31, 120, 212, 304]; 

k=0
for i = 1997:2021
    for j = 1:12
        k=k+1;
        for l = 1:size(mon_d,1)
            if isnan(mon_d(l,k))==0
               ref_date{l,k} = plt_time(1,j) + mon_d(l,k);
               if leapyear(i)==1 & j > 2 & ref_date{l,k} >= 59
                   disp(['leapyear ',num2str(i),'-',num2str(j)])
                   ref_date{l,k} = ref_date{l,k} + leapyear(i);
               end
            else
                ref_date{l,k}=NaN;
            end
        end
    end
end

k=0
for i = 1997:2021
    for j = 1:12
        k=k+1;
        ref_mon{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

cd J:\장기생태_2021\Dynamic\result\3regime\daily\2022_run
% regime1_kdc
% regime1_kdc_zerosw=load('kd_double_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('merg_part_t4_2reg_mix.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_mix_2.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_3.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_t4_1_mix_bio_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_1_mix_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_cla.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('2reg_kd_t4_cir.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime1_kdc_zerosw=load(['J:\장기생태_2021\Dynamic\result\3regime\','2reg_cir_2yr.mat'],'obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('kd_t4w_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_kdc_zerosw=load('part_kd_t6_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime1_riv_chl=load('1regime_part_res_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc=
% regime2_kdc_zerosw=load('kd_double_3reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('kd_t4_2reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('2reg_kd_t4_1_mix_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
%% 2022
% regime2_kdc_zerosw=load('2reg_t1.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw1=load('2reg_t2.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw2=load('2reg_t3.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw3=load('2reg_t4.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw1=load('2reg_t4.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw2=load('2reg_t5.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw3=load('2reg_t6.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_kdc_zerosw=load('2reg_t7.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_kdc_zerosw1=load('2reg_t8.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_kdc_zerosw2=load('2reg_t9.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
regime2_kdc_zerosw3=load('2reg_t10.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');



%% 2022%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END

% regime2_kdc_zerosw=load('3reg_t4_1_mix_bio_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% [num2str(refyear),'reg_t4_1_mix_bio_part'];
% regime2_kdc_zerosw=load('3reg_kd_t4_1_mix_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('2reg_kd_t4_1_mix_part.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('kd_t4w_3reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_kdc_zerosw=load('part_kd_t6_3reg.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_noriv_chl=load('2regime_part_res_2regime_noriv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');
% regime2_riv_chl=load('2regime_part_res_2regime_riv_chl.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');

% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\kd_double\2&3regime\4sigma'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\4sigma'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4w_case\2&3regime\4sigma'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\'

sigma = 2;

% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\kd_double\2&3regime\nocutregime\3sig'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\0815_mix_1'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\0815_mix_2'
% plot_dir ='J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4w_case\2&3regime\nocutregime\1sig'
% plot_dir ='J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t6_case\2&3regime\nocutregime\3sig'
%% 2022
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\220702'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\220702\test4to6'
% plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\220702\test6to8'
plot_dir = 'J:\장기생태_2021\Dynamic\result\3regime\kdc_nosew\t4_case\2&3regime\nocutregime\220702\test7to10'

P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;
DO_CONV=0.7.*44.661; %%% mg/L to mM/m^3;

xlim_set = [1 365];

%% plot obs 1st and 2nd at the same time which is time axis
plt_filter_1{1} = 32:45;
plt_filter_2{1} = 46:59;
plt_filter_1{2} = 121:135.5;
plt_filter_2{2} = 136:151;
plt_filter_1{3} = 213:227.5;
plt_filter_2{3} = 228:243;
plt_filter_1{4} = 305:319;
plt_filter_2{4} = 320:334;


if size(regime1_kdc_zerosw.sp_gy_temp_mid,3) <= size(regime2_kdc_zerosw.sp_gy_temp_mid,3)
   size_tt = size(regime1_kdc_zerosw.sp_gy_temp_mid,3);
else
   size_tt = size(regime2_kdc_zerosw.sp_gy_temp_mid,3);
end


clearvars temp_d sp_gy_*_mm
for i = 1:9 %st
    for j = 1: size_tt %time    
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_chl_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_chl_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_do_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_do_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_no3_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_no3_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_nh4_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_nh4_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_zoo_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_zoo_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime1_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime1_kdc_zerosw.sp_gy_phy_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_phy_mm(i,:,j)=temp_d;
        
        
        % t1
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_chl_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_chl_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_do_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_do_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_no3_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_no3_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_nh4_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_nh4_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_zoo_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_zoo_2reg_mm(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw.sp_gy_phy_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_phy_2reg_mm(i,:,j)=temp_d;
        
        % t2
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw1.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw1.sp_gy_chl_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_chl_2reg_mm_1(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw1.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw1.sp_gy_do_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_do_2reg_mm_1(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw1.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw1.sp_gy_no3_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_no3_2reg_mm_1(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw1.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw1.sp_gy_nh4_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_nh4_2reg_mm_1(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw1.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw1.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_2reg_mm_1(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw1.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw1.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_2reg_mm_1(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw1.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw1.sp_gy_zoo_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_zoo_2reg_mm_1(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw1.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw1.sp_gy_phy_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_phy_2reg_mm_1(i,:,j)=temp_d;
        
        
        % t3
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw2.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw2.sp_gy_chl_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_chl_2reg_mm_2(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw2.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw2.sp_gy_do_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_do_2reg_mm_2(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw2.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw2.sp_gy_no3_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_no3_2reg_mm_2(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw2.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw2.sp_gy_nh4_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_nh4_2reg_mm_2(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw2.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw2.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_2reg_mm_2(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw2.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw2.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_2reg_mm_2(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw2.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw2.sp_gy_zoo_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_zoo_2reg_mm_2(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw2.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw2.sp_gy_phy_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_phy_2reg_mm_2(i,:,j)=temp_d;
        
        
        % t4
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw3.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw3.sp_gy_chl_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_chl_2reg_mm_3(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw3.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw3.sp_gy_do_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_do_2reg_mm_3(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw3.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw3.sp_gy_no3_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_no3_2reg_mm_3(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw3.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw3.sp_gy_nh4_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_nh4_2reg_mm_3(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw3.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw3.sp_gy_temp_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_temp_2reg_mm_3(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw3.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw3.sp_gy_salt_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_salt_2reg_mm_3(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw3.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw3.sp_gy_zoo_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_zoo_2reg_mm_3(i,:,j)=temp_d;
        clearvars temp_d
        temp_d=interp1(squeeze(regime2_kdc_zerosw3.sp_gy_z(i,:,j)),squeeze(regime2_kdc_zerosw3.sp_gy_phy_mid(i,:,j)), [0.1,0.5,1:30].*-1);
        sp_gy_phy_2reg_mm_3(i,:,j)=temp_d;
    end
end


                
% %% make same time (dd) with OBS
clearvars mod_*_pick mod_*_pick_2reg
                mod_chl_pick=NaN(9,32,size(ref_date,2));
                mod_do_pick=NaN(9,32,size(ref_date,2));
                mod_no3_pick=NaN(9,32,size(ref_date,2));
                mod_nh4_pick=NaN(9,32,size(ref_date,2));
                mod_temp_pick=NaN(9,32,size(ref_date,2));
                mod_salt_pick=NaN(9,32,size(ref_date,2));
                mod_zoo_pick=NaN(9,32,size(ref_date,2));
                mod_phy_pick=NaN(9,32,size(ref_date,2));
                mod_chl_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_do_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_no3_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_nh4_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_temp_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_salt_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_zoo_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_phy_pick_2reg=NaN(9,32,size(ref_date,2));
                mod_chl_pick_2reg_1=NaN(9,32,size(ref_date,2));
                mod_do_pick_2reg_1=NaN(9,32,size(ref_date,2));
                mod_no3_pick_2reg_1=NaN(9,32,size(ref_date,2));
                mod_nh4_pick_2reg_1=NaN(9,32,size(ref_date,2));
                mod_temp_pick_2reg_1=NaN(9,32,size(ref_date,2));
                mod_salt_pick_2reg_1=NaN(9,32,size(ref_date,2));
                mod_zoo_pick_2reg_1=NaN(9,32,size(ref_date,2));
                mod_phy_pick_2reg_1=NaN(9,32,size(ref_date,2));
                mod_chl_pick_2reg_2=NaN(9,32,size(ref_date,2));
                mod_do_pick_2reg_2=NaN(9,32,size(ref_date,2));
                mod_no3_pick_2reg_2=NaN(9,32,size(ref_date,2));
                mod_nh4_pick_2reg_2=NaN(9,32,size(ref_date,2));
                mod_temp_pick_2reg_2=NaN(9,32,size(ref_date,2));
                mod_salt_pick_2reg_2=NaN(9,32,size(ref_date,2));
                mod_zoo_pick_2reg_2=NaN(9,32,size(ref_date,2));
                mod_phy_pick_2reg_2=NaN(9,32,size(ref_date,2));
                mod_chl_pick_2reg_3=NaN(9,32,size(ref_date,2));
                mod_do_pick_2reg_3=NaN(9,32,size(ref_date,2));
                mod_no3_pick_2reg_3=NaN(9,32,size(ref_date,2));
                mod_nh4_pick_2reg_3=NaN(9,32,size(ref_date,2));
                mod_temp_pick_2reg_3=NaN(9,32,size(ref_date,2));
                mod_salt_pick_2reg_3=NaN(9,32,size(ref_date,2));
                mod_zoo_pick_2reg_3=NaN(9,32,size(ref_date,2));
                mod_phy_pick_2reg_3=NaN(9,32,size(ref_date,2));


for i = 1:9 %st
    for j = 1:size(ref_date,2) 
        if isnan(ref_date{i,j}) == 0 
            if j <= 228 %time ~ 2015.12
                mod_chl_pick(i,:,j)=sp_gy_chl_mm(i,:,ref_date{i,j});
                mod_do_pick(i,:,j)=sp_gy_do_mm(i,:,ref_date{i,j});
                mod_no3_pick(i,:,j)=sp_gy_no3_mm(i,:,ref_date{i,j});
                mod_nh4_pick(i,:,j)=sp_gy_nh4_mm(i,:,ref_date{i,j});
                mod_temp_pick(i,:,j)=sp_gy_temp_mm(i,:,ref_date{i,j});
                mod_salt_pick(i,:,j)=sp_gy_salt_mm(i,:,ref_date{i,j});
                mod_zoo_pick(i,:,j)=sp_gy_zoo_mm(i,:,ref_date{i,j});
                mod_phy_pick(i,:,j)=sp_gy_phy_mm(i,:,ref_date{i,j});
                 % t1
                mod_chl_pick_2reg(i,:,j)=sp_gy_chl_2reg_mm(i,:,ref_date{i,j});
                mod_do_pick_2reg(i,:,j)=sp_gy_do_2reg_mm(i,:,ref_date{i,j});
                mod_no3_pick_2reg(i,:,j)=sp_gy_no3_2reg_mm(i,:,ref_date{i,j});
                mod_nh4_pick_2reg(i,:,j)=sp_gy_nh4_2reg_mm(i,:,ref_date{i,j});
                mod_temp_pick_2reg(i,:,j)=sp_gy_temp_2reg_mm(i,:,ref_date{i,j});
                mod_salt_pick_2reg(i,:,j)=sp_gy_salt_2reg_mm(i,:,ref_date{i,j});
                mod_zoo_pick_2reg(i,:,j)=sp_gy_zoo_2reg_mm(i,:,ref_date{i,j});
                mod_phy_pick_2reg(i,:,j)=sp_gy_phy_2reg_mm(i,:,ref_date{i,j});
                 % t2
                mod_chl_pick_2reg_1(i,:,j)=sp_gy_chl_2reg_mm_1(i,:,ref_date{i,j});
                mod_do_pick_2reg_1(i,:,j)=sp_gy_do_2reg_mm_1(i,:,ref_date{i,j});
                mod_no3_pick_2reg_1(i,:,j)=sp_gy_no3_2reg_mm_1(i,:,ref_date{i,j});
                mod_nh4_pick_2reg_1(i,:,j)=sp_gy_nh4_2reg_mm_1(i,:,ref_date{i,j});
                mod_temp_pick_2reg_1(i,:,j)=sp_gy_temp_2reg_mm_1(i,:,ref_date{i,j});
                mod_salt_pick_2reg_1(i,:,j)=sp_gy_salt_2reg_mm_1(i,:,ref_date{i,j});
                mod_zoo_pick_2reg_1(i,:,j)=sp_gy_zoo_2reg_mm_1(i,:,ref_date{i,j});
                mod_phy_pick_2reg_1(i,:,j)=sp_gy_phy_2reg_mm_1(i,:,ref_date{i,j});
                % t3
                mod_chl_pick_2reg_2(i,:,j)=sp_gy_chl_2reg_mm_2(i,:,ref_date{i,j});
                mod_do_pick_2reg_2(i,:,j)=sp_gy_do_2reg_mm_2(i,:,ref_date{i,j});
                mod_no3_pick_2reg_2(i,:,j)=sp_gy_no3_2reg_mm_2(i,:,ref_date{i,j});
                mod_nh4_pick_2reg_2(i,:,j)=sp_gy_nh4_2reg_mm_2(i,:,ref_date{i,j});
                mod_temp_pick_2reg_2(i,:,j)=sp_gy_temp_2reg_mm_2(i,:,ref_date{i,j});
                mod_salt_pick_2reg_2(i,:,j)=sp_gy_salt_2reg_mm_2(i,:,ref_date{i,j});
                mod_zoo_pick_2reg_2(i,:,j)=sp_gy_zoo_2reg_mm_2(i,:,ref_date{i,j});
                mod_phy_pick_2reg_2(i,:,j)=sp_gy_phy_2reg_mm_2(i,:,ref_date{i,j});
                % t4
                mod_chl_pick_2reg_3(i,:,j)=sp_gy_chl_2reg_mm_3(i,:,ref_date{i,j});
                mod_do_pick_2reg_3(i,:,j)=sp_gy_do_2reg_mm_3(i,:,ref_date{i,j});
                mod_no3_pick_2reg_3(i,:,j)=sp_gy_no3_2reg_mm_3(i,:,ref_date{i,j});
                mod_nh4_pick_2reg_3(i,:,j)=sp_gy_nh4_2reg_mm_3(i,:,ref_date{i,j});
                mod_temp_pick_2reg_3(i,:,j)=sp_gy_temp_2reg_mm_3(i,:,ref_date{i,j});
                mod_salt_pick_2reg_3(i,:,j)=sp_gy_salt_2reg_mm_3(i,:,ref_date{i,j});
                mod_zoo_pick_2reg_3(i,:,j)=sp_gy_zoo_2reg_mm_3(i,:,ref_date{i,j});
                mod_phy_pick_2reg_3(i,:,j)=sp_gy_phy_2reg_mm_3(i,:,ref_date{i,j});
                
                
            else % 2016.01 ~ 2019.12(2020.12) : ~276 (~288)
               disp('no')
            end
        else
             if j <= 228 %time ~ 2015.12
                mod_chl_pick(i,:,j)=NaN;
                mod_do_pick(i,:,j)=NaN;
                mod_no3_pick(i,:,j)=NaN;
                mod_nh4_pick(i,:,j)=NaN;
                mod_temp_pick(i,:,j)=NaN;
                mod_salt_pick(i,:,j)=NaN;
                mod_zoo_pick(i,:,j)=NaN;
                mod_phy_pick(i,:,j)=NaN;
%                 t1
                 mod_chl_pick_2reg(i,:,j)=NaN;
                mod_do_pick_2reg(i,:,j)=NaN;
                mod_no3_pick_2reg(i,:,j)=NaN;
                mod_nh4_pick_2reg(i,:,j)=NaN;
                mod_temp_pick_2reg(i,:,j)=NaN;
                mod_salt_pick_2reg(i,:,j)=NaN;
                mod_zoo_pick_2reg(i,:,j)=NaN;
                mod_phy_pick_2reg(i,:,j)=NaN;
                %                 t2
                 mod_chl_pick_2reg_1(i,:,j)=NaN;
                mod_do_pick_2reg_1(i,:,j)=NaN;
                mod_no3_pick_2reg_1(i,:,j)=NaN;
                mod_nh4_pick_2reg_1(i,:,j)=NaN;
                mod_temp_pick_2reg_1(i,:,j)=NaN;
                mod_salt_pick_2reg_1(i,:,j)=NaN;
                mod_zoo_pick_2reg_1(i,:,j)=NaN;
                mod_phy_pick_2reg_1(i,:,j)=NaN;
                %                 t3
                 mod_chl_pick_2reg_2(i,:,j)=NaN;
                mod_do_pick_2reg_2(i,:,j)=NaN;
                mod_no3_pick_2reg_2(i,:,j)=NaN;
                mod_nh4_pick_2reg_2(i,:,j)=NaN;
                mod_temp_pick_2reg_2(i,:,j)=NaN;
                mod_salt_pick_2reg_2(i,:,j)=NaN;
                mod_zoo_pick_2reg_2(i,:,j)=NaN;
                mod_phy_pick_2reg_2(i,:,j)=NaN;
                
                %                 t4
                 mod_chl_pick_2reg_3(i,:,j)=NaN;
                mod_do_pick_2reg_3(i,:,j)=NaN;
                mod_no3_pick_2reg_3(i,:,j)=NaN;
                mod_nh4_pick_2reg_3(i,:,j)=NaN;
                mod_temp_pick_2reg_3(i,:,j)=NaN;
                mod_salt_pick_2reg_3(i,:,j)=NaN;
                mod_zoo_pick_2reg_3(i,:,j)=NaN;
                mod_phy_pick_2reg_3(i,:,j)=NaN;
                
            else % 2016.01 ~ 2021.12
                % 2 regime
               
            end
        end
    end
end
cd(plot_dir)



% plot(squeeze(mean(sp_gy_chl_mm(:,3,:),1)))
% squeeze(regime1_kdc_zerosw.sp_gy_z(:,1,1)) % 9st's depth

% clearvars clim2*
% clim2=clim1.obm_gy_po4_ft_2(:,sigma);
% clim2_std = clim1.obs_std_gy_po4_ft_2(:,sigma);
% 
% clim2_2=clim1.obm_gy_po4_ft_3(:,sigma);
% clim2_2_std = clim1.obs_std_gy_po4_ft_3(:,sigma);

% fig = figure; hold on;   
% %             plot(nanmean(mer_po4_1,1),'g','linew',2);
% %             clearvars *_bound_plt nonan_data_plt discon_p temp_d
% %             temp_d=nanmean(mer_po4_1,1);
% %             nonan_data_plt = find(isnan(temp_d)==0);
% %             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
% %             upper_bound_plt = temp_d + regm_std_po4_1;
% %             lower_bound_plt = temp_d - regm_std_po4_1;
% %             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
% 
%             clearvars *_bound_plt nonan_data_plt discon_p
%             nonan_data_plt = find(isnan(clim2)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             discon_p = find(nonan_diff ~= 1); % find discontinuous point
%             discon_p(4) = length(nonan_data_plt); %end point
%             upper_bound_plt = clim2./30.973762 + clim2_std./30.973762;
%             lower_bound_plt = clim2./30.973762 - clim2_std./30.973762;
%             for i = 1:4
%             if i == 1
%             patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'g');
%             else
%             patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'g');
%             end
%             end 
%             plot(clim2./30.973762,'g-','linew',2); 
%             plot(1:length(clim2),upper_bound_plt,'g-','linew',2);
%             plot(1:length(clim2),lower_bound_plt,'g-','linew',2);
% 
%             alpha(0.3) %transparency      
% 
%             title(['GY KOEM OBS vs. daily MODEL po4']);
%             xlabel('time','fontsize',13)
%             ylabel('po4 (umol/m^3)','fontsize',13)
%             grid on
%             ylim([0 5])
%             xlim(xlim_set)
%             set(gca,'fontsize',13)
%     %         %xticks(com_eom(1:12:end))
%     %         %xticklabels(com_t_tic)
%     %         %xtickangle(90)
%     print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4'),'-dpng')

fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_zoo_mm(:,3,:),1)),'g','linew',2);
%             plot((squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)).*N_MW + squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)).*N_MW)./N_MW,'k--','linew',2);
            plot(squeeze(nanmean(sp_gy_zoo_2reg_mm(:,3,:),1)),'b','linew',2);  
            plot(squeeze(nanmean(sp_gy_zoo_2reg_mm_1(:,3,:),1)),'r','linew',2);  
            plot(squeeze(nanmean(sp_gy_zoo_2reg_mm_2(:,3,:),1)),'k','linew',2);
            plot(squeeze(nanmean(sp_gy_zoo_2reg_mm_3(:,3,:),1)),'m','linew',2);
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL ZOO']);
            xlabel('time','fontsize',13)
            ylabel('ZOO (umol/m^3)','fontsize',13)
            grid on
            ylim([0 inf])
            xlim(xlim_set)
            set(gca,'fontsize',13)
             print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_zoo'),'-dpng')
             
fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_phy_mm(:,3,:),1)),'g','linew',2);
%             plot((squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)).*N_MW + squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)).*N_MW)./N_MW,'k--','linew',2);
            plot(squeeze(nanmean(sp_gy_phy_2reg_mm(:,3,:),1)),'b','linew',2);
            plot(squeeze(nanmean(sp_gy_phy_2reg_mm_1(:,3,:),1)),'r','linew',2);  
            plot(squeeze(nanmean(sp_gy_phy_2reg_mm_2(:,3,:),1)),'k','linew',2);
            plot(squeeze(nanmean(sp_gy_phy_2reg_mm_3(:,3,:),1)),'m','linew',2);
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL PHY']);
            xlabel('time','fontsize',13)
            ylabel('Phy (umol/m^3)','fontsize',13)
            grid on
            ylim([0 inf])
            xlim(xlim_set)
            set(gca,'fontsize',13)
             print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_phy'),'-dpng')

   clearvars clim2*
clim2=clim1.obm_gy_din_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_din_ft_2(:,sigma);
clim2_2=clim1.obm_gy_din_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_din_ft_3(:,sigma);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)) + squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)),'g','linew',2);
%             plot((squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)).*N_MW + squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)).*N_MW)./N_MW,'k--','linew',2);
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm(:,3,:),1)) + squeeze(nanmean(sp_gy_nh4_2reg_mm(:,3,:),1)),'b','linew',2);
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm_1(:,3,:),1)) + squeeze(nanmean(sp_gy_nh4_2reg_mm_1(:,3,:),1)),'r','linew',2);  
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm_2(:,3,:),1)) + squeeze(nanmean(sp_gy_nh4_2reg_mm_2(:,3,:),1)),'k','linew',2);
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm_3(:,3,:),1)) + squeeze(nanmean(sp_gy_nh4_2reg_mm_3(:,3,:),1)),'m','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL DIN']);
            xlabel('time','fontsize',13)
            ylabel('DIN (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_DIN'),'-dpng')

    clearvars clim2*
clim2=clim1.obm_gy_no3_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_no3_ft_2(:,sigma);
clim2_2=clim1.obm_gy_no3_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_no3_ft_3(:,sigma);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)),'g','linew',2);
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm(:,3,:),1)),'b','linew',2);
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm_1(:,3,:),1)) ,'r','linew',2);  
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm_2(:,3,:),1)) ,'k','linew',2);
            plot(squeeze(nanmean(sp_gy_no3_2reg_mm_3(:,3,:),1)) ,'m','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_no3_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_no3_1;
%             lower_bound_plt = temp_d - regm_std_no3_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
% 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL NO3']);
            xlabel('time','fontsize',13)
            ylabel('NO3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')

    
        clearvars clim2*
clim2=clim1.obm_gy_nh4_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_nh4_ft_2(:,sigma);
clim2_2=clim1.obm_gy_nh4_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_nh4_ft_3(:,sigma);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)),'g','linew',2);
            plot(squeeze(nanmean(sp_gy_nh4_2reg_mm(:,3,:),1)),'b','linew',2);
            plot(squeeze(nanmean(sp_gy_nh4_2reg_mm_1(:,3,:),1)),'r','linew',2);  
            plot(squeeze(nanmean(sp_gy_nh4_2reg_mm_2(:,3,:),1)),'k','linew',2);
            plot(squeeze(nanmean(sp_gy_nh4_2reg_mm_3(:,3,:),1)),'m','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_nh4_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_nh4_1;
%             lower_bound_plt = temp_d - regm_std_nh4_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
             % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')


    clearvars clim2*
clim2=clim1.obm_gy_chl_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_chl_ft_2(:,sigma);
clim2_2=clim1.obm_gy_chl_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_chl_ft_3(:,sigma);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_chl_mm(:,3,:),1)),'g','linew',2);
            plot(squeeze(nanmean(sp_gy_chl_2reg_mm(:,3,:),1)),'b','linew',2);
            plot(squeeze(nanmean(sp_gy_chl_2reg_mm_1(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_chl_2reg_mm_2(:,3,:),1)),'k','linew',2);
            plot(squeeze(nanmean(sp_gy_chl_2reg_mm_3(:,3,:),1)),'m','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_chl_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_chl_1;
%             lower_bound_plt = temp_d - regm_std_chl_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
             % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 21])
            xlim(xlim_set)
            set(gca,'fontsize',13)

    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


     clearvars clim2*
clim2=clim1.obm_gy_temp_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_temp_ft_2(:,sigma);
clim2_2=clim1.obm_gy_temp_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_temp_ft_3(:,sigma);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_temp_mm(:,3,:),1)),'g','linew',2);
            plot(squeeze(nanmean(sp_gy_temp_2reg_mm(:,3,:),1)),'b','linew',2);
            plot(squeeze(nanmean(sp_gy_temp_2reg_mm_1(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_temp_2reg_mm_2(:,3,:),1)),'k','linew',2);
            plot(squeeze(nanmean(sp_gy_temp_2reg_mm_3(:,3,:),1)),'m','linew',2);
            %             clearvars *_bound_plt nonan_data_plt discon_p temp_d
%             temp_d=nanmean(mer_temp_1,1);
%             nonan_data_plt = find(isnan(temp_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_d + regm_std_temp_1;
%             lower_bound_plt = temp_d - regm_std_temp_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp_1; c1998.gy_temp_1; c1999.gy_temp_1; c2000.gy_temp_1; c2001.gy_temp_1; c2002.gy_temp_1; c2003.gy_temp_1;],1),'g','linew',2);
                         % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')
    
     
             clearvars clim2*
clim2=clim1.obm_gy_salt_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_salt_ft_1(:,sigma);
clim2_2=clim1.obm_gy_salt_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_salt_ft_2(:,sigma);   

    fig = figure; hold on; 
            plot(squeeze(nanmean(sp_gy_salt_mm(:,3,:),1)),'g','linew',2);
            plot(squeeze(nanmean(sp_gy_salt_2reg_mm(:,3,:),1)),'b','linew',2);
            plot(squeeze(nanmean(sp_gy_salt_2reg_mm_1(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_salt_2reg_mm_2(:,3,:),1)),'k','linew',2);
            plot(squeeze(nanmean(sp_gy_salt_2reg_mm_3(:,3,:),1)),'m','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p salt_d
%             salt_d=nanmean(mer_salt_1,1);
%             nonan_data_plt = find(isnan(salt_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = salt_d + regm_std_salt_1;
%             lower_bound_plt = salt_d - regm_std_salt_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_1; c1998.gy_salt_1; c1999.gy_salt_1; c2000.gy_salt_1; c2001.gy_salt_1; c2002.gy_salt_1; c2003.gy_salt_1;],1),'r','linew',2);
                 % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')
    
    
    clearvars clim2*
clim2=clim1.obm_gy_do_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_do_ft_2(:,sigma);
clim2_2=clim1.obm_gy_do_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_do_ft_3(:,sigma);

    fig = figure; hold on;
            plot(squeeze(nanmean(sp_gy_do_mm(:,3,:),1)),'g','linew',2);
            plot(squeeze(nanmean(sp_gy_do_2reg_mm(:,3,:),1)),'b','linew',2);
            plot(squeeze(nanmean(sp_gy_do_2reg_mm_1(:,3,:),1)),'r','linew',2);
            plot(squeeze(nanmean(sp_gy_do_2reg_mm_2(:,3,:),1)),'k','linew',2);
            plot(squeeze(nanmean(sp_gy_do_2reg_mm_3(:,3,:),1)),'m','linew',2);
% 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2.* DO_CONV + clim2_std.* DO_CONV;
            lower_bound_plt = clim2.* DO_CONV - clim2_std.* DO_CONV;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}).* DO_CONV,'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2.* DO_CONV + clim2_2_std.* DO_CONV;
            lower_bound_plt = clim2_2.* DO_CONV - clim2_2_std.* DO_CONV;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}).* DO_CONV,'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL DO']);
            xlabel('time','fontsize',13)
            ylabel('DO (umol/m^3)','fontsize',13)
            grid on
            ylim([0 380])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_do'),'-dpng')
    
    
    clearvars clim2*
clim2=clim1.obm_gy_salt_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_salt_ft_2(:,sigma);
    
%     fig = figure; hold on; 
%             plot(nanmean(regime1_noriv_chl.gy_salt,1),'g','linew',2);
    
     %% bot
% fig = figure; hold on;
%     %         plot(zeros(365,1),'g--','linew',2); % zero lines
% 
%     %data=NaN(1,365);        
%     % data(1:length(c2004.gy_po4_b))=c2004.gy_po4_b_1; 
%             plot(nanmean(mer_po4_b_1,1),'g','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p po4_b_d
%             po4_b_d=nanmean(mer_po4_b_1,1);
%             nonan_data_plt = find(isnan(po4_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = po4_b_d + regm_std_po4_b_1;
%             lower_bound_plt = po4_b_d - regm_std_po4_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
%     %         plot(nanmean([c1997.gy_po4_b_1; c1998.gy_po4_b_1; c1999.gy_po4_b_1; c2000.gy_po4_b_1; c2001.gy_po4_b_1; c2002.gy_po4_b_1; c2003.gy_po4_b_1;],1),'g','linew',2);
%                    clearvars *_bound_plt nonan_data_plt discon_p
%             nonan_data_plt = find(isnan(clim1.obm_gy_po4_b_ft_1)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             discon_p = find(nonan_diff ~= 1); % find discontinuous point
%             discon_p(4) = length(nonan_data_plt); %end point
%             upper_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 + clim1.obs_std_gy_po4_b_ft_1./30.973762;
%             lower_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 - clim1.obs_std_gy_po4_b_ft_1./30.973762;
%             for i = 1:4
%             if i == 1
%             patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'g');
%             else
%             patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'g');
%             end
%             end 
%             plot(clim1.obm_gy_po4_b_ft_1./30.973762,'g-','linew',2); 
%             plot(1:length(clim1.obm_gy_po4_b_ft_1),upper_bound_plt,'g-','linew',2);
%             plot(1:length(clim1.obm_gy_po4_b_ft_1),lower_bound_plt,'g-','linew',2);
% 
% 
%             alpha(0.3) %transparency      
% 
%             title(['GY KOEM OBS vs. daily MODEL po4 bot.']);
%             xlabel('time','fontsize',13)
%             ylabel('po4 (umol/m^3)','fontsize',13)
%             grid on
%  ylim([0 5])
%             xlim(xlim_set)
%             set(gca,'fontsize',13)
%     %         %xticks(com_eom(1:12:end))
%     %         %xticklabels(com_t_tic)
%     %         %xtickangle(90)
%             print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4_bot'),'-dpng')

   clearvars clim2*
clim2=clim1.obm_gy_din_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_din_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_din_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_din_b_ft_3(:,sigma);

    fig = figure; hold on;
            plot(squeeze(nanmean(regime1_kdc_zerosw.sp_gy_no3_b,1)) + squeeze(nanmean(regime1_kdc_zerosw.sp_gy_nh4_b,1)),'g','linew',2);
%             plot((squeeze(nanmean(sp_gy_no3_mm(:,3,:),1)).*N_MW + squeeze(nanmean(sp_gy_nh4_mm(:,3,:),1)).*N_MW)./N_MW,'k--','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw.sp_gy_no3_b,1)) + squeeze(nanmean(regime2_kdc_zerosw.sp_gy_nh4_b,1)),'b','linew',2);
             plot(squeeze(nanmean(regime2_kdc_zerosw1.sp_gy_no3_b,1)) + squeeze(nanmean(regime2_kdc_zerosw1.sp_gy_nh4_b,1)),'r','linew',2);
              plot(squeeze(nanmean(regime2_kdc_zerosw2.sp_gy_no3_b,1)) + squeeze(nanmean(regime2_kdc_zerosw2.sp_gy_nh4_b,1)),'k','linew',2);
              plot(squeeze(nanmean(regime2_kdc_zerosw3.sp_gy_no3_b,1)) + squeeze(nanmean(regime2_kdc_zerosw3.sp_gy_nh4_b,1)),'m','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL DIN']);
            xlabel('time','fontsize',13)
            ylabel('DIN (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_DIN_bot'),'-dpng')

clearvars clim2*
clim2=clim1.obm_gy_no3_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_no3_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_no3_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_no3_b_ft_3(:,sigma);

     fig = figure; hold on;
            plot(squeeze(nanmean(regime1_kdc_zerosw.sp_gy_no3_b,1)),'g','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw.sp_gy_no3_b,1)),'b','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw1.sp_gy_no3_b,1)),'r','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw2.sp_gy_no3_b,1)),'k','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw3.sp_gy_no3_b,1)),'m','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p no3_b_d
%             no3_b_d=nanmean(mer_no3_b_1,1);
%             nonan_data_plt = find(isnan(no3_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = no3_b_d + regm_std_no3_b_1;
%             lower_bound_plt = no3_b_d - regm_std_no3_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_no3_b_1; c1998.gy_no3_b_1; c1999.gy_no3_b_1; c2000.gy_no3_b_1; c2001.gy_no3_b_1; c2002.gy_no3_b_1; c2003.gy_no3_b_1;],1),'g','linew',2);
% 1regime
             clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
            xlabel('time','fontsize',13)
            ylabel('no3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

      clearvars clim2*
clim2=clim1.obm_gy_nh4_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_nh4_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_nh4_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_nh4_b_ft_3(:,sigma);

                
    fig = figure; hold on;
            plot(squeeze(nanmean(regime1_kdc_zerosw.sp_gy_nh4_b,1)),'g','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw.sp_gy_nh4_b,1)),'b','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw1.sp_gy_nh4_b,1)),'r','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw2.sp_gy_nh4_b,1)),'k','linew',2);
             plot(squeeze(nanmean(regime2_kdc_zerosw3.sp_gy_nh4_b,1)),'k','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p nh4_b_d
%             nh4_b_d=nanmean(mer_nh4_b_1,1);
%             nonan_data_plt = find(isnan(nh4_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = nh4_b_d + regm_std_nh4_b_1;
%             lower_bound_plt = nh4_b_d - regm_std_nh4_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
%     %         plot(nanmean([c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b;],1),'g','linew',2);
            % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2./N_MW + clim2_std./N_MW;
            lower_bound_plt = clim2./N_MW - clim2_std./N_MW;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i})./N_MW,'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2./N_MW + clim2_2_std./N_MW;
            lower_bound_plt = clim2_2./N_MW - clim2_2_std./N_MW;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i})./N_MW,'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')

          clearvars clim2*
clim2=clim1.obm_gy_chl_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_chl_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_chl_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_chl_b_ft_3(:,sigma);
                
    fig = figure; hold on;
            plot(squeeze(nanmean(regime1_kdc_zerosw.sp_gy_chl_b,1)),'g','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw.sp_gy_chl_b,1)),'b','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw1.sp_gy_chl_b,1)),'r','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw2.sp_gy_chl_b,1)),'k','linew',2);
            plot(squeeze(nanmean(regime2_kdc_zerosw3.sp_gy_chl_b,1)),'m','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p chl_b_d
%             chl_b_d=nanmean(mer_chl_b_1,1);
%             nonan_data_plt = find(isnan(chl_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = chl_b_d + regm_std_chl_b_1;
%             lower_bound_plt = chl_b_d - regm_std_chl_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b;],1),'g','linew',2);
                        % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL chl bot.']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 20])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')
      


            clearvars clim2*
clim2=clim1.obm_gy_temp_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_temp_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_temp_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_temp_b_ft_3(:,sigma);  
                
    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_temp_b,1),'g','linew',2);
            plot(nanmean(regime2_kdc_zerosw.gy_temp_b,1),'b','linew',2); 
            plot(nanmean(regime2_kdc_zerosw1.gy_temp_b,1),'r','linew',2); 
            plot(nanmean(regime2_kdc_zerosw2.gy_temp_b,1),'k','linew',2); 
             plot(nanmean(regime2_kdc_zerosw3.gy_temp_b,1),'m','linew',2); 
%             clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
%             temp_b_d=nanmean(mer_temp_b_1,1);
%             nonan_data_plt = find(isnan(temp_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = temp_b_d + regm_std_temp_b_1;
%             lower_bound_plt = temp_b_d - regm_std_temp_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'g','linew',2);
                     % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')

            clearvars clim2*
clim2=clim1.obm_gy_salt_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_salt_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_salt_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_salt_b_ft_3(:,sigma);   
                
    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.gy_salt_b,1),'g','linew',2);  
            plot(nanmean(regime2_kdc_zerosw.gy_salt_b,1),'b','linew',2);
            plot(nanmean(regime2_kdc_zerosw1.gy_salt_b,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw2.gy_salt_b,1),'k','linew',2);
            plot(nanmean(regime2_kdc_zerosw2.gy_salt_b,1),'m','linew',2);
%             clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
%             salt_b_d=nanmean(mer_salt_b_1,1);
%             nonan_data_plt = find(isnan(salt_b_d)==0);
%             nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
%             upper_bound_plt = salt_b_d + regm_std_salt_b_1;
%             lower_bound_plt = salt_b_d - regm_std_salt_b_1;
%             patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'g','linew',2);
              % 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2 + clim2_std;
            lower_bound_plt = clim2 - clim2_std;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2 + clim2_2_std;
            lower_bound_plt = clim2_2 - clim2_2_std;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            alpha(0.3) %transparency      
            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng') 
            
                clearvars clim2*
clim2=clim1.obm_gy_do_b_ft_2(:,sigma);
clim2_std = clim1.obs_std_gy_do_b_ft_2(:,sigma);
clim2_2=clim1.obm_gy_do_b_ft_3(:,sigma);
clim2_2_std = clim1.obs_std_gy_do_b_ft_3(:,sigma);

    fig = figure; hold on;
            plot(nanmean(regime1_kdc_zerosw.sp_gy_do_b,1),'g','linew',2);
            plot(nanmean(regime2_kdc_zerosw.sp_gy_do_b,1),'b','linew',2);
            plot(nanmean(regime2_kdc_zerosw1.sp_gy_do_b,1),'r','linew',2);
            plot(nanmean(regime2_kdc_zerosw2.sp_gy_do_b,1),'k','linew',2);
            plot(nanmean(regime2_kdc_zerosw3.sp_gy_do_b,1),'m','linew',2);
% 1regime
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2.* DO_CONV + clim2_std.* DO_CONV;
            lower_bound_plt = clim2.* DO_CONV - clim2_std.* DO_CONV;
            for i = 1:4
            patch([plt_filter_1{i} fliplr(plt_filter_1{i})], [upper_bound_plt(plt_filter_1{i})' fliplr(lower_bound_plt(plt_filter_1{i}))'],'g');
            plot(plt_filter_1{i},clim2(plt_filter_1{i}).* DO_CONV,'g-','linew',2);
            plot(plt_filter_1{i},upper_bound_plt(plt_filter_1{i}),'g-','linew',2);
            plot(plt_filter_1{i},lower_bound_plt(plt_filter_1{i}),'g-','linew',2);
            end             
 % 2regime           
            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim2_2)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim2_2.* DO_CONV + clim2_2_std.* DO_CONV;
            lower_bound_plt = clim2_2.* DO_CONV - clim2_2_std.* DO_CONV;
            for i = 1:4
            patch([plt_filter_2{i} fliplr(plt_filter_2{i})], [upper_bound_plt(plt_filter_2{i})' fliplr(lower_bound_plt(plt_filter_2{i}))'],'b');
            plot(plt_filter_2{i},clim2_2(plt_filter_2{i}).* DO_CONV,'b-','linew',2);
            plot(plt_filter_2{i},upper_bound_plt(plt_filter_2{i}),'b-','linew',2);
            plot(plt_filter_2{i},lower_bound_plt(plt_filter_2{i}),'b-','linew',2);
            end 
            
            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL DO bot']);
            xlabel('time','fontsize',13)
            ylabel('DO (umol/m^3)','fontsize',13)
            grid on
            ylim([0 380])
            xlim(xlim_set)
            set(gca,'fontsize',13)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_do_bot'),'-dpng')
            
            close all;
            
            return
%% mean            
for i = 1:4; sp_yr_chl_2reg_s(i)=mean(sp_gy_chl_2reg_mm(:,3,[plt_filter_1{i},plt_filter_2{i}]),[1 2 3],'omitnan'); end
for i = 1:4; sp_yr_chl_s(i)=mean(sp_gy_chl_mm(:,3,[plt_filter_1{i},plt_filter_2{i}]),[1 2 3],'omitnan'); end
%% monthly max and mean
for i = 1:4; sp_yr_chl_2reg_s2(i)=max(sp_gy_chl_2reg_mm(:,3,[plt_filter_1{i},plt_filter_2{i}]),[],'all','omitnan'); end
for i = 1:4; sp_yr_chl_s2(i)=max(sp_gy_chl_mm(:,3,[plt_filter_1{i},plt_filter_2{i}]),[],'all','omitnan'); end

yline(mean(sp_yr_chl_2reg_s),'r--','linew',2);
yline(mean(sp_yr_chl_s),'k--','linew',2)
yline(mean(sp_yr_chl_2reg_s2),'r-o','linew',2);
yline(mean(sp_yr_chl_s2),'k-o','linew',2)
ylim([0 13])
ylabel('chl (mg/m^3)','fontsize',13)

mean(mod_chl_pick(:,3,1:228),[1 2 3],'omitnan') - mean(mod_chl_pick_2reg(:,3,229:276),[1 2 3],'omitnan')
