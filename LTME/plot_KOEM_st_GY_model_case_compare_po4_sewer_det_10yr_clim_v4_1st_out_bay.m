close all; clear;clc

case_num = [1:3];
sigsig=[1;2;3;];
% 1: spatial & time std.
% 2: time std.
% 3: spatial std.
 dir_name ={'all_std','t_std','sp_std'};


cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
ltme =load('LTME_observation_data_v2.mat'); %LTME obs.

% cd D:\장기생태\Dynamic\KOEM\gy_2001\gy_2001_koem_daily\spmean
% load('koem_result_processed_std_spmean.mat'); % KOEM obs.
% % mpdata= load('mpdata_result.mat');

cd D:\장기생태\Dynamic\result\1997
c1997 = load('1997_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*','name_tag','sp_gy');

cd D:\장기생태\Dynamic\result\1998
c1998 = load('1998_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\1999
c1999 = load('1999_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2000
c2000 = load('2000_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');



% first_10;
c1997_eom = [1 (c1997.eom_d_each(1:end-1) + 1)];
c1998_eom = [(c1997.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(1:end-1) + 1)];
c1999_eom = [(c1997.eom_d_each(end)+c1998.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(1:end-1) + 1)];
c2000_eom = [(c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(end)+c2000.eom_d_each(1:end-1) + 1)];

com_eom1 = [c1997_eom c1998_eom c1999_eom c2000_eom];

clearvars com_t_tic
% com_t_tic = {'2001.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2002.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2003.01','02','03','04','05','06','07','08','09','10','11','12',...}
%             '2004.01','02','03','04','05','06','07','08','09','10','11','12'}
com_t_tic = {'1997','1998','1999','2000','2001', ...
            '2002', ...
            '2003'};
        
c1997.sp_gy_no3(:,60)=[]; c1997.sp_gy_chl(:,60)=[]; c1997.sp_gy_nh4(:,60)=[]; c1997.sp_gy_temp(:,60)=[]; c1997.sp_gy_salt(:,60)=[];
c1998.sp_gy_no3(:,60)=[]; c1998.sp_gy_chl(:,60)=[]; c1998.sp_gy_nh4(:,60)=[]; c1998.sp_gy_temp(:,60)=[]; c1998.sp_gy_salt(:,60)=[];
c1999.sp_gy_no3(:,60)=[]; c1999.sp_gy_chl(:,60)=[]; c1999.sp_gy_nh4(:,60)=[]; c1999.sp_gy_temp(:,60)=[]; c1999.sp_gy_salt(:,60)=[];

c1997.sp_gy_no3_b(:,60)=[]; c1997.sp_gy_chl_b(:,60)=[]; c1997.sp_gy_nh4_b(:,60)=[]; c1997.sp_gy_temp_b(:,60)=[]; c1997.sp_gy_salt_b(:,60)=[];
c1998.sp_gy_no3_b(:,60)=[]; c1998.sp_gy_chl_b(:,60)=[]; c1998.sp_gy_nh4_b(:,60)=[]; c1998.sp_gy_temp_b(:,60)=[]; c1998.sp_gy_salt_b(:,60)=[];
c1999.sp_gy_no3_b(:,60)=[]; c1999.sp_gy_chl_b(:,60)=[]; c1999.sp_gy_nh4_b(:,60)=[]; c1999.sp_gy_temp_b(:,60)=[]; c1999.sp_gy_salt_b(:,60)=[];


c1997.gy_no3(60)=[]; c1997.gy_chl(60)=[]; c1997.gy_nh4(60)=[]; c1997.gy_temp(60)=[]; c1997.gy_salt(60)=[];
c1998.gy_no3(60)=[]; c1998.gy_chl(60)=[]; c1998.gy_nh4(60)=[]; c1998.gy_temp(60)=[]; c1998.gy_salt(60)=[];
c1999.gy_no3(60)=[]; c1999.gy_chl(60)=[]; c1999.gy_nh4(60)=[]; c1999.gy_temp(60)=[]; c1999.gy_salt(60)=[];


c1997.gy_no3_b(60)=[]; c1997.gy_chl_b(60)=[]; c1997.gy_nh4_b(60)=[]; c1997.gy_temp_b(60)=[]; c1997.gy_salt_b(60)=[];
c1998.gy_no3_b(60)=[]; c1998.gy_chl_b(60)=[]; c1998.gy_nh4_b(60)=[]; c1998.gy_temp_b(60)=[]; c1998.gy_salt_b(60)=[];
c1999.gy_no3_b(60)=[]; c1999.gy_chl_b(60)=[]; c1999.gy_nh4_b(60)=[]; c1999.gy_temp_b(60)=[]; c1999.gy_salt_b(60)=[];

c2000.sp_gy_no3(:,60)=[]; c2000.sp_gy_chl(:,60)=[]; c2000.sp_gy_nh4(:,60)=[]; c2000.sp_gy_temp(:,60)=[]; c2000.sp_gy_salt(:,60)=[];    
c2000.sp_gy_no3_b(:,60)=[]; c2000.sp_gy_chl_b(:,60)=[]; c2000.sp_gy_nh4_b(:,60)=[]; c2000.sp_gy_temp_b(:,60)=[]; c2000.sp_gy_salt_b(:,60)=[];
c2000.gy_no3(60)=[]; c2000.gy_chl(60)=[]; c2000.gy_nh4(60)=[]; c2000.gy_temp(60)=[]; c2000.gy_salt(60)=[];
c2000.gy_no3_b(60)=[]; c2000.gy_chl_b(60)=[]; c2000.gy_nh4_b(60)=[]; c2000.gy_temp_b(60)=[]; c2000.gy_salt_b(60)=[];
   
c1997.sp_gy_po4(:,60)=[]; c1998.sp_gy_po4(:,60)=[]; c1999.sp_gy_po4(:,60)=[]; c1997.sp_gy_po4_b(:,60)=[]; 
c1998.sp_gy_po4_b(:,60)=[]; c1999.sp_gy_po4_b(:,60)=[]; c2000.sp_gy_po4(:,60)=[]; c2000.sp_gy_po4_b(:,60)=[];

c1997.gy_po4(60)=[]; c1998.gy_po4(60)=[]; c1999.gy_po4(60)=[]; c1997.gy_po4_b(60)=[]; 
c1998.gy_po4_b(60)=[]; c1999.gy_po4_b(60)=[]; c2000.gy_po4(60)=[]; c2000.gy_po4_b(60)=[];

%% out bay only left on model results
% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 
[res I]=sort([1,5,4,3,2,6,8,9,7]);
% sorting confirm
nn=['광양항'; '광양4'; '광양3'; '광양2'; '광양1'; '광양5'; '여수2'; '여수3'; '여수1']; % order of origin koem
nn(I,:) % order of yoonjas koem
sp_gy=c1997.sp_gy;
name_tag=c1997.name_tag;
name_tag{sp_gy(I)} % means nn(I,:)

sp_gy(I(5:9))
sp_gy(I(1))

yoonja_indx= I;

c1997.sp_gy_no3(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_chl(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_nh4(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_temp(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_salt(yoonja_indx(1:5),:)=NaN;
c1998.sp_gy_no3(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_chl(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_nh4(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_temp(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_salt(yoonja_indx(1:5),:)=NaN;
c1999.sp_gy_no3(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_chl(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_nh4(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_temp(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_salt(yoonja_indx(1:5),:)=NaN;

c1997.sp_gy_no3_b(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_chl_b(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_nh4_b(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_temp_b(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_salt_b(yoonja_indx(1:5),:)=NaN;
c1998.sp_gy_no3_b(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_chl_b(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_nh4_b(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_temp_b(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_salt_b(yoonja_indx(1:5),:)=NaN;
c1999.sp_gy_no3_b(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_chl_b(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_nh4_b(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_temp_b(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_salt_b(yoonja_indx(1:5),:)=NaN;


c1997.gy_no3(yoonja_indx(1:5),:)=NaN; c1997.gy_chl(yoonja_indx(1:5),:)=NaN; c1997.gy_nh4(yoonja_indx(1:5),:)=NaN; c1997.gy_temp(yoonja_indx(1:5),:)=NaN; c1997.gy_salt(yoonja_indx(1:5),:)=NaN;
c1998.gy_no3(yoonja_indx(1:5),:)=NaN; c1998.gy_chl(yoonja_indx(1:5),:)=NaN; c1998.gy_nh4(yoonja_indx(1:5),:)=NaN; c1998.gy_temp(yoonja_indx(1:5),:)=NaN; c1998.gy_salt(yoonja_indx(1:5),:)=NaN;
c1999.gy_no3(yoonja_indx(1:5),:)=NaN; c1999.gy_chl(yoonja_indx(1:5),:)=NaN; c1999.gy_nh4(yoonja_indx(1:5),:)=NaN; c1999.gy_temp(yoonja_indx(1:5),:)=NaN; c1999.gy_salt(yoonja_indx(1:5),:)=NaN;


c1997.gy_no3_b(yoonja_indx(1:5),:)=NaN; c1997.gy_chl_b(yoonja_indx(1:5),:)=NaN; c1997.gy_nh4_b(yoonja_indx(1:5),:)=NaN; c1997.gy_temp_b(yoonja_indx(1:5),:)=NaN; c1997.gy_salt_b(yoonja_indx(1:5),:)=NaN;
c1998.gy_no3_b(yoonja_indx(1:5),:)=NaN; c1998.gy_chl_b(yoonja_indx(1:5),:)=NaN; c1998.gy_nh4_b(yoonja_indx(1:5),:)=NaN; c1998.gy_temp_b(yoonja_indx(1:5),:)=NaN; c1998.gy_salt_b(yoonja_indx(1:5),:)=NaN;
c1999.gy_no3_b(yoonja_indx(1:5),:)=NaN; c1999.gy_chl_b(yoonja_indx(1:5),:)=NaN; c1999.gy_nh4_b(yoonja_indx(1:5),:)=NaN; c1999.gy_temp_b(yoonja_indx(1:5),:)=NaN; c1999.gy_salt_b(yoonja_indx(1:5),:)=NaN;

c2000.sp_gy_no3(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_chl(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_nh4(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_temp(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_salt(yoonja_indx(1:5),:)=NaN;    
c2000.sp_gy_no3_b(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_chl_b(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_nh4_b(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_temp_b(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_salt_b(yoonja_indx(1:5),:)=NaN;
c2000.gy_no3(yoonja_indx(1:5),:)=NaN; c2000.gy_chl(yoonja_indx(1:5),:)=NaN; c2000.gy_nh4(yoonja_indx(1:5),:)=NaN; c2000.gy_temp(yoonja_indx(1:5),:)=NaN; c2000.gy_salt(yoonja_indx(1:5),:)=NaN;
c2000.gy_no3_b(yoonja_indx(1:5),:)=NaN; c2000.gy_chl_b(yoonja_indx(1:5),:)=NaN; c2000.gy_nh4_b(yoonja_indx(1:5),:)=NaN; c2000.gy_temp_b(yoonja_indx(1:5),:)=NaN; c2000.gy_salt_b(yoonja_indx(1:5),:)=NaN;
   
c1997.sp_gy_po4(yoonja_indx(1:5),:)=NaN; c1998.sp_gy_po4(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_po4(yoonja_indx(1:5),:)=NaN; c1997.sp_gy_po4_b(yoonja_indx(1:5),:)=NaN; 
c1998.sp_gy_po4_b(yoonja_indx(1:5),:)=NaN; c1999.sp_gy_po4_b(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_po4(yoonja_indx(1:5),:)=NaN; c2000.sp_gy_po4_b(yoonja_indx(1:5),:)=NaN;

c1997.gy_po4(yoonja_indx(1:5),:)=NaN; c1998.gy_po4(yoonja_indx(1:5),:)=NaN; c1999.gy_po4(yoonja_indx(1:5),:)=NaN; c1997.gy_po4_b(yoonja_indx(1:5),:)=NaN; 
c1998.gy_po4_b(yoonja_indx(1:5),:)=NaN; c1999.gy_po4_b(yoonja_indx(1:5),:)=NaN; c2000.gy_po4(yoonja_indx(1:5),:)=NaN; c2000.gy_po4_b(yoonja_indx(1:5),:)=NaN;

% %
% c1997.sp_gy_no3(yoonja_indx(1),:)=NaN; c1997.sp_gy_chl(yoonja_indx(1),:)=NaN; c1997.sp_gy_nh4(yoonja_indx(1),:)=NaN; c1997.sp_gy_temp(yoonja_indx(1),:)=NaN; c1997.sp_gy_salt(yoonja_indx(1),:)=NaN;
% c1998.sp_gy_no3(yoonja_indx(1),:)=NaN; c1998.sp_gy_chl(yoonja_indx(1),:)=NaN; c1998.sp_gy_nh4(yoonja_indx(1),:)=NaN; c1998.sp_gy_temp(yoonja_indx(1),:)=NaN; c1998.sp_gy_salt(yoonja_indx(1),:)=NaN;
% c1999.sp_gy_no3(yoonja_indx(1),:)=NaN; c1999.sp_gy_chl(yoonja_indx(1),:)=NaN; c1999.sp_gy_nh4(yoonja_indx(1),:)=NaN; c1999.sp_gy_temp(yoonja_indx(1),:)=NaN; c1999.sp_gy_salt(yoonja_indx(1),:)=NaN;
% 
% c1997.sp_gy_no3_b(yoonja_indx(1),:)=NaN; c1997.sp_gy_chl_b(yoonja_indx(1),:)=NaN; c1997.sp_gy_nh4_b(yoonja_indx(1),:)=NaN; c1997.sp_gy_temp_b(yoonja_indx(1),:)=NaN; c1997.sp_gy_salt_b(yoonja_indx(1),:)=NaN;
% c1998.sp_gy_no3_b(yoonja_indx(1),:)=NaN; c1998.sp_gy_chl_b(yoonja_indx(1),:)=NaN; c1998.sp_gy_nh4_b(yoonja_indx(1),:)=NaN; c1998.sp_gy_temp_b(yoonja_indx(1),:)=NaN; c1998.sp_gy_salt_b(yoonja_indx(1),:)=NaN;
% c1999.sp_gy_no3_b(yoonja_indx(1),:)=NaN; c1999.sp_gy_chl_b(yoonja_indx(1),:)=NaN; c1999.sp_gy_nh4_b(yoonja_indx(1),:)=NaN; c1999.sp_gy_temp_b(yoonja_indx(1),:)=NaN; c1999.sp_gy_salt_b(yoonja_indx(1),:)=NaN;
% 
% 
% c1997.gy_no3(yoonja_indx(1),:)=NaN; c1997.gy_chl(yoonja_indx(1),:)=NaN; c1997.gy_nh4(yoonja_indx(1),:)=NaN; c1997.gy_temp(yoonja_indx(1),:)=NaN; c1997.gy_salt(yoonja_indx(1),:)=NaN;
% c1998.gy_no3(yoonja_indx(1),:)=NaN; c1998.gy_chl(yoonja_indx(1),:)=NaN; c1998.gy_nh4(yoonja_indx(1),:)=NaN; c1998.gy_temp(yoonja_indx(1),:)=NaN; c1998.gy_salt(yoonja_indx(1),:)=NaN;
% c1999.gy_no3(yoonja_indx(1),:)=NaN; c1999.gy_chl(yoonja_indx(1),:)=NaN; c1999.gy_nh4(yoonja_indx(1),:)=NaN; c1999.gy_temp(yoonja_indx(1),:)=NaN; c1999.gy_salt(yoonja_indx(1),:)=NaN;
% 
% 
% c1997.gy_no3_b(yoonja_indx(1),:)=NaN; c1997.gy_chl_b(yoonja_indx(1),:)=NaN; c1997.gy_nh4_b(yoonja_indx(1),:)=NaN; c1997.gy_temp_b(yoonja_indx(1),:)=NaN; c1997.gy_salt_b(yoonja_indx(1),:)=NaN;
% c1998.gy_no3_b(yoonja_indx(1),:)=NaN; c1998.gy_chl_b(yoonja_indx(1),:)=NaN; c1998.gy_nh4_b(yoonja_indx(1),:)=NaN; c1998.gy_temp_b(yoonja_indx(1),:)=NaN; c1998.gy_salt_b(yoonja_indx(1),:)=NaN;
% c1999.gy_no3_b(yoonja_indx(1),:)=NaN; c1999.gy_chl_b(yoonja_indx(1),:)=NaN; c1999.gy_nh4_b(yoonja_indx(1),:)=NaN; c1999.gy_temp_b(yoonja_indx(1),:)=NaN; c1999.gy_salt_b(yoonja_indx(1),:)=NaN;
% 
% c2000.sp_gy_no3(yoonja_indx(1),:)=NaN; c2000.sp_gy_chl(yoonja_indx(1),:)=NaN; c2000.sp_gy_nh4(yoonja_indx(1),:)=NaN; c2000.sp_gy_temp(yoonja_indx(1),:)=NaN; c2000.sp_gy_salt(yoonja_indx(1),:)=NaN;    
% c2000.sp_gy_no3_b(yoonja_indx(1),:)=NaN; c2000.sp_gy_chl_b(yoonja_indx(1),:)=NaN; c2000.sp_gy_nh4_b(yoonja_indx(1),:)=NaN; c2000.sp_gy_temp_b(yoonja_indx(1),:)=NaN; c2000.sp_gy_salt_b(yoonja_indx(1),:)=NaN;
% c2000.gy_no3(yoonja_indx(1),:)=NaN; c2000.gy_chl(yoonja_indx(1),:)=NaN; c2000.gy_nh4(yoonja_indx(1),:)=NaN; c2000.gy_temp(yoonja_indx(1),:)=NaN; c2000.gy_salt(yoonja_indx(1),:)=NaN;
% c2000.gy_no3_b(yoonja_indx(1),:)=NaN; c2000.gy_chl_b(yoonja_indx(1),:)=NaN; c2000.gy_nh4_b(yoonja_indx(1),:)=NaN; c2000.gy_temp_b(yoonja_indx(1),:)=NaN; c2000.gy_salt_b(yoonja_indx(1),:)=NaN;
%    
% c1997.sp_gy_po4(yoonja_indx(1),:)=NaN; c1998.sp_gy_po4(yoonja_indx(1),:)=NaN; c1999.sp_gy_po4(yoonja_indx(1),:)=NaN; c1997.sp_gy_po4_b(yoonja_indx(1),:)=NaN; 
% c1998.sp_gy_po4_b(yoonja_indx(1),:)=NaN; c1999.sp_gy_po4_b(yoonja_indx(1),:)=NaN; c2000.sp_gy_po4(yoonja_indx(1),:)=NaN; c2000.sp_gy_po4_b(yoonja_indx(1),:)=NaN;
% 
% c1997.gy_po4(yoonja_indx(1),:)=NaN; c1998.gy_po4(yoonja_indx(1),:)=NaN; c1999.gy_po4(yoonja_indx(1),:)=NaN; c1997.gy_po4_b(yoonja_indx(1),:)=NaN; 
% c1998.gy_po4_b(yoonja_indx(1),:)=NaN; c1999.gy_po4_b(yoonja_indx(1),:)=NaN; c2000.gy_po4(yoonja_indx(1),:)=NaN; c2000.gy_po4_b(yoonja_indx(1),:)=NaN;


%% 1st - 10yr        
mer_no3_1 = [nanmean(c1997.sp_gy_no3,1); nanmean(c1998.sp_gy_no3,1); nanmean(c1999.sp_gy_no3,1); nanmean(c2000.sp_gy_no3,1);];
mer_nh4_1 = [nanmean(c1997.sp_gy_nh4,1); nanmean(c1998.sp_gy_nh4,1); nanmean(c1999.sp_gy_nh4,1); nanmean(c2000.sp_gy_nh4,1);];
mer_chl_1 = [nanmean(c1997.sp_gy_chl,1); nanmean(c1998.sp_gy_chl,1); nanmean(c1999.sp_gy_chl,1); nanmean(c2000.sp_gy_chl,1);];
mer_temp_1 = [nanmean(c1997.sp_gy_temp,1); nanmean(c1998.sp_gy_temp,1); nanmean(c1999.sp_gy_temp,1); nanmean(c2000.sp_gy_temp,1);];
mer_salt_1 = [nanmean(c1997.sp_gy_salt,1); nanmean(c1998.sp_gy_salt,1); nanmean(c1999.sp_gy_salt,1); nanmean(c2000.sp_gy_salt,1);];
mer_po4_1 = [nanmean(c1997.sp_gy_po4,1); nanmean(c1998.sp_gy_po4,1); nanmean(c1999.sp_gy_po4,1); nanmean(c2000.sp_gy_po4,1);];

%bot merged
mer_no3_b_1 = [nanmean(c1997.sp_gy_no3_b,1); nanmean(c1998.sp_gy_no3_b,1); nanmean(c1999.sp_gy_no3_b,1); nanmean(c2000.sp_gy_no3_b,1);];
mer_nh4_b_1 = [nanmean(c1997.sp_gy_nh4_b,1); nanmean(c1998.sp_gy_nh4_b,1); nanmean(c1999.sp_gy_nh4_b,1); nanmean(c2000.sp_gy_nh4_b,1);];
mer_chl_b_1 = [nanmean(c1997.sp_gy_chl_b,1); nanmean(c1998.sp_gy_chl_b,1); nanmean(c1999.sp_gy_chl_b,1); nanmean(c2000.sp_gy_chl_b,1);];
mer_temp_b_1 = [nanmean(c1997.sp_gy_temp_b,1); nanmean(c1998.sp_gy_temp_b,1); nanmean(c1999.sp_gy_temp_b,1); nanmean(c2000.sp_gy_temp_b,1);];
mer_salt_b_1 = [nanmean(c1997.sp_gy_salt_b,1); nanmean(c1998.sp_gy_salt_b,1); nanmean(c1999.sp_gy_salt_b,1); nanmean(c2000.sp_gy_salt_b,1);];
mer_po4_b_1 = [nanmean(c1997.sp_gy_po4_b,1); nanmean(c1998.sp_gy_po4_b,1); nanmean(c1999.sp_gy_po4_b,1); nanmean(c2000.sp_gy_po4_b,1);];

% sptatial value
%surf merged
mer_sp_no3_1 = [c1997.sp_gy_no3 c1998.sp_gy_no3 c1999.sp_gy_no3 c2000.sp_gy_no3];
mer_sp_nh4_1 = [c1997.sp_gy_nh4 c1998.sp_gy_nh4 c1999.sp_gy_nh4 c2000.sp_gy_nh4];
mer_sp_chl_1 = [c1997.sp_gy_chl c1998.sp_gy_chl c1999.sp_gy_chl c2000.sp_gy_chl];
mer_sp_temp_1 = [c1997.sp_gy_temp c1998.sp_gy_temp c1999.sp_gy_temp c2000.sp_gy_temp];
mer_sp_salt_1 = [c1997.sp_gy_salt c1998.sp_gy_salt c1999.sp_gy_salt c2000.sp_gy_salt];
mer_sp_po4_1 = [c1997.sp_gy_po4 c1998.sp_gy_po4 c1999.sp_gy_po4 c2000.sp_gy_po4];

%bot merged
mer_sp_no3_b_1 = [c1997.sp_gy_no3_b c1998.sp_gy_no3_b c1999.sp_gy_no3_b c2000.sp_gy_no3_b];
mer_sp_nh4_b_1 = [c1997.sp_gy_nh4_b c1998.sp_gy_nh4_b c1999.sp_gy_nh4_b c2000.sp_gy_nh4_b];
mer_sp_chl_b_1 = [c1997.sp_gy_chl_b c1998.sp_gy_chl_b c1999.sp_gy_chl_b c2000.sp_gy_chl_b];
mer_sp_temp_b_1 = [c1997.sp_gy_temp_b c1998.sp_gy_temp_b c1999.sp_gy_temp_b c2000.sp_gy_temp_b];
mer_sp_salt_b_1 = [c1997.sp_gy_salt_b c1998.sp_gy_salt_b c1999.sp_gy_salt_b c2000.sp_gy_salt_b];
mer_sp_po4_b_1 = [c1997.sp_gy_po4_b c1998.sp_gy_po4_b c1999.sp_gy_po4_b c2000.sp_gy_po4_b];


for case_num = 1:3
    for ixx=1:3
        clearvars -except ixx case_num sigsig c199* c200* ltme *_eom dir_name mer_*
        
        sig=sigsig(ixx)

cd D:\장기생태\Dynamic\koem
% clim1=load('koem_climate_3regime_v2_2sig_gy_to_2003_12.mat'); % from 'plot_KOEM_st_GY_model_compare_horizontal_plt_climate_3regime_v2_gy_2sig.m'
% clim1=load('koem_climate_3regime_v2_to_2003_12.mat');  % from 'plot_KOEM_st_GY_model_compare_horizontal_plt_seasonal_climate_3regime_v2.m'
% clim1=load('koem_climate_3regime_v2_3sig_gy_to_2003_12.mat');  % from 'plot_KOEM_st_seasonal_climate_3regime_v3_spatial_sigma_extract.m'

clim1=load(['koem_climate_10yr_v5_',num2str(sig),'sig_gy_out_bay.mat']);

% 1997~2018 : 22yr
% 1st regime : 1997 ~ 2003.12 (from koem nh4 shift)
% 2nd regime : 2004.01 ~ 2009.12 (from koem no3 shift)
% 3rd regime : 2010.01 ~ 2018.12 (from koem no3 shift)

% 5yr cutting
% 1997~1999 : first_5;
% 2000~2004 : second_5;
% 2005~2009 : third_5;

% % mpdata.gy_temp


for i =  1:365
%% model 
% time & spatial std
    %% 1st
    clearvars size2
    size2=size(mer_sp_no3_1(:,i:365:end),2);
    regm_std_no3_1(i)=nanstd(reshape(mer_sp_no3_1(:,i:365:end),1,9*size2));
    regm_std_nh4_1(i)=nanstd(reshape(mer_sp_nh4_1(:,i:365:end),1,9*size2));
    regm_std_chl_1(i)=nanstd(reshape(mer_sp_chl_1(:,i:365:end),1,9*size2));
    regm_std_temp_1(i)=nanstd(reshape(mer_sp_temp_1(:,i:365:end),1,9*size2));
    regm_std_salt_1(i)=nanstd(reshape(mer_sp_salt_1(:,i:365:end),1,9*size2));
    regm_std_po4_1(i)=nanstd(reshape(mer_sp_po4_1(:,i:365:end),1,9*size2));
    
    clearvars size2
    size2=size(mer_sp_no3_1(:,i:365:end),2);
    regm_std_no3_b_1(i)=nanstd(reshape(mer_sp_no3_b_1(:,i:365:end),1,9*size2));
    regm_std_nh4_b_1(i)=nanstd(reshape(mer_sp_nh4_b_1(:,i:365:end),1,9*size2));
    regm_std_chl_b_1(i)=nanstd(reshape(mer_sp_chl_b_1(:,i:365:end),1,9*size2));
    regm_std_temp_b_1(i)=nanstd(reshape(mer_sp_temp_b_1(:,i:365:end),1,9*size2));
    regm_std_salt_b_1(i)=nanstd(reshape(mer_sp_salt_b_1(:,i:365:end),1,9*size2));
    regm_std_po4_b_1(i)=nanstd(reshape(mer_sp_po4_b_1(:,i:365:end),1,9*size2));
    
% time std  
    %% 1st
    spm_regm_std_no3_1(i)=nanstd(nanmean(mer_sp_no3_1(:,i:365:end),1));
    spm_regm_std_nh4_1(i)=nanstd(nanmean(mer_sp_nh4_1(:,i:365:end),1));
    spm_regm_std_chl_1(i)=nanstd(nanmean(mer_sp_chl_1(:,i:365:end),1));
    spm_regm_std_temp_1(i)=nanstd(nanmean(mer_sp_temp_1(:,i:365:end),1));
    spm_regm_std_salt_1(i)=nanstd(nanmean(mer_sp_salt_1(:,i:365:end),1));
    spm_regm_std_po4_1(i)=nanstd(nanmean(mer_sp_po4_1(:,i:365:end),1));
    
    spm_regm_std_no3_b_1(i)=nanstd(nanmean(mer_sp_no3_b_1(:,i:365:end),1));
    spm_regm_std_nh4_b_1(i)=nanstd(nanmean(mer_sp_nh4_b_1(:,i:365:end),1));
    spm_regm_std_chl_b_1(i)=nanstd(nanmean(mer_sp_chl_b_1(:,i:365:end),1));
    spm_regm_std_temp_b_1(i)=nanstd(nanmean(mer_sp_temp_b_1(:,i:365:end),1));
    spm_regm_std_salt_b_1(i)=nanstd(nanmean(mer_sp_salt_b_1(:,i:365:end),1));
    spm_regm_std_po4_b_1(i)=nanstd(nanmean(mer_sp_po4_b_1(:,i:365:end),1));
    
%spatial std  
    %% 1st
    tm_regm_std_no3_1(i)=nanstd(nanmean(mer_sp_no3_1(:,i:365:end),2));
    tm_regm_std_nh4_1(i)=nanstd(nanmean(mer_sp_nh4_1(:,i:365:end),2));
    tm_regm_std_chl_1(i)=nanstd(nanmean(mer_sp_chl_1(:,i:365:end),2));
    tm_regm_std_temp_1(i)=nanstd(nanmean(mer_sp_temp_1(:,i:365:end),2));
    tm_regm_std_salt_1(i)=nanstd(nanmean(mer_sp_salt_1(:,i:365:end),2));
    tm_regm_std_po4_1(i)=nanstd(nanmean(mer_sp_po4_1(:,i:365:end),2));
    
    tm_regm_std_no3_b_1(i)=nanstd(nanmean(mer_sp_no3_b_1(:,i:365:end),2));
    tm_regm_std_nh4_b_1(i)=nanstd(nanmean(mer_sp_nh4_b_1(:,i:365:end),2));
    tm_regm_std_chl_b_1(i)=nanstd(nanmean(mer_sp_chl_b_1(:,i:365:end),2));
    tm_regm_std_temp_b_1(i)=nanstd(nanmean(mer_sp_temp_b_1(:,i:365:end),2));
    tm_regm_std_salt_b_1(i)=nanstd(nanmean(mer_sp_salt_b_1(:,i:365:end),2));
    tm_regm_std_po4_b_1(i)=nanstd(nanmean(mer_sp_po4_b_1(:,i:365:end),2));
    
end   

% %  dir_name ={'all_std','t_std','sp_std'};
cd(['D:\장기생태\Dynamic\KOEM\compare_sewer\1st_10yr\v5_out_bay\1st_10yr_',num2str(sig),'sig_',dir_name{case_num}]);

if case_num == 1 
    
    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines
    %data=NaN(1,365);        
    %data(1:length(c2004.gy_po4))=c2004.gy_po4_1; 
    %         plot(nanmean([c1997.gy_po4_1; c1998.gy_po4_1; c1999.gy_po4_1; c2000.gy_po4_1; c2001.gy_po4_1; c2002.gy_po4_1; c2003.gy_po4_1;],1),'b','linew',2);
            plot(nanmean(mer_po4_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_po4_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_po4_1;
            lower_bound_plt = temp_d - regm_std_po4_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_po4_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_po4_ft_1./30.973762 + clim1.obs_std_gy_po4_ft_1./30.973762;
            lower_bound_plt = clim1.obm_gy_po4_ft_1./30.973762 - clim1.obs_std_gy_po4_ft_1./30.973762;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_po4_ft_1./30.973762,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_po4_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_po4_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL po4']);
            xlabel('time','fontsize',13)
            ylabel('po4 (umol/m^3)','fontsize',13)
            grid on
 ylim([0 5])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4'),'-dpng')
   
    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines
    %data=NaN(1,365);        
    %data(1:length(c2004.gy_no3))=c2004.gy_no3_1; 
    %         plot(nanmean([c1997.gy_no3_1; c1998.gy_no3_1; c1999.gy_no3_1; c2000.gy_no3_1; c2001.gy_no3_1; c2002.gy_no3_1; c2003.gy_no3_1;],1),'b','linew',2);
            plot(nanmean(mer_no3_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_no3_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_no3_1;
            lower_bound_plt = temp_d - regm_std_no3_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_ft_1./14 + clim1.obs_std_gy_no3_ft_1./14;
            lower_bound_plt = clim1.obm_gy_no3_ft_1./14 - clim1.obs_std_gy_no3_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL NO3']);
            xlabel('time','fontsize',13)
            ylabel('NO3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4))=c2004.gy_nh4_1; 
            plot(nanmean(mer_nh4_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_nh4_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_nh4_1;
            lower_bound_plt = temp_d - regm_std_nh4_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_1; c1998.gy_nh4_1; c1999.gy_nh4_1; c2000.gy_nh4_1; c2001.gy_nh4_1; c2002.gy_nh4_1; c2003.gy_nh4_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_ft_1./14 + clim1.obs_std_gy_nh4_ft_1./14;
            lower_bound_plt = clim1.obm_gy_nh4_ft_1./14 - clim1.obs_std_gy_nh4_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    % data=NaN(1,366);        
    % data(1:length(c2004.gy_chl))=c2004.gy_chl_1; 
            plot(nanmean(mer_chl_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_chl_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_chl_1;
            lower_bound_plt = temp_d - regm_std_chl_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_1; c1998.gy_chl_1; c1999.gy_chl_1; c2000.gy_chl_1; c2001.gy_chl_1; c2002.gy_chl_1; c2003.gy_chl_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_ft_1 + clim1.obs_std_gy_chl_ft_1;
            lower_bound_plt = clim1.obm_gy_chl_ft_1 - clim1.obs_std_gy_chl_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp))=c2004.gy_temp_1;
            plot(nanmean(mer_temp_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_temp_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + regm_std_temp_1;
            lower_bound_plt = temp_d - regm_std_temp_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp_1; c1998.gy_temp_1; c1999.gy_temp_1; c2000.gy_temp_1; c2001.gy_temp_1; c2002.gy_temp_1; c2003.gy_temp_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_ft_1 + clim1.obs_std_gy_temp_ft_1;
            lower_bound_plt = clim1.obm_gy_temp_ft_1 - clim1.obs_std_gy_temp_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt))=c2004.gy_salt_1; 
            plot(nanmean(mer_salt_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_d
            salt_d=nanmean(mer_salt_1,1);
            nonan_data_plt = find(isnan(salt_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_d + regm_std_salt_1;
            lower_bound_plt = salt_d - regm_std_salt_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_1; c1998.gy_salt_1; c1999.gy_salt_1; c2000.gy_salt_1; c2001.gy_salt_1; c2002.gy_salt_1; c2003.gy_salt_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_ft_1 + clim1.obs_std_gy_salt_ft_1;
            lower_bound_plt = clim1.obm_gy_salt_ft_1 - clim1.obs_std_gy_salt_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')            
     %% bot
fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_po4_b))=c2004.gy_po4_b_1; 
            plot(nanmean(mer_po4_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p po4_b_d
            po4_b_d=nanmean(mer_po4_b_1,1);
            nonan_data_plt = find(isnan(po4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = po4_b_d + regm_std_po4_b_1;
            lower_bound_plt = po4_b_d - regm_std_po4_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_po4_b_1; c1998.gy_po4_b_1; c1999.gy_po4_b_1; c2000.gy_po4_b_1; c2001.gy_po4_b_1; c2002.gy_po4_b_1; c2003.gy_po4_b_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_po4_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 + clim1.obs_std_gy_po4_b_ft_1./30.973762;
            lower_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 - clim1.obs_std_gy_po4_b_ft_1./30.973762;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_po4_b_ft_1./30.973762,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_po4_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_po4_b_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL po4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('po4 (umol/m^3)','fontsize',13)
            grid on
 ylim([0 5])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4_bot'),'-dpng')

     fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_no3_b))=c2004.gy_no3_b_1; 
            plot(nanmean(mer_no3_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p no3_b_d
            no3_b_d=nanmean(mer_no3_b_1,1);
            nonan_data_plt = find(isnan(no3_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = no3_b_d + regm_std_no3_b_1;
            lower_bound_plt = no3_b_d - regm_std_no3_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_no3_b_1; c1998.gy_no3_b_1; c1999.gy_no3_b_1; c2000.gy_no3_b_1; c2001.gy_no3_b_1; c2002.gy_no3_b_1; c2003.gy_no3_b_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_b_ft_1./14 + clim1.obs_std_gy_no3_b_ft_1./14;
            lower_bound_plt = clim1.obm_gy_no3_b_ft_1./14 - clim1.obs_std_gy_no3_b_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_b_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_b_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
            xlabel('time','fontsize',13)
            ylabel('no3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4_b))=c2004.gy_nh4_b; 
            plot(nanmean(mer_nh4_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p nh4_b_d
            nh4_b_d=nanmean(mer_nh4_b_1,1);
            nonan_data_plt = find(isnan(nh4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = nh4_b_d + regm_std_nh4_b_1;
            lower_bound_plt = nh4_b_d - regm_std_nh4_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_b_ft_1./14 + clim1.obs_std_gy_nh4_b_ft_1./14;
            lower_bound_plt = clim1.obm_gy_nh4_b_ft_1./14 - clim1.obs_std_gy_nh4_b_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_b_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_b_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_chl_b))=c2004.gy_chl_b; 
            plot(nanmean(mer_chl_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p chl_b_d
            chl_b_d=nanmean(mer_chl_b_1,1);
            nonan_data_plt = find(isnan(chl_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = chl_b_d + regm_std_chl_b_1;
            lower_bound_plt = chl_b_d - regm_std_chl_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_b_ft_1 + clim1.obs_std_gy_chl_b_ft_1;
            lower_bound_plt = clim1.obm_gy_chl_b_ft_1 - clim1.obs_std_gy_chl_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl bot.']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp_b))=c2004.gy_temp_b; 
            plot(nanmean(mer_temp_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
            temp_b_d=nanmean(mer_temp_b_1,1);
            nonan_data_plt = find(isnan(temp_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_b_d + regm_std_temp_b_1;
            lower_bound_plt = temp_b_d - regm_std_temp_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_b_ft_1 + clim1.obs_std_gy_temp_b_ft_1;
            lower_bound_plt = clim1.obm_gy_temp_b_ft_1 - clim1.obs_std_gy_temp_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt_b))=c2004.gy_salt_b; 
            plot(nanmean(mer_salt_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
            salt_b_d=nanmean(mer_salt_b_1,1);
            nonan_data_plt = find(isnan(salt_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_b_d + regm_std_salt_b_1;
            lower_bound_plt = salt_b_d - regm_std_salt_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_b_ft_1 + clim1.obs_std_gy_salt_b_ft_1;
            lower_bound_plt = clim1.obm_gy_salt_b_ft_1 - clim1.obs_std_gy_salt_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng') 
            
            close all;
elseif case_num==2
         fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    %data(1:length(c2004.gy_po4))=c2004.gy_po4_1; 
    %         plot(nanmean([c1997.gy_po4_1; c1998.gy_po4_1; c1999.gy_po4_1; c2000.gy_po4_1; c2001.gy_po4_1; c2002.gy_po4_1; c2003.gy_po4_1;],1),'b','linew',2);
            plot(nanmean(mer_po4_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_po4_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_po4_1;
            lower_bound_plt = temp_d - spm_regm_std_po4_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_po4_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_po4_ft_1./30.973762 + clim1.spm_obs_std_gy_po4_ft_1./30.973762;
            lower_bound_plt = clim1.obm_gy_po4_ft_1./30.973762 - clim1.spm_obs_std_gy_po4_ft_1./30.973762;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_po4_ft_1./30.973762,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_po4_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_po4_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL po4']);
            xlabel('time','fontsize',13)
            ylabel('po4 (umol/m^3)','fontsize',13)
            grid on
 ylim([0 5])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4'),'-dpng')
    
    
        fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    %data(1:length(c2004.gy_no3))=c2004.gy_no3_1; 
    %         plot(nanmean([c1997.gy_no3_1; c1998.gy_no3_1; c1999.gy_no3_1; c2000.gy_no3_1; c2001.gy_no3_1; c2002.gy_no3_1; c2003.gy_no3_1;],1),'b','linew',2);
            plot(nanmean(mer_no3_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_no3_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_no3_1;
            lower_bound_plt = temp_d - spm_regm_std_no3_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_ft_1./14 + clim1.spm_obs_std_gy_no3_ft_1./14;
            lower_bound_plt = clim1.obm_gy_no3_ft_1./14 - clim1.spm_obs_std_gy_no3_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL NO3']);
            xlabel('time','fontsize',13)
            ylabel('NO3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4))=c2004.gy_nh4_1; 
            plot(nanmean(mer_nh4_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_nh4_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_nh4_1;
            lower_bound_plt = temp_d - spm_regm_std_nh4_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_1; c1998.gy_nh4_1; c1999.gy_nh4_1; c2000.gy_nh4_1; c2001.gy_nh4_1; c2002.gy_nh4_1; c2003.gy_nh4_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_ft_1./14 + clim1.spm_obs_std_gy_nh4_ft_1./14;
            lower_bound_plt = clim1.obm_gy_nh4_ft_1./14 - clim1.spm_obs_std_gy_nh4_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    % data=NaN(1,366);        
    % data(1:length(c2004.gy_chl))=c2004.gy_chl_1; 
            plot(nanmean(mer_chl_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_chl_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_chl_1;
            lower_bound_plt = temp_d - spm_regm_std_chl_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_1; c1998.gy_chl_1; c1999.gy_chl_1; c2000.gy_chl_1; c2001.gy_chl_1; c2002.gy_chl_1; c2003.gy_chl_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_ft_1 + clim1.spm_obs_std_gy_chl_ft_1;
            lower_bound_plt = clim1.obm_gy_chl_ft_1 - clim1.spm_obs_std_gy_chl_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp))=c2004.gy_temp_1;
            plot(nanmean(mer_temp_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_temp_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + spm_regm_std_temp_1;
            lower_bound_plt = temp_d - spm_regm_std_temp_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp_1; c1998.gy_temp_1; c1999.gy_temp_1; c2000.gy_temp_1; c2001.gy_temp_1; c2002.gy_temp_1; c2003.gy_temp_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_ft_1 + clim1.spm_obs_std_gy_temp_ft_1;
            lower_bound_plt = clim1.obm_gy_temp_ft_1 - clim1.spm_obs_std_gy_temp_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt))=c2004.gy_salt_1; 
            plot(nanmean(mer_salt_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_d
            salt_d=nanmean(mer_salt_1,1);
            nonan_data_plt = find(isnan(salt_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_d + spm_regm_std_salt_1;
            lower_bound_plt = salt_d - spm_regm_std_salt_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_1; c1998.gy_salt_1; c1999.gy_salt_1; c2000.gy_salt_1; c2001.gy_salt_1; c2002.gy_salt_1; c2003.gy_salt_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_ft_1 + clim1.spm_obs_std_gy_salt_ft_1;
            lower_bound_plt = clim1.obm_gy_salt_ft_1 - clim1.spm_obs_std_gy_salt_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')            
     %% bot
 fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_po4_b))=c2004.gy_po4_b_1; 
            plot(nanmean(mer_po4_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p po4_b_d
            po4_b_d=nanmean(mer_po4_b_1,1);
            nonan_data_plt = find(isnan(po4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = po4_b_d + spm_regm_std_po4_b_1;
            lower_bound_plt = po4_b_d - spm_regm_std_po4_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_po4_b_1; c1998.gy_po4_b_1; c1999.gy_po4_b_1; c2000.gy_po4_b_1; c2001.gy_po4_b_1; c2002.gy_po4_b_1; c2003.gy_po4_b_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_po4_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 + clim1.spm_obs_std_gy_po4_b_ft_1./30.973762;
            lower_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 - clim1.spm_obs_std_gy_po4_b_ft_1./30.973762;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_po4_b_ft_1./30.973762,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_po4_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_po4_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL po4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('po4 (umol/m^3)','fontsize',13)
            grid on
 ylim([0 5])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4_bot'),'-dpng')

     
     fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_no3_b))=c2004.gy_no3_b_1; 
            plot(nanmean(mer_no3_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p no3_b_d
            no3_b_d=nanmean(mer_no3_b_1,1);
            nonan_data_plt = find(isnan(no3_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = no3_b_d + spm_regm_std_no3_b_1;
            lower_bound_plt = no3_b_d - spm_regm_std_no3_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_no3_b_1; c1998.gy_no3_b_1; c1999.gy_no3_b_1; c2000.gy_no3_b_1; c2001.gy_no3_b_1; c2002.gy_no3_b_1; c2003.gy_no3_b_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_b_ft_1./14 + clim1.spm_obs_std_gy_no3_b_ft_1./14;
            lower_bound_plt = clim1.obm_gy_no3_b_ft_1./14 - clim1.spm_obs_std_gy_no3_b_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_b_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
            xlabel('time','fontsize',13)
            ylabel('no3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4_b))=c2004.gy_nh4_b; 
            plot(nanmean(mer_nh4_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p nh4_b_d
            nh4_b_d=nanmean(mer_nh4_b_1,1);
            nonan_data_plt = find(isnan(nh4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = nh4_b_d + spm_regm_std_nh4_b_1;
            lower_bound_plt = nh4_b_d - spm_regm_std_nh4_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_b_ft_1./14 + clim1.spm_obs_std_gy_nh4_b_ft_1./14;
            lower_bound_plt = clim1.obm_gy_nh4_b_ft_1./14 - clim1.spm_obs_std_gy_nh4_b_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_b_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_b_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_chl_b))=c2004.gy_chl_b; 
            plot(nanmean(mer_chl_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p chl_b_d
            chl_b_d=nanmean(mer_chl_b_1,1);
            nonan_data_plt = find(isnan(chl_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = chl_b_d + spm_regm_std_chl_b_1;
            lower_bound_plt = chl_b_d - spm_regm_std_chl_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_b_ft_1 + clim1.spm_obs_std_gy_chl_b_ft_1;
            lower_bound_plt = clim1.obm_gy_chl_b_ft_1 - clim1.spm_obs_std_gy_chl_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl bot.']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp_b))=c2004.gy_temp_b; 
            plot(nanmean(mer_temp_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
            temp_b_d=nanmean(mer_temp_b_1,1);
            nonan_data_plt = find(isnan(temp_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_b_d + spm_regm_std_temp_b_1;
            lower_bound_plt = temp_b_d - spm_regm_std_temp_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_b_ft_1 + clim1.spm_obs_std_gy_temp_b_ft_1;
            lower_bound_plt = clim1.obm_gy_temp_b_ft_1 - clim1.spm_obs_std_gy_temp_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt_b))=c2004.gy_salt_b; 
            plot(nanmean(mer_salt_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
            salt_b_d=nanmean(mer_salt_b_1,1);
            nonan_data_plt = find(isnan(salt_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_b_d + spm_regm_std_salt_b_1;
            lower_bound_plt = salt_b_d - spm_regm_std_salt_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_b_ft_1 + clim1.spm_obs_std_gy_salt_b_ft_1;
            lower_bound_plt = clim1.obm_gy_salt_b_ft_1 - clim1.spm_obs_std_gy_salt_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng')
            
            close all;
elseif case_num == 3

fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    %data(1:length(c2004.gy_po4))=c2004.gy_po4_1; 
    %         plot(nanmean([c1997.gy_po4_1; c1998.gy_po4_1; c1999.gy_po4_1; c2000.gy_po4_1; c2001.gy_po4_1; c2002.gy_po4_1; c2003.gy_po4_1;],1),'b','linew',2);
            plot(nanmean(mer_po4_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_po4_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_po4_1;
            lower_bound_plt = temp_d - tm_regm_std_po4_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_po4_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_po4_ft_1./30.973762 + clim1.tm_obs_std_gy_po4_ft_1./30.973762;
            lower_bound_plt = clim1.obm_gy_po4_ft_1./30.973762 - clim1.tm_obs_std_gy_po4_ft_1./30.973762;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_po4_ft_1./30.973762,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_po4_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_po4_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL po4']);
            xlabel('time','fontsize',13)
            ylabel('po4 (umol/m^3)','fontsize',13)
            grid on
 ylim([0 5])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4'),'-dpng')

        fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    %data(1:length(c2004.gy_no3))=c2004.gy_no3_1; 
    %         plot(nanmean([c1997.gy_no3_1; c1998.gy_no3_1; c1999.gy_no3_1; c2000.gy_no3_1; c2001.gy_no3_1; c2002.gy_no3_1; c2003.gy_no3_1;],1),'b','linew',2);
            plot(nanmean(mer_no3_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_no3_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_no3_1;
            lower_bound_plt = temp_d - tm_regm_std_no3_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);

            clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_ft_1./14 + clim1.tm_obs_std_gy_no3_ft_1./14;
            lower_bound_plt = clim1.obm_gy_no3_ft_1./14 - clim1.tm_obs_std_gy_no3_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL NO3']);
            xlabel('time','fontsize',13)
            ylabel('NO3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4))=c2004.gy_nh4_1; 
            plot(nanmean(mer_nh4_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_nh4_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_nh4_1;
            lower_bound_plt = temp_d - tm_regm_std_nh4_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_1; c1998.gy_nh4_1; c1999.gy_nh4_1; c2000.gy_nh4_1; c2001.gy_nh4_1; c2002.gy_nh4_1; c2003.gy_nh4_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_ft_1./14 + clim1.tm_obs_std_gy_nh4_ft_1./14;
            lower_bound_plt = clim1.obm_gy_nh4_ft_1./14 - clim1.tm_obs_std_gy_nh4_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    % data=NaN(1,366);        
    % data(1:length(c2004.gy_chl))=c2004.gy_chl_1; 
            plot(nanmean(mer_chl_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_chl_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_chl_1;
            lower_bound_plt = temp_d - tm_regm_std_chl_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_1; c1998.gy_chl_1; c1999.gy_chl_1; c2000.gy_chl_1; c2001.gy_chl_1; c2002.gy_chl_1; c2003.gy_chl_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_ft_1 + clim1.tm_obs_std_gy_chl_ft_1;
            lower_bound_plt = clim1.obm_gy_chl_ft_1 - clim1.tm_obs_std_gy_chl_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp))=c2004.gy_temp_1;
            plot(nanmean(mer_temp_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_d
            temp_d=nanmean(mer_temp_1,1);
            nonan_data_plt = find(isnan(temp_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_d + tm_regm_std_temp_1;
            lower_bound_plt = temp_d - tm_regm_std_temp_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_temp_1; c1998.gy_temp_1; c1999.gy_temp_1; c2000.gy_temp_1; c2001.gy_temp_1; c2002.gy_temp_1; c2003.gy_temp_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_ft_1 + clim1.tm_obs_std_gy_temp_ft_1;
            lower_bound_plt = clim1.obm_gy_temp_ft_1 - clim1.tm_obs_std_gy_temp_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt))=c2004.gy_salt_1; 
            plot(nanmean(mer_salt_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_d
            salt_d=nanmean(mer_salt_1,1);
            nonan_data_plt = find(isnan(salt_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_d + tm_regm_std_salt_1;
            lower_bound_plt = salt_d - tm_regm_std_salt_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_1; c1998.gy_salt_1; c1999.gy_salt_1; c2000.gy_salt_1; c2001.gy_salt_1; c2002.gy_salt_1; c2003.gy_salt_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_ft_1 + clim1.tm_obs_std_gy_salt_ft_1;
            lower_bound_plt = clim1.obm_gy_salt_ft_1 - clim1.tm_obs_std_gy_salt_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
    print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt'),'-dpng')            
     %% bot
 fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_po4_b))=c2004.gy_po4_b_1; 
            plot(nanmean(mer_po4_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p po4_b_d
            po4_b_d=nanmean(mer_po4_b_1,1);
            nonan_data_plt = find(isnan(po4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = po4_b_d + tm_regm_std_po4_b_1;
            lower_bound_plt = po4_b_d - tm_regm_std_po4_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_po4_b_1; c1998.gy_po4_b_1; c1999.gy_po4_b_1; c2000.gy_po4_b_1; c2001.gy_po4_b_1; c2002.gy_po4_b_1; c2003.gy_po4_b_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_po4_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 + clim1.tm_obs_std_gy_po4_b_ft_1./30.973762;
            lower_bound_plt = clim1.obm_gy_po4_b_ft_1./30.973762 - clim1.tm_obs_std_gy_po4_b_ft_1./30.973762;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_po4_b_ft_1./30.973762,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_po4_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_po4_b_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL po4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('po4 (umol/m^3)','fontsize',13)
            grid on
 ylim([0 5])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_po4_bot'),'-dpng')
     
     fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_no3_b))=c2004.gy_no3_b_1; 
            plot(nanmean(mer_no3_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p no3_b_d
            no3_b_d=nanmean(mer_no3_b_1,1);
            nonan_data_plt = find(isnan(no3_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = no3_b_d + tm_regm_std_no3_b_1;
            lower_bound_plt = no3_b_d - tm_regm_std_no3_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_no3_b_1; c1998.gy_no3_b_1; c1999.gy_no3_b_1; c2000.gy_no3_b_1; c2001.gy_no3_b_1; c2002.gy_no3_b_1; c2003.gy_no3_b_1;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_no3_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_no3_b_ft_1./14 + clim1.tm_obs_std_gy_no3_b_ft_1./14;
            lower_bound_plt = clim1.obm_gy_no3_b_ft_1./14 - clim1.tm_obs_std_gy_no3_b_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_no3_b_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_no3_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_no3_b_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL no3 bot.']);
            xlabel('time','fontsize',13)
            ylabel('no3 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 70])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_no3_bot'),'-dpng')

    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_nh4_b))=c2004.gy_nh4_b; 
            plot(nanmean(mer_nh4_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p nh4_b_d
            nh4_b_d=nanmean(mer_nh4_b_1,1);
            nonan_data_plt = find(isnan(nh4_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = nh4_b_d + tm_regm_std_nh4_b_1;
            lower_bound_plt = nh4_b_d - tm_regm_std_nh4_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b; c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_nh4_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_nh4_b_ft_1./14 + clim1.tm_obs_std_gy_nh4_b_ft_1./14;
            lower_bound_plt = clim1.obm_gy_nh4_b_ft_1./14 - clim1.tm_obs_std_gy_nh4_b_ft_1./14;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_nh4_b_ft_1./14,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_nh4_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_nh4_b_ft_1),lower_bound_plt,'r-','linew',2);


            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL nh4 bot.']);
            xlabel('time','fontsize',13)
            ylabel('nh4 (umol/m^3)','fontsize',13)
            grid on
            ylim([0 14])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_nh4_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_chl_b))=c2004.gy_chl_b; 
            plot(nanmean(mer_chl_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p chl_b_d
            chl_b_d=nanmean(mer_chl_b_1,1);
            nonan_data_plt = find(isnan(chl_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = chl_b_d + tm_regm_std_chl_b_1;
            lower_bound_plt = chl_b_d - tm_regm_std_chl_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b; c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_chl_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_chl_b_ft_1 + clim1.tm_obs_std_gy_chl_b_ft_1;
            lower_bound_plt = clim1.obm_gy_chl_b_ft_1 - clim1.tm_obs_std_gy_chl_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_chl_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_chl_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_chl_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL chl bot.']);
            xlabel('time','fontsize',13)
            ylabel('chl (umol/m^3)','fontsize',13)
            grid on
            ylim([0 15])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_chl_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_temp_b))=c2004.gy_temp_b; 
            plot(nanmean(mer_temp_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p temp_b_d
            temp_b_d=nanmean(mer_temp_b_1,1);
            nonan_data_plt = find(isnan(temp_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = temp_b_d + tm_regm_std_temp_b_1;
            lower_bound_plt = temp_b_d - tm_regm_std_temp_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    % plot(nanmean([c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b; c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_temp_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_temp_b_ft_1 + clim1.tm_obs_std_gy_temp_b_ft_1;
            lower_bound_plt = clim1.obm_gy_temp_b_ft_1 - clim1.tm_obs_std_gy_temp_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_temp_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_temp_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_temp_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL temp bot.']);
            xlabel('time','fontsize',13)
            ylabel('temp (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)

            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_temp_bot'),'-dpng')


    fig = figure; hold on;
    %         plot(zeros(365,1),'g--','linew',2); % zero lines

    %data=NaN(1,365);        
    % data(1:length(c2004.gy_salt_b))=c2004.gy_salt_b; 
            plot(nanmean(mer_salt_b_1,1),'b','linew',2);
            clearvars *_bound_plt nonan_data_plt discon_p salt_b_d
            salt_b_d=nanmean(mer_salt_b_1,1);
            nonan_data_plt = find(isnan(salt_b_d)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            upper_bound_plt = salt_b_d + tm_regm_std_salt_b_1;
            lower_bound_plt = salt_b_d - tm_regm_std_salt_b_1;
            patch([nonan_data_plt fliplr(nonan_data_plt)], [upper_bound_plt fliplr(lower_bound_plt)],[0.5 0.5 0.5]);
    %         plot(nanmean([c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b; c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b;],1),'b','linew',2);
                   clearvars *_bound_plt nonan_data_plt discon_p
            nonan_data_plt = find(isnan(clim1.obm_gy_salt_b_ft_1)==0);
            nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
            discon_p = find(nonan_diff ~= 1); % find discontinuous point
            discon_p(4) = length(nonan_data_plt); %end point
            upper_bound_plt = clim1.obm_gy_salt_b_ft_1 + clim1.tm_obs_std_gy_salt_b_ft_1;
            lower_bound_plt = clim1.obm_gy_salt_b_ft_1 - clim1.tm_obs_std_gy_salt_b_ft_1;
            for i = 1:4
            if i == 1
            patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],'r');
            else
            patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'r');
            end
            end 
            plot(clim1.obm_gy_salt_b_ft_1,'r-','linew',2); 
            plot(1:length(clim1.obm_gy_salt_b_ft_1),upper_bound_plt,'r-','linew',2);
            plot(1:length(clim1.obm_gy_salt_b_ft_1),lower_bound_plt,'r-','linew',2);

            alpha(0.3) %transparency      

            title(['GY KOEM OBS vs. daily MODEL salt bot.']);
            xlabel('time','fontsize',13)
            ylabel('salt (umol/m^3)','fontsize',13)
            grid on
            ylim([0 35])
            xlim([1 365])
            set(gca,'fontsize',13)
    %         %xticks(com_eom(1:12:end))
    %         %xticklabels(com_t_tic)
    %         %xtickangle(90)
            print(fig,strcat('multi_yr_KOEM_OBS_vs_daily_MODEL_salt_bot'),'-dpng') 
            
            close all;
end

    end
end
       